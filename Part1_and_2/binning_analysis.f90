! -------------------------------------------------------------------

!   Program for binning analysis

!   Takes as an input file the list of energy values for each temperature and each step "energy_ouput.txt"

!   Creates an output file with the temperature, mean energy, error "mean_e_and_error.txt"

! -------------------------------------------------------------------
PROGRAM binning_analysis

  implicit none
  
  integer, parameter :: N_T = 400, L = 20
  real, parameter :: T_min = 1.0, T_max = 500.0
  integer :: k, n, n_bins, min_bin_size = 10, max_bin_size, bin_size
  integer :: i, j
  real, dimension(N_T) :: T, E_mean, E_error
  real, dimension(N_T, 1000) :: E
  real,dimension(:, :), allocatable :: E_bins
  character(100) :: filename = "energy_output.txt"
  real, dimension(:, :), allocatable :: data


  ! Generate T array
  T = [(T_min + (i-1) * (T_max - T_min) / (N_T-1), i = 1, N_T)]

  ! Read data from file
  open(unit=10, file=filename, status='old')
  allocate(data(400000,1))
  read(10,*) data
  close(10)
  
  ! Reshape data to 2D array
  do i = 1, N_T
    do j = 1, 1000
      E(i,j) = data((i-1)*1000+j, 1)
    end do
  end do

  open(5, file="mean_e_and_error.txt", status = "replace")

  ! Calculate E_mean and E_error
  max_bin_size = size(E, 2) / 2
  do k = 1, N_T
    do n = 1, max_bin_size - min_bin_size + 1
      bin_size = min_bin_size + n - 1
      n_bins = size(E(k,:)) / bin_size
      allocate(E_bins(n_bins, bin_size))
      do i = 1, n_bins
        do j = 1, bin_size
          E_bins(i,j) = E(k,(i-1)*bin_size+j)
        end do
      end do
      E_mean(k) = sum(E_bins(:,1))/n_bins
      do i = 2, bin_size
        E_mean(k) = E_mean(k) + sum(E_bins(:,i))/n_bins
        E_error(k) = E_error(k) + sum((E_bins(:,i) - E_mean(k))**2)/n_bins
      end do
      !E_error(k) = sqrt(E_error(k)/n_bins) / sqrt(n_bins)
      E_error(k) = sqrt(real(E_error(k), kind=8) / n_bins) / sqrt(real(n_bins, kind=8))

      deallocate(E_bins)
      if (n_bins == 1) then
        exit
      end if
    end do
  end do
  
  do k = 1, N_T
    write(5, *) T(k), E_mean(k)/N_T, E_error(k)/(N_T*1000)
  end do

  close(5)

END