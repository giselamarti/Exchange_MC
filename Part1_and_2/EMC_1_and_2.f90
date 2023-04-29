PROGRAM EMC_Ising_spin_glass

    implicit none

    ! Define constants and variables
    integer :: i, j, k, t
    integer :: status
    integer, parameter :: dp = kind(0.0d0)
    integer, parameter :: L = 20 ! linear size of the system
    integer, parameter :: N_T = L*L 
    integer, parameter :: n_sw = 10
    integer, parameter :: n_meas = 10
    integer, parameter :: n_MCS = 10000 ! number of MC steps
    integer, parameter :: seed = 1234
    real(dp), parameter :: T_max = 3.0_dp ! maximum temperature (changed to 500 to compare with Ferdinand)
    real(dp), parameter :: T_min = 1.0_dp ! minimum temperature
    real(dp) :: Temp(N_T)
    real(dp) :: S_a_tmp(L,L)
    real(dp) :: S_b_tmp(L,L)
    real(dp) :: M(L,L) ! Random varibales with Gaussian distribution
    real(dp), dimension(:,:,:), allocatable :: S_a, S_b ! spin configuration
    real(dp), dimension(:,:), allocatable :: Q, E ! Energy and Overlaps
    real(dp) :: dE, dT
    logical :: accept
    character*20 :: filename
    
    filename = "gaussian_numbers.txt"
    open(2, file=filename, status='old', action='read')
    do i = 1, L
        do j = 1, L
            read(2, *) M(i,j)
        end do
    end do
    close(2)

    ! allocate spin configuration and initialize to random values
    allocate(S_a(N_T,L,L))
    allocate(S_b(N_T,L,L))    
    call random_number(S_a)
    call random_number(S_b)
    S_a = sign(1.0_dp, S_a-0.5_dp)
    S_b = sign(1.0_dp, S_b-0.5_dp)

    ! Initialize temperatures
    Temp = [(T_min + (i-1) * (T_max - T_min) / (N_T-1), i = 1, N_T)]

    ! Initialize arrays to store overlaps and energies
    allocate(Q(N_T, n_MCS/n_meas))
    allocate(E(N_T, n_MCS/n_meas))

    ! Open output file
    open(1, file = "output_all_data.txt", status="replace") ! temperature, overlap and energy values
    open(2, file = "energy_output.txt", status = "replace") ! only energies

    ! Perform Monte Carlo simulation
    do t = 1, n_MCS
        ! Perform Metropolis updates
        do k = 1, N_T
            ! Choose random spin to flip
            i = int(rand() * L) + 1
            j = int(rand() * L) + 1

            ! Compute energy difference
            dE = 2.0 * M(i,j) * S_a(k, i, j) * (S_a(k, mod(i, L) + 1, j) + S_a(k, i, mod(j, L) + 1))
        
            ! Accept or reject update
            accept = rand() < exp(-dE / Temp(k))
            if (accept) S_a(k, i, j) = -S_a(k, i, j)

            ! Repeat for second set of replicas
            i = int(rand() * L) + 1
            j = int(rand() * L) + 1
            dE = 2.0 * M(i, j) * S_b(k, i, j) * (S_b(k, mod(i, L) + 1, j) + S_b(k, i, mod(j, L) + 1))
            accept = rand() < exp(-dE / Temp(k))
            if (accept) S_b(k, i, j) = -S_b(k, i, j)
        end do ! end of metropolis updates

        ! Attempt replica exchanges
        if (mod(t, n_sw) == 0) then
            do k = 1, N_T-1, 2
                dE = (Hamiltonian(S_a(k+1,:,:), M) + Hamiltonian(S_b(k,:,:), M)) - &
                     (Hamiltonian(S_a(k,:,:), M) + Hamiltonian(S_b(k+1,:,:), M))
                dT = Temp(k+1) - Temp(k)
                if (rand() < exp(-dE / dT)) then
                    S_a_tmp = S_a(k+1,:,:)
                    S_b_tmp = S_b(k+1,:,:)
                    S_a(k+1,:,:) = S_a(k,:,:)
                    S_b(k+1,:,:) = S_b(k,:,:)
                    S_a(k,:,:) = S_a_tmp
                    S_b(k,:,:) = S_b_tmp
                end if
            end do

            do k = 2, N_T-1, 2
                dE = (Hamiltonian(S_a(k-1,:,:), M) + Hamiltonian(S_b(k,:,:), M)) - &
                     (Hamiltonian(S_a(k,:,:), M) + Hamiltonian(S_b(k-1,:,:), M))
                dT = Temp(k-1) - Temp(k)
                if (rand() < exp(-dE / dT)) then
                    S_a_tmp = S_a(k-1,:,:)
                    S_b_tmp = S_b(k-1,:,:)
                    S_a(k-1,:,:) = S_a(k,:,:)
                    S_b(k-1,:,:) = S_b(k,:,:)
                    S_a(k,:,:) = S_a_tmp
                    S_b(k,:,:) = S_b_tmp
                end if
            end do           
        end if ! end attempt replica exchanges

        ! Measure overlaps and energies
        if (mod(t, n_meas) == 0) then
            do k = 1, N_T
                Q(k, t/n_meas) = sum(S_a(k,:,:) * S_b(k,:,:)) / L**2
                E(k, t/n_meas) = -(0.5*(Hamiltonian(S_a(k,:,:), M) + Hamiltonian(S_b(k,:,:), M)) / L**2)
            end do

            ! Write results to output file
            write(1, '(a,i0,a)') "Step ", t, ":"
            do k = 1, N_T
                write(1, '(a,f6.3,a,f6.3,a,f6.3)', advance='yes') "Temperature ", Temp(k), ": Q = ", Q(k, t/n_meas), &
                 ", E = ", E(k, t/n_meas)
                write(2,*) E(k, t/n_meas)
            end do
            write(1, '(a)', advance='yes') " " 
        end if ! end measure overlaps and energies

        ! Print progress
        if (mod(t, 1000) == 0) write(*,*) "Step = ", t, "/ 10000"

    end do ! end of Monte Carlo simulation

    close(1)
    close(2)

    contains

    ! Hamiltonian function
    real function Hamiltonian(S, M) !result(H)
        integer, parameter :: dp = kind(0.0d0)
        integer :: i, j, L
        real(dp), dimension(:,:) :: S, M
        real(dp) :: H

        ! Input: S (spin configuration) and M (random variables) matrices
        ! Output: H matrix (Hamiltonian for a determinate configuration)

        L = size(S,1)
        H = 0.0d0
        do i = 1, L
            do j = 1, L
                H = H - M(i,j) * S(i,j) * (S(mod(i,L)+1,j) + S(i,mod(j,L)+1))
            end do
        end do
        Hamiltonian = H
    end function Hamiltonian

END