PROGRAM EMC_Ising_spin_glass

    implicit none

    ! Define constants and varibales
    integer :: i, j, k, t, p, sample, a, b
    integer, parameter :: dp = kind(0.0d0)
    integer, parameter :: L = 4 ! linear size of the system
    integer, parameter :: N_T = L*L*L
    integer, parameter :: n_sw = 10
    integer, parameter :: n_meas = 10
    integer, parameter :: n_MCS = 10000 ! number of MC steps
    integer, parameter :: n_samples = 1000
    integer, parameter :: last_steps = 7000
    real(dp), parameter :: T_max = 2.0_dp ! maximum temperature
    real(dp), parameter :: T_min = 0.2_dp ! minimum temperature
    real(dp) :: Temp(N_T)
    real(dp) :: S_a_tmp(L,L,L)
    real(dp) :: S_b_tmp(L,L,L)
    real, dimension(700000) :: q02, q05   ! new array to store Q values
    real(dp) :: M(L,L,L)
    real(dp) :: M_samples(n_samples,L,L,L) ! Random varibales with Gaussian distribution
    real(dp), dimension(:,:,:,:), allocatable :: S_a, S_b ! spin configuration
    real(dp), dimension(:,:), allocatable :: Q, E ! Energy and Overlaps
    real(dp) :: dE, dT
    logical :: accept
    logical :: found_temp
    integer :: status
    character*20 :: filename

    filename = "gaussian_numbers.txt"
    open(2, file=filename, status='old', action='read')
    do i = 1, L
        do j = 1, L
            do p = 1, L
                read(2, *) M_samples(n_samples,i,j,p)
            end do
        end do
    end do
    close(2)

    ! allocate spin configuration and initialize to random values
    allocate(S_a(N_T,L,L,L))
    allocate(S_b(N_T,L,L,L))    
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
    open(3, file = "histogram_data_0.2.txt", status="replace")
    open(4, file = "histogram_data_0.5.txt", status="replace")

    ! Loop over all samples
    do sample = 1, n_samples
        M = M_samples(sample,:,:,:)

        ! Perform Monte Carlo simulation
        do t = 1, n_MCS
            ! Perform Metropolis updates
            do k = 1, N_T
                ! Choose random spin to flip
                i = int(rand() * L) + 1
                j = int(rand() * L) + 1
                p = int(rand() * L) + 1

                ! Compute energy difference
                dE = 2.0 * M(i,j,p) * S_a(k,i,j,p) * (S_a(k,mod(i,L)+1,j,p) + S_a(k,mod(i,L)-1,j,p) &
                                                    + S_a(k,i,mod(j,L)+1,p) + S_a(k,i,mod(j,L)-1,p) &
                                                    + S_a(k,i,j,mod(p,L)+1) + S_a(k,i,j,mod(p,L)-1))
        
                ! Accept or reject update
                accept = rand() < exp(-dE / Temp(k))
                if (accept) S_a(k,i,j,p) = -S_a(k,i,j,p)

                ! Repeat for second set of replicas
                i = int(rand() * L) + 1
                j = int(rand() * L) + 1
                p = int(rand() * L) + 1

                dE = 2.0 * M(i,j,p) * S_b(k,i,j,p) * (S_b(k,mod(i,L)+1,j,p) + S_b(k,mod(i,L)-1,j,p) &
                                                    + S_b(k,i,mod(j,L)+1,p) + S_b(k,i,mod(j,L)-1,p) &
                                                    + S_b(k,i,j,mod(p,L)+1) + S_b(k,i,j,mod(p,L)-1))

                accept = rand() < exp(-dE / Temp(k))
                if (accept) S_b(k,i,j,p) = -S_b(k,i,j,p)
            end do ! end of metropolis updates

            ! Attempt replica exchanges
            if (mod(t, n_sw) == 0) then
                do k = 1, N_T-1, 2
                    dE = (Hamiltonian(S_a(k+1,:,:,:), M) + Hamiltonian(S_b(k,:,:,:), M)) - &
                        (Hamiltonian(S_a(k,:,:,:), M) + Hamiltonian(S_b(k+1,:,:,:), M))
                    dT = Temp(k+1) - Temp(k)
                    if (rand() < exp(-dE / dT)) then
                        S_a_tmp = S_a(k+1,:,:,:)
                        S_b_tmp = S_b(k+1,:,:,:)
                        S_a(k+1,:,:,:) = S_a(k,:,:,:)
                        S_b(k+1,:,:,:) = S_b(k,:,:,:)
                        S_a(k,:,:,:) = S_a_tmp
                        S_b(k,:,:,:) = S_b_tmp
                    end if
                end do

                do k = 2, N_T-1, 2
                    dE = (Hamiltonian(S_a(k-1,:,:,:), M) + Hamiltonian(S_b(k,:,:,:), M)) - &
                        (Hamiltonian(S_a(k,:,:,:), M) + Hamiltonian(S_b(k-1,:,:,:), M))
                    dT = Temp(k-1) - Temp(k)
                    if (rand() < exp(-dE / dT)) then
                        S_a_tmp = S_a(k-1,:,:,:)
                        S_b_tmp = S_b(k-1,:,:,:)
                        S_a(k-1,:,:,:) = S_a(k,:,:,:)
                        S_b(k-1,:,:,:) = S_b(k,:,:,:)
                        S_a(k,:,:,:) = S_a_tmp
                        S_b(k,:,:,:) = S_b_tmp
                    end if
                end do           
            end if ! end attempt replica exchanges

            ! Measure overlaps and energies
            if (mod(t, n_meas) == 0) then
                found_temp = .false.
                do k = 1, N_T
                    Q(k, t/n_meas) = sum(S_a(k,:,:,:) * S_b(k,:,:,:)) / L**3
                    E(k, t/n_meas) = -(0.5*(Hamiltonian(S_a(k,:,:,:), M) + Hamiltonian(S_b(k,:,:,:), M)) / L**3)
                    
                    ! Store Q values in new array q02 for T=0.2
                    if (t > n_MCS - last_steps .and. Temp(k) == 0.2_dp) then
                        i = (sample-1)*last_steps/N_T + (t-n_MCS+last_steps)/n_meas
                        q02(i) = Q(k, t/n_meas)/N_T
                        write(3,*) q02(i)
                    end if

                    ! Write Q values in new array q05 only for the last 7000 steps and first T=0.5
                    if (t > n_MCS - last_steps .and. abs(Temp(k) - 0.5) < 0.015 .and. .not. found_temp) then
                        found_temp = .true.
                        j = (sample-1)*last_steps/N_T + (t-n_MCS+last_steps)/n_meas
                        q05(j) = Q(k, t/n_meas)/N_T
                        write(4,*) q05(j)
                    end if
                    ! Exit loop over k if temperature already found
                    if (found_temp) exit
                end do
            end if ! end measure overlaps and energies
        end do ! end of Monte Carlo simulation
    end do ! end of loop over samples

    close(3)
    close(4)

    contains

    ! Hamiltonian function
    real function Hamiltonian(S, M) !result(H)
        integer, parameter :: dp = kind(0.0d0)
        integer :: i, j, p, L
        real(dp), dimension(:,:,:) :: S, M
        real(dp) :: H

        ! Input: S (spin configuration) and M (random variables) matrices
        ! Output: H matrix (Hamiltonian for a determinate configuration)

        L = size(S,1)
        H = 0.0d0
        do i = 1, L
            do j = 1, L
                do p = 1, L
                    H = H - M(i,j,p) * S(i,j,p) * (S(mod(i,L)+1,j,p) + S(mod(i,L)-1,j,p) + &
                                                   S(i,mod(j,L)+1,p) + S(i,mod(j,L)-1,p) + &
                                                   S(i,j,mod(p,L)+1) + S(i,j,mod(p,L)-1))
                end do
            end do
        end do
        Hamiltonian = H
    end function Hamiltonian

END