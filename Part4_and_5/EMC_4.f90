PROGRAM EMC_Ising_spin_glass_part4

    implicit none

    ! Define constants and variables
    integer :: i, j, k, t, p, sample
    integer, parameter :: dp = kind(0.0d0)
    integer, parameter :: L = 8 ! linear size of the system
    integer, parameter :: N_T = L*L*L
    integer, parameter :: n_sw = 10
    integer, parameter :: n_meas = 10
    integer, parameter :: n_MCS = 1000000 ! number of MC steps
    real(dp), parameter :: T_max = 2.0_dp ! maximum temperature
    real(dp), parameter :: T_min = 0.2_dp ! minimum temperature
    real(dp) :: Temp(N_T)
    real(dp) :: S_a_tmp(L,L,L)
    real(dp) :: S_b_tmp(L,L,L)    
    real(dp) :: M(L,L,L)
    real(dp), dimension(:,:,:,:), allocatable :: S_a, S_b ! spin configuration
    real(dp), dimension(:,:), allocatable :: Q, E ! Energy and Overlaps
    real(dp) :: sum_q1(N_T), sum_e1(N_T), sum_q2(N_T), sum_e2(N_T), sum_q3(N_T), sum_e3(N_T), sum_q4(N_T), sum_e4(N_T)
    real(dp) :: mean_q1(N_T), mean_e1(N_T), mean_q2(N_T), mean_e2(N_T), mean_q3(N_T), mean_e3(N_T), mean_q4(N_T), mean_e4(N_T)
    integer :: count1, count2, count3, count4
    real(dp) :: dE, dT
    logical :: accept
    integer :: status
    character*20 :: filename

    filename = "gaussian_numbers.txt"
    open(2, file=filename, status='old', action='read')
    do i = 1, L
        do j = 1, L
            do p = 1, L
                read(2, *) M(i,j,p)! (8,8,8)
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
    open(3, file = "4_mean_values.txt", status="replace")

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
            do k = 1, N_T
                Q(k, t/n_meas) = sum(S_a(k,:,:,:) * S_b(k,:,:,:)) / L**3
                E(k, t/n_meas) = -(0.5*(Hamiltonian(S_a(k,:,:,:), M) + Hamiltonian(S_b(k,:,:,:), M)) / L**3)
            end do

            ! Check if the current time interval is within the specified range
            if (t/n_meas >= 1000 .and. t/n_meas <= 3000) then !10.000-30.000 MCS
                do k = 1, N_T
                    sum_q1(k) = sum_q1(k) + (Q(k, t/n_meas)/N_T)**2
                    sum_e1(k) = sum_e1(k) + E(k, t/n_meas)
                end do
                count1 = count1 + 1
            else if (t/n_meas >= 3000 .and. t/n_meas <= 10000) then !30.000-100.000 MCS
                do k = 1, N_T
                    sum_q2(k) = sum_q2(k) + (Q(k, t/n_meas)/N_T)**2
                    sum_e2(k) = sum_e2(k) + E(k, t/n_meas)
                end do
                count2 = count2 + 1
            else if (t/n_meas >= 10001 .and. t/n_meas <= 30000) then !100.000-300.000 MCS
                do k = 1, N_T
                    sum_q3(k) = sum_q3(k) + (Q(k, t/n_meas)/N_T)**2
                    sum_e3(k) = sum_e3(k) + E(k, t/n_meas)
                end do
                count3 = count3 + 1
            else if (t/n_meas >= 30001 .and. t/n_meas <= 100000) then !300.000-1.000.000 MCS
                do k = 1, N_T
                    sum_q4(k) = sum_q4(k) + (Q(k, t/n_meas)/N_T)**2
                    sum_e4(k) = sum_e4(k) + E(k, t/n_meas)
                end do
                count4 = count4 + 1
            end if
        end if ! end measure overlaps and energies

        ! Print progress
        if (mod(t, 100000) == 0) write(*,*) "Step = ", t, "/ 1000000"

    end do ! end of Monte Carlo simulation

    ! Calculate mean values
    mean_q1 = sum_q1 / count1
    mean_e1 = sum_e1 / count1
    mean_q2 = sum_q2 / count2
    mean_e2 = sum_e2 / count2
    mean_q3 = sum_q3 / count3
    mean_e3 = sum_e3 / count3
    mean_q4 = sum_q4 / count4
    mean_e4 = sum_e4 / count4

    write(3,*) "Temperature (K)", "<q^2>, 10^4 : 3·10^4", "<E>/N, 10^4 : 3·10^4", "<q^2>, 3·10^4 : 10^5", &
    "<E>/N, 3·10^4 : 10^5", "<q^2>, 10^5 : 3·10^5", "<E>/N, 10^5 : 3·10^5", "<q^2>, 3·10^5 : 10^6", "<E>/N, 3·10^5 : 10^6"
    do k = 1, N_T
        write(3,*) Temp(k), mean_q1(k), mean_e1(k)/N_T, mean_q2(k), mean_e2(k)/N_T, mean_q3(k), mean_e3(k)/N_T, mean_q4(k), &
                mean_e4(k)/N_T
    end do

    close(3)

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