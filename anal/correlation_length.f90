program jackv5_coratio_spin
    use m_variance
    implicit none

    integer, parameter :: &
        long = selected_int_kind(9),            &
        single = kind(1.0), 			&
        double = kind(1.0d0)
    real(kind=double), dimension(:), allocatable  :: kx, ky
    real(kind=double), dimension(:,:), allocatable  :: obs
    real(kind=double), dimension(:), allocatable  :: coratio
    real(kind=double) :: c1,c2,c3,c4,c5,c6, xm, xerr, dq
    real(kind=double) :: pi = dacos(-1.d0)
    integer :: nbins, i, j, nb, l, lq, idx_gamma
    character(len=64) :: fname, lx

    call get_command_argument(1, fname)
    call get_command_argument(2, lx)

    read(lx, *) l
    lq = l*l

    open(unit=5, file=fname, status='old')
    nbins = 0
    ! count the number of bins
    do
        do j =1, lq
            read(5,*,iostat=i) c1
        enddo
        if(i /= 0) exit
        nbins = nbins + 1
    end do

    !write(*,'(3x,a,i0,a)') 'there are ',nbins,' valid bin(s)'
    allocate(obs(lq, nbins))
    allocate(kx(lq))
    allocate(ky(lq))
    allocate(coratio(nbins))
    obs = 0.0d0

    ! read in data and calculate the correlation ratio for each bin
    rewind(5)
    do nb = 1, nbins
        do j = 1, lq
            read(5,*) kx(j), ky(j), obs(j, nb)
            ! for afm, the ordering vector is (-\pi, -\pi)
            if ( (kx(j)+pi) .lt. 0.00001d0 .and. (ky(j)+pi) .lt. 0.00001d0 ) then
                idx_gamma = j
            endif
        enddo

        ! for test
        !print*, idx_gamma, dq, obs(idx_gamma, nb), obs(idx_gamma+1,nb)

        ! set dq
        dq = abs(kx(2) - kx(1))
        coratio(nb) = dsqrt(1.d0/dq**2 * ( obs(idx_gamma, nb)/obs(idx_gamma+1,nb) - 1.d0)) ! the first point the gamma, the second is the nearest

        ! output each bin
        write(*, '(2e16.8)') coratio(nb)
    end do
    close(5)

    ! use jackknife to estimate the error
    !call errcalcj(coratio, xm, xerr)
    ! just the standard deviation
    !call errcalc(coratio(bin_start:nbins), xm, xerr)

    ! output avg and error bar
    !write(*, '(2e16.8)') xm, xerr

    ! free up space
    deallocate(obs, kx, ky, coratio)
end program jackv5_coratio_spin
