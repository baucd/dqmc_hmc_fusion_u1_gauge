PROGRAM rebin
    implicit none
    Integer, Parameter :: &
        long = selected_int_kind(9),            &
        single = kind(1.0), 			&
        double = kind(1.0D0)
    REAL(kind=double), DIMENSION(:), ALLOCATABLE :: obs
    REAL(kind=double) :: c1, binvalue
    ! nbins: original bins  rebinsize: expected bins  binsize: number of samples of each rebin
    INTEGER :: n, nbins, rebinsize, binsize, i, start, end_
    CHARACTER(len=64) :: fname

    call get_command_argument(1, fname)

    ! read in the rebinsize parameter
    OPEN(UNIT=6, FILE='rebinsize', STATUS='OLD')
    read(6, *) rebinsize
    close(6)

    ! read the data file and count the number of bins
    OPEN(UNIT=5, FILE=fname, STATUS='OLD')
    nbins = 0
    ! Count the number of bins
    DO
        READ(5,*,iostat=i) c1
        IF(i /= 0) EXIT
        nbins = nbins + 1
    END DO

    ! allocate the vector
    allocate( obs(nbins) )
    obs = 0.d0

    ! read in data
    rewind(5)
    do n = 1, nbins
        read(5, *) obs(n)
    enddo
    close(5)

    ! decide the number of sample of each bin given the binsize
    binsize = int(dble(nbins)/dble(rebinsize))
    !write(*,*) 'test, binsize', binsize

    ! rebinning
    do i = 1, rebinsize
        start = binsize * (i-1) + 1
        end_ = binsize * (i-1) + binsize
        binvalue = sum( obs( start:end_ )) / dble(binsize)
        ! write out to the screen
        write(*, '(e16.8)') binvalue
    enddo
END PROGRAM rebin
