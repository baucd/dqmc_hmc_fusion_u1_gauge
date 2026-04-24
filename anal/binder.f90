PROGRAM binder
    use m_variance
    IMPLICIT NONE

    Integer, Parameter :: &
        long = selected_int_kind(9),            &
        single = kind(1.0), 			&
        double = kind(1.0D0)
    REAL(kind=double), DIMENSION(:), ALLOCATABLE  :: obs
    REAL(kind=double) :: c1,c2,c3,c4,c5,c6, XM, XERR
    INTEGER :: nbins, i, nb
    CHARACTER(len=64) :: fname

    call get_command_argument(1, fname)

    OPEN(UNIT=5, FILE=fname, STATUS='OLD')
    nbins = 0
    ! Count the number of bins
    DO
        READ(5,*,iostat=i) c1
        IF(i /= 0) EXIT
        nbins = nbins + 1
    END DO

    !WRITE(*,'(3X,A,I0,A)') 'There are ',nbins,' valid bin(s)'
    ALLOCATE(obs(Nbins))
    obs = 0.0D0

    ! Read in data
    REWIND(5)
    DO nb = 1, Nbins
        READ(5,*) obs(nb)
    END DO
    CLOSE(5)

    ! Use Jackknife to estimate the error
    call errcalc_binder(obs, nbins, XM, XERR)
    write(*, '(2e16.8)') XM, XERR
END PROGRAM binder
