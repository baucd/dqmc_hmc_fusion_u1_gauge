subroutine errcalc_binder(en,nbin, xm,xerr)
    !   calculates jacknife error on the input vector en.  mean and  variance.
    !   the input are the bins.
    use m_variance
    implicit none

    integer :: nbin
    !real (kind=8), dimension(:) ::  en
    real (kind=8) ::  en(nbin)
    real (kind=8)               ::  xm, xerr, x2, x4
    real (kind=8), dimension(:), allocatable ::  en1
    integer     :: np, n, n1

    np = size(en)
    allocate (en1(np))

    ! build the jackknife averages and send to errcalc.
    do n = 1,np
        x2 = 0.d0
        x4 = 0.d0
        do n1 = 1,np
            if (n1.ne.n) then
                x2 = x2 + en(n1)**2
                x4 = x4 + en(n1)**4
            endif
        enddo
        en1(n) = 1.d0 - x4/dble(np-1) / (3.d0*(x2/dble(np-1))**2)
    enddo
    call errcalc(en1,xm,xerr)
    deallocate  ( en1 )

    return
end subroutine errcalc_binder
