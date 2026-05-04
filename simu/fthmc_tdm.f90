module fthmc_tdm
    use ftdqmc_hamilt
    use ftdqmc_latt_class
    use fthmc_gfun

    implicit none
    ! index for the array variables
    integer, parameter :: IGFUN  = 1
    integer, parameter :: ICHIZ  = 2

type tdm
    ! measurement part
    integer :: nclass       ! number of distance pair
    integer :: narrays = 2  ! GFUN and CHIZ
    integer :: nobs

    !integer :: IARR(0: narrays + 1)
    integer, allocatable :: IARR(:)

    complex(dp), pointer :: AllProp(:,:)    ! vector of all physical properties

    ! pointer to all correlation function
    complex(dp), pointer :: GFUN(:,:)     ! G(tau)
    complex(dp), pointer :: CHIZ(:,:)     ! chiszsz

    ! flag
    logical :: init
endtype tdm

contains

    subroutine fthmc_tdm_alloc(T0)
        implicit none
        type(tdm), intent(inout) :: T0 ! tdm to be initalized

        ! local variables
        integer :: i, n

        allocate(T0%IARR(T0%narrays + 1))

        T0%nclass = lq

        ! size of the vector
        n = T0%nclass*T0%narrays

        ! allocate storage for all properties
        allocate(T0%AllProp(n, ltrot))

        ! pointer to beginning of each array
        T0%IARR(:) = 1
        ! pointer for correlation function
        do i = 1, T0%narrays + 1
            T0%IARR(i) = 1 + (i-1)*T0%nclass
        enddo

        ! set all pointers
        T0%GFUN    =>  T0%AllProp( T0%IARR(IGFUN) : T0%IARR(IGFUN+1) - 1, :)
        T0%CHIZ    =>  T0%AllProp( T0%IARR(ICHIZ) : T0%IARR(ICHIZ+1) - 1, :)

        T0%init = .true.

    endsubroutine fthmc_tdm_alloc


    subroutine fthmc_tdm_init(T0)
        implicit none
        type(tdm), intent(inout) :: T0 ! tdm to be initalized

        ! initalization
        T0%GFUN   = czero
        T0%CHIZ   = czero

        ! reset count
        T0%nobs   = 0

    endsubroutine fthmc_tdm_init

    subroutine fthmc_tdm_free(T0)
        implicit none
        type(tdm), intent(inout) :: T0 ! tdm to be freed

        ! executable
        if ( T0%init ) then
            nullify(T0%GFUN, T0%CHIZ)
            deallocate(T0%AllProp)
            deallocate(T0%IARR)
        endif
    endsubroutine fthmc_tdm_free

    subroutine fthmc_tdm_meas(nt, gfun0, T0, latt)
#IFDEF _OPENMP
        use OMP_LIB
#ENDIF
        implicit none
        type(tdm), intent(inout) :: T0
        type(gfun), intent(in), target :: gfun0
        integer, intent(in) :: nt
        class(ftdqmc_latt), intent(in) :: latt

        ! local variables
        complex(dp) :: ztmp
        complex(dp), dimension(:,:), pointer :: grt0_up, gr0t_up, grtt_up, gr00_up
        complex(dp), dimension(:,:), allocatable :: gfun_pnt, chiz_pnt

        ! local
        integer :: i, j, imj

        ! set the count, just for nt = 1 is enough
        if (nt .eq. 1) T0%nobs = T0%nobs + 1

        grt0_up => gfun0%gt0up(:,:,nt)
        gr0t_up => gfun0%g0tup(:,:,nt)
        grtt_up => gfun0%gttup(:,:,nt)
        gr00_up => gfun0%g00up(:,:,nt)

        ! for the observables
        allocate(gfun_pnt(lq,ltrot))
        allocate(chiz_pnt(lq,ltrot))
        gfun_pnt = T0%GFUN
        chiz_pnt = T0%CHIZ

!!$OMP PARALLEL &
!!$OMP PRIVATE ( j, i, imj, ztmp)
!!$OMP DO REDUCTION ( + : gfun_pnt, chiz_pnt)
        do j = 1, ndim
            do i = 1, ndim
                imj  = latt%latt_imj(i,j)
                gfun_pnt(imj, nt) = gfun_pnt(imj, nt) + grt0_up(i,j)
                ! szsz
                ztmp = - gr0t_up(j,i)*grt0_up(i,j)
                chiz_pnt(imj,nt) = chiz_pnt(imj,nt) + ztmp
            end do
        end do
!!$OMP END DO
!!$OMP END PARALLEL

        ! copy back the results
        T0%GFUN = gfun_pnt
        T0%CHIZ = chiz_pnt

        nullify(grt0_up, gr0t_up, grtt_up, gr00_up)
        deallocate(gfun_pnt, chiz_pnt)
    endsubroutine fthmc_tdm_meas

    subroutine fthmc_tdm_corFT(T0, latt)
        use mpi
        type(tdm), intent(in) :: T0
        class(ftdqmc_latt), intent(inout) :: latt

        ! local variables
        integer :: i
        complex(dp) :: mpi_cor_bin(lq, ltrot)
        character(40) :: file_root

        do i =1, T0%narrays
            if ( i .eq. IGFUN ) then
                call mpi_reduce( T0%GFUN, mpi_cor_bin, lq*ltrot, mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
                file_root='_gk.bin'
            elseif ( i .eq. ICHIZ ) then
                call mpi_reduce( T0%CHIZ, mpi_cor_bin, lq*ltrot, mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
                file_root='_sz.bin'
            endif

            ! Fourier transformation
            if( irank .eq. 0 ) then
                mpi_cor_bin = mpi_cor_bin / dcmplx(dble(isize*T0%nobs), 0.d0)
                call fthmc_tdm_ft(mpi_cor_bin, file_root, latt)
            endif

        enddo
        call mpi_barrier( mpi_comm_world, ierr )
    endsubroutine fthmc_tdm_corFT

    subroutine fthmc_tdm_ft(gr, file_root, latt)
#IFDEF _OPENMP
        use OMP_LIB
#ENDIF
        implicit none
        complex(dp), intent(in), dimension(:,:) :: gr
        character (40), intent(in) :: file_root
        class(ftdqmc_latt), intent(inout) :: latt

        ! local variables
        integer :: imj, iq, itau
        real(dp) :: qvec(2)
        complex(dp) :: gk
        complex(dp), external :: zdotu
        character (40) :: outname, cTemp

        !!! loop over k points
        do iq = 1, lq
            ! prepare the output file
            write(cTemp, '(i3)') iq
            outname = trim(adjustl(cTemp))//file_root
            open(unit=177,file=outname,status='unknown', action="write", position="append")

            ! retrive the k-point
            qvec = dble( latt%listk(iq,1))*latt%b1_p + dble( latt%listk(iq,2))*latt%b2_p

            ! The first line of data is blank
            write(177, *)

            !!! for itau = ltrot, which is just tau=0 which is the equal time correlation
            gk = czero
            gk = zdotu(lq, gr(:, ltrot), 1, cone/latt%zexpiqr(:, iq), 1) ! use blas level 1 function: zdotu
            write(177, '(2e16.8)') 0.d0, dble(gk*dcmplx(1.d0/dble(lq),0.d0))

            !!! output the tgrid to file for only one time
            if ( iq .eq. 1) then
                open(unit=100,file='tgrid.dat',status='unknown', action="write")
                write(100, *) ltrot
                write(100, '(f8.3)') 0.d0
            endif

            !!! for other tau
            do itau = 1, ltrot-1
                gk = czero
                gk = zdotu(lq, gr(:, itau), 1, cone/latt%zexpiqr(:, iq), 1) ! use blas level 1 function: zdotu
                write(177, '(2e16.8)') dtau*dble(itau), dble(gk*dcmplx(1.d0/dble(lq),0.d0))

                !!! output to tgrid.dat
                if ( iq .eq. 1) then
                    write(100, '(f8.3)') dtau*dble(itau)
                endif
            enddo

            ! close file
            close(177)
            close(100)
        end do
    endsubroutine fthmc_tdm_ft
endmodule fthmc_tdm
