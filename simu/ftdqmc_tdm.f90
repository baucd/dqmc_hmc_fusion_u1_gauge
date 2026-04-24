#ifdef __GFORTRAN__
#define GET_FILENAME_DYN(x) generate_filename_dyn("x")
#else
#define GET_FILENAME_DYN(x) generate_filename_dyn(#x)
#endif

module ftdqmc_tdm
    use ftdqmc_hamilt
    use ftdqmc_latt
    use ftdqmc_gfun

    implicit none
    ! set number of dyn correlators
    integer :: ndcorrelators = 3

    ! index for the array variables
    integer, parameter :: ICT_GRUP    = 1
    integer, parameter :: ICT_GRDN    = 2
    integer, parameter :: ICT_SPSM    = 3
    !integer, parameter :: ICT_SWAVE   = 4
    !integer, parameter :: IISING   = 3
    !integer, parameter :: IPHYGFUN = 4

    ! pointer array to store pointers
    type :: pnt1d_t
    	complex(dp), pointer :: pnt(:,:)
	end type pnt1d_t

    type tdm
        ! measurement part
        integer :: nclass       ! number of distance pair
        integer :: nobs

        integer, allocatable :: IARR(:)
        complex(dp), pointer :: AllProp(:,:)    ! vector of all physical properties

        ! pointer to dyn correlators
        type(pnt1d_t), allocatable :: pnt_dcr(:)

        ! flag
        logical :: init
    endtype tdm

contains

    subroutine ftdqmc_tdm_alloc(T0)
        type(tdm), intent(inout) :: T0 ! tdm to be initalized

        ! local variables
        integer :: i, n

        allocate(T0%IARR(ndcorrelators + 1))

        T0%nclass = lq

        ! size of the vector
        n = T0%nclass*ndcorrelators

        ! allocate storage for all properties
        allocate(T0%AllProp(n, ltrot))

        ! pointer to beginning of each array
        T0%IARR(:) = 1
        ! pointer for correlation function
        do i = 1, ndcorrelators + 1
            T0%IARR(i) = 1 + (i-1)*T0%nclass
        enddo

        ! set all pointers
        allocate(T0%pnt_dcr(ndcorrelators))
        do i = 1, ndcorrelators
            T0%pnt_dcr(i)%pnt => T0%AllProp( T0%IARR(i) : T0%IARR(i+1) - 1, :)
        enddo
        ! old code
        !T0%GFUN       =>  T0%AllProp( T0%IARR(IGFUN) : T0%IARR(IGFUN+1) - 1, :)
        !T0%CHIZ       =>  T0%AllProp( T0%IARR(ICHIZ) : T0%IARR(ICHIZ+1) - 1, :)
        !T0%ISING_GFUN =>  T0%AllProp( T0%IARR(IISING) : T0%IARR(IISING+1) - 1, :)
        !T0%PHYGFUN    =>  T0%AllProp( T0%IARR(IPHYGFUN) : T0%IARR(IPHYGFUN+1) - 1, :)
        !T0%SWAVE      =>  T0%AllProp( T0%IARR(ISWAVE) : T0%IARR(ISWAVE+1) - 1, :)

        T0%init = .true.

    endsubroutine ftdqmc_tdm_alloc


    subroutine ftdqmc_tdm_free(T0)
        type(tdm), intent(inout) :: T0 ! tdm to be freed

        integer :: i

        ! executable
        if ( T0%init ) then
            do i = 1, ndcorrelators
                nullify(T0%pnt_dcr(i)%pnt)
            enddo

            deallocate(T0%AllProp)
            deallocate(T0%IARR)
        endif
    endsubroutine ftdqmc_tdm_free


    subroutine ftdqmc_tdm_init(T0)
        type(tdm), intent(inout) :: T0 ! tdm to be initalized

        integer :: i

        ! initalization
        !T0%GFUN   = czero
        !T0%CHIZ   = czero
        !T0%ISING_GFUN   = 0
        !T0%PHYGFUN   = czero
        !T0%SWAVE   = czero
        do i = 1, ndcorrelators
            T0%pnt_dcr(i)%pnt(:,:) = czero
        enddo

        ! reset count
        T0%nobs   = 0

    endsubroutine ftdqmc_tdm_init


    subroutine ftdqmc_tdm_meas(nt_remote, gfun0, T0, phi)
#ifdef _OPENMP
        use OMP_LIB
#endif
        use ftdqmc_auxfield_f4_class

        implicit none
        type(tdm), intent(inout) :: T0
        type(gfun), intent(in), target :: gfun0
        integer, intent(in) :: nt_remote
        type(ftdqmc_auxfieldHolder), intent(inout) :: phi(num_fields)

        ! local variables
        complex(dp) :: ztmp, ztmp1, ztmp2, ztmp3, ztmp0
        complex(dp), dimension(:,:), allocatable :: grt0_up, gr0t_up, grtt_up, gr00_up
        complex(dp), dimension(:,:), allocatable :: grt0_dnc, gr0t_dnc, grtt_dnc, gr00_dnc
        complex(dp), dimension(:,:), allocatable :: fupt_fdn0, fup0_fdnt
        complex(dp), dimension(:,:), allocatable :: fddnt_fdup0, fddn0_fdupt

        ! local
        integer :: i, j, imj, nt, nx, ny, idx1, idx2, nt_real
        complex(dp), pointer :: chargon(:,:,:)
        integer, external :: npbc

        ! set nt:  nt_remote start from 0 to ltrot-1
        nt = nt_remote + 1
        nt_real = npbc(nt_remote, ltrot)
        ! meas count
        if (nt .eq. 1) T0%nobs = T0%nobs + 1

        ! initial pointer
        chargon => phi(idx_chargon)%pnt%chargon

        allocate(grt0_up(lq, lq))
        allocate(gr0t_up(lq, lq))
        !allocate(grtt_up(lq, lq))
        !allocate(gr00_up(lq, lq))
        allocate(grt0_dnc(lq, lq))
        allocate(gr0t_dnc(lq, lq))
        !allocate(grtt_dnc(lq, lq))
        !allocate(gr00_dnc(lq, lq))
        allocate(fupt_fdn0(lq, lq))
        allocate(fup0_fdnt(lq, lq))
        allocate(fddnt_fdup0(lq, lq))
        allocate(fddn0_fdupt(lq, lq))

        ! setup dyn green's function: convert to original basis f_up and f_down
        do i =1, lq
            do j =1, lq
                ! get dyn grup
                nx = list_prim(i, 1)
                ny = list_prim(i, 2)
                idx1 = invlist(nx, ny, 1)
                nx = list_prim(j, 1)
                ny = list_prim(j, 2)
                idx2 = invlist(nx, ny, 1)
                grt0_up(i,j) = gfun0%gt0up(idx1, idx2)
                gr0t_up(i,j) = gfun0%g0tup(idx1, idx2)

                ! get dyn grdn
                nx = list_prim(i, 1)
                ny = list_prim(i, 2)
                idx1 = invlist(nx, ny, 2)
                nx = list_prim(j, 1)
                ny = list_prim(j, 2)
                idx2 = invlist(nx, ny, 2)
                grt0_dnc(i,j) = gfun0%gt0up(idx1, idx2) / dcmplx(epsj(i), 0.d0) / dcmplx(epsj(j), 0.d0)
                gr0t_dnc(i,j) = gfun0%g0tup(idx1, idx2) / dcmplx(epsj(i), 0.d0) / dcmplx(epsj(j), 0.d0)

                ! get fupfdn
                nx = list_prim(i, 1)
                ny = list_prim(i, 2)
                idx1 = invlist(nx, ny, 1)
                nx = list_prim(j, 1)
                ny = list_prim(j, 2)
                idx2 = invlist(nx, ny, 2)
                fupt_fdn0(i,j) = gfun0%gt0up(idx1, idx2) / dcmplx(epsj(j), 0.d0)
                fup0_fdnt(i,j) = gfun0%g0tup(idx1, idx2) / dcmplx(epsj(j), 0.d0)

                ! get fddnfdup
                nx = list_prim(i, 1)
                ny = list_prim(i, 2)
                idx1 = invlist(nx, ny, 2)
                nx = list_prim(j, 1)
                ny = list_prim(j, 2)
                idx2 = invlist(nx, ny, 1)
                fddnt_fdup0(i,j) = gfun0%gt0up(idx1, idx2) / dcmplx(epsj(i), 0.d0)
                fddn0_fdupt(i,j) = gfun0%g0tup(idx1, idx2) / dcmplx(epsj(i), 0.d0)
            enddo
        enddo

        ! measurements
        do j = 1, lq
            do i = 1, lq
                ! set imj
                imj  = latt_imj(i,j)

                ! reset tmp variables
                ztmp1 = czero
                ztmp2 = czero
                ztmp3 = czero

                ! up electron
                ztmp1 = ztmp1 + dconjg(chargon(1,i,nt_real))*chargon(1,j,nt_real) * grt0_up(i,j)
                ztmp1 = ztmp1 + dconjg(chargon(1,i,nt_real))*chargon(2,j,nt_real) * fupt_fdn0(i,j)
                ztmp1 = ztmp1 + dconjg(chargon(2,i,nt_real))*chargon(1,j,nt_real) * fddnt_fdup0(i,j)
                ztmp1 = ztmp1 + dconjg(chargon(2,i,nt_real))*chargon(2,j,nt_real) * grt0_dnc(i,j)

                ! dn electron
                ztmp2 = ztmp2 - dconjg(chargon(1,i,nt_real))*chargon(1,j,nt_real) * gr0t_dnc(j,i)
                ztmp2 = ztmp2 + dconjg(chargon(1,i,nt_real))*chargon(2,j,nt_real) * fup0_fdnt(j,i)
                ztmp2 = ztmp2 + dconjg(chargon(2,i,nt_real))*chargon(1,j,nt_real) * fddn0_fdupt(j,i)
                ztmp2 = ztmp2 - dconjg(chargon(2,i,nt_real))*chargon(2,j,nt_real) * gr0t_up(j,i)

                ! spsm
                ztmp3 = ztmp3 + gr0t_up(j, i) * gr0t_dnc(j,i) *dcmplx(1.5d0, 0.d0)
                ztmp3 = ztmp3 - fddn0_fdupt(j,i) * fup0_fdnt(j, i)*dcmplx(1.5d0, 0.d0)

                ! record the meas
                T0%pnt_dcr(ICT_GRUP)%pnt(imj, nt) = T0%pnt_dcr(ICT_GRUP)%pnt(imj, nt) + ztmp1
                T0%pnt_dcr(ICT_GRDN)%pnt(imj, nt) = T0%pnt_dcr(ICT_GRDN)%pnt(imj, nt) + ztmp2
                T0%pnt_dcr(ICT_SPSM)%pnt(imj, nt) = T0%pnt_dcr(ICT_SPSM)%pnt(imj, nt) + ztmp3
            end do
        end do

        if (allocated(grt0_up)) deallocate(grt0_up)
        if (allocated(gr0t_up)) deallocate(gr0t_up)
        if (allocated(grtt_up)) deallocate(grtt_up)
        if (allocated(gr00_up)) deallocate(gr00_up)
        if (allocated(grt0_dnc)) deallocate(grt0_dnc)
        if (allocated(gr0t_dnc)) deallocate(gr0t_dnc)
        if (allocated(grtt_dnc)) deallocate(grtt_dnc)
        if (allocated(gr00_dnc)) deallocate(gr00_dnc)
        if (allocated(fupt_fdn0)) deallocate(fupt_fdn0)
        if (allocated(fup0_fdnt)) deallocate(fup0_fdnt)
        if (allocated(fddnt_fdup0)) deallocate(fddnt_fdup0)
        if (allocated(fddn0_fdupt)) deallocate(fddn0_fdupt)

    endsubroutine ftdqmc_tdm_meas


    subroutine ftdqmc_tdm_corFT(T0)
        type(tdm), intent(inout) :: T0

        !do i =1, T0%narrays
        !    if ( i .eq. IGFUN ) then
        !        call mpi_reduce( T0%GFUN, mpi_cor_bin, lq*ltrot, mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
        !        file_root='_gk.bin'
        !    elseif ( i .eq. ICHIZ ) then
        !        call mpi_reduce( T0%CHIZ, mpi_cor_bin, lq*ltrot, mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
        !        file_root='_sz.bin'
        !    elseif ( i .eq. ISWAVE) then
        !        call mpi_reduce( T0%SWAVE, mpi_cor_bin, lq*ltrot, mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
        !        file_root='_swave.bin'
        !    endif

        !    ! Fourier transformation
        !    if ( i .eq. IGFUN .or. i .eq. ICHIZ .or. i .eq. ISWAVE) then
        !        if( irank .eq. 0 ) then
        !            mpi_cor_bin = mpi_cor_bin / dcmplx(dble(isize*T0%nobs), 0.d0)
        !            call ftdqmc_tdm_ft(mpi_cor_bin, file_root)
        !        endif
        !    endif
        !enddo
        !call mpi_barrier( mpi_comm_world, ierr )

        call ftdqmc_tdm_perfFT(T0, ICT_GRUP, GET_FILENAME_DYN(ICT_GRUP))
        call ftdqmc_tdm_perfFT(T0, ICT_GRDN, GET_FILENAME_DYN(ICT_GRDN))
        call ftdqmc_tdm_perfFT(T0, ICT_SPSM, GET_FILENAME_DYN(ICT_SPSM))

    endsubroutine ftdqmc_tdm_corFT


    subroutine ftdqmc_tdm_perfFT(T0, ICT_index, filename)
		! perform FT for input correlation function
        use mpi

        type(tdm), intent(inout) :: T0
        integer, intent(in) :: ICT_index
        character(40), intent(in) :: filename

        ! local variables
        complex(dp) :: mpi_cor_bin(lq, ltrot)

        ! get data and perform FT
        call mpi_reduce( T0%pnt_dcr(ICT_index)%pnt, mpi_cor_bin, lq*ltrot, mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
        ! Fourier transformation
        if( irank .eq. 0 ) then
            mpi_cor_bin = mpi_cor_bin / dcmplx(dble(isize*T0%nobs), 0.d0)
            call ftdqmc_tdm_ft(mpi_cor_bin, filename, ICT_index, T0)
        endif
        call mpi_barrier( mpi_comm_world, ierr )
	endsubroutine ftdqmc_tdm_perfFT


    subroutine ftdqmc_tdm_ft(gr, filename, ICT_index, T0)
#ifdef _OPENMP
        use OMP_LIB
#endif
        implicit none
        complex(dp), intent(in), dimension(:,:) :: gr
        character (40), intent(in) :: filename
        integer, intent(in) :: ICT_index
        type(tdm), intent(inout) :: T0

        ! local variables
        integer :: imj, iq, itau
        real(dp) :: qvec(2)
        complex(dp) :: gk
        complex(dp), external :: zdotu
        character (40) :: outname, cTemp

        ! output time grid to file tgrid.dat
        open(unit=100,file='tgrid.dat',status='unknown', action="write")
        write(100, *) ltrot
        do itau = 1, ltrot
            write(100, '(f8.3)') dtau*dble(itau-1) ! first element \Delta tau = 0.d0
        enddo
        close(100)

        ! prepare the dsq output file
        open(unit=177,file=filename,status='unknown', action="write", position="append")
        ! loop over k points
        do iq = 1, lq
            ! retrive the k-point
            qvec = dble( listk(iq,1))*b1_p + dble( listk(iq,2))*b2_p

            ! The first line of data is index of qvec
            write(177, *) iq
            do itau = 1, ltrot
                ! use blas level 1 function: zdotu
                gk = zdotu(lq, gr(:, itau), 1, cone/zexpiqr(:, iq), 1)
                write(177, '(2e16.8)') dtau*dble(itau-1), dble(gk*dcmplx(1.d0/dble(lq),0.d0))
            enddo
        end do
        close(177)

    endsubroutine ftdqmc_tdm_ft


    function generate_filename_dyn(raw_name) result(f_name)
        character(len=*), intent(in) :: raw_name
        character(len=64) :: f_name
        character(len=64) :: temp_part
        integer :: i, offset

        ! Logic is the same: extract the part after ICT_.
        temp_part = raw_name(5:len_trim(raw_name))

        ! uppercase to lowercase
        offset = ichar('a') - ichar('A')
        do i = 1, len_trim(temp_part)
            if (temp_part(i:i) >= 'A' .and. temp_part(i:i) <= 'Z') then
                temp_part(i:i) = char(ichar(temp_part(i:i)) + offset)
            end if
        end do

        f_name = "dsq_" // trim(temp_part) // ".bin"
    end function generate_filename_dyn


endmodule ftdqmc_tdm
