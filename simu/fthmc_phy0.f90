module fthmc_phy0
    use ftdqmc_hamilt
    use ftdqmc_latt_class
    use ftdqmc_auxfield_f5_class
    use fthmc_gfun

    implicit none
    ! index for the array variables
    integer, parameter :: IMEAS  = 0
    integer, parameter :: IDIME  = 1 ! dimer-dimer correlation
    integer, parameter :: ISPSM  = 2 ! S+ S- correlation

    ! index of the scalar variables (IMEAS)
    integer, parameter :: P0_FILL    = 1    ! electron filling
    integer, parameter :: P0_ELKIN   = 2    ! electron kinetic energy
    integer, parameter :: P0_ZEJS    = 3    ! js term
    integer, parameter :: P0_ZEJPI   = 4    ! jpi term
    integer, parameter :: P0_ZETOT   = 5    ! total energy
    integer, parameter :: P0_MZ      = 6    ! mz
    integer, parameter :: P0_NETFLUX = 7    ! netflux
    integer, parameter :: P0_NE_FLUC = 8    ! netflux

type phy0
    ! measurement part
    integer :: nclass       ! number of distance pair
    integer :: nobs         ! number of measurement
    integer :: nmeas = 8    ! number of scalar properties
    integer :: narrays = 2  ! dimer s+s-

    integer, allocatable :: IARR(:)

    complex(dp), pointer :: meas(:)       ! scalar variables
    complex(dp), pointer :: AllProp_scalar(:)    ! vector of all physical properties
    complex(dp), pointer :: AllProp(:)    ! vector of all physical properties

    ! pointer to all correlation function
    complex(dp), pointer :: Dimer(:)    ! dimer-dimer correlation
    complex(dp), pointer :: Spsm(:)    ! dimer-dimer correlation

    ! flag
    logical :: init
endtype phy0

contains

    subroutine fthmc_phy0_alloc(P0)
        implicit none
        type(phy0), intent(inout) :: P0 ! phy0 to be initalized

        ! local variables
        integer :: i, n

        allocate(P0%IARR(P0%narrays+1))

        P0%nclass = lq

        ! size of the vector
        n = P0%nclass*P0%narrays

        ! allocate storage for all properties
        allocate(P0%AllProp_scalar(P0%nmeas))
        allocate(P0%AllProp(n))

        ! pointer to beginning of each array
        P0%IARR(:) = 1
        ! pointer for correlation function
        do i = 1, P0%narrays + 1
            P0%IARR(i) = 1 + (i-1)*P0%nclass
        enddo

        ! set all pointers
        P0%meas    =>  P0%AllProp_scalar( : )
        P0%Dimer   =>  P0%AllProp( P0%IARR(IDIME) : P0%IARR(IDIME+1) - 1)
        P0%Spsm    =>  P0%AllProp( P0%IARR(ISPSM) : P0%IARR(ISPSM+1) - 1)
        P0%init = .true.
    endsubroutine fthmc_phy0_alloc


    subroutine fthmc_phy0_init(P0)
        implicit none
        type(phy0), intent(inout) :: P0 ! phy0 to be initalized
        ! initalization
        P0%meas   = czero
        P0%Dimer  = czero
        P0%Spsm   = czero
        P0%nobs   = 0

    endsubroutine fthmc_phy0_init

    subroutine fthmc_phy0_free(P0)
        implicit none
        type(phy0), intent(inout) :: P0 ! phy0 to be freed

        ! executable
        if ( P0%init ) then
            nullify(P0%meas)
            nullify(P0%Dimer, P0%Spsm)

            deallocate(P0%AllProp)
            deallocate(P0%IARR)
        endif
    endsubroutine fthmc_phy0_free

    subroutine fthmc_phy0_meas(gfun0, P0, phi_u1, latt)
#IFDEF _OPENMP
        use OMP_LIB
#ENDIF
        implicit none
        type(phy0), intent(inout) :: P0
        type(gfun), intent(in) :: gfun0
        class(ftdqmc_auxfield_f5), intent(in) :: phi_u1
        class(ftdqmc_latt), intent(in) :: latt

        ! local
        integer :: i, j, imj, i1, i2, i3, nf, i_plaqA, i_plaqB, ib, inf, ilf, b_plaqA, b_plaqB, jax, jmx, iax, imx
        integer :: nt
        integer :: nx_i, ny_i, nx_j, ny_j
        integer, external :: npbc
        real(dp) :: phi_Aplaq, phi_Bplaq, dphi
        real(dp) :: xi, xj, netflux_tot
        real(dp), allocatable :: netflux(:)
        complex(dp) :: zne, zkint, zejs, zejpi, zetot, zmz, zne_fluc
        complex(dp) :: ztmp1, ztmp2, ztmp3
        complex(dp), pointer :: grup(:,:,:)
        real(dp), pointer :: xfield(:,:,:)
        complex(dp), allocatable :: grupc(:,:,:), grdn(:,:,:), grdnc(:,:,:)
#IFDEF CUDA_MEAS
        complex(dp), allocatable :: grupc_vec(:), grdn_vec(:), grdnc_vec(:)
#ENDIF

        ! for test
        !real(dp) :: tstart2, tend2

        !P0%nobs = P0%nobs + 1
        P0%nobs = P0%nobs + ltrot

        ! convention
        !grup (i,j) = < c_i c^+_j >
        !grupc (i,j) = < c^+_i c_j >

        ! initial pointer, point to the nt equal time green's function
        grup => gfun0%grup(:,:,:)
        xfield => phi_u1%gauge_u1

        ! alloc matrices
        allocate(grupc(ndim, ndim, ltrot))
        allocate(grdn(ndim, ndim, ltrot))
        allocate(grdnc(ndim, ndim, ltrot))
        allocate(netflux(ltrot))

        ! gfun setup and correlation function
        !call CPU_TIME(tstart2)
#IFNDEF CUDA_MEAS
        do nt = 1, ltrot

            ! prepare grupc, grdn, grdnc together
            do i = 1, ndim
                do j = 1, ndim
                    ! first grupc
                    if ( i .eq. j) then
                        grupc(j,i,nt) = -grup(i,j,nt) + cone
                    else
                        grupc(j,i,nt) = -grup(i,j,nt)
                    endif

                    ! second grdn and grdnc
                    if ( mu .eq. 0) then
                        grdn (j,i,nt) = grup (j,i,nt)
                        grdnc(j,i,nt) = grupc(j,i,nt)
                    else
                        nx_j = latt%list(j,1);ny_j = latt%list(j,2)
                        xj = 1.d0; if ( mod(nx_j,2) .ne. mod(ny_j,2) ) xj = -1.d0
                        nx_i = latt%list(i,1);ny_i = latt%list(i,2)
                        xi = 1.d0; if ( mod(nx_i,2) .ne. mod(ny_i,2) ) xi = -1.d0

                        grdn (j,i,nt) = dcmplx(xj*xi, 0.d0)*dconjg ( grupc(j,i,nt) )
                        grdnc(j,i,nt) = dcmplx(xj*xi, 0.d0)*dconjg ( grup (j,i,nt) )
                    endif
                end do
            end do

            ! fermion correlation functions
            do j = 1, ndim
                do i = 1, ndim
                    imj = latt%latt_imj(i,j)

                    ! first dimer-dimer correlation
                    jax = latt%nnlist(j,1) ! j+x
                    jmx = latt%nnlist(j,3) ! j-x
                    iax = latt%nnlist(i,1) ! i+x
                    imx = latt%nnlist(i,3) ! i-x
                    P0%Dimer(imj) = P0%Dimer(imj) + grupc(i,iax,nt)*grup(i,iax,nt)*grupc(j,jax,nt)  *grup(j,jax,nt)  *z4  &
                                                  + grupc(i,j,nt)  *grup(i,j,nt)  *grupc(iax,jax,nt)*grup(iax,jax,nt)*z2  &
                                                  + grupc(i,jax,nt)*grup(i,jax,nt)*grupc(iax,j,nt)  *grup(iax,j,nt)  *z2  &
                                                  + grupc(i,jax,nt)*grup(i,iax,nt)*grup(iax,j,nt)   *grup(j,jax,nt)  *z3  &
                                                  + grupc(i,iax,nt)*grup(i,jax,nt)*grup(j,iax,nt)   *grup(jax,j,nt)  *z3  &
                                                  + grupc(i,j,nt)  *grup(i,iax,nt)*grup(iax,jax,nt) *grup(jax,j,nt)  *z3  &
                                                  + grupc(i,iax,nt)*grup(i,j,nt)  *grup(jax,iax,nt) *grup(j,jax,nt)  *z3  &
                                                  + grupc(i,jax,nt)*grup(i,j,nt)  *grup(iax,jax,nt) *grup(j,iax,nt)  *z1  &
                                                  + grupc(i,j,nt)  *grup(i,jax,nt)*grup(jax,iax,nt) *grup(iax,j,nt)  *z1

                    ! second s+ s- correlation
                    P0%Spsm(imj) = P0%Spsm(imj) + grupc(i,j,nt)*grup(i,j,nt)
                end do
            end do
        enddo
#ELSE
        ! prepare input
        allocate(grupc_vec(ndim*ndim*ltrot))
        allocate(grdn_vec (ndim*ndim*ltrot))
        allocate(grdnc_vec(ndim*ndim*ltrot))
        grupc_vec = reshape(grupc, (/ndim*ndim*ltrot/))
        grdn_vec  = reshape(grdn,  (/ndim*ndim*ltrot/))
        grdnc_vec = reshape(grdnc, (/ndim*ndim*ltrot/))

        ! gpu code
        call ftdqmc_gpu_phy0_meas_gf_cr(ndim, ltrot, reshape(grup, (/ndim*ndim*ltrot/)), grupc_vec, grdn_vec, grdnc_vec, P0%Dimer, P0%Spsm)

        ! make sure grupc, grdn, grdnc arrive host memory
        call ftdqmc_gpu_stream1sync()

        ! restore grdn etc
        grupc = reshape(grupc_vec, (/ndim, ndim, ltrot/))
        grdn  = reshape(grdn_vec,  (/ndim, ndim, ltrot/))
        grdnc = reshape(grdnc_vec, (/ndim, ndim, ltrot/))

#ENDIF
        !call CPU_TIME(tend2); time_vec(8) = time_vec(8) + (tend2 - tstart2)

        ! measurements on CPU
        ! init temp variables
        zkint = czero
        zejs = czero
        zejpi = czero
        netflux(:) = 0.d0
        zne = czero
        zmz = czero
        zne_fluc = czero

        ! electron filling and gauge field part, on CPU always
        do nt = 1, ltrot
            ! electron fillling
            do i = 1, ndim
                zne = zne + grupc(i,i,nt) + grdnc(i,i,nt)
                zne_fluc = zne_fluc + grupc(i,i,nt)*grupc(i,i,nt) + grupc(i,i,nt)*grup(i,i,nt) + &
                    grupc(i,i,nt)*grdnc(i,i,nt) + grdnc(i,i,nt)*grupc(i,i,nt) + &
                    grdnc(i,i,nt)*grdnc(i,i,nt) + grdnc(i,i,nt)*grdn(i,i,nt) - &
                    dcmplx(2.d0, 0.d0) * (grupc(i,i,nt) + grdnc(i,i,nt)) + cone
                zmz = zmz +  (grupc(i,i,nt) - grdnc(i,i,nt))*chalf
            end do

            ! gauge field part
            do nf = 1, nfam
                do i = 1, lfam
	               i1 = latt%l_bonds(1,i,nf)
	               i2 = latt%l_bonds(2,i,nf)

                   ztmp1 = zep_rsigl_k(i,nf, npbc(nt-1, ltrot) ) /cinvsqrt2 * (grupc(i1,i2,nt) + grdnc(i1,i2,nt))
                   zkint = zkint  +  ztmp1 + dconjg(ztmp1)

                   ! A plaq
                   phi_Aplaq = 0.d0
                   i_plaqA = latt%inv_Aplaq_bondcord(1,i,nf)
                   b_plaqA = latt%inv_Aplaq_bondcord(2,i,nf)
                   do ib = 1, 4
                       inf = latt%plaq_bondcord(1,ib,i_plaqA)
                       ilf = latt%plaq_bondcord(2,ib,i_plaqA)
                       phi_Aplaq = phi_Aplaq + xfield(ilf,inf,nt)*latt%sgnA_plaq(ib)
                   end do
                   netflux(nt) = netflux(nt) + floor(phi_Aplaq/(2.d0*pi))

                   ! B plaq
                   phi_Bplaq = 0.d0
                   i_plaqB = latt%inv_Bplaq_bondcord(1,i,nf)
                   b_plaqB = latt%inv_Bplaq_bondcord(2,i,nf)
                   do ib = 1, 4
                       inf = latt%plaq_bondcord(1,ib,i_plaqB)
                       ilf = latt%plaq_bondcord(2,ib,i_plaqB)
                       phi_Bplaq = phi_Bplaq + xfield(ilf,inf,nt)*latt%sgnB_plaq(ib)
                   end do
                   netflux(nt) = netflux(nt) + floor(phi_Bplaq/(2.d0*pi))

                   dphi = xfield(i,nf,npbc(nt+1,ltrot)) - xfield(i,nf,nt)
                   zejs = zejs + dcmplx( 2.d0*( 1.d0 - dcos(dphi) ) , 0.d0 )
                   zejpi = zejpi + dcmplx( dcos(phi_Aplaq) + dcos(phi_Bplaq), 0.d0 )
                end do
            enddo
        enddo

        ! set scalars
        zne = zne * dcmplx(Nflavor/2.d0,0.d0)
        zmz = zmz * dcmplx(Nflavor/2.d0,0.d0)
        zne_fluc = zne_fluc * dcmplx(Nflavor/2.d0,0.d0)
        zkint = zkint*dcmplx(-rt*Nflavor/2.d0,0.d0)
        zejs = zejs/dcmplx( js*dtau*dtau, 0.d0 )
        zejpi = zejpi * dcmplx(0.25d0,0.d0)  ! 0.25 is because of 4 times conuting
        zetot = zkint + zejs + zejpi
        netflux_tot = sum(abs(netflux(:)*0.25d0+ndim/2.d0))
        P0%meas(P0_FILL) = P0%meas(P0_FILL) + zne
        P0%meas(P0_MZ) = P0%meas(P0_MZ) + zmz
        P0%meas(P0_NE_FLUC) = P0%meas(P0_NE_FLUC) + zne_fluc
        P0%meas(P0_ELKIN) = P0%meas(P0_ELKIN) + zkint
        P0%meas(P0_ZEJS) = P0%meas(P0_ZEJS) + zejs
        P0%meas(P0_ZEJPI) = P0%meas(P0_ZEJPI) + zejpi
        P0%meas(P0_ZETOT) = P0%meas(P0_ZETOT) + zetot
        P0%meas(P0_NETFLUX) = P0%meas(P0_NETFLUX) + dcmplx(netflux_tot, 0.d0)

#IFDEF CUDA_MEAS
        DEALLOCATE(grupc_vec)
        DEALLOCATE(grdn_vec)
        DEALLOCATE(grdnc_vec)
#ENDIF
        deallocate(netflux)
        deallocate(grupc)
        deallocate(grdn)
        deallocate(grdnc)
        nullify(grup)
        nullify(xfield)
    endsubroutine fthmc_phy0_meas


    subroutine fthmc_phy0_getavg(P0)
        ! scalar quantities
        use mpi

        implicit none
        type(phy0), intent(inout) :: P0

        ! local variables
        complex(dp) :: mpi_meas_bin(P0%nmeas)

        ! executable
        call mpi_reduce( P0%meas, mpi_meas_bin, P0%nmeas, mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr)
        P0%meas = mpi_meas_bin

        if ( irank .eq. 0) then
            ! average
            P0%meas = P0%meas / dcmplx( dble(isize * P0%nobs ), 0.d0)

            ! output
            open( unit=90, file='ener1.bin', status='unknown', action='write', position='append')

            write(90, '(30(e16.8, 2x))') dble(P0%meas(P0_FILL))/ndim,&  ! #1
                                         dble(P0%meas(P0_ELKIN))/ndim,& ! #2
                                         dble(P0%meas(P0_ZEJS))/ndim,&  ! #3
                                         dble(P0%meas(P0_ZEJPI))/ndim,& ! #4
                                         dble(P0%meas(P0_ZETOT))/ndim,& ! #5
                                         dble(P0%meas(P0_MZ))/ndim,&    ! #6
                                         dble(P0%meas(P0_NE_FLUC))/ndim,& ! #7
                                         dble(P0%meas(P0_NETFLUX))      ! #8
            close(90)
        endif

        call mpi_barrier( mpi_comm_world, ierr )
    endsubroutine fthmc_phy0_getavg

    subroutine fthmc_phy0_corFT(P0, latt)
        use mpi
        implicit none

        type(phy0), intent(in) :: P0
        class(ftdqmc_latt), intent(inout) :: latt

        ! local variables
        integer :: i
        complex(dp) :: mpi_cor_bin_dimer(lq)
        complex(dp) :: mpi_cor_bin_afm(lq)
        character(40) :: filename1, filename2

        ! firstly reduce data
        call mpi_reduce( P0%Dimer, mpi_cor_bin_dimer, lq, mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
        call mpi_reduce( P0%Spsm, mpi_cor_bin_afm, lq, mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )

        mpi_cor_bin_afm = mpi_cor_bin_afm / dcmplx(dble(isize*P0%nobs), 0.d0)
        mpi_cor_bin_dimer = mpi_cor_bin_dimer / dcmplx(dble(isize*P0%nobs), 0.d0) &
            - mpi_cor_bin_afm(lq-l)*mpi_cor_bin_afm(lq-l) * dcmplx((Nflavor*Nflavor-1)**2/dble(lq), 0.d0) ! deduct the background

        do i =1, P0%narrays
            if ( i .eq. IDIME ) then
                filename1='sq_dimer.bin'; filename2='dimer.bin'
            elseif ( i .eq. ISPSM ) then
                filename1='sq_spsm.bin'; filename2='spsm.bin'
            endif

            ! Fourier transformation
            if( irank .eq. 0 ) then
                if ( i .eq. IDIME ) then
                    call fthmc_phy0_ft(mpi_cor_bin_dimer, filename1, filename2, i, P0, latt)
                elseif ( i .eq. ISPSM ) then
                    call fthmc_phy0_ft(mpi_cor_bin_afm, filename1, filename2, i, P0, latt)
                endif
            endif

        enddo
        call mpi_barrier( mpi_comm_world, ierr )
    endsubroutine fthmc_phy0_corFT

    subroutine fthmc_phy0_ft(cor, filename1, filename2, idx, P0, latt)
#IFDEF _OPENMP
        USE OMP_LIB
#ENDIF

        implicit none
        type(phy0), intent(in) :: P0
        integer, intent(in) :: idx
        complex(dp), intent(in), dimension(:) :: cor
        character(40), intent(in) :: filename1, filename2
        class(ftdqmc_latt), intent(in) :: latt

        ! local variables
        integer :: imj, iq, i, j, order_idx
        real(dp) :: qvec(2)
        complex(dp) :: cork
        complex(dp) :: background
        complex(dp), external :: zdotu

        background = czero

        ! determine the order vector index

        if ( idx .eq. ISPSM ) then ! den0
            order_idx = 1
        elseif ( idx .eq. IDIME ) then
            order_idx = lfam + 1
        endif

        ! prepare the output file
        open(unit=177,file=filename1,status='unknown', action="write", position="append")
        open(unit=178,file=filename2,status='unknown', action="write", position="append")

        ! loop over k points
        do iq = 1, lq
            qvec = dble( latt%listk(iq,1))*latt%b1_p + dble( latt%listk(iq,2))*latt%b2_p
            cork = czero

            ! use blas level 1 function: zdotu
            cork = zdotu(lq, cor(:), 1, cone/latt%zexpiqr(:, iq), 1)

            if ( iq .ne. order_idx ) then
                write(177, '(2f8.4, 2x, e16.8)') qvec(1), qvec(2), dble(cork*dcmplx(1.d0/dble(lq),0.d0))
            endif

            ! output the Gamma point
            if ( iq .eq. order_idx ) then
                write(177, '(2f8.4, 2x, e16.8)') qvec(1), qvec(2), dble(cork*dcmplx(1.d0/dble(lq),0.d0))
                write(178, '(e16.8)') dble(cork*dcmplx(1.d0/dble(lq),0.d0))
            endif

        end do
        close(178)
        close(177)
    endsubroutine fthmc_phy0_ft

endmodule fthmc_phy0
