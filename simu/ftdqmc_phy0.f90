#ifdef __GFORTRAN__
#define GET_FILENAME(x) generate_filename("x")
#else
#define GET_FILENAME(x) generate_filename(#x)
#endif

module ftdqmc_phy0
    use ftdqmc_hamilt
    use ftdqmc_auxfield_class
    use ftdqmc_auxfield_f1_class
    use ftdqmc_auxfield_f2_class
    use ftdqmc_latt
    use ftdqmc_gfun

    implicit none

    ! set number of scalars and correlation functions
    integer :: nscalars = 23      ! number of scalar properties, need to change accordingly
    integer :: ncorrelators = 15  ! number of correlators, den0 szsz swave dwave bond etc

    ! index of the scalar variables: IS
    integer, parameter :: IS_FILL   = 1   ! electron filling
    integer, parameter :: IS_SpinZ  = 2   ! SpinZ
    integer, parameter :: IS_EMU    = 3   ! chemical potential term
    integer, parameter :: IS_EK     = 4   ! kappa term
    integer, parameter :: IS_EJ     = 5   ! J term
    integer, parameter :: IS_EKtau  = 6   ! kappa_tau term
    integer, parameter :: IS_EE2    = 7   ! varepsilon2 term
    integer, parameter :: IS_EE4    = 8   ! varepsilon4 term
    integer, parameter :: IS_BI     = 9   ! <B> itself
    integer, parameter :: IS_QIJ2   = 10  ! <Q_ij^2>
    integer, parameter :: IS_JIJ2   = 11  ! <J_ij^2>
    integer, parameter :: IS_DIJ2   = 12  ! <|\Delta_ij|^2>
    integer, parameter :: IS_DtauB2 = 13  ! <|D_\tau B_i|^2>
    integer, parameter :: IS_EW     = 14  ! w term
    integer, parameter :: IS_EL_DEN   = 15  ! <O> for background cal
    integer, parameter :: IS_EL_BOND  = 16  ! <O> for background cal
    integer, parameter :: IS_EL_BONDY = 17  ! <O> for background cal
    integer, parameter :: IS_EL_BONDC = 18  ! <O> for background cal
    integer, parameter :: IS_EL_SWAV  = 19  ! <O> for background cal
    integer, parameter :: IS_EL_DWAV  = 20  ! <O> for background cal
    integer, parameter :: IS_EL_DBOND = 21  ! <O> for background cal
    integer, parameter :: IS_EL_DWAV2 = 22  ! for binder
    integer, parameter :: IS_EL_DWAV4 = 23  ! for binder

    ! index for the correlation variables: IC
    integer, parameter :: IC_DEN0  = 1  ! Density-Density correlation
    integer, parameter :: IC_SPZZ  = 2  ! Sz-Sz correlation
    integer, parameter :: IC_SWAV  = 3  ! Swave pairing-pairing correlation
    integer, parameter :: IC_DWAV  = 4  ! Dwave pairing-pairing correlation
    integer, parameter :: IC_BOND  = 5  ! bond-bond correlation for VBS
    integer, parameter :: IC_DIME  = 6  ! dimer-dimer correlation for VBS
    integer, parameter :: IC_BONDX = 7  ! bond background
    integer, parameter :: IC_EL_DEN  = 8  ! electron (chargon) density
    integer, parameter :: IC_EL_BOND = 9  ! electron (chargon) x-bond density
    integer, parameter :: IC_EL_BONDC= 10 ! electron (chargon) bond current
    integer, parameter :: IC_EL_SWAV = 11 ! electron (chargon) s-wave pairing
    integer, parameter :: IC_EL_DWAV = 12 ! electron (chargon) d-wave pairing
    integer, parameter :: IC_EL_DBOND = 13! electron (chargon) d-wave bond density
    integer, parameter :: IC_EL_BONDY = 14! electron (chargon) y-bond density
    integer, parameter :: IC_SPSM  = 15  ! S+S- correlation (missing the S-S+ part, okay with no Zeeman field)

    ! pointer array to store pointers
    type :: pnt1d
    	complex(dp), pointer :: pnt(:)
	end type pnt1d

    type phy0
        integer :: nobs         ! number of measurement
        integer :: nclass       ! number of distance pair

        ! array to store start and end indices for scalars and correlation functions
        integer, allocatable :: IARR(:)

        ! array for all physical properties including scalars and correlationn functions
        complex(dp), pointer :: AllProp_scalar(:)    ! vector of all physical properties
        complex(dp), pointer :: AllProp_correlator(:)    ! vector of all physical properties

        ! pointer to scalars and correlators, not mixed
        complex(dp), pointer :: pnt_sc(:)
		type(pnt1d), allocatable :: pnt_cr(:)

        ! flag
        logical :: init
    endtype phy0

contains

    subroutine ftdqmc_phy0_alloc(P0)
        type(phy0), intent(inout) :: P0 ! phy0 to be initalized

        ! local variables
        integer :: i, n

        ! allocate space for scalars
        allocate(P0%AllProp_scalar(nscalars))
        ! set pointers pnt_sc
        P0%pnt_sc => P0%AllProp_scalar( : )

        ! allocate space for correlators
        P0%nclass = lq ! number of relative distance pairs
        allocate(P0%AllProp_correlator(P0%nclass * ncorrelators))
        ! allocate IARR, which is useful for pnt_cr setup
        allocate(P0%IARR(ncorrelators+1))
        ! set indices in IARR
        P0%IARR(:) = 1
        do i = 1, ncorrelators+1
            P0%IARR(i) = 1 + (i-1)*P0%nclass
        enddo
        ! set pointers pnt_cr
        allocate(P0%pnt_cr(ncorrelators))
        do i = 1, ncorrelators
            ! new method
            P0%pnt_cr(i)%pnt => P0%AllProp_correlator( P0%IARR(i) : P0%IARR(i+1) - 1 )
        enddo

        P0%init = .true.
    endsubroutine ftdqmc_phy0_alloc


    subroutine ftdqmc_phy0_init(P0)
        type(phy0), intent(inout) :: P0 ! phy0 to be initalized

        ! local variables
        integer :: i

        P0%pnt_sc(:) = czero
        do i = 1, ncorrelators
            P0%pnt_cr(i)%pnt(:) = czero
        enddo

        P0%nobs   = 0
    endsubroutine ftdqmc_phy0_init

    subroutine ftdqmc_phy0_free(P0)
        type(phy0), intent(inout) :: P0 ! phy0 to be freed

        ! local variable
        integer :: i

        ! executable
        if ( P0%init ) then
            nullify(P0%pnt_sc)
            do i = 1, ncorrelators
                nullify(P0%pnt_cr(i)%pnt)
            enddo

            deallocate(P0%AllProp_scalar)
            deallocate(P0%AllProp_correlator)
            deallocate(P0%IARR)
        endif
    endsubroutine ftdqmc_phy0_free

    subroutine ftdqmc_phy0_meas(nt, gfun0, P0, phi)
#ifdef _OPENMP
        use OMP_LIB
#endif
        use ftdqmc_auxfield_f3_class

        implicit none
        integer, intent(in) :: nt
        type(phy0), intent(inout) :: P0
        type(gfun), target, intent(in) :: gfun0
        type(ftdqmc_auxfieldHolder), intent(inout) :: phi(num_fields)

        ! local
        integer :: i, j, imj, gauge_path, dis, disp, tmp_x, lf, nf, i1, i2, j1, it, at, ii, jj
        integer :: ix, iy, tmp_j, jx, jy, nx, ny, nix, niy, njx, njy, jax, jmx
        integer :: i_plaqA, i_plaqB, ib, inf, ilf, b_plaqA, b_plaqB, idx1, idx2, i0
        integer :: iay
        real(dp), pointer :: gauge_su2(:,:,:,:)
        complex(dp), pointer :: chargon(:,:,:)
        real(dp) :: flux_tot, eij_l, eij_l2, epsj_l, rtmp, rtmp1, dfactor
        complex(dp) :: phi_Aplaq(2,2), phi_Bplaq(2,2), tmp_mat(2,2), Mmat_tmp(2,2), Mmat_tmp2(2,2)
        complex(dp) :: zkint, zchop, zne, sz, ztmp1, ztmp2, ztmp3, ztmp4, ztmp, zehx, zej, zek, zektau, ztmp5, ztmp6
        complex(dp) :: zqij2, zjij2, zdij2, ze2, ze4, zew
        complex(dp) :: zbi
        complex(dp) :: z1, z2, z3, z4
        complex(dp), dimension(:,:), allocatable :: grupc(:,:), grdnc(:,:), grup(:,:), grdn(:,:)
        ! pairing green's function
        complex(dp), dimension(:,:), allocatable :: fupfdn(:,:), fdnfup(:,:), fddnfdup(:,:), fdupfddn(:,:)

        ! local for electron (chargon) observables
        integer :: nd, iax, imx, iaxay, iaxmy, ia2xa2y, ia2xm2y, nta
        complex(dp) :: rho(lq), bondx_den(lq), bondx_current(lq), elswave(lq), eldwave(lq), eldbond(lq)
        complex(dp) :: bondy_den(lq), zeij, Bidel(2), zd2bi
        complex(dp) :: zvec(2), BdUB, rhoiax, rhoiay, rhoiaxay, rhoiaxmy, rhoia2xa2y, rhoia2xm2y
        complex(dp) :: dSC_avg
        ! for test
        COMPLEX(dp) :: BduB1, BduB2

        ! external: npbc and zdotc
        integer, external :: npbc
        complex(dp), external :: zdotc

        ! for test
        integer :: i4, i5

        ! meas count
        P0%nobs = P0%nobs + 1
        ! initial pointer
        gauge_su2 => phi(idx_gauge)%pnt%gauge_su2
        chargon => phi(idx_chargon)%pnt%chargon

        ! setup green's function
        if (.not. lclassic .and. .not. lrmfermion) then
            ! allocate memory
            allocate(grup (lq, lq))
            allocate(grdn (lq, lq))
            allocate(grupc(lq, lq))
            allocate(grdnc(lq, lq))
            ! pairing
            allocate(fupfdn(lq, lq))
            allocate(fdnfup(lq, lq))
            allocate(fddnfdup(lq, lq))
            allocate(fdupfddn(lq, lq))

            ! convention
            !grup (i,j) = < c_i c^+_j >
            !grupc (i,j) = < c^+_i c_j >

            ! get up and down from grup (full Green's function)
            ! and fupfdn, fddnfdup
            do i =1, lq
                do j =1, lq
                    ! get grup
                    nx = list_prim(i, 1)
                    ny = list_prim(i, 2)
                    idx1 = invlist(nx, ny, 1)
                    nx = list_prim(j, 1)
                    ny = list_prim(j, 2)
                    idx2 = invlist(nx, ny, 1)
                    grup(i,j) = gfun0%grup(idx1, idx2)

                    ! get grdnc
                    nx = list_prim(i, 1)
                    ny = list_prim(i, 2)
                    idx1 = invlist(nx, ny, 2)
                    nx = list_prim(j, 1)
                    ny = list_prim(j, 2)
                    idx2 = invlist(nx, ny, 2)
                    grdnc(i,j) = gfun0%grup(idx1, idx2) / dcmplx(epsj(i), 0.d0) / dcmplx(epsj(j), 0.d0)

                    ! get fupfdn
                    nx = list_prim(i, 1)
                    ny = list_prim(i, 2)
                    idx1 = invlist(nx, ny, 1)
                    nx = list_prim(j, 1)
                    ny = list_prim(j, 2)
                    idx2 = invlist(nx, ny, 2)
                    fupfdn(i,j) = gfun0%grup(idx1, idx2) / dcmplx(epsj(j), 0.d0)

                    ! get fddnfdup
                    nx = list_prim(i, 1)
                    ny = list_prim(i, 2)
                    idx1 = invlist(nx, ny, 2)
                    nx = list_prim(j, 1)
                    ny = list_prim(j, 2)
                    idx2 = invlist(nx, ny, 1)
                    fddnfdup(i,j) = gfun0%grup(idx1, idx2) / dcmplx(epsj(i), 0.d0)
                enddo
            enddo

            ! get grupc and grdn
            ! and fdnfup, fdupfddn
            do i = 1, lq
                do j = 1, lq
                    ! normal Green's function
                    grupc(j,i) = - grup(i,j)
                    grdn(j,i)  = - grdnc(i,j)

                    ! abnormal part
                    fdnfup(j,i) = - fupfdn(i,j)
                    fdupfddn(j,i) = - fddnfdup(i,j)
                end do
                grupc(i,i) = grupc(i,i) + cone
                grdn(i,i)  = grdn(i,i)  + cone
            end do
        endif

        ! electron fillling, spinz and <B>
        zne = czero
        sz = czero
        zbi = czero
        do i = 1, lq
            if (.not. lclassic .and. .not. lrmfermion) then
                zne = zne + grupc(i,i) + grupc(i,i)
                sz = sz +  (grupc(i,i) - grupc(i,i)) * dcmplx(0.5d0, 0.d0)
            endif
            zbi = zbi + zdotc(2, chargon(:, i, nt), 1, chargon(:, i, nt), 1)
        end do
        P0%pnt_sc(IS_FILL) = P0%pnt_sc(IS_FILL) + zne/dcmplx(dble(lq), 0.d0)
        P0%pnt_sc(IS_EMU)  = P0%pnt_sc(IS_EMU) + dcmplx(-mu, 0.d0)*zne/dcmplx(dble(lq), 0.d0)
        P0%pnt_sc(IS_SpinZ)= P0%pnt_sc(IS_SpinZ) + sz/dcmplx(dble(lq), 0.d0)
        P0%pnt_sc(IS_BI)   = P0%pnt_sc(IS_BI) + zbi/dcmplx(dble(lq), 0.d0)

        ! get J(zkint, the kinetic energy), kappa term and kappa_tau term energy
        zkint = czero
        zek = czero
        zektau = czero
        do nf = 1, nfam
            do lf = 1, lfam
                i1 = l_bonds(1, lf, nf)
                i2 = l_bonds(2, lf, nf)

                ! get the kinetic energy, the J term
                ! get the SU(2) mat first
                call phi(idx_gauge)%pnt%quat_to_SU2mat_sub(gauge_su2(:,lf,nf,nt), Mmat_tmp)
                ! get the 2lq index
                nix = list_prim(i1, 1); niy = list_prim(i1, 2)
                njx = list_prim(i2, 1); njy = list_prim(i2, 2)
                ! set eij and epsj
                eij_l = eij(lf, nf)
                epsj_l = epsj(i2)

                if ( .not. lclassic .and. .not. lrmfermion) then
                    ! ij part
                    ! AA
                    idx1 = invlist(nix, niy, 1); idx2 = invlist(njx, njy, 1)
                    ztmp1 = dcmplx(0.d0, rj*eij_l ) * Mmat_tmp(1,1)                 * gfun0%grup(idx2, idx1) * (-cone)
                    ! AB
                    idx1 = invlist(nix, niy, 1); idx2 = invlist(njx, njy, 2)
                    ztmp2 = dcmplx(0.d0, rj*eij_l*epsj_l ) * Mmat_tmp(1,2)          * gfun0%grup(idx2, idx1) * (-cone)
                    ! BA
                    idx1 = invlist(nix, niy, 2); idx2 = invlist(njx, njy, 1)
                    ztmp3 = dcmplx(0.d0, rj*eij_l*epsj_l ) * dconjg(Mmat_tmp(1,2))  * gfun0%grup(idx2, idx1) * (-cone)
                    ! BB
                    idx1 = invlist(nix, niy, 2); idx2 = invlist(njx, njy, 2)
                    ztmp4 = dcmplx(0.d0, rj*eij_l*(-1.d0) ) * dconjg(Mmat_tmp(1,1)) * gfun0%grup(idx2, idx1) * (-cone)
                    ! add to zkint
                    zkint = zkint + ztmp1 + ztmp2 + ztmp3 + ztmp4
                    ! for test: pi flux case
                    !zkint = zkint + ztmp1
                    !zkint = zkint + ztmp4

                    ! ji part
                    ! AA
                    idx1 = invlist(nix, niy, 1); idx2 = invlist(njx, njy, 1)
                    ztmp1 = dcmplx(0.d0, -rj*eij_l ) * dconjg(Mmat_tmp(1,1))        * gfun0%grup(idx1, idx2) * (-cone)
                    ! AB
                    idx1 = invlist(nix, niy, 1); idx2 = invlist(njx, njy, 2)
                    ztmp2 = dcmplx(0.d0, -rj*eij_l*epsj_l ) * dconjg(Mmat_tmp(1,2)) * gfun0%grup(idx1, idx2) * (-cone)
                    ! BA
                    idx1 = invlist(nix, niy, 2); idx2 = invlist(njx, njy, 1)
                    ztmp3 = dcmplx(0.d0, -rj*eij_l*epsj_l ) * Mmat_tmp(1,2)         * gfun0%grup(idx1, idx2) * (-cone)
                    ! BB
                    idx1 = invlist(nix, niy, 2); idx2 = invlist(njx, njy, 2)
                    ztmp4 = dcmplx(0.d0, -rj*eij_l*(-1.d0) ) * Mmat_tmp(1,1)        * gfun0%grup(idx1, idx2) * (-cone)
                    ! add to zkint
                    zkint = zkint + ztmp1 + ztmp2 + ztmp3 + ztmp4
                    ! for test: pi flux case
                    !zkint = zkint + ztmp1
                    !zkint = zkint + ztmp4
                endif

                ! get the kappa term, Wilson action
                phi_Aplaq = I2mat
                i_plaqA = inv_Aplaq_bondcord(1, lf, nf)
                b_plaqA = inv_Aplaq_bondcord(2, lf, nf)
                do ib = 1, 4
                    inf = plaq_bondcord(1, ib, i_plaqA)
                    ilf = plaq_bondcord(2, ib, i_plaqA)

                    call phi(idx_gauge)%pnt%quat_to_SU2mat_sub(gauge_su2(:,ilf,inf,nt), Mmat_tmp)
                    !if (sgnA_plaq(ib) .eq. 1) then
                    if (inf .eq. 1 .or. inf .eq. 3) then
                        call zgemm('N','N', 2, 2, 2, cone, Mmat_tmp, 2, phi_Aplaq, 2, czero, tmp_mat, 2)
                        phi_Aplaq(:,:) = tmp_mat(:,:)
                    else
                        call zgemm('C','N', 2, 2, 2, cone, Mmat_tmp, 2, phi_Aplaq, 2, czero, tmp_mat, 2)
                        phi_Aplaq(:,:) = tmp_mat(:,:)
                    endif
                enddo

                phi_Bplaq = I2mat
                i_plaqB = inv_Bplaq_bondcord(1, lf, nf)
                b_plaqB = inv_Bplaq_bondcord(2, lf, nf)
                do ib = 1, 4
                    inf = plaq_bondcord(1, ib, i_plaqB)
                    ilf = plaq_bondcord(2, ib, i_plaqB)

                    call phi(idx_gauge)%pnt%quat_to_SU2mat_sub(gauge_su2(:,ilf,inf,nt), Mmat_tmp)
                    !if (sgnB_plaq(ib) .eq. 1) then
                    if (inf .eq. 1 .or. inf .eq. 3) then
                        call zgemm('N','N', 2, 2, 2, cone, Mmat_tmp, 2, phi_Bplaq, 2, czero, tmp_mat, 2)
                        phi_Bplaq(:,:) = tmp_mat(:,:)
                    else
                        call zgemm('C','N', 2, 2, 2, cone, Mmat_tmp, 2, phi_Bplaq, 2, czero, tmp_mat, 2)
                        phi_Bplaq(:,:) = tmp_mat(:,:)
                    endif
                enddo

                zek = zek + rk*(1.d0 - 0.5d0*dble(phi_Aplaq(1,1)+phi_Aplaq(2,2))) &
                    + rk*(1.d0 - 0.5d0*dble(phi_Bplaq(1,1)+phi_Bplaq(2,2)))

                ! get the kappa_tau part, nt and nt+1, with peroidic boundary condition in time direction
                if (.not. lclassic) then
                    call phi(idx_gauge)%pnt%quat_to_SU2mat_sub(gauge_su2(:,lf,nf,nt), Mmat)
                    call phi(idx_gauge)%pnt%quat_to_SU2mat_sub(gauge_su2(:,lf,nf,npbc(nt+1, ltrot)), Mmat_tmp)
                    call zgemm('C','N', 2, 2, 2, cone, Mmat, 2, Mmat_tmp, 2, czero, tmp_mat, 2)
                    zektau = zektau + 1.d0/(rktau*dtau*dtau)*(1.d0 - 0.5d0*dble(tmp_mat(1,1) + tmp_mat(2,2) ))
                endif

            enddo
        enddo
        zek = zek*dcmplx(0.25d0, 0.d0) ! 0.25 because of 4 times counting
        zkint = zkint * dcmplx(Nflavor/2.d0, 0.d0)
        P0%pnt_sc(IS_EJ)    = P0%pnt_sc(IS_EJ) + zkint/dcmplx(dble(lq), 0.d0)
        P0%pnt_sc(IS_EK)    = P0%pnt_sc(IS_EK) + zek  /dcmplx(dble(lq), 0.d0)
        P0%pnt_sc(IS_EKtau) = P0%pnt_sc(IS_EKtau) + zektau /dcmplx(dble(lq), 0.d0)

        ! electron density-density, sz-sz and pairing-pairing correlation, bond-bond and dimer-dimer correlation function
        if (.not. lclassic .and. .not. lrmfermion) then
            z2 = dcmplx( Nflavor*Nflavor-1.d0, 0.d0 )
            z4 = z2*z2
            z3 = dcmplx( Nflavor*Nflavor*Nflavor - 2.d0*Nflavor + 1.d0/Nflavor, 0.d0 )
            z1 = dcmplx( -Nflavor + 1.d0/Nflavor, 0.d0 )
            do j = 1, lq
                jax = nnlist(j,1) ! j+x
                jmx = nnlist(j,3) ! j-x
                do i = 1, lq
                    iax = nnlist(i,1) ! i+x
                    imx = nnlist(i,3) ! i-x
                    imj = latt_imj(i,j)

                    ! density-density
                    ztmp1 =   grupc(i,i)*grupc(j,j) + grupc(i,i)*grdnc(j,j) &
                            + grdnc(i,i)*grupc(j,j) + grdnc(i,i)*grdnc(j,j) &
                            + grupc(i,j)*grup(i,j)  + grdnc(i,j)*grdn(i,j)
                    ! sz-sz, from u1sl_mag code
                    ztmp2 = (grupc(i,i)-grdnc(i,i))*(grupc(j,j)-grdnc(j,j)) + (grupc(i,j)*grup(i,j)+grdnc(i,j)*grdn(i,j))
                    ztmp2 = ztmp2 * dcmplx(0.25d0, 0.d0)
                    ! swave pairing-pairing
                    ztmp3 = grupc(i,j)*grupc(i,j)
                    ! d-wave pairing correlation functions
                    ztmp = czero
                    do i1 = 1, 4
                        do j1 = 1, 4
                            ii = nnlist(i, i1)
                            jj = nnlist(j, j1)
                            ztmp = ztmp + (-1)**(i1+1) * (-1)**(j1+1) * grupc(ii,jj)*grupc(i,j)
                        enddo
                    enddo
                    ! dimer-dimer
                    ztmp5 = grupc(i,iax)*grup(i,iax)*grupc(j,jax)  *grup(j,jax)  *z4  &
                                                      + grupc(i,j)  *grup(i,j)  *grupc(iax,jax)*grup(iax,jax)*z2  &
                                                      + grupc(i,jax)*grup(i,jax)*grupc(iax,j)  *grup(iax,j)  *z2  &
                                                      + grupc(i,jax)*grup(i,iax)*grup(iax,j)   *grup(j,jax)  *z3  &
                                                      + grupc(i,iax)*grup(i,jax)*grup(j,iax)   *grup(jax,j)  *z3  &
                                                      + grupc(i,j)  *grup(i,iax)*grup(iax,jax) *grup(jax,j)  *z3  &
                                                      + grupc(i,iax)*grup(i,j)  *grup(jax,iax) *grup(j,jax)  *z3  &
                                                      + grupc(i,jax)*grup(i,j)  *grup(iax,jax) *grup(j,iax)  *z1  &
                                                      + grupc(i,j)  *grup(i,jax)*grup(jax,iax) *grup(iax,j)  *z1

                    ! spsm, wick decomposition with pairing part
                    ztmp6 = (grupc(i,j)*grdn(i,j) - fdupfddn(i,j)*fdnfup(i,j))*dcmplx(1.5d0, 0.d0)
                    ! for test
                    !write(111,'(4e16.8)') fdupfddn(i,j), fdnfup(i,j)

                    P0%pnt_cr(IC_DEN0)%pnt(imj) = P0%pnt_cr(IC_DEN0)%pnt(imj) + ztmp1
                    P0%pnt_cr(IC_SPZZ)%pnt(imj) = P0%pnt_cr(IC_SPZZ)%pnt(imj) + ztmp2
                    P0%pnt_cr(IC_SWAV)%pnt(imj) = P0%pnt_cr(IC_SWAV)%pnt(imj) + ztmp3
                    P0%pnt_cr(IC_DWAV)%pnt(imj) = P0%pnt_cr(IC_DWAV)%pnt(imj) + ztmp
                    P0%pnt_cr(IC_DIME)%pnt(imj) = P0%pnt_cr(IC_DIME)%pnt(imj) + ztmp5
                    P0%pnt_cr(IC_SPSM)%pnt(imj) = P0%pnt_cr(IC_SPSM)%pnt(imj) + ztmp6
                end do
            end do

            ! calculate the bond-bond correlation function
            do j = 1, lq
                jax = nnlist(j,1) ! j+x
                jmx = nnlist(j,3) ! j-x

                jx = list_prim(j,1); jy = list_prim(j,2)
                if( mod(jx+jy,2) .eq. 0 ) then
                    nf = 1
                    i0 = allsites_bonds(nf,j)

                    call phi(idx_gauge)%pnt%quat_to_SU2mat_sub(gauge_su2(:,i0,nf,nt), Mmat_tmp)
                    eij_l = eij(i0, nf)
                    z1 = eij_l * Mmat_tmp(1,1)
                    !z1 = zep_rsigl_k(i0,nf,nt) /cinvsqrt2

                    ztmp1 = z1 * grupc(j,jax)
                else
                    nf = 3
                    i0 = allsites_bonds(nf,j)

                    call phi(idx_gauge)%pnt%quat_to_SU2mat_sub(gauge_su2(:,i0,nf,nt), Mmat_tmp)
                    eij_l = eij(i0, nf)
                    z1 = eij_l * Mmat_tmp(1,1)
                    !z1 = dconjg( zep_rsigl_k(i0,nf,nt) ) /cinvsqrt2

                    ztmp1 = z1 * grupc(j,jax)
                end if

                do i = 1, lq
                    imj = latt_imj(i,j)

                    iax = nnlist(i,1) ! i+x
                    imx = nnlist(i,3) ! i-x

                    ix = list_prim(i,1); iy = list_prim(i,2)
                    if( mod(ix+iy,2) .eq. 0 ) then
                        nf = 1
                        i0 = allsites_bonds(nf,i)

                        call phi(idx_gauge)%pnt%quat_to_SU2mat_sub(gauge_su2(:,i0,nf,nt), Mmat_tmp)
                        eij_l = eij(i0, nf)
                        z2 = eij_l * Mmat_tmp(1,1)
                        !z2 = zep_rsigl_k(i0,nf,nt) /cinvsqrt2

                        ztmp2 = z2 * grupc(i,iax)
                    else
                        nf = 3
                        i0 = allsites_bonds(nf,i)

                        call phi(idx_gauge)%pnt%quat_to_SU2mat_sub(gauge_su2(:,i0,nf,nt), Mmat_tmp)
                        eij_l = eij(i0, nf)
                        z2 = eij_l * Mmat_tmp(1,1)
                        !z2 = dconjg( zep_rsigl_k(i0,nf,nt) ) /cinvsqrt2

                        ztmp2 = z2 * grupc(i,iax)
                    end if

                    ! get bondx and bond-bond
                    if(j.eq.1) P0%pnt_cr(IC_BONDX)%pnt(i) = ztmp2 + dconjg(ztmp2) ! only need to store once
                    P0%pnt_cr(IC_BOND)%pnt(imj) = P0%pnt_cr(IC_BOND)%pnt(imj) + (ztmp1+dconjg(ztmp1))*(ztmp2+dconjg(ztmp2))*dcmplx(Nflavor*Nflavor,0.d0) &
                                                    + (  z1*z2*grupc(i,jax)*grup(iax,j) + z1/z2*grupc(i,j)*grup(iax,jax) &
                                                       + z2/z1*grupc(iax,jax)*grup(i,j) + grupc(iax,j)*grup(i,jax)/z1/z2 )*dcmplx(Nflavor,0.d0)
                end do
            end do
        endif

        ! calculate electron correlation fucntion from chargon as well as calculate energy in \varepsilon_2(U, B) and \varepsilon_4(U, B)
        zqij2 = czero
        zjij2 = czero
        zdij2 = czero
        ze2 = czero
        ze4 = czero
        zd2bi = czero
        zew = czero
        if (irank .eq. 0 .and. loutput) then
            open(unit=10, file='Qiiax.dat', status='replace')
            open(unit=11, file='Diiax.dat', status='replace')
            open(unit=12, file='Jiiax.dat', status='replace')
        endif
        do i = 1, lq
            ! get (x, y) cord
            ix = list_prim(i,1); iy = list_prim(i,2)

            ! on-site density: r and u terms
            rho(i) = zdotc(2, chargon(:, i, nt), 1, chargon(:, i, nt), 1)
            ze2 = ze2 + (rr+sqrt(8.d0)*rw) * rho(i)
            ze4 = ze4 + ru/2.d0 * rho(i)**2

            ! density-density interaction: v1, v11, v22 terms
            ! v1
            iax = nnlist(i, 1)
            iay = nnlist(i, 2)
            rhoiax = zdotc(2, chargon(:, iax, nt), 1, chargon(:, iax, nt), 1)
            rhoiay = zdotc(2, chargon(:, iay, nt), 1, chargon(:, iay, nt), 1)
            ze4 = ze4 + rv1 * rho(i)*(rhoiax + rhoiay)
            ! v11
            iaxay = nnlist(i, 5)
            iaxmy = nnlist(i, 8)
            rhoiaxay = zdotc(2, chargon(:, iaxay, nt), 1, chargon(:, iaxay, nt), 1)
            rhoiaxmy = zdotc(2, chargon(:, iaxmy, nt), 1, chargon(:, iaxmy, nt), 1)
            ze4 = ze4 + rv11 * rho(i) * (rhoiaxay + rhoiaxmy)
            ! v22
            nx = list_prim(i,1); ny = list_prim(i,2)
            ia2xa2y = invlist_prim( npbc(nx+2,l) , npbc(ny+2,l) )
            ia2xm2y = invlist_prim( npbc(nx+2,l) , npbc(ny-2,l) )
            rhoia2xa2y = zdotc(2, chargon(:, ia2xa2y, nt), 1, chargon(:, ia2xa2y, nt), 1)
            rhoia2xm2y = zdotc(2, chargon(:, ia2xm2y, nt), 1, chargon(:, ia2xm2y, nt), 1)
            ze4 = ze4 + rv22 * rho(i) * (rhoia2xa2y + rhoia2xm2y)

            ! bond: density, current, pairing
            ! x-direction bonds
            ! get lf
            if( mod(ix+iy,2) .eq. 0 ) then
                nf = 1
                ii = i
                jj = nnlist(i,1) ! i+x
            else
                nf = 3
                ii = nnlist(i,1) ! i+x
                jj = i
            end if
            i0 = allsites_bonds(nf,i) ! lf
            ! set eij
            eij_l = eij(i0, nf)
            zeij = dcmplx(eij_l, 0.d0)

            ! get Uij
            call phi(idx_gauge)%pnt%quat_to_SU2mat_sub(gauge_su2(:,i0,nf,nt), Mmat_tmp)

            ! cal BdUB
            call zgemm('N', 'N', 2, 1, 2, cone, Mmat_tmp, 2, chargon(:, jj, nt), 2, czero, zvec, 2)
            BdUB = zdotc(2, chargon(:, ii, nt ), 1, zvec(:), 1)

            ! set correlator: bond_density and bond_current
            bondx_den(i) = aimag(BdUB*dcmplx(eij_l, 0.d0))
            bondx_current(i) = dble(BdUB*dcmplx(eij_l, 0.d0))

            ! cal energy
            ! w term
            ze2 = ze2 + (BduB - dconjg(BdUB)) * dcmplx(eij_l, 0.d0) * dcmplx(0.d0, rw)
            zew = zew + (BduB - dconjg(BdUB)) * dcmplx(eij_l, 0.d0) * dcmplx(0.d0, rw)
            ! j1
            ze4 = ze4 + rj1 * bondx_den(i)**2
            zqij2 = zqij2 + bondx_den(i)**2
            ! k1
            ze4 = ze4 + rk1 * bondx_current(i)**2
            zjij2 = zjij2 + bondx_current(i)**2
            ! g
            ztmp = chargon(1, ii, nt) * Mmat_tmp(2,1) * chargon(1, jj, nt) &
                 + chargon(1, ii, nt) * Mmat_tmp(2,2) * chargon(2, jj, nt) &
                 - chargon(2, ii, nt) * Mmat_tmp(1,1) * chargon(1, jj, nt) &
                 - chargon(2, ii, nt) * Mmat_tmp(1,2) * chargon(2, jj, nt)
            ztmp = ztmp * zeij
            ze4 = ze4 + rg * ztmp*dconjg(ztmp)
            zdij2 = zdij2 + ztmp*dconjg(ztmp)

            ! for test: output Q_{i, i+x} and D_{i, i+x}
            if (irank .eq. 0 .and. loutput) then
                write(10,'(2i5, e16.8)') list_prim(i, 1), list_prim(i, 2), dble(bondx_den(i))
                write(11,'(2i5, 2e16.8)') list_prim(i, 1), list_prim(i, 2), cdabs(ztmp), datan2(dimag(ztmp), dble(ztmp))
                write(12,'(2i5, e16.8)') list_prim(i, 1), list_prim(i, 2), dble(bondx_current(i))
            endif

            ! for test: output Bi1 amp and phase
            if (loutput) then
                write(333,*) list_prim(i, 1), list_prim(i, 2), cdabs(chargon(1, i, nt))

                rtmp =datan2(dimag(chargon(1, i, nt)), dble(chargon(1, i, nt)))
                rtmp = rtmp * 180.0d0 / acos(-1.0d0)
                if ( rtmp .lt. 0.d0) then
                    rtmp = rtmp + 360.0d0
                endif
                rtmp1 =datan2(dimag(chargon(2, i, nt)), dble(chargon(2, i, nt)))
                rtmp1 = rtmp1 * 180.0d0 / acos(-1.0d0)
                if ( rtmp1 .lt. 0.d0) then
                    rtmp1 = rtmp1 + 360.0d0
                endif
                write(444,*) list_prim(i, 1), list_prim(i, 2), rtmp, rtmp1

                ! output w term in x-direction
                write(555,*) list_prim(i, 1), list_prim(i, 2), dble((BduB - dconjg(BdUB)) * dcmplx(eij_l, 0.d0) * dcmplx(0.d0, rw))
                write(999,'(2f8.3, e16.8)') dble(list_prim(i, 1)), dble(list_prim(i, 2)), dble((BduB - dconjg(BdUB)) * dcmplx(eij_l, 0.d0) * dcmplx(0.d0, 1.d0))/2.d0/2.d0

                ! output phase diff of bi bj from BduB
                rtmp = asin(dimag((BduB1 - dconjg(BdUB1)))/2.d0)
                rtmp = rtmp * 180.0d0 / acos(-1.0d0)
                rtmp1 = asin(dimag((BduB2 - dconjg(BdUB2)))/2.d0)
                rtmp1 = rtmp1 * 180.0d0 / acos(-1.0d0)
                write(666,*) list_prim(i, 1), list_prim(i, 2), rtmp, rtmp1
            endif

            ! y-direction bonds
            ! get lf
            if( mod(ix+iy,2) .eq. 0 ) then
                nf = 2
                ii = i
                jj = nnlist(i,2)
            else
                nf = 4
                ii = nnlist(i,2)
                jj = i
            end if
            i0 = allsites_bonds(nf,i) ! lf
            ! set eij
            eij_l = eij(i0, nf)
            zeij = dcmplx(eij_l, 0.d0)

            ! get Uij
            call phi(idx_gauge)%pnt%quat_to_SU2mat_sub(gauge_su2(:,i0,nf,nt), Mmat_tmp)

            ! cal BdUB
            call zgemm('N', 'N', 2, 1, 2, cone, Mmat_tmp, 2, chargon(:, jj, nt), 2, czero, zvec, 2)
            BdUB = zdotc(2, chargon(:, ii, nt ), 1, zvec(:), 1)

            ! set correlator: bond_density
            bondy_den(i) = aimag(BdUB*dcmplx(eij_l, 0.d0))

            ! cal energy
            ! w term
            ze2 = ze2 + (BduB - dconjg(BdUB)) * dcmplx(eij_l, 0.d0) * dcmplx(0.d0, rw)
            zew = zew + (BduB - dconjg(BdUB)) * dcmplx(eij_l, 0.d0) * dcmplx(0.d0, rw)
            ! j1
            ze4 = ze4 + rj1 * bondy_den(i)**2
            zqij2 = zqij2 + bondy_den(i)**2
            ! k1
            ze4 = ze4 + rk1 * (dble(BdUB*dcmplx(eij_l, 0.d0)))**2
            zjij2 = zjij2 + (dble(BdUB*dcmplx(eij_l, 0.d0)))**2
            if (irank .eq. 0 .and. loutput) then
                write(12,'(2f8.3, e16.8)') list_prim(i, 1)-0.5d0, list_prim(i, 2)+0.5d0, dble(BdUB*dcmplx(eij_l, 0.d0))
            endif
            ! g
            ztmp = chargon(1, ii, nt) * Mmat_tmp(2,1) * chargon(1, jj, nt) &
                 + chargon(1, ii, nt) * Mmat_tmp(2,2) * chargon(2, jj, nt) &
                 - chargon(2, ii, nt) * Mmat_tmp(1,1) * chargon(1, jj, nt) &
                 - chargon(2, ii, nt) * Mmat_tmp(1,2) * chargon(2, jj, nt)
            ztmp = ztmp * zeij
            ze4 = ze4 + rg * ztmp*dconjg(ztmp)
            zdij2 = zdij2 + ztmp*dconjg(ztmp)

            if (loutput) then
                ! output iw term in the y-direction
                write(555,*) list_prim(i, 1), list_prim(i, 2), dble((BduB - dconjg(BdUB)) * dcmplx(eij_l, 0.d0) * dcmplx(0.d0, rw))
                write(999,'(2f8.3, e16.8)') list_prim(i, 1)-0.5d0, list_prim(i, 2)+0.5d0, dble((BduB - dconjg(BdUB)) * dcmplx(eij_l, 0.d0) * dcmplx(0.d0, 1.d0))/2.d0/2.d0
                write(*,'(2i5, 2f16.8)') list_prim(i, 1), list_prim(i, 2), (BduB - dconjg(BdUB)) * dcmplx(eij_l, 0.d0)

                ! output phase diff of bi bj from BduB
                rtmp = asin(dimag((BduB - dconjg(BdUB)))/4.d0)
                rtmp = rtmp * 180.0d0 / acos(-1.0d0)
                rtmp1 = asin(dimag((BduB2 - dconjg(BdUB2)))/2.d0)
                rtmp1 = rtmp1 * 180.0d0 / acos(-1.0d0)
                write(666,*) list_prim(i, 1), list_prim(i, 2), rtmp, rtmp1
            endif

            ! kinetic term
            if (.not. lclassic) then
                nta = npbc(nt+1, ltrot)
                Bidel(:) = chargon(:, i, nt) - chargon(:, i, nta)
                ztmp = zdotc(2, Bidel(:), 1, Bidel(:), 1)
                zd2bi = zd2bi + 1.d0/rm * 1.d0/(dtau)**2 * ztmp
            endif

            ! prepare d-wave pairing and bond density with d form factor
            ztmp1 = czero
            ztmp2 = czero
            ztmp3 = czero
            do nd = 1, 4
                ! get the ij bond index
                i0 = allsites_bonds(nd, i)

                ! get Uij
                call phi(idx_gauge)%pnt%quat_to_SU2mat_sub(gauge_su2(:,i0,nd,nt), Mmat_tmp)
                ! set eij
                eij_l = eij(i0, nd)

                ! get the j index and set eij_l
                if( mod(ix+iy,2) .eq. 0 ) then
                    ii = i
                    jj = l_bonds(2, i0, nd)
                else
                    ii = l_bonds(1, i0, nd)
                    jj = i
                endif

                !dfactor = (-1)**(i+1) * (-1)**(j+1)
                if (nd .eq. 1 .or. nd .eq. 3) then
                    dfactor = 1.d0
                else
                    dfactor = -1.d0
                endif

                ! d-wave pairing
                ztmp1 = ztmp1 + dfactor * dcmplx(eij_l, 0.d0) * ( chargon(1, ii, nt) * Mmat_tmp(2,1) * chargon(1, jj, nt) &
                                                                + chargon(1, ii, nt) * Mmat_tmp(2,2) * chargon(2, jj, nt) &
                                                                - chargon(2, ii, nt) * Mmat_tmp(1,1) * chargon(1, jj, nt) &
                                                                - chargon(2, ii, nt) * Mmat_tmp(1,2) * chargon(2, jj, nt) &
                                                                )
                ! s-wave pairing
                ztmp3 = ztmp3 + dcmplx(eij_l, 0.d0) * (chargon(1, ii, nt) * Mmat_tmp(2,1) * chargon(1, jj, nt) &
                                                     + chargon(1, ii, nt) * Mmat_tmp(2,2) * chargon(2, jj, nt) &
                                                     - chargon(2, ii, nt) * Mmat_tmp(1,1) * chargon(1, jj, nt) &
                                                     - chargon(2, ii, nt) * Mmat_tmp(1,2) * chargon(2, jj, nt) &
                                                     )
                ! bond current, new meas
                call zgemm('N', 'N', 2, 1, 2, cone, Mmat_tmp, 2, chargon(:, jj, nt), 2, czero, zvec, 2)
                BdUB = zdotc(2, chargon(:, ii, nt ), 1, zvec(:), 1)
                ! extra factor not necessary anymore
                !ztmp2 = ztmp2 + dfactor * dble(BdUB*dcmplx(eij_l, 0.d0)) * (-1.d0)**(ix+iy)
                ztmp2 = ztmp2 + dfactor * dble(BdUB*dcmplx(eij_l, 0.d0))
            enddo
            eldwave(i) = ztmp1
            eldbond(i) = ztmp2
            elswave(i) = ztmp3

            ! for test
            !write(111, "(2i5, 2f16.8)") list_prim(i, 1), list_prim(i, 2), eldwave(i)
        enddo
        ! for test
        !stop

        P0%pnt_sc(IS_QIJ2) = P0%pnt_sc(IS_QIJ2) + zqij2/dcmplx(dble(lq)*2.d0, 0.d0)
        P0%pnt_sc(IS_JIJ2) = P0%pnt_sc(IS_JIJ2) + zjij2/dcmplx(dble(lq)*2.d0, 0.d0)
        P0%pnt_sc(IS_DIJ2) = P0%pnt_sc(IS_DIJ2) + zdij2/dcmplx(dble(lq)*2.d0, 0.d0)
        P0%pnt_sc(IS_EE2) = P0%pnt_sc(IS_EE2) + ze2/dcmplx(dble(lq), 0.d0)
        P0%pnt_sc(IS_EE4) = P0%pnt_sc(IS_EE4) + ze4/dcmplx(dble(lq), 0.d0)
        P0%pnt_sc(IS_DtauB2) = P0%pnt_sc(IS_DtauB2) + zd2bi/dcmplx(dble(lq), 0.d0)
        P0%pnt_sc(IS_EW) = P0%pnt_sc(IS_EW) + zew/dcmplx(dble(lq), 0.d0)
        ! get <O>
        P0%pnt_sc(IS_EL_DEN) = P0%pnt_sc(IS_EL_DEN) + sum(rho(:))/dble(lq)
        P0%pnt_sc(IS_EL_BOND) = P0%pnt_sc(IS_EL_BOND) + sum(bondx_den(:))/dble(lq)
        P0%pnt_sc(IS_EL_BONDY) = P0%pnt_sc(IS_EL_BONDY) + sum(bondy_den(:))/dble(lq)
        P0%pnt_sc(IS_EL_BONDC) = P0%pnt_sc(IS_EL_BONDC) + sum(bondx_current(:))/dble(lq)
        P0%pnt_sc(IS_EL_DBOND) = P0%pnt_sc(IS_EL_DBOND) + sum(eldbond(:))/dble(lq)
        ! electron pairing
        P0%pnt_sc(IS_EL_SWAV)  = P0%pnt_sc(IS_EL_SWAV)  + sum(elswave(:))/dble(lq)
        dSC_avg = sum(eldwave(:))/dble(lq)
        P0%pnt_sc(IS_EL_DWAV)  = P0%pnt_sc(IS_EL_DWAV)  + dSC_avg
        ! for dSC binder ratio
        !P0%pnt_sc(IS_EL_DWAV2) = P0%pnt_sc(IS_EL_DWAV2) + sum((cdabs(eldwave(:)))**2)/dble(lq)
        !P0%pnt_sc(IS_EL_DWAV4) = P0%pnt_sc(IS_EL_DWAV4) + sum((cdabs(eldwave(:)))**4)/dble(lq)
        P0%pnt_sc(IS_EL_DWAV2) = P0%pnt_sc(IS_EL_DWAV2) + cdabs(dSC_avg)**2
        P0%pnt_sc(IS_EL_DWAV4) = P0%pnt_sc(IS_EL_DWAV4) + cdabs(dSC_avg)**4

        if (irank .eq. 0 .and. loutput) then
            close(10)
            close(11)
            close(12)
            stop
        endif


        ! get the correlation functions
        do j = 1, lq
            do i = 1, lq
                imj = latt_imj(i,j)

                P0%pnt_cr(IC_EL_DEN)%pnt(imj) = P0%pnt_cr(IC_EL_DEN)%pnt(imj) + rho(i) * rho(j)
                P0%pnt_cr(IC_EL_BOND)%pnt(imj) = P0%pnt_cr(IC_EL_BOND)%pnt(imj) + bondx_den(i) * bondx_den(j)
                P0%pnt_cr(IC_EL_BONDY)%pnt(imj) = P0%pnt_cr(IC_EL_BONDY)%pnt(imj) + bondy_den(i) * bondy_den(j)
                P0%pnt_cr(IC_EL_BONDC)%pnt(imj) = P0%pnt_cr(IC_EL_BONDC)%pnt(imj) + bondx_current(i) * bondx_current(j)
                P0%pnt_cr(IC_EL_DBOND)%pnt(imj) = P0%pnt_cr(IC_EL_DBOND)%pnt(imj) + eldbond(i) * eldbond(j)
                P0%pnt_cr(IC_EL_SWAV)%pnt(imj) = P0%pnt_cr(IC_EL_SWAV)%pnt(imj) + (elswave(i)) * dconjg(elswave(j))
                P0%pnt_cr(IC_EL_DWAV)%pnt(imj) = P0%pnt_cr(IC_EL_DWAV)%pnt(imj) + (eldwave(i)) * dconjg(eldwave(j))
            end do
        end do


        if (.not. lclassic .and. .not. lrmfermion) then
            deallocate(grupc)
            deallocate(grdnc)
            deallocate(grup)
            deallocate(grdn)
            ! pairing
            deallocate(fupfdn)
            deallocate(fdnfup)
            deallocate(fddnfdup)
            deallocate(fdupfddn)
        endif

        nullify(gauge_su2)
        nullify(chargon)
    endsubroutine ftdqmc_phy0_meas

    subroutine ftdqmc_phy0_getavg(P0)
        ! scalar quantities
        use mpi
        type(phy0), intent(inout) :: P0

        ! local variables
        integer :: i
        complex(dp) :: mpi_meas_bin(nscalars)

        ! executable
        call mpi_reduce( P0%pnt_sc, mpi_meas_bin, nscalars, mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr)
        P0%pnt_sc = mpi_meas_bin

        if ( irank .eq. 0) then
            ! average
            P0%pnt_sc = P0%pnt_sc / dcmplx( dble(isize * P0%nobs ), 0.d0)

            ! output
            open( unit=90, file='ener1.bin', status='unknown', action='write', position='append')
            write(90, '(20(e16.8, 2x))') dble(P0%pnt_sc(IS_FILL))  , & ! col 1
                                         dble(P0%pnt_sc(IS_SpinZ)) , & ! col 2
                                         dble(P0%pnt_sc(IS_EMU))   , & ! col 3
                                         dble(P0%pnt_sc(IS_EK))    , & ! col 4
                                         dble(P0%pnt_sc(IS_EJ))    , & ! col 5
                                         dble(P0%pnt_sc(IS_EKtau)) , & ! col 6
                                         dble(P0%pnt_sc(IS_BI))    , & ! col 7
                                         dble(P0%pnt_sc(IS_QIJ2))  , & ! col 8
                                         dble(P0%pnt_sc(IS_JIJ2))  , & ! col 9
                                         dble(P0%pnt_sc(IS_DIJ2))  , & ! col 10
                                         dble(P0%pnt_sc(IS_EE2))   , & ! col 11
                                         dble(P0%pnt_sc(IS_EE4))   , & ! col 12
                                         dble(P0%pnt_sc(IS_DtauB2)), & ! col 13
                                         dble(P0%pnt_sc(IS_EW))    , & ! col 14
                                         1.d0 - dble(P0%pnt_sc(IS_EL_DWAV4))/(2.d0*dble(P0%pnt_sc(IS_EL_DWAV2))**2) ! col 15
            close(90)
        endif

        call mpi_barrier( mpi_comm_world, ierr )
    endsubroutine ftdqmc_phy0_getavg

    subroutine ftdqmc_phy0_corFT(P0)
        type(phy0), intent(inout) :: P0

        ! one by one, ft as needed
        if (.not. lclassic .and. .not. lrmfermion) then
            call ftdqmc_phy0_perfFT(P0, IC_DEN0, GET_FILENAME(IC_DEN0))
            call ftdqmc_phy0_perfFT(P0, IC_SPZZ, GET_FILENAME(IC_SPZZ))
            call ftdqmc_phy0_perfFT(P0, IC_SWAV, GET_FILENAME(IC_SWAV))
            call ftdqmc_phy0_perfFT(P0, IC_DWAV, GET_FILENAME(IC_DWAV))
            call ftdqmc_phy0_perfFT(P0, IC_BOND, GET_FILENAME(IC_BOND))
            !call ftdqmc_phy0_perfFT(P0, IC_BONDX, GET_FILENAME(IC_BONDX))
            call ftdqmc_phy0_perfFT(P0, IC_DIME, GET_FILENAME(IC_DIME))
            call ftdqmc_phy0_perfFT(P0, IC_SPSM, GET_FILENAME(IC_SPSM))
        endif
        call ftdqmc_phy0_perfFT(P0, IC_EL_DEN, GET_FILENAME(IC_EL_DEN))
        call ftdqmc_phy0_perfFT(P0, IC_EL_BOND, GET_FILENAME(IC_EL_BOND))
        call ftdqmc_phy0_perfFT(P0, IC_EL_BONDY, GET_FILENAME(IC_EL_BONDY))
        call ftdqmc_phy0_perfFT(P0, IC_EL_BONDC, GET_FILENAME(IC_EL_BONDC))
        call ftdqmc_phy0_perfFT(P0, IC_EL_SWAV, GET_FILENAME(IC_EL_SWAV))
        call ftdqmc_phy0_perfFT(P0, IC_EL_DWAV, GET_FILENAME(IC_EL_DWAV))
        call ftdqmc_phy0_perfFT(P0, IC_EL_DBOND, GET_FILENAME(IC_EL_DBOND))
    endsubroutine ftdqmc_phy0_corFT


    subroutine ftdqmc_phy0_perfFT(P0, IC_index, filename)
		! perform FT for input correlation function
        use mpi

        type(phy0), intent(inout) :: P0
        integer, intent(in) :: IC_index
        character(40), intent(in) :: filename

        ! local variables
        complex(dp) :: mpi_cor_bin(lq)

        ! get data and perform FT
        call mpi_reduce( P0%pnt_cr(IC_index)%pnt, mpi_cor_bin, lq, mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
        ! Fourier transformation
        if( irank .eq. 0 ) then
            mpi_cor_bin = mpi_cor_bin / dcmplx(dble(isize*P0%nobs), 0.d0)
            call ftdqmc_phy0_ft(mpi_cor_bin, filename, IC_index, P0)
        endif
        call mpi_barrier( mpi_comm_world, ierr )
	endsubroutine ftdqmc_phy0_perfFT


    subroutine ftdqmc_phy0_ft(cor, filename, idx, P0)
#ifdef _OPENMP
        USE OMP_LIB
#endif
        use ftdqmc_latt

        type(phy0), intent(in) :: P0
        integer, intent(in) :: idx
        complex(dp), intent(in), dimension(:) :: cor
        character(40), intent(in) :: filename

        ! local variables
        integer :: imj, iq
        real(dp) :: qvec(2)
        complex(dp) :: cork
        complex(dp) :: background


        ! calculate the background because of the Gamma point
        ! background = <O>*<O>*lq, <O> is intensive
        if ( idx .eq. IC_DEN0 ) then ! den0
            background = dcmplx( dble(P0%pnt_sc(IS_FILL)*P0%pnt_sc(IS_FILL))*dble(lq), 0.d0 )
        elseif ( idx .eq. IC_SPZZ ) then
            background = dcmplx( dble(P0%pnt_sc(IS_SpinZ)*P0%pnt_sc(IS_SpinZ))*dble(lq), 0.d0 )
        !elseif ( idx .eq. IC_EL_DEN ) then
        !    background = dcmplx( dble(P0%pnt_sc(IS_EL_DEN)*P0%pnt_sc(IS_EL_DEN))*dble(lq), 0.d0 )
        !elseif ( idx .eq. IC_EL_BOND ) then
        !    background = dcmplx( dble(P0%pnt_sc(IS_EL_BOND)*P0%pnt_sc(IS_EL_BOND))*dble(lq), 0.d0 )
        !elseif ( idx .eq. IC_EL_BONDY ) then
        !    background = dcmplx( dble(P0%pnt_sc(IS_EL_BONDY)*P0%pnt_sc(IS_EL_BONDY))*dble(lq), 0.d0 )
        !elseif ( idx .eq. IC_EL_BONDC ) then
        !    background = dcmplx( dble(P0%pnt_sc(IS_EL_BONDC)*P0%pnt_sc(IS_EL_BONDC))*dble(lq), 0.d0 )
        !elseif ( idx .eq. IC_EL_DBOND ) then
        !    background = dcmplx( dble(P0%pnt_sc(IS_EL_DBOND)*P0%pnt_sc(IS_EL_DBOND))*dble(lq), 0.d0 )
        !elseif ( idx .eq. IC_EL_DWAV ) then
            !background = dcmplx( dble(P0%pnt_sc(IS_EL_DWAV)*P0%pnt_sc(IS_EL_DWAV))*dble(lq), 0.d0 )
            !background = dcmplx( dble(P0%pnt_sc(IS_EL_DWAV)*dconjg(P0%pnt_sc(IS_EL_DWAV)))*dble(lq), 0.d0 )
        !elseif ( idx .eq. IC_EL_SWAV ) then
            !background = dcmplx( dble(P0%pnt_sc(IS_EL_SWAV)*P0%pnt_sc(IS_EL_SWAV))*dble(lq), 0.d0 )
            !background = dcmplx( dble(P0%pnt_sc(IS_EL_SWAV)*dconjg(P0%pnt_sc(IS_EL_SWAV)))*dble(lq), 0.d0 )
        else
            background = czero
        endif

        ! prepare the output file
        open(unit=177,file=filename,status='unknown', action="write", position="append")
        ! loop over k points
        do iq = 1, lq
            qvec = dble( listk(iq,1))*b1_p + dble( listk(iq,2))*b2_p
            cork = czero

            do imj = 1, lq
                cork = cork + cor(imj) /  zexpiqr(imj, iq)
            end do

            ! every kpoint to sq_*.bin file
            if ( qvec(1) .eq. 0.d0 .and. qvec(2) .eq. 0.d0 ) then
                ! remove the background for Gamma point
                write(177, '(2f8.4, 2x, e16.8)') qvec(1), qvec(2), dble(cork*dcmplx(1.d0/dble(lq),0.d0) - background)
            else
                write(177, '(2f8.4, 2x, e16.8)') qvec(1), qvec(2), dble(cork*dcmplx(1.d0/dble(lq),0.d0))
            endif
        enddo
        close(177)
    endsubroutine ftdqmc_phy0_ft


    function generate_filename(raw_name) result(f_name)
        character(len=*), intent(in) :: raw_name
        character(len=64) :: f_name
        character(len=64) :: temp_part
        integer :: i, offset

        ! Logic is the same: extract the part after IC_.
        temp_part = raw_name(4:len_trim(raw_name))

        ! uppercase to lowercase
        offset = ichar('a') - ichar('A')
        do i = 1, len_trim(temp_part)
            if (temp_part(i:i) >= 'A' .and. temp_part(i:i) <= 'Z') then
                temp_part(i:i) = char(ichar(temp_part(i:i)) + offset)
            end if
        end do

        f_name = "sq_" // trim(temp_part) // ".bin"
    end function generate_filename

endmodule ftdqmc_phy0
