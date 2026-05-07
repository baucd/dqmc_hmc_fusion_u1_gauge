module fthmc_core
    use spring
    use ftdqmc_hamilt
    use ftdqmc_latt_class
    use ftdqmc_auxfield_f5_class
    use fthmc_phi_class
    use fthmc_io
    use fthmc_gfun
    use fthmc_phy0
    use fthmc_tdm
    implicit none

contains
    subroutine fthmc_core_initfield(phi, phi_u1, latt)
        ! initialize the momentum field and the R field at the start of MD, draw from Gaussian distribution
        implicit none
        class(fthmc_phi), intent(inout) :: phi(int(Nflavor/2.d0))
        class(ftdqmc_auxfield_f5), intent(inout) :: phi_u1
        class(ftdqmc_latt), intent(inout) :: latt

        ! local variables
        integer :: i, nn, nt, ntau, i1, i2, j1, j2

        ! for the momentum field
        do nt =1, ltrot
            do nn = 1, nfam
                do i = 1, lfam
                    pfield(i, nn, nt) = spring_sfmt_gaussian() / dsqrt(pm)
                enddo
            enddo
        enddo

        ! for the rfield
        do i = 1, int(Nflavor/2.d0)
            call phi(i)%fthmc_phi_rfield_init()
            ! also set the x_vec_old to zero
            phi(i)%x_vec_old = czero
            phi(i)%x_vec = chalf
        enddo

        ! calculate the phifield, M matrix (actually the B matrix) is always set
        call cpu_time(tstart)
        !call fthmc_matrix_calphi_sparse(phi)
        call phi_u1%Md_mv_to_phi_sparse(phi, latt)
        call cpu_time(tend); time_vec(3) = time_vec(3) + (tend-tstart)

        ! obtain the x_vec
        call cpu_time(tstart)
        call fthmc_core_cg(phi, phi_u1, latt)
        call cpu_time(tend); time_vec(1) = time_vec(1) + (tend-tstart)
    endsubroutine fthmc_core_initfield


    subroutine fthmc_core_cg(phi, phi_u1, latt)
#ifdef CUDA_CG
    use fthmc_conjugate_cuda
#else
    use fthmc_conjugate
#endif

        implicit none
        class(fthmc_phi), intent(inout) :: phi(int(Nflavor/2.d0))
        class(ftdqmc_auxfield_f5), intent(inout) :: phi_u1
        class(ftdqmc_latt), intent(inout) :: latt

        ! conjugate gradient to calculate M^T * M for force as well as the Green's function
        ! input: MM output: (M^T*M)^-1 phi = x_vec, x_vec is enough

        ! local
        integer  :: iter, i, j
        real(dp) :: error
        complex(dp) :: x_vec_tmp(ndim*ltrot)

        do i = 1, int(Nflavor/2.d0)
#ifdef CGOPT
            ! setting the x_0 according to this extrapolation formula from Richard's paper
            x_vec_tmp = phi(i)%x_vec
            phi(i)%x_vec = dcmplx(2.d0, 0.d0) * phi(i)%x_vec - phi(i)%x_vec_old
            ! store the old solution to x_vec_old
            phi(i)%x_vec_old = x_vec_tmp
#else
            ! setting the x_0 to be cone
            phi(i)%x_vec = cone
#endif

            ! double precision cg
#ifdef CUDA_CG
            ! CUDA Fortran ver
            !call fthmc_conjugate_cg_cuda(ndim*ltrot, lfam*nfam*ltrot, reshape(zep_rsigl_k, (/lfam*nfam*ltrot/)), reshape(zem_rsigl_k, (/lfam*nfam*ltrot/)), phi(i)%phifield, phi(i)%x_vec, iter, error) ! sparse version
            call fthmc_conjugate_cg_cuda_graph(reshape(zep_rsigl_k, (/lfam*nfam*ltrot/)), reshape(zem_rsigl_k, (/lfam*nfam*ltrot/)), phi(i)%phifield, phi(i)%x_vec, iter, error) ! sparse version
#else
            call fthmc_conjugate_cg(itermax, errate, phi(i)%phifield, phi(i)%x_vec, phi_u1, iter, error, latt) ! sparse version
#endif
            if ( irank .eq. 0) then
                write(fout2, *) iter, error
            endif

        enddo
    endsubroutine fthmc_core_cg


    subroutine fthmc_core_cg_gfun(phi, phi_u1, latt)
#ifdef CUDA_CG
    use fthmc_conjugate_cuda
#else
    use fthmc_conjugate
#endif
        implicit none
        class(fthmc_phi), intent(inout) :: phi
        class(ftdqmc_auxfield_f5), intent(inout) :: phi_u1
        class(ftdqmc_latt), intent(inout) :: latt

        ! conjugate gradient to calculate M^T * M for force as well as the Green's function
        ! input: MM output: (M^T*M)^-1 phi = x_vec, x_vec is enough

        ! local
        integer  :: iter, iter_sp, i, j
        real(dp) :: error
        complex :: x_vec_sp(ndim*ltrot)

        phi%x_vec = cone
        ! double precision cg
#ifdef CUDA_CG
        ! CUDA Fortran ver
        !call fthmc_conjugate_cg_cuda(ndim*ltrot, lfam*nfam*ltrot, reshape(zep_rsigl_k, (/lfam*nfam*ltrot/)), reshape(zem_rsigl_k, (/lfam*nfam*ltrot/)), phi%phifield, phi%x_vec, iter, error) ! sparse version
        call fthmc_conjugate_cg_cuda_graph(reshape(zep_rsigl_k, (/lfam*nfam*ltrot/)), reshape(zem_rsigl_k, (/lfam*nfam*ltrot/)), phi%phifield, phi%x_vec, iter, error) ! sparse version
#else
        call fthmc_conjugate_cg(itermax, errate, phi%phifield, phi%x_vec, phi_u1, iter, error, latt) ! sparse version
#endif

        if ( irank .eq. 0) then
            write(fout2, *) iter, error
        endif
    endsubroutine fthmc_core_cg_gfun


    subroutine fthmc_core_cal_ener(energy, phi, phi_u1, latt)
        implicit none
        ! calculate energy of the effective Hamiltonian
        real(dp), intent(out) :: energy
        class(fthmc_phi), intent(inout) :: phi(int(Nflavor/2.d0))
        class(ftdqmc_auxfield_f5), intent(inout) :: phi_u1
        class(ftdqmc_latt), intent(inout) :: latt

        ! local variables
        integer :: nt, nf, i, i1, i2, nn, ntau, j1, j2
        integer :: i_plaqA, b_plaqA, i_plaqB, b_plaqB, ib, inf, ilf, ijs
        real(dp) :: dphi
        real(dp) :: phi_Aplaq, phi_Bplaq, ejpi, ejs, ener_fer
        real(dp), pointer :: xfield(:,:,:)

        ! external functions
        integer, external :: npbc
        real(dp), external :: ddot
        complex(dp), external :: zdotc

        ! use xfield to associate phi_u1%gauge_u1
        xfield => phi_u1%gauge_u1

        energy = 0.d0
        ejpi = 0.d0; ejs = 0.d0
        do nt = 1, ltrot
            do nf = 1, nfam
                do i = 1, lfam
                    i1 = latt%l_bonds(1,i,nf)
                    i2 = latt%l_bonds(2,i,nf)

                    ! Jpi term
                    IF ( abs(jpi) .ne. 0.d0 ) THEN
                        ! A plaq
                        phi_Aplaq = 0.d0
                        i_plaqA = latt%inv_Aplaq_bondcord(1,i,nf)
                        b_plaqA = latt%inv_Aplaq_bondcord(2,i,nf)
                        do ib = 1, 4
                            inf = latt%plaq_bondcord(1,ib,i_plaqA)
                            ilf = latt%plaq_bondcord(2,ib,i_plaqA)
                            phi_Aplaq = phi_Aplaq + xfield(ilf,inf,nt)*latt%sgnA_plaq(ib)
                        end do

                        ! B plaq
                        phi_Bplaq = 0.d0
                        i_plaqB = latt%inv_Bplaq_bondcord(1,i,nf)
                        b_plaqB = latt%inv_Bplaq_bondcord(2,i,nf)
                        do ib = 1, 4
                            inf = latt%plaq_bondcord(1,ib,i_plaqB)
                            ilf = latt%plaq_bondcord(2,ib,i_plaqB)
                            phi_Bplaq = phi_Bplaq + xfield(ilf,inf,nt)*latt%sgnB_plaq(ib)
                        end do
#IFDEF COMPACT
                        ejpi = ejpi + cos(phi_Aplaq) + cos(phi_Bplaq)
#ELSE
                        ejpi = ejpi + ((phi_Aplaq-pi)**2 + (phi_Bplaq+pi)**2)/2.d0
#ENDIF
                    END IF

                    ! Js term
                    dphi = xfield(i,nf,npbc(nt+1,ltrot)) - xfield(i,nf,nt)
#IFDEF COMPACT
                    ejs = ejs + 2.d0*( 1.d0 - cos(dphi) )
#ELSE
                    ejs = ejs + dphi * dphi
#ENDIF
                end do
            end do
        end do
        ejs = ejs / ( js*dtau*dtau ) * dtau
        ejpi = ejpi * 0.25d0 * jpi * dtau  ! 0.25 is because of 4 times conuting
        energy = ejs + ejpi

        ! save energy
        ener_js = ejs
        ener_jpi = ejpi

        ! set potential part of energy, exclude the momentum part
        ener_pot = energy

        ! energy from the determinants
        ener_fer = 0.d0
        do i = 1, int(Nflavor/2.d0)
            ener_fer = ener_fer + dble(zdotc(ndim*ltrot, phi(i)%phifield, 1, phi(i)%x_vec,1))
        enddo
        ener_ek = ener_fer

        ! add the momentum field  p^2 and the contribution from fermion determinant
        energy = energy + ener_fer

        ! set ener_q and ener_p
        ener_q = energy
        ener_p = 0.5d0*pm*ddot(2*ndim*ltrot, reshape(pfield, (/lfam*nfam*ltrot/)), 1, reshape(pfield,(/lfam*nfam*ltrot/)),1)

        energy = energy + ener_p

        ! nullify
        NULLIFY(xfield)
    endsubroutine fthmc_core_cal_ener


    subroutine fthmc_core_cal_gfun_meas(gfun0, lmeasure_equaltime, lmeasure_dyn, P0, T0, phi_u1, latt)
        ! Input: rfield, M matrix to compute the new phifield, MM and the new phifield to compute the new x_vec
        ! so rfield, M matrix and MM matrix are needed

        ! calculate the inverse of M matrix which is the Green's function
        use fthmc_phy0
        use fthmc_tdm
        implicit none

        type(gfun), intent(inout) :: gfun0
        logical, intent(in) :: lmeasure_equaltime, lmeasure_dyn
        type(phy0), intent(inout) :: P0
        type(tdm), intent(inout) :: T0
        class(ftdqmc_auxfield_f5), intent(inout) :: phi_u1
        class(ftdqmc_latt), intent(inout) :: latt

#IFDEF DIRINV
        ! DIRECT INVERSION
        !!! benchmark with direct inversion of matrix
        integer :: nt, i, j, i1, i2, i3, i4, j1, j2, j3, j4, nt1, nt2
        integer, external :: npbc
        complex(dp), dimension(:,:), allocatable :: Mmat
#IFDEF DETSGN
        ! test the determinant sign
        complex(dp) :: det_z
#ENDIF
#ELSE
        ! local
        type(fthmc_phi) :: phi
        integer  :: nt, i, j, i1, i2, i3, i4, j1, j2, j3, j4, nt1, nt2, ncount
        integer, external :: npbc
#ENDIF


        !!! allocate memory for array
#IFDEF DIRINV
        ! Mmat
        allocate(Mmat(ndim*ltrot, ndim*ltrot))
#ELSE
        ! allocate space for object phi
        call phi%fthmc_phi_alloc
#ENDIF

        !!! EXECUTABLE
#IFDEF DIRINV
        ! DIRECT INVERSION
        ! set up M matrix
        Mmat = czero
        ! the diagonal blocks
        do i =1, ndim*ltrot
            Mmat(i,i) = cone
        enddo
        ! set up nt from 1 to ltrot - 1
        do nt = 1, ltrot-1
            i1 = (nt-0)*ndim+1; i2 = i1 + ndim - 1
            j1 = (nt-1)*ndim+1; j2 = j1 + ndim - 1
            Mmat(i1:i2, j1:j2) = -Bmat_up(:,:,nt)
        enddo
        ! just ltrot
        i1 = 1; i2 = ndim
        j1 = (ltrot-1)*ndim + 1; j2 = j1+ndim-1
        Mmat(i1:i2,j1:j2) = Bmat_up(:,:,ltrot)

        M_inv = Mmat
        ! s_inv_z will destory M_inv and replace it with the result
        call s_inv_z(ndim*ltrot, M_inv)
        ! END DIRECT INVERSION

#IFDEF DETSGN
        ! test: the determinant sign
        call s_logdet_z(ndim*ltrot, Mmat, det_z)
        write(*, '(2e16.8)') dble(det_z), dcos(aimag(det_z))
#ENDIF

        ! set up the equaltime Green's function
        do nt = 1, ltrot
            do i = 1, ndim
                do j = 1, ndim
                    i1 = i + (nt-1) * ndim
                    j1 = j + (nt-1) * ndim
                    gfun0%grup(i,j,nt) = M_inv(i1, j1)
                enddo
            enddo
        enddo
        ! equaltime measurements
        if (lmeasure_equaltime) call fthmc_phy0_meas(gfun0, P0, latt)


        if ( lmeasure_dyn ) then
            ! set up the equaltime and dynamical Green's function
            do nt = 1, ltrot ! ltrot is the beta point
                do nt1 = 1, ltrot
                    nt2 = npbc(nt1+nt, ltrot)
                    do i = 1, ndim
                        do j = 1, ndim
                            i1 = i + (nt1-1) * ndim
                            j1 = j + (nt1-1) * ndim

                            i2 = i + (nt1-1) * ndim
                            j2 = j + (nt2-1) * ndim

                            i3 = i + (nt2-1) * ndim
                            j3 = j + (nt1-1) * ndim

                            i4 = i + (nt2-1) * ndim
                            j4 = j + (nt2-1) * ndim
                            gfun0%g00up(i,j,nt) = M_inv(i1, j1)
                            gfun0%g0tup(i,j,nt) = M_inv(i2, j2)
                            gfun0%gt0up(i,j,nt) = M_inv(i3, j3)
                            gfun0%gttup(i,j,nt) = M_inv(i4, j4)
                        enddo
                    enddo
                    ! dealing with nt = ltrot, equaltime green's function
                    if (nt .eq. ltrot) gfun0%g0tup(:,:,ltrot) = gfun0%g0tup(:,:,ltrot) - Imat(:,:)

                    ! dynamical measurements
                    call fthmc_tdm_meas(nt, gfun0, T0, latt)
                enddo
            enddo

        endif

#ELSE

        ! STOCHASTIC ESTIMATOR
        gfun0%grup(:,:,:) = czero
        if ( lmeasure_dyn ) then
            gfun0%g00up(:,:,:)= czero
            gfun0%g0tup(:,:,:)= czero
            gfun0%gt0up(:,:,:)= czero
            gfun0%gttup(:,:,:)= czero
        endif
        do ncount = 1, nsamples

            ! assign a new rfield
            call phi%fthmc_phi_rfield_init()

            ! calculate the phifield
            call cpu_time(tstart)
            !call fthmc_matrix_calphi_gfun_sparse(phi)
            call phi_u1%Md_mv_to_phi_gfun_sparse(phi, latt)
            call cpu_time(tend); time_vec(3) = time_vec(3) + (tend-tstart)

            ! do another conjugate residual method to compute the new x_vec with the new phifield
            call cpu_time(tstart)
            call fthmc_core_cg_gfun(phi, phi_u1, latt)
            call cpu_time(tend); time_vec(5) = time_vec(5) + (tend-tstart)

            ! only equaltime measurements
            ! prepare Green's function
            do nt = 1, ltrot
                do i = 1, ndim
                    do j = 1, ndim
                        i1 = i + (nt-1) * ndim
                        j1 = j + (nt-1) * ndim
                        gfun0%grup(i,j,nt) = gfun0%grup(i,j,nt) + phi%x_vec(i1) * dconjg(phi%rfield(j1)) / dble(nsamples)
                    enddo
                enddo
            enddo

            ! dynamical measurements
            if ( lmeasure_dyn ) then
                do nt = 1, ltrot ! ltrot is the beta point

                    do nt1 = 1, ltrot
                       nt2 = npbc(nt1+nt, ltrot)
                        do i = 1, ndim
                            do j = 1, ndim
                                i1 = i + (nt1-1) * ndim
                                j1 = j + (nt1-1) * ndim

                                i2 = i + (nt1-1) * ndim
                                j2 = j + (nt2-1) * ndim

                                i3 = i + (nt2-1) * ndim
                                j3 = j + (nt1-1) * ndim

                                i4 = i + (nt2-1) * ndim
                                j4 = j + (nt2-1) * ndim
                                gfun0%g00up(i,j,nt) = gfun0%g00up(i,j,nt) + phi%x_vec(i1) * dconjg(phi%rfield(j1)) / dble(nsamples)
                                gfun0%g0tup(i,j,nt) = gfun0%g0tup(i,j,nt) + phi%x_vec(i2) * dconjg(phi%rfield(j2)) / dble(nsamples)
                                gfun0%gt0up(i,j,nt) = gfun0%gt0up(i,j,nt) + phi%x_vec(i3) * dconjg(phi%rfield(j3)) / dble(nsamples)
                                gfun0%gttup(i,j,nt) = gfun0%gttup(i,j,nt) + phi%x_vec(i4) * dconjg(phi%rfield(j4)) / dble(nsamples)
                            enddo
                        enddo
                    enddo

                enddo
            endif

        enddo

        ! perform the measurements
        call cpu_time(tstart)

        ! equaltime
        if (lmeasure_equaltime) call fthmc_phy0_meas(gfun0, P0, phi_u1, latt)
        ! dynamical
        if ( lmeasure_dyn ) then
            ! dealing with nt = ltrot, equaltime green's function
            gfun0%g0tup(:,:,ltrot) = gfun0%g0tup(:,:,ltrot) - Imat(:,:)
            call fthmc_tdm_meas(nt, gfun0, T0, latt)
        endif

        call cpu_time(tend); time_vec(6) = time_vec(6) + (tend-tstart)

        ! END STOCHASTIC ESTIMATOR
#ENDIF

        ! free up space
#IFDEF DIRINV
        deallocate(Mmat)
#ELSE
        call phi%fthmc_phi_free
#ENDIF
    endsubroutine fthmc_core_cal_gfun_meas


    subroutine fthmc_core_cal_gfun_meas_singlesource(gfun0, lmeasure_equaltime, lmeasure_dyn, P0, T0, phi_u1, latt)
        ! calculate the inverse of M matrix which is the Green's function
        !!! with only on source
        use fthmc_phy0
        use fthmc_tdm

        implicit none
        type(gfun), intent(inout) :: gfun0
        logical, intent(in) :: lmeasure_equaltime, lmeasure_dyn
        type(phy0), intent(inout) :: P0
        type(tdm), intent(inout) :: T0
        class(ftdqmc_auxfield_f5), intent(inout) :: phi_u1
        class(ftdqmc_latt), intent(inout) :: latt

        ! local variables
        type(fthmc_phi) :: phi
        integer  :: nt, i, j, i1, nt1, nt2, imj
        integer, external :: npbc
        complex(dp) :: tmp_vec(ndim)

        ! initialize tmp_vec
        tmp_vec(:)=czero
        tmp_vec(1)=cone

        ! get the first column of B^{\dagger}_{ltrot}
        call phi_u1%Bmat_mv(tmp_vec, ltrot, latt)

        ! initialize the phi vector
        call phi%fthmc_phi_alloc()
        phi%phifield(:) = czero
        phi%phifield(1) = cone
        phi%phifield( ((ltrot-1)*ndim+1) : (ltrot*ndim) ) = tmp_vec(:)

        ! do another conjugate residual method to compute the new x_vec with the new phifield
        call fthmc_core_cg_gfun(phi, phi_u1, latt)

        ! only equaltime measurements
        do nt = 1, ltrot
            ! set the one source correlation function
            do i = 1, ndim
                j=1 ! vector, many rows and one column
                imj = latt%latt_imj(j,i)
                gfun0%grup_onesource(imj, nt) = phi%x_vec(i)
            enddo

            ! set up the whole equaltime Green's function making use of translational invariance
            do i = 1, ndim
                do j = 1, ndim
                    imj = latt%latt_imj(j,i)
                    gfun0%grup(i,j,nt) = gfun0%grup_onesource(imj,nt)
                enddo
            enddo
        enddo
        ! equaltime measurements
        if (lmeasure_equaltime) call fthmc_phy0_meas(gfun0, P0, phi_u1, latt)


        ! dynamical measurements
        if ( lmeasure_dyn ) then
            do nt = 1, ltrot ! ltrot is the beta point

                ! get the dynamical Green's function
                nt1 = 1
                nt2 = npbc(nt1+nt, ltrot)
                do i = 1, ndim
                    j = 1 ! nt=1, site index=1
                    imj = latt%latt_imj(j,i)
                    i1 = i + (nt2-1) * ndim

                    gfun0%g00up_onesource(imj,nt) = phi%x_vec(i)
                    gfun0%g0tup_onesource(imj,nt) = phi%x_vec(i1)
                    gfun0%gt0up_onesource(imj,nt) = phi%x_vec(i1)
                    gfun0%gttup_onesource(imj,nt) = phi%x_vec(i)
                enddo

                ! from one source to construct the whole dyn Green's function
                do i = 1,ndim
                    do j = 1, ndim
                        imj = latt%latt_imj(j,i)
                        gfun0%g00up(i,j,nt) = gfun0%g00up_onesource(imj,nt)
                        gfun0%g0tup(i,j,nt) = gfun0%g0tup_onesource(imj,nt)
                        gfun0%gt0up(i,j,nt) = gfun0%gt0up_onesource(imj,nt)
                        gfun0%gttup(i,j,nt) = gfun0%gttup_onesource(imj,nt)
                    enddo
                enddo
                    ! dealing with nt = ltrot, equaltime green's function
                if (nt .eq. ltrot) gfun0%g0tup(:,:,ltrot) = gfun0%g0tup(:,:,ltrot) - Imat(:,:)

                ! dynamical measurements
                call fthmc_tdm_meas(nt, gfun0, T0, latt)
            enddo
        endif

        ! free up space
        call phi%fthmc_phi_free()
    endsubroutine fthmc_core_cal_gfun_meas_singlesource


    subroutine fthmc_core_cal_force(phi, phi_u1, latt)
#IFDEF _OPENMP
        use OMP_LIB
#ENDIF
        implicit none
        class(fthmc_phi), intent(inout) :: phi(int(Nflavor/2.d0))
        class(ftdqmc_auxfield_f5), intent(inout) :: phi_u1
        class(ftdqmc_latt), intent(inout) :: latt

        ! calculate the force to update the momentum field, ps: no dt multiply here
        ! input: xfield, rfield, phifield, M matrix  output: force matrix
        ! derivative comes from pure gauge action and the pseudo-fermion action

        integer, external :: npbc
        complex(dp), external :: zdotu, zdotc

        ! local
        integer :: i, ib, j, nt, nf, i_plaq, b_plaq, i4, i5, inf, ilf, icount
        integer :: s1, e1, s2, e2
        real(dp) :: force_bos, force_fer, force_newcode, phi_Aplaq_old, phi_Bplaq_old
        complex(dp) :: vec_tau1(ndim), vec_tau2(ndim)
        real(dp), pointer :: xfield(:,:,:)

        ! important
        integer :: nt_in, nf_in, lf_in

        xfield => phi_u1%gauge_u1

        ! in order to reach O(N) scaling, initialize Btau_vtau vectors
        call cpu_time(tstart)
        call phi_u1%ftdqmc_auxfield_force_p1_preBvec(phi, latt) ! now all e^vi's * phi%(i)vtau are set, total 16 position
        call cpu_time(tend); time_vec(2) = time_vec(2) + (tend-tstart)

        call cpu_time(tstart)
        hybrid_force = 0.d0
        hybrid_force1 = 0.d0
        hybrid_force2 = 0.d0
        ! the order: nfam =4, lfam, and for each time slices with the support of openmp acceleration
        do nf = nfam, 1, -1
            do i = 1, lfam
                do nt = 1, ltrot
                    nf_in = nf
                    lf_in = i
                    nt_in = nt

                    ! force from the pseudofermion part, new code scale with O(N) in total, 1*2 * 2*1 inner product
                    force_fer = 0.d0
                    !call fthmc_matrix_force_Btau_vtau(lf_in, nf_in, nt_in, force_fer) ! return the force
                    call phi_u1%ftdqmc_auxfield_force_p2_Btau_vtau(lf_in, nf_in, nt_in, force_fer, latt) ! return the force

                    ! force from the gauge action part
#ifdef COMPACT
                    force_bos = ( dsin(xfield(i,nf,nt)-xfield(i,nf,npbc(nt-1,ltrot))) + &
                                  dsin(xfield(i,nf,nt)-xfield(i,nf,npbc(nt+1,ltrot))) ) / (0.5d0*js*dtau)
#else
                    force_bos = ( xfield(i,nf,nt)-xfield(i,nf,npbc(nt-1,ltrot)) +&
                                  xfield(i,nf,nt)-xfield(i,nf,npbc(nt+1,ltrot)) ) / (0.5d0*js*dtau)
#endif

                    ! boson_ratio for pi-flux pinning field
                    if ( abs(jpi) .ne. 0.d0 ) then
                        ! A plaq
                        phi_Aplaq_old = 0.d0
                        i_plaq = latt%inv_Aplaq_bondcord(1,i,nf)
                        b_plaq = latt%inv_Aplaq_bondcord(2,i,nf)
                        do ib = 1, 4
                            if( ib .eq. b_plaq ) then
                                phi_Aplaq_old = phi_Aplaq_old + xfield(i,nf,nt)*latt%sgnA_plaq(ib)
                            else
                                inf = latt%plaq_bondcord(1,ib,i_plaq)
                                ilf = latt%plaq_bondcord(2,ib,i_plaq)
                                phi_Aplaq_old = phi_Aplaq_old + xfield(ilf,inf,nt)*latt%sgnA_plaq(ib)
                            end if
                        end do
                        ! Jpi part force, A plaq
#ifdef COMPACT
                        force_bos = force_bos - dtau*jpi*latt%sgnA_plaq(b_plaq)*dsin(phi_Aplaq_old)
#else
                        force_bos = force_bos + dtau*jpi*latt%sgnA_plaq(b_plaq)*(phi_Aplaq_old-pi)
#endif

                        ! B plaq
                        phi_Bplaq_old = 0.d0
                        i_plaq = latt%inv_Bplaq_bondcord(1,i,nf)
                        b_plaq = latt%inv_Bplaq_bondcord(2,i,nf)

                        do ib = 1, 4
                            if( ib .eq. b_plaq ) then
                                phi_Bplaq_old = phi_Bplaq_old + xfield(i,nf,nt)*latt%sgnB_plaq(ib)
                            else
                                inf = latt%plaq_bondcord(1,ib,i_plaq)
                                ilf = latt%plaq_bondcord(2,ib,i_plaq)
                                phi_Bplaq_old = phi_Bplaq_old + xfield(ilf,inf,nt)*latt%sgnB_plaq(ib)
                            end if
                        end do
#ifdef COMPACT
                        force_bos = force_bos - dtau*jpi*latt%sgnB_plaq(b_plaq)*dsin(phi_Bplaq_old)
#else
                        force_bos = force_bos + dtau*jpi*latt%sgnB_plaq(b_plaq)*(phi_Bplaq_old+pi)
#endif
                    endif

                    ! combine the force from gauge boson and fermion part together
                    hybrid_force(i, nf, nt) = force_fer + force_bos
                    hybrid_force1(i, nf, nt) = force_fer
                    hybrid_force2(i, nf, nt) = force_bos

                    ! for test: output force
                    !write(111, '(2e16.8)') force_fer, force_bos

                enddo
            enddo
        enddo
        call cpu_time(tend); time_vec(7) = time_vec(7) + (tend-tstart)

        ! output the force log
        !if ( irank .eq. 0) then
        !    write(fout4, '(3f10.6)') sum(abs(hybrid_force))/(dble(nfam*lfam*ltrot)), minval(abs(reshape(hybrid_force,(/lfam*nfam*ltrot/)))), maxval(abs(reshape(hybrid_force,(/lfam*nfam*ltrot/))))
        !endif

        nullify(xfield)
    endsubroutine fthmc_core_cal_force


    subroutine fthmc_core_cal_force_fer(phi, phi_u1, latt)
#ifdef _OPENMP
        use OMP_LIB
#endif
        implicit none
        class(fthmc_phi), intent(inout) :: phi(int(Nflavor/2.d0))
        class(ftdqmc_auxfield_f5), intent(inout) :: phi_u1
        class(ftdqmc_latt), intent(inout) :: latt

        ! calculate the force to update the momentum field, ps: no dt multiply here
        ! input: xfield, rfield, phifield, M matrix  output: force matrix
        ! derivative comes from pure gauge action and the pseudo-fermion action

        integer, external :: npbc
        complex(dp), external :: zdotu, zdotc

        ! local
        integer :: i, ib, j, nt, nf, i_plaq, b_plaq, i4, i5, inf, ilf, icount
        integer :: s1, e1, s2, e2
        real(dp) :: force_bos, force_fer, force_newcode, phi_Aplaq_old, phi_Bplaq_old
        complex(dp) :: vec_tau1(ndim), vec_tau2(ndim)

        ! important
        integer :: nt_in, nf_in, lf_in

        ! in order to reach O(N) scaling, initialize Btau_vtau vectors
        call cpu_time(tstart)
        !call fthmc_matrix_preBvec(phi) ! now all e^vi's * phi%(i)vtau are set, total 16 position
        call phi_u1%ftdqmc_auxfield_force_p1_preBvec(phi, latt) ! now all e^vi's * phi%(i)vtau are set, total 16 position
        call cpu_time(tend); time_vec(2) = time_vec(2) + (tend-tstart)

        call cpu_time(tstart)
        ! the order: nfam =4, lfam, and for each time slices with the support of openmp acceleration
        hybrid_force1 = 0.d0
        do nf = nfam, 1, -1
            do i = 1, lfam
                do nt = 1, ltrot
                    nf_in = nf
                    lf_in = i
                    nt_in = nt

                    ! new code scale with O(N) in total, 1*2 * 2*1 inner product
                    force_fer = 0.d0
                    !call fthmc_matrix_force_Btau_vtau(lf_in, nf_in, nt_in, force_fer) ! return the force
                    call phi_u1%ftdqmc_auxfield_force_p2_Btau_vtau(lf_in, nf_in, nt_in, force_fer, latt) ! return the force

                    ! force from the fermion part
                    hybrid_force1(i, nf, nt) = force_fer
                enddo
            enddo
        enddo
        call cpu_time(tend); time_vec(7) = time_vec(7) + (tend-tstart)
    endsubroutine fthmc_core_cal_force_fer


    subroutine fthmc_core_cal_force_bos(phi, phi_u1, latt)
#ifdef _OPENMP
        use OMP_LIB
#endif
        implicit none
        class(fthmc_phi), intent(inout) :: phi(int(Nflavor/2.d0))
        class(ftdqmc_auxfield_f5), intent(inout) :: phi_u1
        class(ftdqmc_latt), intent(inout) :: latt

        ! calculate the force to update the momentum field, ps: no dt multiply here
        ! input: xfield, rfield, phifield, M matrix  output: force matrix
        ! derivative comes from pure gauge action and the pseudo-fermion action

        integer, external :: npbc
        complex(dp), external :: zdotu, zdotc

        ! local
        integer :: i, ib, j, nt, nf, i_plaq, b_plaq, i4, i5, inf, ilf, icount
        integer :: s1, e1, s2, e2
        real(dp) :: force_bos, force_fer, force_newcode, phi_Aplaq_old, phi_Bplaq_old
        complex(dp) :: vec_tau1(ndim), vec_tau2(ndim)
        real(dp), POINTER :: xfield(:,:,:)

        ! important
        integer :: nt_in, nf_in, lf_in

        xfield => phi_u1%gauge_u1

        hybrid_force2 = 0.d0
        ! the order: nfam =4, lfam, and for each time slices with the support of openmp acceleration
        do nf = nfam, 1, -1
            do i = 1, lfam
                do nt = 1, ltrot
                    nf_in = nf
                    lf_in = i
                    nt_in = nt

                    !!! force from the gauge action part
#ifdef COMPACT
                    force_bos = ( dsin(xfield(i,nf,nt)-xfield(i,nf,npbc(nt-1,ltrot))) + &
                                  dsin(xfield(i,nf,nt)-xfield(i,nf,npbc(nt+1,ltrot))) ) / (0.5d0*js*dtau)
#else
                    force_bos = ( xfield(i,nf,nt)-xfield(i,nf,npbc(nt-1,ltrot)) +&
                                  xfield(i,nf,nt)-xfield(i,nf,npbc(nt+1,ltrot)) ) / (0.5d0*js*dtau)
#endif

                    ! boson_ratio for pi-flux pinning field
                    if ( abs(jpi) .ne. 0.d0 ) then
                        ! A plaq
                        phi_Aplaq_old = 0.d0
                        i_plaq = latt%inv_Aplaq_bondcord(1,i,nf)
                        b_plaq = latt%inv_Aplaq_bondcord(2,i,nf)
                        do ib = 1, 4
                            if( ib .eq. b_plaq ) then
                                phi_Aplaq_old = phi_Aplaq_old + xfield(i,nf,nt)*latt%sgnA_plaq(ib)
                            else
                                inf = latt%plaq_bondcord(1,ib,i_plaq)
                                ilf = latt%plaq_bondcord(2,ib,i_plaq)
                                phi_Aplaq_old = phi_Aplaq_old + xfield(ilf,inf,nt)*latt%sgnA_plaq(ib)
                            end if
                        end do
                        ! Jpi part force, A plaq
#ifdef COMPACT
                        force_bos = force_bos - dtau*jpi*latt%sgnA_plaq(b_plaq)*dsin(phi_Aplaq_old)
#else
                        force_bos = force_bos + dtau*jpi*latt%sgnA_plaq(b_plaq)*(phi_Aplaq_old-pi)
#endif

                        ! B plaq
                        phi_Bplaq_old = 0.d0
                        i_plaq = latt%inv_Bplaq_bondcord(1,i,nf)
                        b_plaq = latt%inv_Bplaq_bondcord(2,i,nf)

                        do ib = 1, 4
                            if( ib .eq. b_plaq ) then
                                phi_Bplaq_old = phi_Bplaq_old + xfield(i,nf,nt)*latt%sgnB_plaq(ib)
                            else
                                inf = latt%plaq_bondcord(1,ib,i_plaq)
                                ilf = latt%plaq_bondcord(2,ib,i_plaq)
                                phi_Bplaq_old = phi_Bplaq_old + xfield(ilf,inf,nt)*latt%sgnB_plaq(ib)
                            end if
                        end do
#ifdef COMPACT
                        force_bos = force_bos - dtau*jpi*latt%sgnB_plaq(b_plaq)*dsin(phi_Bplaq_old)
#else
                        force_bos = force_bos + dtau*jpi*latt%sgnB_plaq(b_plaq)*(phi_Bplaq_old+pi)
#endif
                    endif

                    ! force from the gauge boson part
                    hybrid_force2(i, nf, nt) = force_bos
                enddo
            enddo
        enddo

        NULLIFY(xfield)
    endsubroutine fthmc_core_cal_force_bos


    subroutine fthmc_core_md(phi, phi_u1, latt)
        implicit none
        class(fthmc_phi), intent(inout) :: phi(int(Nflavor/2.d0))
        class(ftdqmc_auxfield_f5), intent(inout) :: phi_u1
        class(ftdqmc_latt), intent(inout) :: latt

        ! local variables
        integer :: nf, i, ncount
        real(dp), pointer :: xfield(:,:,:)
        ! Fourier acceleration
        complex(dp), allocatable :: z_tmp(:), z_tmp_x(:), z_tmp_p(:)

        xfield => phi_u1%gauge_u1

        ! move the momentum field to t + dt/2, prepare the first force
        call fthmc_core_cal_force(phi, phi_u1, latt)
        ! fourier acceleration, fft of force matrix
        if ( lfourier ) then
            do nf = 1, nfam
                do i = 1, lfam
                    z_tmp(:) = dcmplx(hybrid_force(i, nf, :), 0.d0)
                    call onedimension_fft(ltrot, z_tmp)
                    hybrid_force_k_z(i, nf, :) = z_tmp(:)
                enddo
            enddo
            ! allocate tmp arrays
            allocate(z_tmp(ltrot))
            allocate(z_tmp_x(ltrot))
            allocate(z_tmp_p(ltrot))
        endif

        do ncount = 1, mdstep
            if ( lfourier ) then
                !!! molecular dynamics with Fourier acceleration
                ! update gauge field
                do nf = 1, nfam
                    do i = 1, lfam
                        !! pfield
                        z_tmp_p(:) = dcmplx(pfield(i, nf, :), 0.d0)
                        call onedimension_fft(ltrot, z_tmp_p)

                        !! update the xfield in momentum space
                        z_tmp_x(:) = dcmplx(dt*pm , 0.d0) * a_k(:) * z_tmp_p(:) - dcmplx(pm*dt**2/2.d0, 0.d0) * a_k(:) * a_k(:) * hybrid_force_k_z(i, nf, :)

                        !! inversed Fourier transformation
                        call onedimension_invfft(ltrot, z_tmp_x)
                        xfield(i, nf, :) = xfield(i, nf, :) + 1.d0/dble(ltrot) * dble(z_tmp_x(:))
                    enddo
                enddo

                ! calculation the force at (n+1)
                hybrid_force_n_k_z(:,:,:) = hybrid_force_k_z(:,:,:) ! preserve the F(n)

                ! update the B matrix
                !call fthmc_core_x2B
                call phi_u1%vi_to_expvi(latt)
                ! calculate the inverse of MM matrix: (M^\dagger *M)^-1 phi = x_vec
                call fthmc_core_cg(phi, phi_u1, latt)
                ! calculate the force with the updated xfield and M matrix
                call fthmc_core_cal_force(phi, phi_u1, latt)

                ! update momentum field
                do nf = 1, nfam
                    do i = 1, lfam
                        ! force vector
                        z_tmp(:) = dcmplx(hybrid_force(i, nf, :), 0.d0)
                        call onedimension_fft(ltrot, z_tmp)
                        hybrid_force_k_z(i, nf, :) = z_tmp(:)

                        !! pfield
                        z_tmp_p(:) = dcmplx(pfield(i, nf, :), 0.d0)
                        call onedimension_fft(ltrot, z_tmp_p)

                        !! update the conjugate momentum field in momentum space
                        z_tmp_p(:) = (-1.d0) * dcmplx(dt/2.d0, 0.d0) * a_k(:) * ( hybrid_force_n_k_z(i,nf,:) + hybrid_force_k_z(i,nf,:) )

                        !! inversed Fourier transformation
                        call onedimension_invfft(ltrot, z_tmp_p)
                        pfield(i, nf, :) = pfield(i, nf, :) + 1.d0/dble(ltrot) * dble(z_tmp_p(:))
                    enddo
                enddo

            else
                ! md, leapfrog version 2
                call cpu_time(tstart)
                xfield(:,:,:) = xfield(:,:,:) + dt*pm*pfield(:,:,:) - pm*(dt**2/2.d0) * hybrid_force(:,:,:)
                call cpu_time(tend); time_vec(4) = time_vec(4) + (tend-tstart)

                ! update the B matrix
                !call fthmc_core_x2B
                call phi_u1%vi_to_expvi(latt)

                ! calculate the inverse of MM matrix: (M^\dagger *M)^-1 phi = x_vec
                call cpu_time(tstart)
                call fthmc_core_cg(phi, phi_u1, latt)
                call cpu_time(tend); time_vec(1) = time_vec(1) + (tend-tstart)

                ! store n force
                hybrid_force_n(:,:,:) = hybrid_force(:,:,:)
                ! calculate the force with the updated xfield and M matrix
                call fthmc_core_cal_force(phi, phi_u1, latt)

                call cpu_time(tstart)
                pfield(:,:,:) = pfield(:,:,:) - dt/2.d0 * (hybrid_force(:,:,:) + hybrid_force_n(:,:,:))
                call cpu_time(tend); time_vec(4) = time_vec(4) + (tend-tstart)
            endif
        enddo

        if ( lfourier ) then
            deallocate(z_tmp, z_tmp_x, z_tmp_p)
        endif

        NULLIFY(xfield)
    endsubroutine fthmc_core_md


    subroutine fthmc_core_md_splitting(phi, phi_u1, latt)
        implicit none
        class(fthmc_phi), intent(inout) :: phi(int(Nflavor/2.d0))
        class(ftdqmc_auxfield_f5), intent(inout) :: phi_u1
        class(ftdqmc_latt), intent(inout) :: latt

        ! local variables
        integer :: nf, i, ncount
        real(dp), POINTER :: xfield(:,:,:)
        ! for test
        integer :: j, lf, nt

        xfield => phi_u1%gauge_u1

        ! move the momentum field to t + dt/2, prepare the first force
        call fthmc_core_cal_force(phi, phi_u1, latt)
        ! for test: output gauge field after inner dynamics
        !do nt =1, ltrot
        !do nf =1, nfam
        !do lf =1, lfam
        !    write(*,'(2e16.8)') hybrid_force1(lf, nf, nt), hybrid_force2(lf, nf, nt)
        !enddo
        !enddo
        !enddo
        !stop

        do i = 1, mdstep
            ! firstly update the momentum using force fermion
            call cpu_time(tstart)
            pfield(:,:,:) = pfield(:,:,:) - dt/2.d0*hybrid_force1(:,:,:)
            call cpu_time(tend); time_vec(4) = time_vec(4) + (tend-tstart)

            ! M inner leapfrog dynamics with effective stepsize \epsilon /M
            ! boson force is 10 times larger that fermion force for J=1.0, K=0.0
            do ncount =1, splitting_M
                call cpu_time(tstart)
                pfield(:,:,:) = pfield(:,:,:) - (dt/dble(splitting_M))/2.d0*hybrid_force2(:,:,:)
                xfield(:,:,:) = xfield(:,:,:) + (dt/dble(splitting_M))*pm*pfield(:,:,:)
                call cpu_time(tend); time_vec(4) = time_vec(4) + (tend-tstart)

                ! calculate the force boson
                call cpu_time(tstart)
                call fthmc_core_cal_force_bos(phi, phi_u1, latt)
                call cpu_time(tend); time_vec(7) = time_vec(7) + (tend-tstart)

                call cpu_time(tstart)
                pfield(:,:,:) = pfield(:,:,:) - (dt/dble(splitting_M))/2.d0*hybrid_force2(:,:,:)
                call cpu_time(tend); time_vec(4) = time_vec(4) + (tend-tstart)

                ! for test
                !do j = 1, lfam
                !    write(111, '(e16.8)') hybrid_force2(j,1,1)
                !enddo
            enddo

            ! for test: output gauge field after inner dynamics
            !do nt =1, ltrot
            !do nf =1, nfam
            !do lf =1, lfam
            !    write(*,*) xfield(lf, nf, nt)
            !enddo
            !enddo
            !enddo
            !stop

            ! lastly update the momentum using force fermion (prepare the force fermion before update)
            ! update the B matrix
            !call fthmc_core_x2B
            call phi_u1%vi_to_expvi(latt)

            ! calculate the inverse of MM matrix: (M^\dagger *M)^-1 phi = x_vec
            call cpu_time(tstart)
            call fthmc_core_cg(phi, phi_u1, latt)
            call cpu_time(tend); time_vec(1) = time_vec(1) + (tend-tstart)

            ! calculate the force with the updated xfield and M matrix
            call fthmc_core_cal_force_fer(phi, phi_u1, latt)

            call cpu_time(tstart)
            pfield(:,:,:) = pfield(:,:,:) - dt/2.d0 * hybrid_force1(:,:,:)
            call cpu_time(tend); time_vec(4) = time_vec(4) + (tend-tstart)

            ! for test
            !do j = 1, lfam
            !    write(222, '(e16.8)') hybrid_force1(j,1,1)
            !enddo
        enddo

        NULLIFY(xfield)
    endsubroutine fthmc_core_md_splitting


    subroutine fthmc_sweep_hybrid(lupdate, lmeasure_equaltime, lmeasure_dyn, P0, T0, gfun0, phi, phi_u1, latt0)
        implicit none
        logical, intent(in) :: lupdate, lmeasure_equaltime, lmeasure_dyn
        type(phy0), intent(inout) :: P0
        type(tdm),  intent(inout) :: T0
        type(gfun),  intent(inout) :: gfun0
        class(fthmc_phi), intent(inout) :: phi(int(Nflavor/2.d0))
        class(ftdqmc_auxfield_f5), intent(inout) :: phi_u1
        class(ftdqmc_latt), intent(inout) :: latt0

        ! local variables
        integer :: nt, i, nf
        real(dp) :: tmp, ratio, ener_old, ener_new, random

        if ( lupdate ) then
            ! initial momentum field and R field drawn from Gaussian distribution, also prepare the x_vec
            call fthmc_core_initfield(phi, phi_u1, latt0)

            ! calculate the energy before the MD
            call fthmc_core_cal_ener(ener_old, phi, phi_u1, latt0)
            ener_pot_old = ener_pot
            ! save the gauge field configuration in case of rejection
            !call fthmc_core_save
            call phi_u1%ftdqmc_auxfield_save()

            ! molecular dynamics
            !call fthmc_core_md(phi)
            call fthmc_core_md_splitting(phi, phi_u1, latt0)

            ! calculate the energy after the MD
            call fthmc_core_cal_ener(ener_new, phi, phi_u1, latt0)
            ener_pot_new = ener_pot
            ! output the energy log just the main process
            if ( irank .eq. 0) then
                write(fout3,'(3e16.8)') ener_new, ener_old, ener_new - ener_old
            endif

            ! calculate the ratio and accept/reject
            ratio = ener_new - ener_old
            if ( ratio .lt. 0.d0 ) then
                ratio = 1.01d0
            else
                ratio = dexp ( -ratio )
            endif

            ! accept/reject
            random = spring_sfmt_stream()
            !if ( ratio .gt. random) then
            ! for test
            if ( random .gt. 0.5d0) then
            !if ( .true. ) then
            !if ( .false. ) then
                main_obs(1) = main_obs(1) + dcmplx( 1.d0, 1.0d0)
            else
                main_obs(1) = main_obs(1) + dcmplx( 0.d0, 1.0d0)
                ! restore old configuration
                !call fthmc_core_restore
                call phi_u1%ftdqmc_auxfield_restore()
                ener_pot_new = ener_pot_old
            endif
        endif

        ! equaltime measurement and possible dynamical measurements
        if ( lmeasure_equaltime ) then
            ! initialize the green's function matrices
            call fthmc_gfun_initgt(gfun0)

            ! calculate green's function and do measurements
            call fthmc_core_cal_gfun_meas(gfun0, lmeasure_equaltime, lmeasure_dyn, P0, T0, phi_u1, latt0)
            !call fthmc_hybrid_cal_gfun_meas_singlesource(gfun0, lmeasure_equaltime, lmeasure_dyn, P0, T0, phi_u1, latt0)
        endif
    end subroutine fthmc_sweep_hybrid

end module fthmc_core
