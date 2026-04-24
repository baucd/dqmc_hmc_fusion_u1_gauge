module ftdqmc_core
    use spring
    use ftdqmc_hamilt
    use ftdqmc_asvqrd
    use ftdqmc_latt
    use ftdqmc_auxfield_class
    use ftdqmc_gfun
    use ftdqmc_phy0
    use ftdqmc_tdm
    implicit none

contains
    subroutine Bmat_tau_R(phi, nt1, nt2, bmat_up)
        ! B(tau1,tau2) *
        ! make sure nt1 > nt2
        implicit none
        type(ftdqmc_auxfieldHolder), intent(inout) :: phi(num_fields)
        integer, intent(in) :: nt1, nt2
        complex(dp), dimension(ndim,ndim), intent(inout) :: bmat_up

        ! local
        integer :: nt, nflag, nf, i
        complex(dp) :: phaseu

        phaseu = cone
        do nt = nt2, nt1
            if( llocal ) then
                ! Ising spin
                !do i = 1, ndim
                !    if ( mod( list(i,1) + list(i,2), 2 ) .eq. phi(idx_spin)%pnt%idx_sub_ab ) then
                !        phi(idx_spin)%pnt%idx_site = i
                !        call phi(idx_spin)%pnt%ftdqmc_auxfield_mmr(nt, bmat_up)
                !    endif
                !enddo

                ! gauge
                do nf = 1, nfam
                    phi(idx_gauge)%pnt%idx_family = nf
                    call phi(idx_gauge)%pnt%ftdqmc_auxfield_mmr(nt, bmat_up)
                enddo

                ! chemical potential
                call phi(idx_hopping)%pnt%ftdqmc_auxfield_mmr(nt, bmat_up)

            end if
        end do
    end subroutine Bmat_tau_R

    subroutine Bmat_tau_RH( phi, nt1, nt2, bmat_up )
        ! B(tau1,tau2) *
        ! make sure nt1 > nt2
        implicit none
        type(ftdqmc_auxfieldHolder), intent(inout) :: phi(num_fields)
        integer, intent(in) :: nt1, nt2
        complex(dp), dimension(ndim,ndim), intent(inout) :: bmat_up

        ! local
        integer :: nt, nflag, nf, i
        complex(dp) :: phaseu

        phaseu = cone
        do nt = nt1, nt2, -1
            if( llocal ) then

                ! chemical potential
                call phi(idx_hopping)%pnt%ftdqmc_auxfield_mmrH(nt, bmat_up)

                ! gauge field
                do nf = nfam, 1, -1
                    phi(idx_gauge)%pnt%idx_family = nf
                    call phi(idx_gauge)%pnt%ftdqmc_auxfield_mmrH(nt, bmat_up)
                enddo

                ! Ising spin field
                !do i = ndim, 1, -1
                !    if ( mod( list(i,1) + list(i,2), 2 ) .eq. phi(idx_spin)%pnt%idx_sub_ab ) then
                !        phi(idx_spin)%pnt%idx_site = i
                !        call phi(idx_spin)%pnt%ftdqmc_auxfield_mmrH(nt, bmat_up)
                !    endif
                !enddo

            end if
        end do
    end subroutine Bmat_tau_RH

    subroutine Bmat_tau_L( phi, nt1, nt2, bmat_up )
        ! * B(tau1,tau2)
        ! make sure nt1 > nt2
        implicit none
        type(ftdqmc_auxfieldHolder), intent(inout) :: phi(num_fields)
        integer, intent(in) :: nt1, nt2
        complex(dp), dimension(ndim,ndim), intent(inout) :: bmat_up

        ! local
        integer :: nt, nflag, nf, i
        complex(dp) :: phaseu

        phaseu = cone
        do nt = nt1, nt2, -1
            if( llocal ) then

                ! chemical potential
                call phi(idx_hopping)%pnt%ftdqmc_auxfield_mml(nt, bmat_up)

                ! gauge field
                do nf = nfam, 1, -1
                    phi(idx_gauge)%pnt%idx_family = nf
                    call phi(idx_gauge)%pnt%ftdqmc_auxfield_mml(nt, bmat_up)
                enddo

                ! Ising spin field
                !do i = ndim, 1, -1
                !    if ( mod( list(i,1) + list(i,2), 2 ) .eq. phi(idx_spin)%pnt%idx_sub_ab ) then
                !        phi(idx_spin)%pnt%idx_site = i
                !        call phi(idx_spin)%pnt%ftdqmc_auxfield_mml(nt, bmat_up)
                !    endif
                !enddo

            end if
        end do
    end subroutine Bmat_tau_L

    subroutine Bmatinv_tau_L( phi, nt1, nt2, bmat_up )
        ! *B(tau1,tau2)^-1
        ! make sure nt1 > nt2
        implicit none
        type(ftdqmc_auxfieldHolder), intent(inout) :: phi(num_fields)
        integer, intent(in) :: nt1, nt2
        complex(dp), dimension(ndim,ndim), intent(inout) :: bmat_up

        ! local
        integer :: nt, nflag, nf, i
        complex(dp) :: phaseu

        phaseu = cone
        do nt = nt2, nt1
            if( llocal ) then
                ! Ising spin field
                !do i = 1, ndim
                !    if ( mod( list(i,1) + list(i,2), 2 ) .eq. phi(idx_spin)%pnt%idx_sub_ab ) then
                !        phi(idx_spin)%pnt%idx_site = i
                !        call phi(idx_spin)%pnt%ftdqmc_auxfield_mmlm1(nt, bmat_up)
                !    endif
                !enddo

                ! gauge field
                do nf = 1, nfam
                    phi(idx_gauge)%pnt%idx_family = nf
                    call phi(idx_gauge)%pnt%ftdqmc_auxfield_mmlm1(nt, bmat_up)
                enddo

                ! chemical potential
                call phi(idx_hopping)%pnt%ftdqmc_auxfield_mmlm1(nt, bmat_up)
            end if
        end do
    end subroutine Bmatinv_tau_L


    subroutine ftdqmc_sweep_start_0b(qr, gfun0, phi)
        implicit none
        type(asvqrd), intent(inout) :: qr
        type(gfun), intent(inout) :: gfun0
        type(ftdqmc_auxfieldHolder), intent(inout) :: phi(num_fields)

        ! for test
        INTEGER :: j,nx, ny, idx1, idx2
        COMPLEX(dp) :: grup(lq, lq)

        integer :: n, i, info
        real(dp) :: tmp

        if ( nst .gt. 0 ) THEN
            ! at tau = 0
            call ftdqmc_asvqrd_initst(qr, 0)
            gfun0%grup(:,:) = Imat(:,:)

            do n = 1, nst
                ! at tau = n * tau1
                qr%Bdtau1_up(:,:) = qr%Ust_up(:,:,n-1)
                call Bmat_tau_R( phi, wrap_step(2,n), wrap_step(1,n), qr%Bdtau1_up)
                call ftdqmc_stablize_0b_svd(n, qr, phi)
            end do

            ! at tau = beta
            call ftdqmc_asvqrd_initR(qr, nst)
            call green_equaltimebb( nst, ndim, qr, gfun0, logweightf_up, info )

            ! for test: benchmark equaltime green's function
            !do i =1, lq
            !    do j =1, lq
            !        ! get grup
            !        nx = list_prim(i, 1)
            !        ny = list_prim(i, 2)
            !        idx1 = invlist(nx, ny, 1)
            !        nx = list_prim(j, 1)
            !        ny = list_prim(j, 2)
            !        idx2 = invlist(nx, ny, 1)
            !        grup(i,j) = gfun0%grup(idx1, idx2)
            !        write(*, '(2e16.8)') grup(i,j)
            !    enddo
            !enddo
            !stop

            if( info .eq. -1 ) then
                write(fout,'(a)') ' WRONG in sweep_start, exit '
            end if

        else
            if (lstglobal) then
                write(fout, *) 'lstglobal is true. gauge block update does not work with nst = 0. stopping...'
                write(fout, *) 'it actually depends on whether gauge block is turned on, chargon block is okay.'
                stop
            endif

            if ( llocal ) then
                write(fout,*) 'sweep_b0, llocal = true'
                qr%Bdtau1_up(:,:) = Imat(:,:)
                call Bmat_tau_R( phi, ltrot, 1, qr%Bdtau1_up)
                do  i = 1, ndim
                    qr%Bdtau1_up(i,i) = qr%Bdtau1_up(i,i) + cone
                end do
                call s_invlu_z( ndim, qr%Bdtau1_up )
                gfun0%grup(:,:) = qr%Bdtau1_up
            else
                gfun0%grup(:,:) = Imat(:,:)
                call Bmat_tau_R( phi, ltrot, 1, gfun0%grup)
                do  i = 1, ndim
                    gfun0%grup(i,i) = gfun0%grup(i,i) + cone
                end do
            end if
        END if
    end subroutine ftdqmc_sweep_start_0b


    subroutine ftdqmc_sweep_start_b0(qr, gfun0, phi)
        implicit none
        type(asvqrd), intent(inout) :: qr
        type(gfun), intent(inout) :: gfun0
        type(ftdqmc_auxfieldHolder), intent(inout) :: phi(num_fields)

        integer :: n, i, info
        real(dp) :: tmp

        if ( nst .gt. 0 ) THEN

            ! at tau = beta
            call ftdqmc_asvqrd_initst(qr, nst)

            do n = nst, 1, -1
                ! at tau = (n-1) * tau1
                ! calculate B(n*tau1,(n-1)*tau1), and set Vst(:,:,n-1), Dst(:,:,n-1), Ust(:,:,n-1)
                qr%Bdtau1_up(:,:) = qr%Ust_up(:,:,n)
                ! Bdtau1_up = U_up*Bup(tau+dtau,tau)
                call Bmat_tau_RH( phi, wrap_step(2,n), wrap_step(1,n), qr%Bdtau1_up) ! we need to get ( U*B(n*tau1,(n-1)*tau1) )^H = B^H * U^H
                call ftdqmc_stablize_b0_svd(n, qr,phi)
            end do

            ! at tau = 0
            call ftdqmc_asvqrd_initL(qr, 0)
            call green_equaltime00( nst, ndim, qr, gfun0, logweightf_up, info )

            if( info .eq. -1 ) then
                write(fout,'(a)') ' WRONG in sweep_start, exit '
                stop
            end if

        else

            if ( llocal ) then
                qr%Bdtau1_up(:,:) = Imat(:,:)
                call Bmat_tau_R( phi, ltrot, 1, qr%Bdtau1_up)
                do  i = 1, ndim
                    qr%Bdtau1_up(i,i) = qr%Bdtau1_up(i,i) + cone
                end do
                call s_invlu_z( ndim, qr%Bdtau1_up )
                gfun0%grup(:,:) = qr%Bdtau1_up
            else
                gfun0%grup(:,:) = Imat(:,:)
                call Bmat_tau_R( phi, ltrot, 1, gfun0%grup)
                do  i = 1, ndim
                    gfun0%grup(i,i) = gfun0%grup(i,i) + cone
                end do
            end if
        END if
    end subroutine ftdqmc_sweep_start_b0


    subroutine ftdqmc_sweep_transgf_b0(phi, qr, gfun0, nt)
        type(ftdqmc_auxfieldHolder), intent(inout) :: phi(num_fields)
        type(asvqrd), intent(inout) :: qr
        type(gfun),   intent(inout) :: gfun0
        integer, intent(in) :: nt

        integer :: n, i, j, nt_ob, ilq, it, nn_ilq, nn_it, inn_st, info, nt1, nt2
        if ( (iwrap_nt(nt-1) .gt. 0 ) .or. ( nt.eq.1 .and. nst .gt. 0 ) ) then
            n = iwrap_nt(nt-1)
            ! at tau = n * tau1
            call ftdqmc_asvqrd_initR(qr, n)

            qr%Bdtau1_up(:,:) = qr%Ust_up(:,:,n+1)
            ! Bdtau1_up = U_up*Bup(tau+dtau,tau)
            call Bmat_tau_RH( phi, wrap_step(2,n+1), wrap_step(1,n+1), qr%Bdtau1_up) ! we need to get ( U*B(n*tau1,(n-1)*tau1) )^H = B^H * U^H
            call ftdqmc_stablize_b0_svd(n+1, qr, phi)

            call ftdqmc_asvqrd_initL(qr, n)

            call green_equaltime( n, ndim, qr, gfun0, logweightf_up, info )

            call s_compare_max_z( ndim, gfun0%grup_tmp, gfun0%grup, max_wrap_error_tmp )
            if( max_wrap_error_tmp .gt. max_wrap_error ) max_wrap_error = max_wrap_error_tmp
            if( info .eq. 0 ) gfun0%grup(:,:) = gfun0%grup_tmp(:,:)
        end if
    endsubroutine ftdqmc_sweep_transgf_b0


    subroutine ftdqmc_sweep_transgf_0b(phi, qr, gfun0, nt, lmeasure_dyn)
        type(ftdqmc_auxfieldHolder), intent(inout) :: phi(num_fields)
        type(asvqrd), intent(inout) :: qr
        type(gfun),   intent(inout) :: gfun0
        integer, intent(in) :: nt
        logical, intent(in) :: lmeasure_dyn

        integer :: n, i, j, nt_ob, ilq, it, nn_ilq, nn_it, inn_st, info, nt1, nt2
        if ( iwrap_nt(nt) .gt. 0 ) then
            n = iwrap_nt(nt)
            ! at tau = n * tau1
            call ftdqmc_asvqrd_initL(qr, n)

            qr%Bdtau1_up(:,:) = qr%Ust_up(:,:,n-1)
            call Bmat_tau_R( phi, wrap_step(2,n), wrap_step(1,n), qr%Bdtau1_up)
            call ftdqmc_stablize_0b_svd(n, qr, phi)

            call ftdqmc_asvqrd_initR(qr, n)

            if( .not. lmeasure_dyn ) then
                call green_equaltime( n, ndim, qr, gfun0, logweightf_up, info )
            else
            ! only when we need measure dynamical quantities, we will call green_tau
#ifdef DYNERROR
                ! B(nt1,nt2) with nt1 >= nt2
                nt1 = nt
                nt2 = nt

                ! G(t',0) = B(t',t) * G(t,0)
                call Bmat_tau_R( phi, nt1, nt2, gfun0%gt0up)

                ! G(0,t') = G(0,t) * B(t',t)^-1
                call Bmatinv_tau_L( phi, nt1, nt2, gfun0%g0tup)
#endif

                ! calculate the dynamical green's function for spin up and spin down
                call  green_tau(n, ndim, qr, gfun0, logweightf_up, info )

                ! gt0
#ifdef DYNERROR
                call s_compare_max_z( ndim, gfun0%gt0up, gfun0%gt0up_tmp, xmax_dyn_tmp )
                if( xmax_dyn_tmp .gt. xmax_dyn ) xmax_dyn = xmax_dyn_tmp
#endif
                gfun0%gt0up = gfun0%gt0up_tmp

                ! g0t
#ifdef DYNERROR
                call s_compare_max_z( ndim, gfun0%g0tup, gfun0%g0tup_tmp, xmax_dyn_tmp )
                if( xmax_dyn_tmp .gt. xmax_dyn ) xmax_dyn = xmax_dyn_tmp
#endif
                gfun0%g0tup = gfun0%g0tup_tmp
            end if

            ! equaltime green's function
            call s_compare_max_z( ndim, gfun0%grup_tmp, gfun0%grup, max_wrap_error_tmp )
            if( max_wrap_error_tmp .gt. max_wrap_error ) max_wrap_error = max_wrap_error_tmp
            if( info .eq. 0 ) gfun0%grup(:,:) = gfun0%grup_tmp(:,:)
         end if
    endsubroutine ftdqmc_sweep_transgf_0b


    subroutine ftdqmc_sweep_b0(lupdate, lmeasure_equaltime, P0, T0, gfun0, qr, phi)
        implicit none
        logical, intent(in) :: lupdate, lmeasure_equaltime
        type(phy0), intent(inout) :: P0
        type(tdm),  intent(inout) :: T0
        type(gfun),  intent(inout) :: gfun0
        type(asvqrd),  intent(inout) :: qr
        type(ftdqmc_auxfieldHolder), intent(inout) :: phi(num_fields)

        ! local variables
        integer :: nf, nt, nflag, i, j, nt_ob, ilq, it, nn_ilq, nn_it, inn_st, info, nt1, nt2
        real(dp) :: tmp, ratiof, ratiofi

        ! at tau = beta
        if ( .not. lclassic .and. .not. lrmfermion) then
            if (nst.gt.0) call ftdqmc_asvqrd_initst(qr, nst)
        endif
        nt_ob = ceiling( spring_sfmt_stream() * ltrot )

        do nt = ltrot, 1, -1
            if (.not. lnoupdate) then
                ! handle chemicial potential
                if (.not. lclassic .and. .not. lrmfermion) then
                    call phi(idx_hopping)%pnt%ftdqmc_auxfield_mml(nt, gfun0%grup)
                    call phi(idx_hopping)%pnt%ftdqmc_auxfield_mmrm1(nt, gfun0%grup)
                endif

                ! update gauge field
                do nf = nfam, 1, -1
                    phi(idx_gauge)%pnt%idx_family = nf
                    if(lupdate .and. .not. lfreezeU) then
                        call phi(idx_gauge)%pnt%ftdqmc_auxfield_upgrade(phi, nt, gfun0%grup)
                        ! for test
                        !write(*,*) 'in over'
                        if (loverrelaxation) then
                            call phi(idx_gauge)%pnt%ftdqmc_auxfield_upgrade_overrelaxation(phi, nt, gfun0%grup)
                        endif
                        ! for test
                        !write(*,*) 'out over'
                    endif

                    if (.not. lclassic .and. .not. lrmfermion ) then
                        ! propagate the green's function
                        call phi(idx_gauge)%pnt%ftdqmc_auxfield_mml(nt, gfun0%grup)
                        call phi(idx_gauge)%pnt%ftdqmc_auxfield_mmrm1(nt, gfun0%grup)
                    endif
                enddo
                if (.not. lclassic .and. .not. lrmfermion ) then
                    ! stablization of green's function
                    call ftdqmc_sweep_transgf_b0(phi, qr, gfun0, nt)
                endif

                ! update chargon field
                if (lupdate .and. .not. lrmchargon) call phi(idx_chargon)%pnt%ftdqmc_auxfield_upgrade(phi, nt, gfun0%grup)
            endif
        end do
    end subroutine ftdqmc_sweep_b0


    subroutine ftdqmc_sweep_0b(lupdate, lmeasure_equaltime, lmeasure_dyn, P0, T0, gfun0, qr, phi)
        implicit none
        logical, intent(in) :: lupdate, lmeasure_equaltime, lmeasure_dyn
        type(phy0), intent(inout) :: P0
        type(tdm),  intent(inout) :: T0
        type(gfun),  intent(inout) :: gfun0
        type(asvqrd),  intent(inout) :: qr
        type(ftdqmc_auxfieldHolder), intent(inout) :: phi(num_fields)

        ! local variables
        integer :: nf, nt, nflag, i, j, nt_ob, ilq, it, nn_ilq, nn_it, inn_st, info, nt1, nt2
        real(dp) :: tmp, ratiof, ratiofi
        logical :: lmeasure_dyn_local
        ! for test
        integer :: ncount

        lmeasure_dyn_local = lmeasure_dyn

        if ( .not. lclassic .and. .not. lrmfermion ) then
            ! at tau = 0
            if (nst.gt.0) call ftdqmc_asvqrd_initst(qr, 0)

            ! and measure the \Delta tau = 0 dyn correlation function
            if( lmeasure_dyn_local .and. (.not. lupdate)) then
                call ftdqmc_gfun_initgt(gfun0)

                ! perform the dynamical measurement
                nt = 0
                call ftdqmc_tdm_meas(nt,gfun0,T0,phi)
            end if
        endif

        nt_ob = ceiling( spring_sfmt_stream() * ltrot )
        do nt = 1, ltrot, 1
            if (.not. lnoupdate) then
                ! update chargon field
                if (lupdate .and. .not. lrmchargon) call phi(idx_chargon)%pnt%ftdqmc_auxfield_upgrade(phi, nt, gfun0%grup)

                ! update the gauge field
                do nf = 1, nfam
                    ! propagate the green's function
                    phi(idx_gauge)%pnt%idx_family = nf
                    if ( .not. lclassic .and. .not. lrmfermion ) then
                        call phi(idx_gauge)%pnt%ftdqmc_auxfield_mmr(nt, gfun0%grup)
                        call phi(idx_gauge)%pnt%ftdqmc_auxfield_mmlm1(nt, gfun0%grup)
                    endif

                    ! update the gauge field
                    if( lupdate .and. .not. lfreezeU) then
                        call phi(idx_gauge)%pnt%ftdqmc_auxfield_upgrade(phi, nt, gfun0%grup)
                        ! for test
                        !write(*,*) 'in over'
                        if ( loverrelaxation) then
                            call phi(idx_gauge)%pnt%ftdqmc_auxfield_upgrade_overrelaxation(phi, nt, gfun0%grup)
                        endif
                        ! for test
                        !write(*,*) 'out over'
                    endif

                enddo

                ! handle chemicial potential
                if ( .not. lclassic .and. .not. lrmfermion ) then
                    call phi(idx_hopping)%pnt%ftdqmc_auxfield_mmr(nt, gfun0%grup)
                    call phi(idx_hopping)%pnt%ftdqmc_auxfield_mmlm1(nt, gfun0%grup)

                    ! stablization of green's function
                    call ftdqmc_sweep_transgf_0b(phi, qr, gfun0, nt, lmeasure_dyn_local)
                endif
            endif

            ! equaltime measurement
            if( lmeasure_equaltime ) then
                call ftdqmc_phy0_meas(nt, gfun0, P0, phi)
            end if

            ! dynamical measurement
            if ( .not. lclassic .and. .not. lrmfermion ) then
                if ( (lmeasure_dyn_local .and. (.not.lupdate)) .and. ( nt .ne. ltrot) ) then
                    if( iwrap_nt(nt) .gt. 0 ) then
                        ! at stablization point, already have gt0,g0t from green_tau
                    else
                        ! B(nt1,nt2) with nt1 >= nt2
                        nt1 = nt
                        nt2 = nt
                        ! G(t',0) = B(t',t) * G(t,0)
                        call Bmat_tau_R( phi, nt1, nt2, gfun0%gt0up)

                        ! G(0,t') = G(0,t) * B(t',t)^-1
                        call Bmatinv_tau_L( phi, nt1, nt2, gfun0%g0tup)
                    end if

                    ! perform the dynamical measurement
                    call ftdqmc_tdm_meas(nt,gfun0,T0,phi)
                endif
            endif

        end do
    end subroutine ftdqmc_sweep_0b


    subroutine ftdqmc_stglobal(gfun0, qr, phi)
        implicit none
        type(gfun),  intent(inout) :: gfun0
        type(asvqrd),  intent(inout) :: qr
        type(ftdqmc_auxfieldHolder), intent(inout) :: phi(num_fields)

        ! local variabels
        integer :: i, nblock_gauge, lf, nf, nt
        real(dp) :: boson_ratio, ratio_re
        complex(dp) :: fermion_ratio

        ! chargon global update, all space-time confs are involved.
        if ( .not. lrmchargon) then
            call phi(idx_chargon)%pnt%ftdqmc_auxfield_upgrade_global(gfun0%grup, phi)
        endif

        ! gauge global update, only one chain is involved.
        if ( .not. lfreezeU) then
            if ( .not. lrmfermion) then
                logweightf_old = weight_track
                nblock_gauge = nblock
            else
                ! number of bonds: lq*2
                nblock_gauge = lq*2
            endif

            do i = 1, nblock_gauge
                ! set lf and nf
                if ( lrmfermion ) then
                    ! without fermion, nblock_gauge is number of bonds
                    nf = mod(i-1, nfam) + 1
                    lf = (i - 1) / nfam + 1
                else
                    ! with fermion, nblock_gauge is fixed
                    lf = int(spring_sfmt_stream()*lfam) + 1
                    nf = int(spring_sfmt_stream()*nfam) + 1
                endif
                ! sync indices
                phi(idx_gauge)%pnt%idx_lf = lf
                phi(idx_gauge)%pnt%idx_family = nf

                ! propose update and get the gauge boson part ratio
                call phi(idx_gauge)%pnt%ftdqmc_auxfield_upgrade_global(gfun0%grup, phi)
                boson_ratio = exp(-dtau*phi(idx_gauge)%pnt%delt_E_global)

                ! for test: output the ratio
                !write(222, '(2i5, f16.8)') phi(idx_gauge)%pnt%idx_lf, phi(idx_gauge)%pnt%idx_family, boson_ratio
                !write(222, '(2i5, f16.8)') phi(idx_gauge)%pnt%idx_lf, phi(idx_gauge)%pnt%idx_family, phi(idx_gauge)%pnt%delt_E_global*dtau

                if (.not. lrmfermion) then
                    ! store old UDV
                    call push_stage(qr, gfun0)
                    ! get fermion ratio
                    call ftdqmc_sweep_start_0b(qr, gfun0, phi)
                    logweightf_new = logweightf_up * dble(Nflavor/2.d0)
                    fermion_ratio = exp(logweightf_new - logweightf_old)

                    ! set ratio_re
                    ratio_re = dble(fermion_ratio) * boson_ratio
                else
                    ratio_re = boson_ratio
                endif

                ! accept or reject
                if (ratio_re .gt. spring_sfmt_stream()) then
                !if ( .false. ) then
                !if ( .true. ) then
                    if (.not. lrmfermion) then
                        weight_track = logweightf_new
                        logweightf_old = logweightf_new
                    endif
                    main_obs(4) = main_obs(4) + dcmplx(1.d0, 1.d0)
                else
                    if (.not. lrmfermion) then
                        ! restore UDV matrices
                        call pop_stage(qr, gfun0)

                        ! restore gauge field
                        do nt = 1, ltrot
                            phi(idx_gauge)%pnt%gauge_su2(:, lf, nf, nt) = phi(idx_gauge)%pnt%gauge_su2_chain(:, nt)
                        enddo
                    endif

                    main_obs(4) = main_obs(4) + dcmplx(0.d0, 1.d0)
                endif
            enddo
        endif
    endsubroutine ftdqmc_stglobal


end module ftdqmc_core
