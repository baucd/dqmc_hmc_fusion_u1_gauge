module fthmc_core
    use spring
    use ftdqmc_hamilt
    use fthmc_phi_class
    use fthmc_hybrid_mat
    use fthmc_hybrid
    use ftdqmc_latt
    use ftdqmc_gfun
    use fthmc_phy0
    use fthmc_tdm
    implicit none

contains

    subroutine ftdqmc_sweep_hybrid(lupdate, lmeasure_equaltime, lmeasure_dyn, lfourier, P0, T0, gfun0, phi)
        implicit none
        logical, intent(in) :: lupdate, lmeasure_equaltime, lmeasure_dyn, lfourier
        type(phy0), intent(inout) :: P0
        type(tdm),  intent(inout) :: T0
        type(gfun),  intent(inout) :: gfun0
        class(ftdqmc_phi), intent(inout) :: phi(int(Nflavor/2.d0))

        ! local variables
        integer :: nt, i, nf
        real(dp) :: tmp, ratio, ener_old, ener_new, random
        logical :: isfourier

        ! Fourier acceleration
        complex(dp) :: z_tmp(ltrot)

        ! set isfourier
        isfourier = lfourier
        if ( lupdate ) then
            ! initialize the Green's function matrices
            call ftdqmc_gfun_initgt(gfun0)

            ! initial momentum field and R field drawn from Gaussian distribution, also prepare the x_vec
            call ftdqmc_hybrid_initfield(phi)
            ! for test: cg gpu
            !stop

            ! calculate the energy before the MD
            call ftdqmc_hybrid_cal_ener(ener_old, phi)
            ener_pot_old = ener_pot
            ! for test
            !ener_q_old = ener_q
            !ener_p_old = ener_p
            !ener_js_old = ener_js
            !ener_jpi_old = ener_jpi
            !ener_ek_old = ener_ek

            ! save the gauge field configuration in case of rejection
            call ftdqmc_hybrid_save

            ! firstly move the momentum field to t + dt/2
            call ftdqmc_hybrid_cal_force(phi)

            ! for Fourier acceleration
            if ( lfourier ) then
                do nf = 1, nfam
                    do i = 1, lfam
                        z_tmp(:) = dcmplx(hybrid_force(i, nf, :), 0.d0)
                        call onedimension_fft(ltrot, z_tmp)
                        hybrid_force_k_z(i, nf, :) = z_tmp(:)
                    enddo
                enddo
            endif

            ! Molecular Dynamics
            do i = 1, mdstep
                !call ftdqmc_hybrid_md(phi, isfourier)
                call ftdqmc_hybrid_md_splitting(phi, isfourier)
            enddo

            ! for test: cg
            !call ftdqmc_hybrid_cg_test(phi)
            !write(*,*) 'end of cg test and stop the program'
            !stop

            ! calculate the energy after the MD
            call ftdqmc_hybrid_cal_ener(ener_new, phi)
            ener_pot_new = ener_pot

            ! calculate the ratio
            ratio = ener_new - ener_old

            ! output the energy log just the main process
            if ( irank .eq. 0) then
                write(fout3,'(3e16.8)') ener_new, ener_old, ratio
            endif

            if ( ratio .lt. 0.d0 ) then
                ratio = 1.01d0
            else
                ratio = dexp ( -ratio )
            endif

            ! Metropolis to ensure detailed balanced
            random = spring_sfmt_stream()
            if ( ratio .gt. random) then
                ! nothing need to be done
                main_obs(1) = main_obs(1) + dcmplx( 1.d0, 1.0d0)

                ! for test
                !write(*,'(5e16.8)') ener_q, ener_p, ener_js, ener_jpi, ener_ek
            else
                main_obs(1) = main_obs(1) + dcmplx( 0.d0, 1.0d0)
                call ftdqmc_hybrid_restore
                ener_pot_new = ener_pot_old

                ! for test
                !write(*,'(5e16.8)') ener_q_old, ener_p_old, ener_js_old, ener_jpi_old, ener_ek_old
            endif
        endif

        ! equaltime measurement and possible dynamical measurements
        if ( lmeasure_equaltime ) then
            call ftdqmc_hybrid_calgfun_meas(gfun0, lmeasure_equaltime, lmeasure_dyn, P0, T0)
            !call ftdqmc_hybrid_calgfun_meas_singlesource(gfun0, lmeasure_equaltime, lmeasure_dyn, P0, T0)
        endif
    end subroutine ftdqmc_sweep_hybrid

end module ftdqmc_core
