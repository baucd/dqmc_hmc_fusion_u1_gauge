module fthmc_core
    use spring
    use ftdqmc_hamilt
    use fthmc_phi_class
    use fthmc_hybrid
    use fthmc_latt
    use fthmc_gfun
    use fthmc_phy0
    use fthmc_tdm
    implicit none

contains

    subroutine fthmc_sweep_hybrid(lupdate, lmeasure_equaltime, lmeasure_dyn, lfourier, P0, T0, gfun0, phi)
        implicit none
        logical, intent(in) :: lupdate, lmeasure_equaltime, lmeasure_dyn, lfourier
        type(phy0), intent(inout) :: P0
        type(tdm),  intent(inout) :: T0
        type(gfun),  intent(inout) :: gfun0
        class(fthmc_phi), intent(inout) :: phi(int(Nflavor/2.d0))

        ! local variables
        integer :: nt, i, nf
        real(dp) :: tmp, ratio, ener_old, ener_new, random
        logical :: isfourier

        ! Fourier acceleration
        complex(dp) :: z_tmp(ltrot)

        ! set isfourier
        isfourier = lfourier
        if ( lupdate ) then
            ! initial momentum field and R field drawn from Gaussian distribution, also prepare the x_vec
            call fthmc_hybrid_initfield(phi)

            ! calculate the energy before the MD
            call fthmc_hybrid_cal_ener(ener_old, phi)
            ener_pot_old = ener_pot
            ! save the gauge field configuration in case of rejection
            call fthmc_hybrid_save

            ! molecular dynamics
            !call fthmc_hybrid_md(phi, isfourier)
            call fthmc_hybrid_md_splitting(phi, isfourier)

            ! calculate the energy after the MD
            call fthmc_hybrid_cal_ener(ener_new, phi)
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
            random = spring_sfmt_stream()
            if ( ratio .gt. random) then
                main_obs(1) = main_obs(1) + dcmplx( 1.d0, 1.0d0)
            else
                main_obs(1) = main_obs(1) + dcmplx( 0.d0, 1.0d0)
                ! restore old configuration
                call fthmc_hybrid_restore
                ener_pot_new = ener_pot_old
            endif
        endif

        ! equaltime measurement and possible dynamical measurements
        if ( lmeasure_equaltime ) then
            ! initialize the green's function matrices
            call fthmc_gfun_initgt(gfun0)

            ! calculate green's function and do measurements
            call fthmc_hybrid_calgfun_meas(gfun0, lmeasure_equaltime, lmeasure_dyn, P0, T0)
            !call fthmc_hybrid_calgfun_meas_singlesource(gfun0, lmeasure_equaltime, lmeasure_dyn, P0, T0)
        endif
    end subroutine fthmc_sweep_hybrid

end module fthmc_core
