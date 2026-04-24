program fthmc_main
    use mpi
#IFDEF _OPENMP
    use OMP_LIB
#ENDIF
    use ftdqmc_hamilt
    use fthmc_phi_class
    use ftdqmc_latt
    use fthmc_hybrid_mat
    use fthmc_hybrid
    use fthmc_core
    use fthmc_gfun
    use fthmc_phy0
    use fthmc_tdm
    use mmpi, only : mp_bcast, mp_barrier

    implicit none

    ! local
    type(ftdqmc_phi), dimension(:), pointer :: phi
    type(gfun)   :: gfun0
    type(phy0)   :: P0
    type(tdm)    :: T0

    integer  :: nsw, icount, nsw2, i, j, nstat, nequi
    real(dp) :: start_time, end_time, time0, time1, time2, ratio
    real(dp) :: ener_pot_old_local, ener_pot_new_local
    character (len = 20) :: date_time_string

    ! timing
    real(dp) :: tstart_main, tend_main

    ! MPI initialization
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,irank,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,isize,ierr)

#IFDEF _OPENMP
    start_time = omp_get_wtime()
#ELSE
    call cpu_time(start_time)
#ENDIF

    if( irank.eq.0 ) then
        open( unit=fout, file='ftdqmc.out', status='unknown' )
        open( unit=fout2, file='cg.log', status='unknown' )
        open( unit=fout3, file='ener.log', status='unknown' )
        open( unit=fout4, file='force.log', status='unknown' )
    end if

    ! array that record the accept ratio
    main_obs(:) = czero
    ! time profiling
    time_vec(:) = zero

    ! print head
    call ftdqmc_initial_head

    ! hamilt init
    call ftdqmc_hamilt_init

    ! latt
    call ftdqmc_latt_alloc
    call ftdqmc_latt_sli
    call sltpf

    ! matrix operation tmp
    call ftdqmc_matrix_alloc

    ! pseudo fermion and initialization
    allocate(phi(int(Nflavor/2.d0)))
    do icount = 1, int(Nflavor/2.d0)
        call phi(icount)%ftdqmc_phi_alloc()
    enddo

    ! print hamilt and lattice info
    call ftdqmc_initial_print

    ! alloc space for hybrid and initialize gauge field
    call ftdqmc_hybrid_alloc
    call ftdqmc_hybrid_inconfc

    ! gfun and observables: P0 and T0
    call ftdqmc_gfun_alloc(gfun0)
    call ftdqmc_phy0_alloc(P0)
    if ( ltau) call ftdqmc_tdm_alloc(T0)

#IFDEF CUDA
    ! use default gpu card 0
    dev = 0
    call ftdqmc_gpu_conjgate_init(ndim, ltrot, nfam, lfam, reshape(l_bonds, (/2*lfam*nfam/)), itermax, errate, coshrtdtau, sinhrtdtau, csqrt2, mu_complex, dev)
#IFDEF CUDA_MEAS
    call ftdqmc_gpu_phy0_init(ndim, ltrot, reshape(list, (/ndim*2/)), reshape(nnlist,(/ndim*9/)), reshape(inv_latt_imj,(/ndim*ndim*2/)), z1, z2, z3, z4)
#ENDIF
#ENDIF

    ! decide whether to warmup
    if ( ltunedt ) lwarmup = .true.
    if ( luseinputdt ) then
        lwarmup = .false.
        dt = inputdt
    endif

    ! warmup
    if( lwarmup ) then
#IFDEF _OPENMP
    time0 = omp_get_wtime()
#ELSE
    call cpu_time(time0)
#ENDIF
        ! prepare warmup log output
        if ( irank .eq. 0) then
            open( unit=fout5, file='warmup_equi.log', status='unknown' )
            write(fout5,'(a)') ' road to equilibrium ...'
            write(fout,'(a,i8)') ' nwarmup = ', nwarmup
        endif

        ! part 1: try to reach convergence
        mdstep = 20
        nstat = 21 ! 20 is too special
        nsw = 0
        nsw2 = 0
        ratio = 0.d0
        main_obs(1) = czero
        lequi = .false.
        lexit = .false.
        do while ( .true. )
            nsw2 = nsw2 + 1

            ! test runs
            do i = 1, nstat
                nsw = nsw + 1
                call ftdqmc_sweep_hybrid(lupdate=.true., lmeasure_equaltime=.false.,lmeasure_dyn=.false., lfourier=.false., P0=P0, T0=T0, gfun0=gfun0, phi=phi)
                if ( i .eq. 1)     ener_pot_old_local = ener_pot_old
                if ( i .eq. nstat) ener_pot_new_local = ener_pot_new
            enddo

            ! calculate weight ratio
            if ( ener_pot_new_local .eq. 0) then
                weight_ratio = 1.d0
            else
                weight_ratio = abs((ener_pot_new_local - ener_pot_old_local)/ener_pot_new_local)
            endif
            ! calculate acceptance ratio
            ratio = dble(main_obs(1))/aimag(main_obs(1))
            ! reset ratio
            main_obs(1) = czero

            ! check the covergence and dt at master process
            if ( irank .eq. 0) then
                ! output log
                write(fout5, '(i6, 3e16.8)') nsw2, dt, ratio, weight_ratio
                if ( weight_ratio .lt. converge_ratio .and. ratio .gt. 0.d0) then
                    lequi = .true.
                    lexit = .true.

                    ! set mdstep for equilibration runs, when accept ratio is not to small, md is stable such that we could increase mdstep if needed
                    !mdstep = nint(mdtime/dt)
                    ! reset ratio
                    main_obs(1) = czero
                    write(fout5, '(a, i5)') ' during warmup, model is equilibrated (might be meta-stable). number of sweeps used: ', nsw
                endif
                ! will not continue to production runs
                if ( nsw .gt. nwarmup) then
                    lequi = .false.
                    lexit = .true.
                    write(fout5, '(a)') ' during warmup, model is not equilibrated. will not continue to meas. number of sweeps used: ', nsw
                endif
            endif

            ! tune step size
            if ( ratio .lt. 0.7d0 ) then
                dt = dt*0.8d0
            elseif (ratio .gt. 0.8d0) then
                dt = dt*1.2d0
            endif

            ! whether to exit
            call mp_bcast( lexit,   0 )
            call mp_bcast( lequi,   0 )
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            if ( lexit ) exit
        enddo

        ! part 2: equilibration runs
        if ( lequi ) then
            lexit = .false.
            ! just 5 equilibration runs
            write(fout5,'(a)') ' runs at equilibrium ...'
            nequi = 200
            ! do 200 sweeps, if ratio is between 0.7 to 0.8 in the middle, then exit in advance.
            do j = 1, nequi
                ! test runs
                do i = 1, nstat
                    call ftdqmc_sweep_hybrid(lupdate=.true., lmeasure_equaltime=.false.,lmeasure_dyn=.false., lfourier=.false., P0=P0, T0=T0, gfun0=gfun0, phi=phi)
                    if ( i .eq. 1)     ener_pot_old_local = ener_pot_old
                    if ( i .eq. nstat) ener_pot_old_local = ener_pot_new
                enddo

                ! calculate weight ratio
                weight_ratio = abs((ener_pot_new_local - ener_pot_old_local)/ener_pot_new_local)
                ! calculate acceptance ratio
                ratio = dble(main_obs(1))/aimag(main_obs(1))
                ! reset ratio
                main_obs(1) = czero

                if ( irank .eq. 0) then
                    write(fout5, '(i6, 3e16.8)') j, dt, ratio, weight_ratio
                    ! check acceptance ratio
                    if ( ratio .gt. 0.7d0 .and. ratio .lt. 0.8d0) then
                        lexit = .true.
                        write(fout5, '(a)') ' ending equilibration runs in advance.'
                    endif

                    if ( j .eq. nequi) then
                        lexit = .true.
                        write(fout5, '(a)') ' ending equilibration runs on time.'
                    endif
                endif

                ! whether to exit
                call mp_bcast( lexit,   0 )
                call MPI_BARRIER(MPI_COMM_WORLD,ierr)
                if ( lexit ) exit

                ! tune step size
                if ( ratio .lt. 0.7d0 ) then
                    dt = dt*0.8d0
                elseif (ratio .gt. 0.8d0) then
                    dt = dt*1.2d0
                endif
            enddo

            ! broadcast new step size
            call mp_bcast( dt,   0 )
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        endif

#IFDEF _OPENMP
    time1 = omp_get_wtime()
#ELSE
    call cpu_time(time1)
#ENDIF
        ! output log
        if ( irank .eq. 0 ) then
            ! output new md parameters to main output
            write(fout,'(a,f8.3)') ' new dt: ', dt
            write(fout,*) ' lequi: ', lequi
            write(fout,'(a,f14.6)') ' >>> accept_ratio  = ', ratio
            write(fout,'(a,f14.6,a)') ' >>> Time spent for warmup:', (time1-time0)/dble(60), 'm'
            write(fout, *)
            ! close warmup log output
            close(fout5)
        endif

        ! if not equi, stop the program.
        if ( .not. lequi ) then
            ! output configurations if not equi
            call ftdqmc_hybrid_outconfc
            stop
        endif
    else
        ! no warmup, directly into production runs
        !mdstep = nint(mdtime/dt)
        !mdstep = 20
        !mdstep = 5
        if ( irank .eq. 0 ) then
            ! output new md parameters to main output
            write(fout,'(a)') ' use dt from confin instead of ftdqmc.in'
            write(fout,'(a,f8.3)') ' new dt: ', dt
            write(fout,*) ' lequi: ', lequi
        endif

        ! time profiling
#IFDEF _OPENMP
    time1 = omp_get_wtime()
#ELSE
    call cpu_time(time1)
#ENDIF

    end if


    ! HMC starts, production runs
    do nbc =  1, nbin
        ! initialize observables
        call ftdqmc_phy0_init(P0)
        if (ltau) call ftdqmc_tdm_init(T0)

        do nsw = 1, nsweep
            call ftdqmc_sweep_hybrid(lupdate=.true., lmeasure_equaltime=.true., lmeasure_dyn=ltau, lfourier=.false., P0=P0, T0=T0, gfun0=gfun0, phi=phi)
        end do

        ! avg, ft and output configurations
        call ftdqmc_phy0_getavg(P0) ! scalar properties
        call ftdqmc_phy0_corFT(P0)  ! sq
        if(ltau) call ftdqmc_tdm_corFT(T0) ! sq(tau)

        call ftdqmc_hybrid_outconfc

        ! --- Timming
        if( nbc .eq. 1 )  then
#IFDEF _OPENMP
            time2 = omp_get_wtime()
#ELSE
            call cpu_time(time2)
#ENDIF
            if(irank.eq.0) then
				if (nbc .eq. 1) write(fout,'(a)') ' proceed to production runs.'
                write(fout,'(a,f8.3,a)') ' time for 1 bin: ', (time2-time1)/dble(60), ' m'
            end if
        end if

        if( irank.eq.0 .and. mod(nbc,max(nbin/10,1) ).eq.0 ) then
            write( fout, '(i5,a,i5,a)' ) nbc, '  /', nbin, '   finished'
        end if
        ! --- END Timming and outconfc
    end do

    if( irank.eq.0 ) then
#IFDEF _OPENMP
        end_time = omp_get_wtime()
#ELSE
        call cpu_time(end_time)
#ENDIF
    endif

    ! calculate the accept ratio
    call mpi_reduce(main_obs, mpi_main_obs, size(main_obs), mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
    if(irank.eq.0) then
        write(fout,'(a,e16.8)') ' >>> accept_ratio  = ', dble(mpi_main_obs(1))/aimag(mpi_main_obs(1))
    end if

    do icount = 1, int(Nflavor/2.d0)
        call phi(icount)%ftdqmc_phi_free()
    enddo
    call ftdqmc_latt_free
    call ftdqmc_matrix_free
    call ftdqmc_hybrid_free
    if ( ltau ) call ftdqmc_tdm_free(T0)
    call ftdqmc_phy0_free(P0)
    call ftdqmc_gfun_free(gfun0)
    call ftdqmc_hamilt_free
#IFDEF CUDA
#IFDEF CUDA_MEAS
    call ftdqmc_gpu_phy0_free()
#ENDIF
    call ftdqmc_gpu_conjgate_free()
#ENDIF

    if( irank.eq.0 ) then
        call s_time_builder(date_time_string)
        write(fout,*)
        write(fout,'(a,f14.6,a)') ' >>> Total time spent:', (end_time-start_time)/dble(60), 'm'
        write(fout,'(a,e16.8,a)') ' >>> Time for one meas bin:', (end_time-time1)/dble(nbin), 's'
        write(fout,'(a)') ' >>> Happy ending at '//date_time_string
        write(fout,*)
        write(fout,'(a)') ' The simulation done !!! '
        write(fout,*)
        write(fout,'(a)') '        o         o    '
        write(fout,'(a)') '       o o       o o   '
        write(fout,'(a)') '       o o       o o   '
        write(fout,'(a)') '        o         o    '
        write(fout,'(a)') '       o o       o o   '
        write(fout,'(a)') '       o o       o o   '
        write(fout,'(a)') '        o         o    '
		! output time profile
        write(fout, *)
        write(fout, '(a20, 2f8.3, a)') 'force_fer_cg: ', time_vec(1), time_vec(1)/(end_time-start_time)*100.d0, '%'
        write(fout, '(a20, 2f8.3, a)') 'force_fer_deri: ', time_vec(2), time_vec(2)/(end_time-start_time)*100.d0, '%'
        write(fout, '(a20, 2f8.3, a)') 'force_loop: ', time_vec(7), time_vec(7)/(end_time-start_time)*100.d0, '%'
        write(fout, '(a20, 2f8.3, a)') 'md_leapfrog: ', time_vec(4), time_vec(4)/(end_time-start_time)*100.d0, '%'
        write(fout, '(a20, 2f8.3, a)') 'meas_estimator_cg: ', time_vec(5), time_vec(5)/(end_time-start_time)*100.d0, '%'
        write(fout, '(a20, 2f8.3, a)') 'meas_others: ', time_vec(6), time_vec(6)/(end_time-start_time)*100.d0, '%'
        write(fout, '(a20, 2f8.3, a)') 'Mr_to_phi: ', time_vec(3), time_vec(3)/(end_time-start_time)*100.d0, '%'
        if ( time_vec(8) .ne. 0.d0) then
            write(fout, '(a20, 2f8.3, a)') 'test: ', time_vec(8), time_vec(8)/(end_time-start_time)*100.d0, '%'
        endif
    end if
    close(fout)
    close(fout2)
    close(fout3)
    close(fout4)

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_FINALIZE(ierr)
    stop
end program fthmc_main
