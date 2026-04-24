program ftdqmc_main
#ifdef _OPENMP
    use OMP_LIB
#endif

#ifdef CUDA
    use magma
#endif
    use mpi
    use ftdqmc_hamilt
    use ftdqmc_io
    use ftdqmc_latt
    use ftdqmc_auxfield_class
    use ftdqmc_auxfield_f0_class
    use ftdqmc_auxfield_f1_class
    use ftdqmc_auxfield_f2_class
    use ftdqmc_auxfield_f3_class
    use ftdqmc_auxfield_f4_class
    use ftdqmc_core
    use ftdqmc_gfun
    use ftdqmc_asvqrd
    use ftdqmc_phy0
    use ftdqmc_tdm

    implicit none

    ! local
    type(ftdqmc_auxfieldHolder), allocatable :: phi(:)
    type(gfun)   :: gfun0
    type(asvqrd) :: qr
    type(phy0)   :: P0
    type(tdm)    :: T0

    integer  :: nsw
    real(dp) :: start_time, end_time, time1, time2
    character (len = 20) :: date_time_string
    character (len = 40) :: filename

    ! MPI initialization
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,irank,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,isize,ierr)

    ! CUDA initialization
#ifdef CUDA
    call magmaf_init()
#endif

#ifdef _OPENMP
    start_time = omp_get_wtime()
#else
    call cpu_time(start_time)
#endif

    if ( irank.eq.0 ) then
        open( unit=fout, file='ftdqmc.out', status='unknown' )
    end if

    ! array that record the accept ratio
    main_obs(:) = czero

    ! hamilt init
    call ftdqmc_hamilt_init

    ! latt
    call ftdqmc_latt_alloc
    call ftdqmc_latt_sli

    ! print head
    call ftdqmc_initial_head
    ! print hamilt and lattice info
    call ftdqmc_initial_print

    ! use holder
    allocate(phi(num_fields))
    allocate(ftdqmc_auxfield_f0::phi(idx_hopping)%pnt)
    allocate(ftdqmc_auxfield_f3::phi(idx_gauge)%pnt) ! gauge refers to SU(2) gauge
    allocate(ftdqmc_auxfield_f4::phi(idx_chargon)%pnt) ! for chargon field

    ! hopping field
    phi(idx_hopping)%pnt%field_index = idx_hopping
    call phi(idx_hopping)%pnt%ftdqmc_auxfield_alloc()

    ! gauge field
    phi(idx_gauge)%pnt%field_index = idx_gauge
    call phi(idx_gauge)%pnt%ftdqmc_auxfield_alloc()
    filename='confin_su2'
    call phi(idx_gauge)%pnt%ftdqmc_auxfield_inconfc(filename)

    ! chargon field
    phi(idx_chargon)%pnt%field_index = idx_chargon
    call phi(idx_chargon)%pnt%ftdqmc_auxfield_alloc()
    filename='confin_cg'
    call phi(idx_chargon)%pnt%ftdqmc_auxfield_inconfc(filename)

    ! gfun, qr and observables: P0 and T0
    if (.not. lclassic .and. .not. lrmfermion) then
        call ftdqmc_asvqrd_alloc(qr)
        call ftdqmc_gfun_alloc(gfun0)
    endif
    call ftdqmc_phy0_alloc(P0)
    if ( ltau) call ftdqmc_tdm_alloc(T0)

    ! matrix operation tmp
    if ( .not. lclassic .and. .not. lrmfermion) then
        call ftdqmc_matrix_alloc
    endif

    ! init numerical stability monitor variables
    max_wrap_error = 0.d0
    if(ltau) xmax_dyn = 0.d0

    ! get the initial UDV matrices and Green's function
    if ( .not. lclassic .and. .not. lrmfermion) then
        call ftdqmc_sweep_start_0b(qr, gfun0, phi)
        ! set weight_track and logweightf_old, fermion det weight
        logweightf_old = logweightf_up * dble(Nflavor/2.d0)
        weight_track = logweightf_old
    endif

    ! warmup
    if( lwarmup ) then
        if( irank.eq.0 ) then
            write(fout,'(a,i6)') ' nwarmup = ', nwarmup
        end if

        do nsw = 1, nwarmup
            ! local update
            if(llocal) then
                call ftdqmc_sweep_b0(lupdate=.true., lmeasure_equaltime=.false., P0=P0, T0=T0, gfun0=gfun0, qr=qr, phi=phi)
                call ftdqmc_sweep_0b(lupdate=.true., lmeasure_equaltime=.false., lmeasure_dyn=.false., P0=P0, T0=T0, gfun0=gfun0, qr=qr, phi=phi)
            end if
            ! cluster update
            if (lstglobal) then
                call ftdqmc_stglobal(gfun0=gfun0, qr=qr, phi=phi)
            endif
        end do

        if(irank.eq.0) write(fout, '(a,e16.8)') ' after warmup, max_wrap_error = ', max_wrap_error
        !if(irank.eq.0 .and. ltau) write(fout,'(a,e16.8)')' after warmup, xmax_dyn = ', xmax_dyn
        ! in warnup, xmax_dyn is not right, reset it here
        xmax_dyn = 0.d0
    end if

#ifdef _OPENMP
    time1 = omp_get_wtime()
#else
    call cpu_time(time1)
#endif

    ! DQMC starts
    if(irank.eq.0) write(fout, *)
    do nbc =  1, nbin
        ! initialize observables
        call ftdqmc_phy0_init(P0)
        if (ltau) call ftdqmc_tdm_init(T0)

        do nsw = 1, nsweep
            ! local update
            if ( llocal ) then
                call ftdqmc_sweep_b0(lupdate=.true., lmeasure_equaltime=.false., P0=P0, T0=T0, gfun0=gfun0, qr=qr, phi=phi)
                ! for dynamical measurements
                if(ltau) then
                    call push_stage(qr, gfun0)
                    call ftdqmc_sweep_0b(lupdate=.false., lmeasure_equaltime=.false., lmeasure_dyn=ltau, P0=P0, T0=T0, gfun0=gfun0, qr=qr, phi=phi)
                    call pop_stage(qr, gfun0)
                end if
                call ftdqmc_sweep_0b(lupdate=.true., lmeasure_equaltime=.true., lmeasure_dyn=.false., P0=P0, T0=T0, gfun0=gfun0, qr=qr, phi=phi)
            else
                stop 'llocal is false!'
            end if
            ! cluster update
            if (lstglobal) then
                call ftdqmc_stglobal(gfun0=gfun0, qr=qr, phi=phi)
            endif
        end do

        ! avg, ft and output configurations
        call ftdqmc_phy0_getavg(P0) ! scalar properties
        call ftdqmc_phy0_corFT(P0)  ! sq
        if(ltau) call ftdqmc_tdm_corFT(T0) ! sq(tau)

        ! save configurations
        filename='confout_su2'
        call phi(idx_gauge)%pnt%ftdqmc_auxfield_outconfc(filename) ! output configuration after every bin
        filename='confout_cg'
        call phi(idx_chargon)%pnt%ftdqmc_auxfield_outconfc(filename) ! output configuration after every bin

        ! --- Timming and outconfc
        if( nbc .eq. 1 )  then
#ifdef _OPENMP
            time2 = omp_get_wtime()
#else
            call cpu_time(time2)
#endif
            if(irank.eq.0) then
                write(fout,'(a,f8.3,a)') ' time for 1 bin: ', (time2-time1)/dble(60), ' m'
            end if
        end if

        if( irank.eq.0 .and. mod(nbc,max(nbin/10,1) ).eq.0 ) then
            write( fout, '(i5,a,i5,a)' ) nbc, '  /', nbin, '   finished'
        end if
        ! --- END Timming and outconfc
    end do

    if(irank.eq.0) write(fout, '(a,e16.8)') ' max_wrap_error = ', max_wrap_error
    if(irank.eq.0 .and. ltau) write(fout,'(a,e16.8)')' xmax_dyn = ', xmax_dyn

    ! save configurations
    filename='confout_su2'
    call phi(idx_gauge)%pnt%ftdqmc_auxfield_outconfc(filename)
    filename='confout_cg'
    call phi(idx_chargon)%pnt%ftdqmc_auxfield_outconfc(filename)

    if( irank.eq.0 ) then
#ifdef _OPENMP
        end_time = omp_get_wtime()
#else
        call cpu_time(end_time)
#endif
    endif

    ! calculate the accept ratio
    call mpi_reduce(main_obs, mpi_main_obs, size(main_obs), mpi_complex16, mpi_sum, 0, mpi_comm_world, ierr )
    if(irank.eq.0) then
        write(fout,*)
        if (aimag(main_obs(1)) .gt. 0.d0) write(fout,'(a,e16.8)') ' >>> accep_gauge  = ', dble(main_obs(1))/aimag(main_obs(1))
        if (aimag(main_obs(4)) .gt. 0.d0) write(fout,'(a,e16.8)') ' >>> accep_gauge_global  = ', dble(main_obs(4))/aimag(main_obs(4))
        if (aimag(main_obs(5)) .gt. 0.d0) write(fout,'(a,e16.8)') ' >>> accep_gauge_overrelaxation  = ', dble(main_obs(5))/aimag(main_obs(5))
        if (aimag(main_obs(2)) .gt. 0.d0) write(fout,'(a,e16.8)') ' >>> accep_chargon  = ', dble(main_obs(2))/aimag(main_obs(2))
        if (aimag(main_obs(3)) .gt. 0.d0) write(fout,'(a,e16.8)') ' >>> accep_chargon_global  = ', dble(main_obs(3))/aimag(main_obs(3))
    end if

    call ftdqmc_latt_free
    if (.not. lclassic .and. .not. lrmfermion) then
        call ftdqmc_matrix_free
    endif
    call phi(idx_gauge)%pnt%ftdqmc_auxfield_free()
    call phi(idx_chargon)%pnt%ftdqmc_auxfield_free()
    if ( ltau ) call ftdqmc_tdm_free(T0)
    call ftdqmc_phy0_free(P0)
    if (.not. lclassic .and. .not. lrmfermion) then
        call ftdqmc_gfun_free(gfun0)
        call ftdqmc_asvqrd_free(qr)
    endif
    call ftdqmc_hamilt_free

    if( irank.eq.0 ) then
        call s_time_builder(date_time_string)
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
    end if
    close(fout)

#ifdef CUDA
    call magmaf_finalize()
#endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_FINALIZE(ierr)

    stop

end program ftdqmc_main
