module ftdqmc_io
      use ftdqmc_hamilt

contains

    subroutine fthmc_hamilt_init
        use mpi
        use parser, only : p_create, p_parse, p_get, p_get_vec, p_destroy
        use mmpi, only : mp_bcast, mp_barrier
        implicit none

        integer :: i, nwrap_mid
        real(dp) :: r1, r2, eta
        logical :: exists

        ! default parameters
        l    = 2
        beta = 20
        dtau = 0.1d0
        mu   = 0.d0 ! default is half filling
        muA   = 0.d0
        muB   = 0.d0

        ltau = .false.
        nuse = 0

        xmag = 0.d0
        flux_x = 0.d0
        flux_y = 0.d0

        Nflavor = 2.d0
        rt = 1.d0
        js = 1.d0
        jpi = 1.d0
        init_xmag = 0.d0 ! default to 0-flux zflux

        nsweep = 100
        nwarmup = 100
        nbin = 20

        ! default hmc parameters
        ! set default values
        dt = 0.001d0      ! stepsize
        pm = 10.d0
        mdstep = 20
        mdtime = 0.5
        nsamples = 1
        itermax = 1500
        splitting_M = 12 ! md splitting_M value
        !errate = 0.0001d0 ! too large
        errate = 0.000001d0

        ! hmc control
        ltunedt = .false.
        luseinputdt = .false.
        lfourier = .false.

        ! read parameters from ftdqmc.in input file
        if ( irank.eq.0 ) then
            exists = .false.
            inquire (file = 'fthmc.in', exist = exists)
            if ( exists .eqv. .true. ) then
                call p_create()
                call p_parse('fthmc.in')
                call p_get( 'L'        , l       )
                call p_get( 'Nflavor'  , Nflavor )
                call p_get( 'beta'     , beta    )
                call p_get( 'dtau'     , dtau    )
                call p_get( 'mu'       , mu      )
                call p_get( 'muA'      , muA     )
                call p_get( 'muB'      , muB     )
                call p_get( 'xmag'     , xmag    )
                call p_get( 'flux_x'   , flux_x  )
                call p_get( 'flux_y'   , flux_y  )
                call p_get( 'nsweep'   , nsweep  )
                call p_get( 'nwarmup'  , nwarmup )
                call p_get( 'nbin'     , nbin    )
                call p_get( 'ltau'     , ltau    )
                call p_get( 'nuse'     , nuse    )
                call p_get( 'rt'       , rt      )
                call p_get( 'js'       , js      )
                call p_get( 'jpi'      , jpi     )
                call p_get( 'init_xmag', init_xmag     )
                call p_get( 'dt'       , dt      )
                call p_get( 'pm'       , pm      )
                call p_get( 'mdtime'   , mdtime  )
                call p_get( 'mdstep'   , mdstep  )
                call p_get( 'nsamples' , nsamples)
                call p_get( 'itermax'  , itermax )
                call p_get( 'errate'   , errate  )
                call p_get( 'ltunedt'  , ltunedt    )
                call p_get( 'lfourier' , lfourier
                call p_get( 'luseinputdt'  , luseinputdt    )
                call p_destroy()
            end if
        end if
        call mp_bcast( l,    0 )
        call mp_bcast( Nflavor,0 )
        call mp_bcast( beta, 0 )
        call mp_bcast( dtau, 0 )
        call mp_bcast( mu,   0 )
        call mp_bcast( muA,  0 )
        call mp_bcast( muB,  0 )
        call mp_bcast( xmag, 0 )
        call mp_bcast( flux_x, 0 )
        call mp_bcast( flux_y, 0 )
        call mp_bcast( nsweep, 0 )
        call mp_bcast( nwarmup, 0 )
        call mp_bcast( nbin, 0 )
        call mp_bcast( ltau, 0 )
        call mp_bcast( nuse, 0 )
        call mp_bcast( rt, 0 )
        call mp_bcast( js, 0 )
        call mp_bcast( jpi, 0 )
        call mp_bcast( init_xmag, 0 )
        call mp_bcast( dt,   0 )
        call mp_bcast( pm,   0 )
        call mp_bcast( mdtime,   0 )
        call mp_bcast( mdstep,   0 )
        call mp_bcast( nsamples,   0 )
        call mp_bcast( itermax,   0 )
        call mp_bcast( errate,   0 )
        call mp_bcast( ltunedt, 0 )
        call mp_bcast( lfourier, 0 )
        call mp_bcast( luseinputdt, 0 )
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        ! tune parameters
        lq = l*l
        lfam = max(lq/2,1)
        nfam = 4       ! checkerboard decomposition
        ndim = lq*norb ! dimension of matrix of determinant
        ltrot = nint( beta / dtau )

        ! mass re-scaling
        pm = pm * dtau
        inputdt = dt

#ifdef CUDA
        ! set mu for gpu, directly use this value, as constant
        mu_complex = dcmplx(exp(-dtau*(-mu/2.d0)), 0.d0)
#endif

        ! warmup parameters
        !converge_ratio = 0.1d0
        converge_ratio = 0.05d0 ! use smaller value
        !mdstep = nint(mdtime/dt)
        !mdstep = 20

        ! symmetrized checkerboard to ensure the Hermiticity of B matrices
        ! if dqmc, remove the /2.0 factor
        zeig_rsigl_k(1) = dcmplx( dexp( rt*dtau/2.d0), 0.d0 )
        zeig_rsigl_k(2) = cone / zeig_rsigl_k(1)
        sinhrtdtau = dcmplx( dsinh(rt*dtau/2.d0), 0.d0 )
        coshrtdtau = dcmplx( dcosh(rt*dtau/2.d0), 0.d0 )

        if ( Nflavor .ne. 0.d0 ) then
            js = js * Nflavor * 0.5d0
            jpi = jpi * Nflavor * 0.5d0
        endif

        ! set constants for dimer-dimer correlation function
        z2 = dcmplx( Nflavor*Nflavor-1.d0, 0.d0 )
        z4 = z2*z2
        z3 = dcmplx( Nflavor*Nflavor*Nflavor - 2.d0*Nflavor + 1.d0/Nflavor, 0.d0 )
        z1 = dcmplx( -Nflavor + 1.d0/Nflavor, 0.d0 )

        ! allocate matrices
        allocate( Imat(ndim,ndim) )
        call s_identity_z(ndim,Imat)
        allocate( Ivec(ndim) )
        do i = 1, ndim
            Ivec(i) = 1.d0
        end do

        ! u1 model expvi
        allocate( zep_rsigl_k(lfam, nfam, ltrot) )
        allocate( zem_rsigl_k(lfam, nfam, ltrot) )

        ! force related
        allocate(hybrid_force1(lfam, nfam, ltrot))
        allocate(hybrid_force2(lfam, nfam, ltrot))

        allocate(pfield(lfam, nfam, ltrot))
        !allocate(xfield(lfam, nfam, ltrot))
        !allocate(xfield_tmp(lfam, nfam, ltrot))
        allocate(hybrid_force(lfam, nfam, ltrot))
        allocate(hybrid_force_n(lfam, nfam, ltrot))
        allocate(Btau_vtau(ndim, 16, ltrot, int(Nflavor/2.d0)))
        allocate(Btau2_vtau(ndim, 8, ltrot, int(Nflavor/2.d0)))
!#ifdef DIRINV
        allocate(M_inv(ndim*ltrot, ndim*ltrot))
        allocate(Bmat_up(ndim, ndim, ltrot))
        allocate(Bmat_up_tmp(ndim, ndim, ltrot))
!#endif
        if ( lfourier) then
            allocate(hybrid_force_k_z(lfam, nfam, ltrot))
            allocate(hybrid_force_n_k_z(lfam, nfam, ltrot))
            allocate(a_k(ltrot))

            !eta =0.005d0
            eta =0.05d0
            r2 =-4.d0+eta
            do i = 0, ltrot-1
                r1 = -1.d0 / (2.d0-2.d0*dcos(2.d0*pi*dble(i)/dble(ltrot)) + eta)
                a_k(i+1) = dcmplx( dsqrt(r2*r1), 0.d0 )
            enddo
        endif


    endsubroutine fthmc_hamilt_init


    subroutine fthmc_hamilt_free
        deallocate( Ivec, Imat )
        deallocate( zep_rsigl_k, zem_rsigl_k)

        deallocate(pfield)
        !deallocate(xfield)
        !deallocate(xfield_tmp)
        deallocate(hybrid_force)
        deallocate(hybrid_force_n)
        deallocate(Btau_vtau)
        deallocate(Btau2_vtau)
#IFDEF DIRINV
        deallocate(M_inv)
        deallocate(Bmat_up)
        deallocate(Bmat_up_tmp)
#ENDIF
        deallocate(hybrid_force1, hybrid_force2)
        if (lfourier) then
            deallocate(hybrid_force_k_z, hybrid_force_n_k_z, a_k)
        endif
    endsubroutine fthmc_hamilt_free



    subroutine fthmc_initial_head
        use spring
        use ftdqmc_hamilt
        use fthmc_hamilt
#ifdef _OPENMP
            use omp_lib
#endif

        integer :: system_time
        integer :: stream_seed
        integer :: nthreads, myid
        integer :: omp_get_num_threads, omp_get_thread_num
        character (len = 20) :: date_time_string

        !================================================
        !%% inital the pseudo random number generator   $
        !------------------------------------------------
        if ( .not. lfixseed) then
            call system_clock(system_time)
            stream_seed = abs( system_time - ( irank * 1981 + 2008 ) * 951049 )
        else
            stream_seed = 1371612732 !!! same seed for all processes
        endif

        call spring_sfmt_init(stream_seed)
        call s_time_builder(date_time_string)

        ! print head
        if(irank.eq.0) then
        write(fout,'(a)') ' ===================================================================================='
        write(fout,*)
        write(fout,'(a)') '        The finite temperature hybrid monte carlo (HMC) package '
        write(fout,*)
        write(fout,'(a)') '            FFFF   TTTTT   H    H     M   M      CCCC                    '
        write(fout,'(a)') '            F        T     H    H    M M M M    C                        '
        write(fout,'(a)') '            FFFF     T     HHHHHH    M M M M    C                        '
        write(fout,'(a)') '            F        T     H    H    M M M M    C                        '
        write(fout,'(a)') '            F        T     H    H   M   M   M    CCCC                    '
        write(fout,*)
        write(fout,*)
        write(fout,'(a)') ' written by Chuang Chen ( chenchuang@iphy.ac.cn )                                '
        write(fout,*)
        write(fout,'(a)') ' history: '
        write(fout,*)
        write(fout,'(a)') '     20/06/2018,  version 1.0  '
        write(fout,*)
        write(fout,'(a)') ' ------------------------------------------------------------------------------------'
        write(fout,*)
        write(fout,'(a)') ' >>> The simulation start running at '//date_time_string
        if( isize .gt. 1 ) then
        write(fout,'(a,i6,a)') ' >>> Parallelism running with', isize, '  processes'
    else
        write(fout,'(a)') ' >>> Serial running '
    end if

        ! output openmp info
#ifdef _OPENMP
!$OMP PARALLEL PRIVATE(myid) &
!$OMP shared(nthreads)
        nthreads = omp_get_num_threads()
        myid = omp_get_thread_num()
!$OMP BARRIER
        if ( myid .eq. 0 ) then
        write(fout,'(a,i6,a)') ' >>> OpenMP running with', nthreads, '  threads'
    endif
!$OMP END PARALLEL
#endif
        end if
    end subroutine fthmc_initial_head

    subroutine fthmc_initial_print(latt)
        use ftdqmc_hamilt
        use ftdqmc_latt_class
        implicit none
        class(ftdqmc_latt), intent(in) :: latt

        ! local variables
        integer :: i, j, iq, imj
        real(dp) :: qvec(2)

        if(irank.eq.0) then
            write(fout,*)
            write(fout,'(a)')' --------------------- '
            write(fout,'(a)')' The input parameters  '
            write(fout,'(a)')' --------------------- '
            write(fout,*)
            write(fout,'(a,f6.2)')    ' t      = ', rt
            write(fout,'(a,f6.2)')    ' js     = ', js
            write(fout,'(a,f6.2)')    ' jpi    = ', jpi
            write(fout,'(a,f6.2)')    ' init_xmag    = ', init_xmag
            write(fout,'(a,f7.3)')    ' mu     = ', mu
            write(fout,'(a,f6.2)')    ' B      = ', xmag
            write(fout,'(a,f6.2)')    ' dimer  = ', dimer
            write(fout,'(a,f8.5)')    ' flux_x = ', flux_x
            write(fout,'(a,f8.5)')    ' flux_y = ', flux_y
            write(fout,'(a,i4)')      ' L      = ', l
            write(fout,'(a,i4)')      ' LQ     = ', lq
            write(fout,'(a,f6.2)')    ' Nflavor= ', Nflavor
            write(fout,'(a,f6.2)')    ' beta   = ', beta
            write(fout,'(a,f7.3)')    ' dtau   = ', dtau
            write(fout,'(a,i6)')      ' ltrot  = ', ltrot
            write(fout,'(a,i6)')      ' nwarmup = ', nwarmup
            write(fout,'(a,i6)')      ' nsweep = ', nsweep
            write(fout,'(a,i6)')      ' nbin   = ', nbin
            write(fout,'(a,f7.3)')    ' dt     = ', dt
            write(fout,'(a,i6)')      ' mdstep = ', mdstep
            write(fout,'(a,f7.3)')    ' pm     = ', pm
            write(fout,'(a,i6)')      ' nsamples = ', nsamples
            write(fout,'(a,i6)')      ' itermax= ', itermax
            write(fout,'(a,f7.3)')    ' errate = ', errate
            write(fout,*)  'ltau = ', ltau
            write(fout,*)  'ltunedt = ', ltunedt
            write(fout,*)  'lfourier = ', lfourier
            write(fout,*)  'luseinputdt = ', luseinputdt
            write(fout,*)

            ! output the momentum vector list
            write(fout,*)
            write(fout,'(a)')' --------------------- '
            write(fout,'(a)')' momentum vector list  '
            write(fout,'(a)')' --------------------- '
            write(fout,'(a)')'       listk(i, 1), listk(i, 2), i, kx, ky'
            do iq = 1, lq
                qvec = dble(latt%listk(iq,1))*latt%b1_p + dble(latt%listk(iq,2))*latt%b2_p
                write(fout, '(3i8, 2f16.8)') latt%listk(iq,1), latt%listk(iq,2), iq, qvec(:)
            end do

            ! simulation banner
            write(fout,*)
            write(fout,'(a)')' --------------------- '
            write(fout,'(a)')' simulation started '
            write(fout,'(a)')' --------------------- '
        endif

    end subroutine fthmc_initial_print

end module ftdqmc_io
