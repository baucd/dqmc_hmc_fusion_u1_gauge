module ftdqmc_io
      use ftdqmc_hamilt

contains

    subroutine ftdqmc_hamilt_init
        use mpi
        use parser, only : p_create, p_parse, p_get, p_get_vec, p_destroy
        use mmpi, only : mp_bcast, mp_barrier
        implicit none

        !include 'mpif.h'

        integer :: i
        logical :: exists

        ! default parameters
        l      = 4
        beta   = 20
        dtau   = 0.1d0
        mu     = 0.d0 ! default is half filling
        muA    = 0.d0
        muB    = 0.d0
        nwrap  = 10
        n_outconf_pace = 1

        ! model parameters
        ! for su2 gauge field
        rk = 1.d0
        rktau = 1.d0
        rj = 1.d0
        Nflavor = 2.d0
        infktau = .false.
        vari = 0.2d0

        ! for chargon field
        rr = 1.d0
        rw = 1.d0
        ru = 1.d0
        rg = 1.d0
        rv1 = 1.d0
        rj1 = 1.d0
        rk1 = 1.d0
        rv11 = 1.d0
        rv22 = 1.d0
        rm = 1.d0
        vari_cg = 0.2d0

        ! model control
        lclassic = .false.
        lfreezeU = .false.
        lnoupdate = .false.
        loutput = .false.
        lrmchargon = .false.
        lrmfermion = .false.
        loverrelaxation = .false.

        ltau = .false.
        nuse = 0

        xmag   = 0.d0
        flux_x = 0.d0
        flux_y = 0.d0

        nsweep  = 20
        nwarmup = 10
        nbin    = 10
        obs_segment_len = 10

        num_st_nn = 6
        nsw_stglobal = -1
        nblock = 4 ! usual value
        icount_nsw_stglobal = 0
        lstglobal = .false.
        llocal = .true.
        lfixseed = .false.

        nublock = 16

        ! get parameters from file: ftdqmc.in
        if ( irank.eq.0 ) then
            exists = .false.
            inquire (file = 'ftdqmc.in', exist = exists)
            if ( exists .eqv. .true. ) then
                call p_create()
                call p_parse('ftdqmc.in')
                call p_get( 'L'        , l       )            ! 1
                call p_get( 'beta'     , beta    )            ! 2
                call p_get( 'dtau'     , dtau    )            ! 3
                call p_get( 'rk'       , rk    )              ! 3
                call p_get( 'rktau'    , rktau    )           ! 3
                call p_get( 'rj'       , rj    )              ! 3
                call p_get( 'Nflavor'  , Nflavor  )           ! 3
                call p_get( 'vari'     , vari    )            ! 3
                call p_get( 'mu'       , mu      )            ! 3.5
                call p_get( 'muA'      , muA     )            ! 3.5
                call p_get( 'muB'      , muB     )            ! 3.5
                call p_get( 'xmag'     , xmag    )            ! 7
                call p_get( 'flux_x'   , flux_x  )            ! 7
                call p_get( 'flux_y'   , flux_y  )            ! 7
                call p_get( 'rr'       , rr    )
                call p_get( 'rw'       , rw    )
                call p_get( 'ru'       , ru    )
                call p_get( 'rg'       , rg    )
                call p_get( 'rv1'      , rv1    )
                call p_get( 'rj1'      , rj1    )
                call p_get( 'rk1'      , rk1    )
                call p_get( 'rv11'     , rv11    )
                call p_get( 'rv22'     , rv22    )
                call p_get( 'rm'       , rm    )
                call p_get( 'nwrap'    , nwrap   )            ! 8
                call p_get( 'nsweep'   , nsweep  )            ! 9
                call p_get( 'nwarmup'  , nwarmup )            ! 9
                call p_get( 'nbin'     , nbin    )            ! 10
                call p_get( 'llocal'   , llocal  )            ! 11
                call p_get( 'nsw_stglobal', nsw_stglobal )    ! 11
                call p_get( 'nblock', nblock )    ! 11
                call p_get( 'ltau'     , ltau    )
                call p_get( 'nuse'     , nuse    )
                call p_get( 'infktau'     , infktau    )
                call p_get( 'lclassic'     , lclassic    )
                call p_get( 'lfreezeU'     , lfreezeU    )
                call p_get( 'lnoupdate'     , lnoupdate    )
                call p_get( 'loutput'     , loutput    )
                call p_get( 'lrmchargon'     , lrmchargon    )
                call p_get( 'lrmfermion'     , lrmfermion    )
                call p_get( 'loverrelaxation'     , loverrelaxation    )
                call p_get( 'lfixseed'     , lfixseed    )
                call p_destroy()
            end if
        end if

        ! broadcast to other mpi nodes
        call mp_bcast( l,    0 )
        call mp_bcast( beta, 0 )
        call mp_bcast( dtau, 0 )
	    call mp_bcast( rk,  0 )
	    call mp_bcast( rktau,  0 )
	    call mp_bcast( rj, 0)
	    call mp_bcast( Nflavor, 0)
	    call mp_bcast( vari, 0)
        call mp_bcast( mu,   0 )
        call mp_bcast( muA,  0 )
        call mp_bcast( muB,  0 )
        call mp_bcast( xmag, 0 )
        call mp_bcast( flux_x, 0 )
        call mp_bcast( flux_y, 0 )
        call mp_bcast( rr, 0 )
        call mp_bcast( rw, 0 )
        call mp_bcast( ru, 0 )
        call mp_bcast( rg, 0 )
        call mp_bcast( rv1, 0 )
        call mp_bcast( rj1, 0 )
        call mp_bcast( rk1, 0 )
        call mp_bcast( rv11, 0 )
        call mp_bcast( rv22, 0 )
        call mp_bcast( rm, 0 )
        call mp_bcast( nwrap, 0 )
        call mp_bcast( nsweep, 0 )
        call mp_bcast( nwarmup, 0 )
        call mp_bcast( nbin, 0 )
        call mp_bcast( llocal, 0 )
        call mp_bcast( nsw_stglobal, 0 )
        call mp_bcast( nblock, 0 )
        call mp_bcast( ltau, 0 )
        call mp_bcast( nuse, 0 )
        call mp_bcast( nublock, 0 )
        call mp_bcast( infktau, 0 )
        call mp_bcast( lclassic, 0 )
        call mp_bcast( lfreezeU, 0 )
        call mp_bcast( lnoupdate, 0 )
        call mp_bcast( loutput, 0 )
        call mp_bcast( lrmchargon, 0 )
        call mp_bcast( lrmfermion, 0 )
        call mp_bcast( loverrelaxation, 0 )
        call mp_bcast( lfixseed, 0 )
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        ! tune parameters
        if (lclassic) nsw_stglobal = 0
        if (lclassic) loverrelaxation = .false.
        if (lclassic) nblock = 0
        if( nsw_stglobal .gt. 0 ) lstglobal = .true.
        if( .not. llocal .and. lstglobal ) nsw_stglobal = 1
        lq = l*l
        ndim = lq * norb
        nfam = 4
        lfam = max(lq/2, 1)
        ltrot = nint( beta / dtau )
        if (lclassic) ltrot = 1
        ! force lrmfermion to be true when simulating the classical model
        if (lclassic) lrmfermion = .true.
        ! disable dyn corr meas when no fermion
        if (lrmfermion) ltau = .false.
        expmmu = dcmplx( dexp(-mu), 0.d0) ! exp(-mu) the chemical potential

        ! tune rk and rktau according to Nflavor
        rk = rk * Nflavor/2.d0
        rktau = rktau * Nflavor/2.d0

        ! set sinhrjdtau and coshrjdtau
        sinhrjdtau = dcmplx( dsinh(rj*dtau), 0.d0)
        coshrjdtau = dcmplx( dcosh(rj*dtau), 0.d0)

        ! set zeig for SU2 gauge field, e^{-\Delta \tau J}
        zeig(1) = dcmplx(exp(-dtau *  rj) , 0.d0)
        zeig(2) = dcmplx(exp(-dtau * (-rj)) , 0.d0)
        zeig(3) = dcmplx(exp(-dtau *  rj) , 0.d0)
        zeig(4) = dcmplx(exp(-dtau * (-rj)) , 0.d0)

        ! tune para for delay update
        if( ndim/5 .lt. 16) then
            nublock = 4
        else if( ndim/5 .lt. 32 ) then
            nublock = 8
        else if( ndim/5 .lt. 64 ) then
            nublock = 16
        else if( ndim/5 .lt. 256 ) then
            nublock = 32
        else ! equal to or greater than 256
            nublock = 64
        end if

        ! numerical stablization
        allocate( iwrap_nt(0:ltrot) )
        iwrap_nt(0:ltrot) = 0
        ! set nst, and wrap_step
        if( ltrot .lt. nwrap ) then
            if(irank.eq.0) write(fout,'(a,i3,a)')  " WARNNING, ltrot is less than nwrap ", ltrot, ', do not need stablization '
            nst = 0
        else
            if( mod(ltrot,nwrap) .eq. 0 ) then
                nst = ltrot/nwrap
                allocate( wrap_step(2,nst) )
                do i = 1, nst
                    iwrap_nt(i*nwrap) = i
                    wrap_step(1,i) = 1+(i-1)*nwrap
                    wrap_step(2,i) = i*nwrap
                end do
            else
                nst = ltrot/nwrap + 1
                allocate( wrap_step(2,nst) )
                do i = 1, nst-1
                    iwrap_nt(i*nwrap) = i
                    wrap_step(1,i) = 1+(i-1)*nwrap
                    wrap_step(2,i) = i*nwrap
                end do
                i = nst
                iwrap_nt(ltrot) = i
                wrap_step(1,i) = (i-1)*nwrap+1
                wrap_step(2,i) = ltrot
            end if
        end if

        nmeas_bin = 2*(2*obs_segment_len+1)*nsweep*isize
        weight_track = czero

        ! flipping table
        nflipl(-1, 1) =  1
        nflipl( 1, 1) = -1

        if ( .not. lclassic) then
            allocate( Imat(ndim,ndim) )
            call s_identity_z(ndim,Imat)
            allocate( Ivec(ndim) )
            do i = 1, ndim
                Ivec(i) = 1.d0
            end do
        endif
        ! for SU(2) matrix
        allocate( I2mat(2,2) )
        call s_identity_z(2,I2mat)

        if( lstglobal ) then
            allocate( tentacle(2,lq*ltrot), tentacle_old(2,lq*ltrot) )
            allocate( stcluster(lq,ltrot) )
            allocate( stbonds_neib(2,6,lq,ltrot) )
        end if

        ! set lnoupdate
        lnoupdate = loutput
    endsubroutine ftdqmc_hamilt_init


    subroutine ftdqmc_hamilt_free
        if( lstglobal ) then
            deallocate( stbonds_neib )
            deallocate( stcluster )
            deallocate( tentacle_old, tentacle )
        end if
        if (.not. lclassic) then
            deallocate( Ivec, Imat )
        endif
        if(allocated(wrap_step) ) deallocate( wrap_step )
        deallocate( iwrap_nt )
    endsubroutine ftdqmc_hamilt_free


    subroutine ftdqmc_initial_head
        use spring
        use ftdqmc_hamilt
#ifdef _OPENMP
        use omp_lib
#endif

        integer :: system_time
        integer :: stream_seed
        integer :: nthreads, myid
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
            write(fout,'(a)') '        The finite temperature determinant quantum monte carlo (DQMC) package '
            write(fout,*)
            write(fout,'(a)') '            FFFF   TTTTT   DDD      QQQ     M   M      CCCC                    '
            write(fout,'(a)') '            F        T     D  D    Q   Q   M M M M    C                        '
            write(fout,'(a)') '            FFFF     T     D   D   Q   Q   M M M M    C                        '
            write(fout,'(a)') '            F        T     D  D    Q   Q   M M M M    C                        '
            write(fout,'(a)') '            F        T     DDD      QQQ   M   M   M    CCCC                    '
            write(fout,'(a)') '                                      \                                       '
            write(fout,*)
            write(fout,*)
            write(fout,'(a)') ' written by Chuang Chen ( chenchuang2020@icloud.com )                                '
            write(fout,*)
            write(fout,'(a)') ' history: '
            write(fout,*)
            write(fout,'(a)') '     21/05/2018,  version 2.0  '
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
!$OMP PARALLEL PRIVATE(nthreads)
              nthreads = omp_get_num_threads()
              myid = omp_get_thread_num()
              if ( myid .eq. 0 ) then
                  write(fout,'(a,i6,a)') ' >>> OpenMP running with', nthreads, '  threads'
              endif
!$OMP END PARALLEL
#endif
        end if
    end subroutine ftdqmc_initial_head


    subroutine ftdqmc_initial_print
        use ftdqmc_hamilt
        use ftdqmc_latt
        implicit none

        ! local variables
        integer :: i, j, iq
        real(dp) :: qvec(2)

        if(irank.eq.0) THEN
            write(fout,*)
            write(fout,'(a)')' --------------------- '
            write(fout,'(a)')' The input parameters  '
            write(fout,'(a)')' --------------------- '
            write(fout,*)
            write(fout,'(a,i4)')      ' l       = ', l
            write(fout,'(a,i4)')      ' lq      = ', lq
            write(fout,'(a,i4)')      ' norb    = ', norb
            write(fout,'(a,i4)')      ' ndim    = ', ndim
            write(fout,'(a,f6.2)')    ' beta    = ', beta
            write(fout,'(a,f6.2)')    ' dtau    = ', dtau
            write(fout,'(a,i6)')      ' ltrot   = ', ltrot
            write(fout,'(a,f6.2)')    ' kappa   = ', rk
            write(fout,'(a,f6.2)')    ' kappa_tau   = ', rktau
            write(fout,'(a,f6.2)')    ' J       = ', rj
            write(fout,'(a,f6.2)')    ' Nflavor = ', Nflavor
            write(fout,'(a,f6.2)')    ' vari    = ', vari
            write(fout,'(a,f6.2)')    ' mu      = ', mu
            write(fout,'(a,f8.5)')    ' flux_x  = ', flux_x
            write(fout,'(a,f8.5)')    ' flux_y  = ', flux_y
            write(fout,'(a,f6.2)')    ' rr  = ', rr
            write(fout,'(a,f6.2)')    ' rw  = ', rw
            write(fout,'(a,f6.2)')    ' ru  = ', ru
            write(fout,'(a,f6.2)')    ' rg  = ', rg
            write(fout,'(a,f6.2)')    ' rv1  = ', rv1
            write(fout,'(a,f6.2)')    ' rj1  = ', rj1
            write(fout,'(a,f6.2)')    ' rk1  = ', rk1
            write(fout,'(a,f6.2)')    ' rv11  = ', rv11
            write(fout,'(a,f6.2)')    ' rv22  = ', rv22
            write(fout,'(a,f8.2)')    ' rm  = ', rm
            write(fout,'(a,i6)')      ' nwarmup = ', nwarmup
            write(fout,'(a,i6)')      ' nsweep  = ', nsweep
            write(fout,'(a,i6)')      ' nbin    = ', nbin
            write(fout,'(a,i6)')      ' nwrap   = ', nwrap
            write(fout,'(a,i6)')      ' nst     = ', nst
            write(fout,'(a,i6)')      ' nsw_stglobal = ', nsw_stglobal
            write(fout,'(a,i6)')      ' nblock = ', nblock
#ifdef DELEY
            write(fout,'(a,i6)')      ' nublock = ', nublock
#endif
            write(fout,'(a,l1)') ' llocal = ', llocal
            write(fout,'(a,l1)') ' lstglobal = ', lstglobal
            write(fout,'(a,l1)') ' ltau = ', ltau
            write(fout,'(a,l1)') ' infktau = ', infktau
            write(fout,'(a,l1)') ' lclassic = ', lclassic
            write(fout,'(a,l1)') ' lfreezeU = ', lfreezeU
            write(fout,'(a,l1)') ' lrmchargon = ', lrmchargon
            write(fout,'(a,l1)') ' lrmfermion = ', lrmfermion
            write(fout,'(a,l1)') ' loverrelaxation = ', loverrelaxation
            write(fout,'(a,l1)') ' lfixseed = ', lfixseed
            write(fout,'(a,i4)')      ' lq      = ', lq
            write(fout,'(a,i6)') ' number of meas for 1 bin = ', nmeas_bin
            write(fout,'(a,i6)') ' total number of measurements = ', nbin*nmeas_bin

            write(fout,*)
            write(fout,'(a)')' --------------------- '
            write(fout,'(a)')' wrapping coordinates  '
            write(fout,'(a)')' --------------------- '
            write(fout,'(a)')'       wrap_step(1,i)   wrap_step(2,i)   iwrap_nt(nt) '
            do i = 1, nst
                write( fout, '(3i16)') wrap_step(1,i), wrap_step(2,i), iwrap_nt( wrap_step(2,i) )
            end do

            ! output the momentum vector list
            write(fout,*)
            write(fout,'(a)')' --------------------- '
            write(fout,'(a)')' momentum vector list  '
            write(fout,'(a)')' --------------------- '
            write(fout,'(a)')'       listk(i, 1), listk(i, 2), i, kx, ky'
            do iq = 1, lq
                qvec = dble(listk(iq,1))*b1_p + dble(listk(iq,2))*b2_p
                write(fout, '(3i8, 2f16.8)') listk(iq,1), listk(iq,2), iq, qvec(:)
            end do

            write(fout,*)
            write(fout,'(a)')' --------------------- '
            write(fout,'(a)')' simulation started '
            write(fout,'(a)')' --------------------- '

        END if
    end subroutine ftdqmc_initial_print

end module ftdqmc_io
