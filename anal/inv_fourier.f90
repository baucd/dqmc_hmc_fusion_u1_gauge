program inv_fourier
    use mpi
    use ftdqmc_hamilt
    use ftdqmc_latt
    use ftdqmc_io

    implicit none

    ! local
    integer :: i, j, nbins, imj, nk, nb
    real(dp) :: aimj_p(2)
    complex(dp) :: c1
    character(len=64) :: fname, outname

    ! array
    real(dp), ALLOCATABLE, DIMENSION(:) :: kx, ky
    complex(dp), ALLOCATABLE, DIMENSION(:) :: obs, gr

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,irank,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,isize,ierr)

    ! set lattice
    !call ftdqmc_initial
    !call make_tables
    !call sli

    ! hamilt init
    call ftdqmc_hamilt_init
    ! latt
    call ftdqmc_latt_alloc
    call ftdqmc_latt_sli

    !call get_command_argument(1, fname)
    fname = "sq_spzz.bin"
    outname = "szz_r.bin"

    ! count the number of bins
    open(unit=5, file=fname, status='old')
    nbins = 0
    do
        do j =1, lq
            read(5,*,iostat=i) c1
        enddo
        if(i /= 0) exit
        nbins = nbins + 1
    end do


    ! one bin by one bin
    allocate(obs(lq))
    allocate(kx(lq))
    allocate(ky(lq))
    allocate(gr(lq))

    open(unit=10,file=outname,status='unknown', action="write")
    rewind(5)
    do nb = 1, nbins
        ! load data
        do j = 1, lq
            read(5,*) kx(j), ky(j), obs(j)
        enddo

        ! perform inv fourier transform
        gr = dcmplx(0.d0,0.d0)
        ! loop over real space
        do imj = 1,lq
           aimj_p = dble(list_prim(imj,1)*a1_p) + dble(list_prim(imj,2)*a2_p)
           ! loop over k space
           do nk = 1,lq
              gr(imj) = gr(imj) +  obs(nk)*zexpiqr(imj,nk)
           enddo
        enddo
        gr = gr/dcmplx(dble(lq),0.d0)

        ! output real space data
        do imj = 1,lq
           aimj_p = dble(list_prim(imj,1)*a1_p) + dble(list_prim(imj,2)*a2_p)
           !! convert ffa convention to correct convention
           write(10,'(2f8.3, 2e16.8)') aimj_p(1), aimj_p(2), gr(imj)
        enddo

    end do
    close(5)
    close(10)

    ! free memory
    deallocate(obs)
    deallocate(kx)
    deallocate(ky)
    deallocate(gr)

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_FINALIZE(ierr)

end program inv_fourier
