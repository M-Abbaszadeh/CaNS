!
!        CCCCCCCCCCCCC                    NNNNNNNN        NNNNNNNN    SSSSSSSSSSSSSSS
!     CCC::::::::::::C                    N:::::::N       N::::::N  SS:::::::::::::::S
!   CC:::::::::::::::C                    N::::::::N      N::::::N S:::::SSSSSS::::::S
!  C:::::CCCCCCCC::::C                    N:::::::::N     N::::::N S:::::S     SSSSSSS
! C:::::C       CCCCCC   aaaaaaaaaaaaa    N::::::::::N    N::::::N S:::::S
!C:::::C                 a::::::::::::a   N:::::::::::N   N::::::N S:::::S
!C:::::C                 aaaaaaaaa:::::a  N:::::::N::::N  N::::::N  S::::SSSS
!C:::::C                          a::::a  N::::::N N::::N N::::::N   SS::::::SSSSS
!C:::::C                   aaaaaaa:::::a  N::::::N  N::::N:::::::N     SSS::::::::SS
!C:::::C                 aa::::::::::::a  N::::::N   N:::::::::::N        SSSSSS::::S
!C:::::C                a::::aaaa::::::a  N::::::N    N::::::::::N             S:::::S
! C:::::C       CCCCCC a::::a    a:::::a  N::::::N     N:::::::::N             S:::::S
!  C:::::CCCCCCCC::::C a::::a    a:::::a  N::::::N      N::::::::N SSSSSSS     S:::::S
!   CC:::::::::::::::C a:::::aaaa::::::a  N::::::N       N:::::::N S::::::SSSSSS:::::S
!     CCC::::::::::::C  a::::::::::aa:::a N::::::N        N::::::N S:::::::::::::::SS
!        CCCCCCCCCCCCC   aaaaaaaaaa  aaaa NNNNNNNN         NNNNNNN  SSSSSSSSSSSSSSS
!-------------------------------------------------------------------------------------
! CaNS -- Canonical Navier-Stokes Solver
! Pedro Costa (p.simoes.costa@gmail.com)
!-------------------------------------------------------------------------------------
program cans
  use iso_c_binding  , only: C_PTR
  use mpi
  use decomp_2d
  use mod_bound      , only: boundp,bounduvw,updt_rhs_b,updthalo
  use mod_chkdiv     , only: chkdiv
  use mod_chkdt      , only: chkdt
  use mod_common_mpi , only: myid,ierr,coord
  use mod_correc     , only: correc
  use mod_debug      , only: chkmean
  use mod_fft        , only: fftini,fftend
  use mod_fillps     , only: fillps
  use mod_initflow   , only: initflow
  use mod_initgrid   , only: initgrid
  use mod_initmpi    , only: initmpi
  use mod_initsolver , only: initsolver
  use mod_load       , only: load
  use mod_rk         , only: rk,rk_id
  use mod_output     , only: out0d,out1d,out1d_2,out2d,out3d,write_log_output,write_visu_2d,write_visu_3d
  use mod_param      , only: itot,jtot,ktot,lx,ly,lz,dx,dy,dz,dxi,dyi,dzi,uref,lref,rey,visc,small, &
                             cbcvel,bcvel,cbcpre,bcpre, &
                             icheck,iout0d,iout1d,iout2d,iout3d,isave, &
                             nstep,time_max,tw_max,stop_type,restart,is_overwrite_save, &
                             rkcoeff,   &
                             datadir,   &
                             cfl,dtmin, &
                             inivel,    &
                             imax,jmax,dims, &
                             nthreadsmax, &
                             gr, &
                             is_forced,bforce, &
                             n,ng,l,dl,dli, &
                             read_input
  use mod_post       , only: rotation_rate,strain_rate,q_criterion
  use mod_sanity     , only: test_sanity
  use mod_solver     , only: solver
  use mod_types
#ifdef USE_CUDA
  use mod_common_mpi, only: mydev
#endif
#ifdef USE_NVTX
  use nvtx
#endif 
#ifdef USE_CATALYST
  use catalyst_interfaces
#endif
  !$ use omp_lib
  implicit none
  real(rp), allocatable, dimension(:,:,:) :: u,v,w,p,up,vp,wp,pp
#ifdef USE_CUDA
  attributes(managed) :: u,v,w,p,up,vp,wp,pp
#endif
  real(rp), allocatable, dimension(:,:,:) :: str,ens,qcr
#ifdef USE_CUDA
  attributes(managed) :: str,ens,qcr
#endif
  real(rp), dimension(:,:,:), allocatable, target  :: dudtrk_A,dvdtrk_A,dwdtrk_A
  real(rp), dimension(:,:,:), allocatable, target  :: dudtrk_B,dvdtrk_B,dwdtrk_B
  real(rp), dimension(:,:,:), pointer :: dudtrko,dvdtrko,dwdtrko
  real(rp), dimension(:,:,:), pointer :: dudtrk,dvdtrk,dwdtrk, rk_tmp
  real(rp), dimension(3) :: tauxo,tauyo,tauzo
  real(rp), dimension(3) :: f
  type(C_PTR), dimension(2,2) :: arrplanp
  real(rp), allocatable, dimension(:,:) :: lambdaxyp
  real(rp), allocatable, dimension(:) :: ap,bp,cp
  real(rp) :: normfftp
  integer  :: i,j,k,im,ip,jm,jp,km,kp
  type rhs_bound
#ifdef USE_CUDA
    real(rp), allocatable, dimension(:,:,:), managed :: x,y,z
#else
    real(rp), allocatable, dimension(:,:,:) :: x,y,z
#endif
  end type rhs_bound 
  type(rhs_bound), allocatable :: rhsbp
#ifdef IMPDIFF
  type(C_PTR), dimension(2,2) :: arrplanu,arrplanv,arrplanw
  real(rp), dimension(:,:), allocatable :: lambdaxyu,lambdaxyv,lambdaxyw
  real(rp), dimension(:)  , allocatable :: au,av,aw,bu,bv,bw,bb,cu,cv,cw
  real(rp) :: normfftu,normfftv,normfftw
  real(rp) :: alpha,alphai
  type(rhs_bound), allocatable :: rhsbu,rhsbv,rhsbw
  real(rp), dimension(:,:,:), allocatable    :: dudtrkd,dvdtrkd,dwdtrkd
  #ifdef USE_CUDA
  attributes(managed) :: dudtrkd,dvdtrkd,dwdtrkd,lambdaxyu,lambdaxyv,lambdaxyw
  attributes(managed) :: au,av,aw,bu,bv,bw,bb,cu,cv,cw
  attributes(managed) :: rhsbu,rhsbv,rhsbw
  #endif
#endif
  real(rp) :: ristep
  real(rp) :: dt,dti,dtmax,time,dtrk,dtrki,divtot,divmax
  integer  :: irk,istep
  real(rp), allocatable, dimension(:) :: dzc,dzf,zc,zf,dzci,dzfi,dzflzi,dzclzi,zclzi
  real(rp), allocatable, dimension(:) :: xc,yc ! for catalyst
#ifdef USE_CUDA
  integer :: istat
  integer(kind=cuda_count_kind) :: freeMem,totMem
  integer(8) :: totEle
  attributes(managed) :: dzc,dzf,zc,zf,dzci,dzfi,dzflzi,dzclzi,zclzi,lambdaxyp,ap,bp,cp,rhsbp
  attributes(managed) :: dudtrko,dvdtrko,dwdtrko,dudtrk,dvdtrk,dwdtrk,dudtrk_A,dvdtrk_A,dwdtrk_A,dudtrk_B,dvdtrk_B,dwdtrk_B
  attributes(managed) :: xc,yc
#endif
  real(rp) :: meanvel,meanvelu,meanvelv,meanvelw
  real(rp), dimension(3) :: dpdl
  !real(rp), allocatable, dimension(:) :: var
  real(rp), dimension(10) :: var
#ifdef TIMING
  real(rp) :: dt12,dt12av,dt12min,dt12max
#endif
  real(rp) :: f1d,f2d,f3d
  real(rp) :: twi,tw
  character(len=7  ) :: fldnum
  character(len=100) :: filename
  integer :: kk
  logical :: is_done,kill
  integer :: rlen
#ifdef USE_CATALYST
  logical(c_bool) :: catalyst_active
  real(rp) :: t0, t1
#endif
  !
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
  !
  ! read parameter file
  !
  call read_input(myid)
  !
  ! initialize MPI/OpenMP
  !
  !$call omp_set_num_threads(nthreadsmax)
  call initmpi(ng,cbcpre)
  n  = (/imax,jmax,ktot/) ! now set in initmpi
  twi = MPI_WTIME()
  !
  ! allocate memory
  !
  dudtrko => dudtrk_A
  dvdtrko => dvdtrk_A
  dwdtrko => dwdtrk_A
  dudtrk  => dudtrk_B
  dvdtrk  => dudtrk_B
  dwdtrk  => dudtrk_B
  !
  allocate(u( 0:n(1)+1,0:n(2)+1,0:n(3)+1) , &
           v( 0:n(1)+1,0:n(2)+1,0:n(3)+1) , &
           w( 0:n(1)+1,0:n(2)+1,0:n(3)+1) , &
           p( 0:n(1)+1,0:n(2)+1,0:n(3)+1) , &
           up(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           vp(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           wp(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           pp(0:n(1)+1,0:n(2)+1,0:n(3)+1))
  allocate(dudtrko(n(1),n(2),n(3)), &
           dvdtrko(n(1),n(2),n(3)), &
           dwdtrko(n(1),n(2),n(3)), &
           dudtrk( n(1),n(2),n(3)), &
           dvdtrk( n(1),n(2),n(3)), &
           dwdtrk( n(1),n(2),n(3)))
  allocate(str(n(1),n(2),n(3)), &
           ens(n(1),n(2),n(3)), &
           qcr(0:n(1)+1,0:n(2)+1,0:n(3)+1))

#ifdef USE_CUDA
  allocate(lambdaxyp(ng(1)/dims(2),ng(2)/dims(1)))
#else
  allocate(lambdaxyp(n(1),n(2)))
#endif
  allocate(ap(n(3)),bp(n(3)),cp(n(3)))
  allocate(dzc(   0:n(3)+1), &
           dzf(   0:n(3)+1), &
           zc(    0:n(3)+1), &
           zf(    0:n(3)+1), &
           dzci(  0:n(3)+1), &
           dzfi(  0:n(3)+1), &
           dzclzi(0:n(3)+1), &
           dzflzi(0:n(3)+1), &
           zclzi( 0:n(3)+1))
  allocate(xc(0:n(1)+1),yc(0:n(2)+1))
  allocate(rhsbp)
  allocate(rhsbp%x(n(2),n(3),0:1), &
           rhsbp%y(n(1),n(3),0:1), &
           rhsbp%z(n(1),n(2),0:1))
#ifdef IMPDIFF
  allocate(dudtrkd(n(1),n(2),n(3)), &
           dvdtrkd(n(1),n(2),n(3)), & 
           dwdtrkd(n(1),n(2),n(3)))
#ifdef USE_CUDA
  allocate(lambdaxyu(ng(1)/dims(2),ng(2)/dims(1)), &
           lambdaxyv(ng(1)/dims(2),ng(2)/dims(1)), &
           lambdaxyw(ng(1)/dims(2),ng(2)/dims(1)))
#else
  allocate(lambdaxyu(n(1),n(2)), &
           lambdaxyv(n(1),n(2)), &
           lambdaxyw(n(1),n(2)))
#endif
  allocate(au(n(3)),bu(n(3)),cu(n(3)), &
           av(n(3)),bv(n(3)),cv(n(3)), &
           aw(n(3)),bw(n(3)),cw(n(3)), &
           bb(n(3)))
  allocate(rhsbu,rhsbv,rhsbw)
  allocate(rhsbu%x(n(2),n(3),0:1), &
           rhsbu%y(n(1),n(3),0:1), &
           rhsbu%z(n(1),n(2),0:1), &
           rhsbv%x(n(2),n(3),0:1), &
           rhsbv%y(n(1),n(3),0:1), &
           rhsbv%z(n(1),n(2),0:1), &
           rhsbw%x(n(2),n(3),0:1), &
           rhsbw%y(n(1),n(3),0:1), &
           rhsbw%z(n(1),n(2),0:1))
#endif

#ifdef USE_CUDA
  istat = cudaMemAdvise( u, size(u), cudaMemAdviseSetPreferredLocation, mydev )
  istat = cudaMemAdvise( v, size(v), cudaMemAdviseSetPreferredLocation, mydev )
  istat = cudaMemAdvise( w, size(w), cudaMemAdviseSetPreferredLocation, mydev )
  istat = cudaMemAdvise( p, size(p), cudaMemAdviseSetPreferredLocation, mydev )
  istat = cudaMemAdvise( up, size(up), cudaMemAdviseSetPreferredLocation, mydev )
  istat = cudaMemAdvise( vp, size(vp), cudaMemAdviseSetPreferredLocation, mydev )
  istat = cudaMemAdvise( wp, size(wp), cudaMemAdviseSetPreferredLocation, mydev )
  istat = cudaMemAdvise( pp, size(pp), cudaMemAdviseSetPreferredLocation, mydev )
  istat = cudaMemAdvise( dudtrko, size(dudtrko), cudaMemAdviseSetPreferredLocation, mydev )
  istat = cudaMemAdvise( dvdtrko, size(dvdtrko), cudaMemAdviseSetPreferredLocation, mydev )
  istat = cudaMemAdvise( dwdtrko, size(dwdtrko), cudaMemAdviseSetPreferredLocation, mydev )
  istat = cudaMemAdvise( dudtrk, size(dudtrk), cudaMemAdviseSetPreferredLocation, mydev )
  istat = cudaMemAdvise( dvdtrk, size(dvdtrk), cudaMemAdviseSetPreferredLocation, mydev )
  istat = cudaMemAdvise( dwdtrk, size(dwdtrk), cudaMemAdviseSetPreferredLocation, mydev )
  istat = cudaMemPrefetchAsync( u, size(u), cudaCpuDeviceId, 0 )
  istat = cudaMemPrefetchAsync( v, size(v), cudaCpuDeviceId, 0 )
  istat = cudaMemPrefetchAsync( w, size(w), cudaCpuDeviceId, 0 )
  istat = cudaMemPrefetchAsync( p, size(p), cudaCpuDeviceId, 0 )
  istat = cudaMemPrefetchAsync( up, size(up), mydev, 0 )
  istat = cudaMemPrefetchAsync( vp, size(vp), mydev, 0 )
  istat = cudaMemPrefetchAsync( wp, size(wp), mydev, 0 )
  istat = cudaMemPrefetchAsync( pp, size(pp), mydev, 0 )
  istat = cudaMemPrefetchAsync( dudtrko, size(dudtrko), mydev, 0 )
  istat = cudaMemPrefetchAsync( dvdtrko, size(dvdtrko), mydev, 0 )
  istat = cudaMemPrefetchAsync( dwdtrko, size(dwdtrko), mydev, 0 )
  istat = cudaMemPrefetchAsync( dudtrk, size(dudtrk), mydev, 0 )
  istat = cudaMemPrefetchAsync( dvdtrk, size(dvdtrk), mydev, 0 )
  istat = cudaMemPrefetchAsync( dwdtrk, size(dwdtrk), mydev, 0 )
  !
  istat = cudaMemAdvise( ens, size(ens), cudaMemAdviseSetPreferredLocation, mydev )
  istat = cudaMemAdvise( str, size(str), cudaMemAdviseSetPreferredLocation, mydev )
  istat = cudaMemAdvise( qcr, size(qcr), cudaMemAdviseSetPreferredLocation, mydev )
  istat = cudaMemPrefetchAsync( ens, size(ens), cudaCpuDeviceId, 0 )
  istat = cudaMemPrefetchAsync( str, size(str), cudaCpuDeviceId, 0 )
  istat = cudaMemPrefetchAsync( qcr, size(qcr), cudaCpuDeviceId, 0 )
  !
  istat = cudaMemAdvise( rhsbp%x, size(rhsbp%x), cudaMemAdviseSetPreferredLocation, mydev )
  istat = cudaMemAdvise( rhsbp%y, size(rhsbp%y), cudaMemAdviseSetPreferredLocation, mydev )
  istat = cudaMemAdvise( rhsbp%z, size(rhsbp%z), cudaMemAdviseSetPreferredLocation, mydev )
  istat = cudaMemPrefetchAsync( rhsbp%x, size(rhsbp%x), cudaCpuDeviceId, 0 )
  istat = cudaMemPrefetchAsync( rhsbp%y, size(rhsbp%y), cudaCpuDeviceId, 0 )
  istat = cudaMemPrefetchAsync( rhsbp%z, size(rhsbp%z), cudaCpuDeviceId, 0 )
  !
  istat = cudaMemAdvise( zc, size(zc), cudaMemAdviseSetReadMostly, 0 )
  istat = cudaMemAdvise( zf, size(zf), cudaMemAdviseSetReadMostly, 0 )
  istat = cudaMemAdvise( dzc, size(dzc), cudaMemAdviseSetReadMostly, 0 )
  istat = cudaMemAdvise( dzf, size(dzf), cudaMemAdviseSetReadMostly, 0 )
  istat = cudaMemAdvise( dzci, size(dzci), cudaMemAdviseSetReadMostly, 0 )
  istat = cudaMemAdvise( dzfi, size(dzfi), cudaMemAdviseSetReadMostly, 0 )
  istat = cudaMemAdvise( dzclzi, size(dzclzi), cudaMemAdviseSetReadMostly, 0 )
  istat = cudaMemAdvise( dzflzi, size(dzflzi), cudaMemAdviseSetReadMostly, 0 )
  istat = cudaMemAdvise( zclzi, size(zclzi), cudaMemAdviseSetReadMostly, 0 )
#endif
  if(myid.eq.0) print*, '*******************************'
  if(myid.eq.0) print*, '*** Beginning of simulation ***'
  if(myid.eq.0) print*, '*******************************'
  if(myid.eq.0) print*, ''
#ifdef USE_CUDA
  if(myid.eq.0) then 
    print*, ' GPU accelerated version, grid size:',n(1)*dims(1),n(2)*dims(2),n(3)
    #ifndef IMPDIFF
    totEle=(8*int((n(1)+2)*(n(2)+2)*(n(3)+2),8) + 6*int(n(1)*n(2)*n(3),8) + 6*(n(1)*n(2)+n(2)*n(3)+n(1)*n(3)) )
    #else
    totEle=(8*int((n(1)+2)*(n(2)+2)*(n(3)+2),8) + 9*int(n(1)*n(2)*n(3),8) +24*(n(1)*n(2)+n(2)*n(3)+n(1)*n(3)) )
    #endif
    totEle = totEle+ int((n(1)*dims(1)*n(2)*dims(2)/dims(1)*n(3)/dims(2)),8) + &
                     int((n(1)*n(2)*dims(2)*n(3)/dims(2)),8) + int(n(1)*n(2)*ng(3),8) + &
                     int((n(1)*dims(1)+2)* n(2)*dims(2)/dims(1)* n(3)/dims(2),8) +  &
                     int((n(2)*dims(2)+2)* n(1)* n(3)/dims(2) , 8 ) + &
                     int(n(2)*dims(2)* n(1)* n(3)/dims(2), 8 ) 
    print*, ' Estimated Memory usage (GB) per MPI task: ',totEle*sizeof(1._rp)/(1024.**3)
    ! missing fft plan
    istat=cudaMemGetInfo(freeMem,totMem)
    print *," Memory on GPU (GB)                      : ",totMem/(1024.**3)
    print *," "
  endif
#endif 
  !
#ifdef USE_NVTX
  call nvtxStartRange("initgrid", 2)
#endif
#ifdef USE_CUDA
  istat = cudaMemPrefetchAsync( zc, size(zc), cudaCpuDeviceId, 0 )
  istat = cudaMemPrefetchAsync( zf, size(zf), cudaCpuDeviceId, 0 )
  istat = cudaMemPrefetchAsync( dzc, size(dzc), cudaCpuDeviceId, 0 )
  istat = cudaMemPrefetchAsync( dzf, size(dzf), cudaCpuDeviceId, 0 )
#endif
  call initgrid(inivel,n(3),gr,lz,dzc,dzf,zc,zf)
  do kk=0,n(1)+1
    xc(kk) = (kk-0.5)*dl(1)
  enddo
  do kk=0,n(2)+1
    yc(kk) = (kk-0.5)*dl(2)
  enddo
#ifdef USE_CUDA
  istat = cudaMemPrefetchAsync(  zc, size( zc), mydev, 0 )
  istat = cudaMemPrefetchAsync(  zf, size( zf), mydev, 0 )
  istat = cudaMemPrefetchAsync( dzc, size(dzc), mydev, 0 )
  istat = cudaMemPrefetchAsync( dzf, size(dzf), mydev, 0 )
#endif
#ifdef USE_NVTX
  call nvtxEndRange
#endif
  if(myid.eq.0) then
    inquire(iolength=rlen) 1._rp
    open(99,file=trim(datadir)//'grid.bin',access='direct',recl=4*n(3)*rlen)
    write(99,rec=1) dzc(1:n(3)),dzf(1:n(3)),zc(1:n(3)),zf(1:n(3))
    close(99)
    open(99,file=trim(datadir)//'grid.out')
    do kk=0,ktot+1
      write(99,'(5E15.7)') 0.,zf(kk),zc(kk),dzf(kk),dzc(kk)
    enddo
    close(99)
    open(99,file=trim(datadir)//'geometry.out')
      write(99,*) ng(1),ng(2),ng(3) 
      write(99,*) l(1),l(2),l(3) 
    close(99)
  endif

#ifdef USE_CUDA
  istat = cudaMemPrefetchAsync( dzci, size(dzci), mydev, 0 )
  istat = cudaMemPrefetchAsync( dzfi, size(dzfi), mydev, 0 )
  istat = cudaMemPrefetchAsync( dzflzi, size(dzflzi), mydev, 0 )
  istat = cudaMemPrefetchAsync( dzclzi, size(dzclzi), mydev, 0 )
  istat = cudaMemPrefetchAsync( zclzi, size(zclzi), mydev, 0 )
  istat = cudaMemPrefetchAsync( zc, size(zc), mydev, 0 )
  !$cuf kernel do(1) <<<*,*>>>
#endif
  do k=0,ktot+1
     dzci(k)=1./dzc(k)
     dzfi(k)=1./dzf(k)
     dzflzi(k)=dzf(k)/lz
     dzclzi(k)=dzc(k)/lz
     zclzi(k)=zc(k)/lz
  enddo
  !@cuf istat=cudaDeviceSynchronize()
  !
  ! test input files before proceeding with the calculation
  !
#ifdef USE_NVTX
  call nvtxStartRange("sanity", 3)
#endif
  call test_sanity(ng,n,dims,stop_type,cbcvel,cbcpre,bcvel,bcpre,is_forced, &
                   dli,dzci,dzfi)
#ifdef USE_NVTX
  call nvtxEndRange
#endif
  !
  if(.not.restart) then
    istep = 0
    time = 0.
#ifdef USE_NVTX
    call nvtxStartRange("initflow", 4)
#endif
#ifdef USE_CUDA
  istat = cudaMemPrefetchAsync( dzci, size(dzci), cudaCpuDeviceId, 0 )
  istat = cudaMemPrefetchAsync( dzfi, size(dzfi), cudaCpuDeviceId, 0 )
  istat = cudaMemPrefetchAsync( dzflzi, size(dzflzi), cudaCpuDeviceId, 0 )
  istat = cudaMemPrefetchAsync( dzclzi, size(dzclzi), cudaCpuDeviceId, 0 )
  istat = cudaMemPrefetchAsync( zclzi, size(zclzi), cudaCpuDeviceId, 0 )
#endif
  call initflow(inivel,n,zclzi,dzclzi,dzflzi,visc,u,v,w,p)
#ifdef USE_NVTX
  call nvtxEndRange
#endif
    if(myid.eq.0) print*, '*** Initial condition succesfully set ***'
  else
    call load('r',trim(datadir)//'fld.bin',n,u(1:n(1),1:n(2),1:n(3)), &
                                             v(1:n(1),1:n(2),1:n(3)), &
                                             w(1:n(1),1:n(2),1:n(3)), &
                                             p(1:n(1),1:n(2),1:n(3)), &
                                             time,ristep)
    istep = nint(ristep)
    if(myid.eq.0) print*, '*** Checkpoint loaded at time = ', time, 'time step = ', istep, '. ***'
  endif
#ifdef USE_CUDA
  istat = cudaMemPrefetchAsync( u, size(u), mydev, 0 )
  istat = cudaMemPrefetchAsync( v, size(v), mydev, 0 )
  istat = cudaMemPrefetchAsync( w, size(w), mydev, 0 )
  istat = cudaMemPrefetchAsync( p, size(p), mydev, 0 )
#endif
  call bounduvw(cbcvel,n,bcvel,.false.,dl,dzc,dzf,u,v,w)
  call boundp(cbcpre,n,bcpre,dl,dzc,dzf,p)
  !
  ! post-process and write initial condition
  !
  write(fldnum,'(i7.7)') istep
#ifndef USE_CATALYST
  call strain_rate(  n,dli,dzci,u,v,w,str)
  call rotation_rate(n,dli,dzci,u,v,w,ens)
  call q_criterion(n,ens,str,qcr)
  call write_visu_3d(datadir,'qcr_fld_'//fldnum//'.bin','log_visu_3d.out','Q_criterion', &
                     (/1,1,1/),(/ng(1),ng(2),ng(3)/),(/1,1,1/),time,istep,qcr(1:n(1),1:n(2),1:n(3)))
  include 'out1d.h90'
  include 'out2d.h90'
  include 'out3d.h90'
#else
  call boundp(cbcpre,n,bcpre,dl,dzc,dzf,qcr)
  catalyst_active = .TRUE.
  call CatalystInitialize(catalyst_active)
  call InitializeFlowGrid(n+[2,2,0],ng+[2,2,0],xc,yc,zc(1), &
                          u(0,0,1), v(0,0,1), w(0,0,1), p(0,0,1), qcr(0,0,1))
  if (myid .eq. 0) print*, "Running Catalyst pipeline..."
  t0 = MPI_WTIME()
  call CatalystCoProcess(0.d0, 0)
  t1 = MPI_WTIME()
  if (myid .eq. 0) print*, "Done. Catalyst time:", t1-t0
#endif

  !$cuf kernel do(3) <<<*,*>>> 
  do k=1,n(3)
    do j=1,n(2)
      do i=1,n(1)
        dudtrko(i,j,k) = 0.
        dvdtrko(i,j,k) = 0.
        dwdtrko(i,j,k) = 0.
      enddo
    enddo
  enddo
  !@cuf istat=cudaDeviceSynchronize()
#ifdef USE_NVTX
  call nvtxStartRange("chkdt", 1)
#endif
  call chkdt(n,dl,dzci,dzfi,visc,u,v,w,dtmax)
#ifdef USE_NVTX
  call nvtxEndRange
#endif
  dt = min(cfl*dtmax,dtmin)
  if(myid.eq.0) print*, 'dtmax = ', dtmax, 'dt = ',dt
  dti = 1./dt
  kill = .false.
  !
  ! initialize Poisson solver
  !
  call initsolver(n,dli,dzci,dzfi,cbcpre,bcpre(:,:),lambdaxyp,(/'c','c','c'/),ap,bp,cp,arrplanp,normfftp,rhsbp%x,rhsbp%y,rhsbp%z)
#ifdef USE_CUDA
  istat = cudaMemPrefetchAsync( rhsbp%x, size(rhsbp%x), mydev, 0 )
  istat = cudaMemPrefetchAsync( rhsbp%y, size(rhsbp%y), mydev, 0 )
  istat = cudaMemPrefetchAsync( rhsbp%z, size(rhsbp%z), mydev, 0 )
#endif
#ifdef IMPDIFF
  call initsolver(n,dli,dzci,dzfi,cbcvel(:,:,1),bcvel(:,:,1),lambdaxyu,(/'f','c','c'/),au,bu,cu,arrplanu,normfftu, &
                  rhsbu%x,rhsbu%y,rhsbu%z)
  call initsolver(n,dli,dzci,dzfi,cbcvel(:,:,2),bcvel(:,:,2),lambdaxyv,(/'c','f','c'/),av,bv,cv,arrplanv,normfftv, &
                  rhsbv%x,rhsbv%y,rhsbv%z)
  call initsolver(n,dli,dzci,dzfi,cbcvel(:,:,3),bcvel(:,:,3),lambdaxyw,(/'c','c','f'/),aw,bw,cw,arrplanw,normfftw, &
                  rhsbw%x,rhsbw%y,rhsbw%z)
#endif
  !
  ! main loop
  !
  if(myid.eq.0) print*, '*** Calculation loop starts now ***'
  is_done = .false.
  do while(.not.is_done)
#ifdef USE_NVTX
    call nvtxStartRange("timestep", istep)
#endif
#ifdef TIMING
    dt12 = MPI_WTIME()
#endif
    istep = istep + 1
    time = time + dt
    if(myid.eq.0) print*, 'Timestep #', istep, 'Time = ', time
    dpdl(:)  = 0.
    tauxo(:) = 0.
    tauyo(:) = 0.
    tauzo(:) = 0.
    do irk=1,3
      dtrk = sum(rkcoeff(:,irk))*dt
      dtrki = 1./dtrk
#ifndef IMPDIFF
      #ifdef USE_NVTX
      call nvtxStartRange("rk", irk)
      #endif
      call rk(rkcoeff(:,irk),n,dli,dzci,dzfi,dzflzi,dzclzi,visc,dt,l, &
              u,v,w,p,dudtrk,dvdtrk,dwdtrk, &
              dudtrko,dvdtrko,dwdtrko,tauxo,tauyo,tauzo,up,vp,wp,f)
      !swap d*dtrk <=> d*dtrko (saves data movement)
      rk_tmp  => dudtrk
      dudtrk  => dudtrko
      dudtrko => rk_tmp
      rk_tmp  => dvdtrk
      dvdtrk  => dvdtrko
      dvdtrko => rk_tmp
      rk_tmp  => dwdtrk
      dwdtrk  => dwdtrko
      dwdtrko => rk_tmp
      #ifdef USE_NVTX
      call nvtxEndRange
      #endif
#else
      #ifdef USE_NVTX
      call nvtxStartRange("rk_id", irk)
      #endif
      call rk_id(rkcoeff(:,irk),n,dli,dzci,dzfi,dzflzi,dzclzi,visc,dt,l, &
                 u,v,w,p,dudtrk,dvdtrk,dwdtrk,dudtrkd,dvdtrkd,dwdtrkd, &
                 dudtrko,dvdtrko,dwdtrko,tauxo,tauyo,tauzo,up,vp,wp,f)
      #ifdef USE_NVTX
      call nvtxEndRange
      #endif
#endif
      #ifdef USE_NVTX
      call nvtxStartRange("is_forced", irk+1)
      #endif
      if(is_forced(1)) then
        f1d=f(1)
        !$cuf kernel do(3) <<<*,*>>> 
        do k=1,n(3)
          do j=1,n(2)
            do i=1,n(1)
              up(i,j,k) = up(i,j,k) + f1d
            enddo
          enddo
        enddo
        !@cuf istat=cudaDeviceSynchronize()
      endif
      if(is_forced(2)) then
        f2d=f(2)
        !$cuf kernel do(3) <<<*,*>>> 
        do k=1,n(3)
          do j=1,n(2)
            do i=1,n(1)
              vp(i,j,k) = vp(i,j,k) + f2d
            enddo
          enddo
        enddo
        !@cuf istat=cudaDeviceSynchronize()
      endif
      if(is_forced(3)) then
        !$cuf kernel do(3) <<<*,*>>> 
        f3d=f(3)
        do k=1,n(3)
          do j=1,n(2)
            do i=1,n(1)
              wp(i,j,k) = wp(i,j,k) + f3d
            enddo
          enddo
        enddo
        !@cuf istat=cudaDeviceSynchronize()
      endif
      #ifdef USE_NVTX
      call nvtxEndRange
      #endif
#ifdef IMPDIFF
      alpha = -1./(.5*visc*dtrk)
      #ifdef USE_CUDA
      !$cuf kernel do(3) <<<*,*>>> 
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            up(i,j,k) = up(i,j,k) *alpha
          enddo
        enddo
      enddo
      !@cuf istat=cudaDeviceSynchronize()
      #else
      !$OMP WORKSHARE
      up(1:n(1),1:n(2),1:n(3)) = up(1:n(1),1:n(2),1:n(3))*alpha
      !$OMP END WORKSHARE
      #endif
      !$cuf kernel do(1) <<<*,*>>> 
      do k=1,ktot
        bb(k) = bu(k) + alpha
      enddo
      !@cuf istat=cudaDeviceSynchronize()
      #ifdef USE_NVTX
      call nvtxStartRange("solver_u", irk+3)
      #endif
      call updt_rhs_b((/'f','c','c'/),cbcvel(:,:,1),n,rhsbu%x,rhsbu%y,rhsbu%z,up)
      call solver(n,arrplanu,normfftu,lambdaxyu,au,bb,cu,cbcvel(:,:,1),(/'f','c','c'/),up)
      #ifdef USE_NVTX
      call nvtxEndRange
      #endif
      #ifdef USE_CUDA
      !$cuf kernel do(3) <<<*,*>>> 
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            vp(i,j,k) = vp(i,j,k) *alpha
          enddo
        enddo
      enddo
      !@cuf istat=cudaDeviceSynchronize()
      #else
      !$OMP WORKSHARE
      vp(1:n(1),1:n(2),1:n(3)) = vp(1:n(1),1:n(2),1:n(3))*alpha
      !$OMP END WORKSHARE
      #endif
      !$cuf kernel do(1) <<<*,*>>> 
      do k=1,ktot
        bb(k) = bv(k) + alpha
      enddo
      !@cuf istat=cudaDeviceSynchronize()
      #ifdef USE_NVTX
      call nvtxStartRange("solver_v", irk+4)
      #endif
      call updt_rhs_b((/'c','f','c'/),cbcvel(:,:,2),n,rhsbv%x,rhsbv%y,rhsbv%z,vp)
      call solver(n,arrplanv,normfftv,lambdaxyv,av,bb,cv,cbcvel(:,:,2),(/'c','f','c'/),vp)
      #ifdef USE_NVTX
      call nvtxEndRange
      #endif
      #ifdef USE_CUDA
      !$cuf kernel do(3) <<<*,*>>> 
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            wp(i,j,k) = wp(i,j,k) *alpha
          enddo
        enddo
      enddo
      !@cuf istat=cudaDeviceSynchronize()
      #else
      !$OMP WORKSHARE
      wp(1:n(1),1:n(2),1:n(3)) = wp(1:n(1),1:n(2),1:n(3))*alpha
      !$OMP END WORKSHARE
      #endif
      !$cuf kernel do(1) <<<*,*>>> 
      do k=1,ktot
        bb(k) = bw(k) + alpha
      enddo
      !@cuf istat=cudaDeviceSynchronize()
      #ifdef USE_NVTX
      call nvtxStartRange("solver_w", irk+5)
      #endif
      call updt_rhs_b((/'c','c','f'/),cbcvel(:,:,3),n,rhsbw%x,rhsbw%y,rhsbw%z,wp)
      call solver(n,arrplanw,normfftw,lambdaxyw,aw,bb,cw,cbcvel(:,:,3),(/'c','c','f'/),wp)
      #ifdef USE_NVTX
      call nvtxEndRange
      #endif
#endif
      dpdl(:) = dpdl(:) + f(:)
#ifdef DEBUG
      #ifdef USE_NVTX
      call nvtxStartRange("chkmean", irk)
      #endif
      if(is_forced(1)) then
        call chkmean(n,dzflzi,up,meanvel)
        if(myid.eq.0) print*,'Mean u = ', meanvel
      if(is_forced(2)) then
        call chkmean(n,dzflzi,vp,meanvel)
        if(myid.eq.0) print*,'Mean v = ', meanvel
      endif
      if(is_forced(3)) then
        call chkmean(n,dzclzi,wp,meanvel)
        if(myid.eq.0) print*,'Mean w = ', meanvel
      endif
      #ifdef USE_NVTX
      call nvtxEndRange
      #endif
#endif
      #ifdef USE_NVTX
      call nvtxStartRange("bounduvw", irk+5)
      #endif
      call bounduvw(cbcvel,n,bcvel,.false.,dl,dzc,dzf,up,vp,wp) ! outflow BC only at final velocity
      #ifdef USE_NVTX
      call nvtxEndRange
      call nvtxStartRange("fillps", irk+6)
      #endif
      call fillps(n,dli,dzfi,dtrki,up,vp,wp,pp)
      #ifdef USE_NVTX
      call nvtxEndRange
      call nvtxStartRange("updt_rhs_b", irk+7)
      #endif
      call updt_rhs_b((/'c','c','c'/),cbcpre,n,rhsbp%x,rhsbp%y,rhsbp%z,pp)
      #ifdef USE_NVTX
      call nvtxEndRange
      call nvtxStartRange("solver", irk+5)
      #endif
      call solver(n,arrplanp,normfftp,lambdaxyp,ap,bp,cp,cbcpre(:,:),(/'c','c','c'/),pp)
      #ifdef USE_NVTX
      call nvtxEndRange
      call nvtxStartRange("boundp", irk+6)
      #endif
      call boundp(cbcpre,n,bcpre,dl,dzc,dzf,pp)
      #ifdef USE_NVTX
      call nvtxEndRange
      call nvtxStartRange("correc", irk+7)
      #endif
      call correc(n,dli,dzci,dtrk,pp,up,vp,wp,u,v,w)
      #ifdef USE_NVTX
      call nvtxEndRange
      call nvtxStartRange("bounduvw", irk+8)
      #endif
      call bounduvw(cbcvel,n,bcvel,.true.,dl,dzc,dzf,u,v,w)
      #ifdef USE_NVTX
      call nvtxEndRange
      #endif
#ifdef IMPDIFF
      alphai = alpha**(-1)
      #ifdef USE_CUDA
      !$cuf kernel do(3) <<<*,*>>> 
      #else
      !$OMP PARALLEL DO DEFAULT(none) &
      !$OMP PRIVATE(i,j,k,im,jm,km,ip,jp,kp) &
      !$OMP SHARED(n,p,pp,dxi,dyi,dzfi,dzci,alphai)
      #endif
      do k=1,n(3)
        kp = k + 1
        km = k - 1
        do j=1,n(2)
          jp = j + 1
          jm = j - 1
          do i=1,n(1)
            ip = i + 1
            im = i - 1
            p(i,j,k) = p(i,j,k) + pp(i,j,k) + alphai*( &
                        (pp(ip,j,k)-2.*pp(i,j,k)+pp(im,j,k))*(dxi**2) + &
                        (pp(i,jp,k)-2.*pp(i,j,k)+pp(i,jm,k))*(dyi**2) + &
                        ((pp(i,j,kp)-pp(i,j,k ))*dzci(k ) - &
                         (pp(i,j,k )-pp(i,j,km))*dzci(km))*dzfi(k) )
          enddo
        enddo
      enddo
      #ifndef USE_CUDA
      !$OMP END PARALLEL DO
      #endif
      !@cuf istat=cudaDeviceSynchronize()
#else
      #ifdef USE_CUDA
      !$cuf kernel do(3) <<<*,*>>> 
      do k=1,n(3)    
        do j=1,n(2)
          do i=1,n(1)
            p(i,j,k) = p(i,j,k) + pp(i,j,k)  
          enddo
        enddo
      enddo
      !@cuf istat=cudaDeviceSynchronize()
      #else
      !$OMP WORKSHARE
      p(1:n(1),1:n(2),1:n(3)) = p(1:n(1),1:n(2),1:n(3)) + pp(1:n(1),1:n(2),1:n(3))
      !$OMP END WORKSHARE
      #endif 
#endif
      #ifdef USE_NVTX
      call nvtxStartRange("boundp", irk+9)
      #endif
      call boundp(cbcpre,n,bcpre,dl,dzc,dzf,p)
      #ifdef USE_NVTX
      call nvtxEndRange
      #endif
    enddo
    dpdl(:) = -dpdl(:)*dti
    !
    ! check simulation stopping criteria
    !
    if(stop_type(1)) then ! maximum number of time steps reached
      if(istep.ge.nstep   ) is_done = is_done.or..true.
    endif
    if(stop_type(2)) then ! maximum simulation time reached
      if(time .ge.time_max) is_done = is_done.or..true.
    endif
    if(stop_type(3)) then ! maximum wall-clock time reached
      tw = (MPI_WTIME()-twi)/3600.
      if(tw   .ge.tw_max  ) is_done = is_done.or..true.
    endif
    if(mod(istep,icheck).eq.0) then
      if(myid.eq.0) print*, 'Checking stability and divergence...'
      call chkdt(n,dl,dzci,dzfi,visc,u,v,w,dtmax)
      dt  = min(cfl*dtmax,dtmin)
      if(myid.eq.0) print*, 'dtmax = ', dtmax, 'dt = ',dt
      if(dtmax.lt.small) then
        if(myid.eq.0) print*, 'ERROR: timestep is too small.'
        if(myid.eq.0) print*, 'Aborting...'
        is_done = .true.
        kill = .true.
      endif
      dti = 1./dt
      call chkdiv(n,dli,dzfi,u,v,w,divtot,divmax)
      if(divmax.gt.small.or.divtot.ne.divtot) then
        if(myid.eq.0) print*, 'ERROR: maximum divergence is too large.'
        if(myid.eq.0) print*, 'Aborting...'
        is_done = .true.
        kill = .true.
      endif
    endif
    !
    ! output routines below
    !
    if(mod(istep,iout0d).eq.0) then
      !allocate(var(4))
      var(1) = 1.*istep
      var(2) = dt
      var(3) = time
      call out0d(trim(datadir)//'time.out',3,var)
      !
      if(any(is_forced(:)).or.any(abs(bforce(:)).gt.0.)) then
        meanvelu = 0.
        meanvelv = 0.
        meanvelw = 0.
        if(is_forced(1).or.abs(bforce(1)).gt.0.) then
          call chkmean(n,dzf/lz,up,meanvelu)
        endif
        if(is_forced(2).or.abs(bforce(2)).gt.0.) then
          call chkmean(n,dzf/lz,vp,meanvelv)
        endif
        if(is_forced(3).or.abs(bforce(3)).gt.0.) then
          call chkmean(n,dzc/lz,wp,meanvelw)
        endif
        if(.not.any(is_forced(:))) dpdl(:) = -bforce(:) ! constant pressure gradient
        var(1)   = time
        var(2:4) = dpdl(1:3)
        var(5:7) = (/meanvelu,meanvelv,meanvelw/)
        call out0d(trim(datadir)//'forcing.out',7,var)
      endif
      !deallocate(var)
    endif
    write(fldnum,'(i7.7)') istep
    if(mod(istep,iout1d).eq.0) then
      include 'out1d.h90'
    endif
    if(mod(istep,iout2d).eq.0) then
      include 'out2d.h90'
    endif
    if(mod(istep,iout3d).eq.0) then
      call strain_rate(  n,dli,dzci,u,v,w,str)
      call rotation_rate(n,dli,dzci,u,v,w,ens)
      call q_criterion(n,ens,str,qcr)
#ifndef USE_CATALYST
      call write_visu_3d(datadir,'qcr_fld_'//fldnum//'.bin','log_visu_3d.out','Q_criterion', &
                         (/1,1,1/),(/ng(1),ng(2),ng(3)/),(/1,1,1/),time,istep, &
                         qcr(1:n(1),1:n(2),1:n(3)))
      include 'out3d.h90'
#else
      call boundp(cbcpre,n,bcpre,dl,dzc,dzf,qcr)
      if (myid .eq. 0) print*, "Running Catalyst pipeline..."
      t0 = MPI_WTIME()
      call CatalystCoProcess(time, istep)
      t1 = MPI_WTIME()
      if (myid .eq. 0) print*, "Done. Catalyst time:", t1-t0
#endif
    endif
    if(mod(istep,isave ).eq.0.or.(is_done.and..not.kill)) then
      ristep = 1.*istep
      if(is_overwrite_save) then
        filename = 'fld.bin'
      else
        filename = 'fld_'//fldnum//'.bin'
      endif
      call load('w',trim(datadir)//trim(filename),n,u(1:n(1),1:n(2),1:n(3)), &
                                                    v(1:n(1),1:n(2),1:n(3)), &
                                                    w(1:n(1),1:n(2),1:n(3)), &
                                                    p(1:n(1),1:n(2),1:n(3)), &
                                                    time,ristep)
      if(.not.is_overwrite_save) then
        !
        ! fld.bin -> last checkpoint file (symbolic link)
        !
        if(myid.eq.0) call execute_command_line('ln -sf '//trim(filename)//' '//trim(datadir)//'fld.bin')
      endif
      if(myid.eq.0) print*, '*** Checkpoint saved at time = ', time, 'time step = ', istep, '. ***'
    endif
#ifdef TIMING
      dt12 = MPI_WTIME()-dt12
      call MPI_ALLREDUCE(dt12,dt12av ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(dt12,dt12min,1,MPI_REAL_RP,MPI_MIN,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(dt12,dt12max,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
      if(myid.eq.0) print*, 'Avrg, min & max elapsed time: '
      if(myid.eq.0) print*, dt12av/(1.*product(dims)),dt12min,dt12max
#endif
#ifdef USE_NVTX
      call nvtxEndRange
#endif
  enddo
  !
  ! clear ffts
  !
  call fftend(arrplanp)
#ifdef IMPDIFF
  call fftend(arrplanu)
  call fftend(arrplanv)
  call fftend(arrplanw)
#endif
  !
  ! deallocate memory
  !
  deallocate(u,v,w,p,up,vp,wp,pp)
  deallocate(dudtrko,dvdtrko,dwdtrko)
  deallocate(dudtrk,dvdtrk,dwdtrk)
  deallocate(lambdaxyp)
  deallocate(ap,bp,cp)
  deallocate(rhsbp%x,rhsbp%y,rhsbp%z)
#ifdef IMPDIFF
  deallocate(dudtrkd,dvdtrkd,dwdtrkd)
  deallocate(lambdaxyu,lambdaxyv,lambdaxyw)
  deallocate(au,bu,cu,av,bv,cv,aw,bw,cw,bb)
  deallocate(rhsbu%x,rhsbu%y,rhsbu%z, &
             rhsbv%x,rhsbv%y,rhsbv%z, &
             rhsbw%x,rhsbw%y,rhsbw%z)
#endif
  deallocate(dzc,dzf,zc,zf,dzci,dzfi,dzflzi,dzclzi,zclzi)
  !
  if(myid.eq.0.and.(.not.kill)) print*, '*** Fim ***'
  call decomp_2d_finalize
  call MPI_FINALIZE(ierr)
  stop
end program cans
