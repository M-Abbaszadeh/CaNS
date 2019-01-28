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
  use mod_bound      , only: boundp,bounduvw,updt_rhs_b
  use mod_chkdiv     , only: chkdiv
  use mod_chkdt      , only: chkdt
  use mod_common_mpi , only: myid,ierr
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
  use mod_output     , only: out0d,out1d,out1d_2,out2d,out3d
  use mod_param      , only: itot,jtot,ktot,lx,ly,lz,dx,dy,dz,dxi,dyi,dzi,visc,small, &
                             cbcvel,bcvel,cbcpre,bcpre, &
                             icheck,iout0d,iout1d,iout2d,iout3d,isave, &
                             nstep,restart, &
                             rkcoeff, &
                             datadir, &
                             cfl,     &
                             inivel,  &
                             uref,lref, &
                             imax,jmax,dims, &
                             nthreadsmax, &
                             gr, &
                             is_outflow,no_outflow,is_forced
  use mod_sanity     , only: test_sanity
  use mod_solver     , only: solver
#ifdef USE_CUDA
  use mod_common_mpi, only: mydev
#endif
#ifdef USE_NVTX
  use nvtx
#endif 
  !$ use omp_lib
  implicit none
  integer, parameter, dimension(3) :: ng = (/itot,jtot,ktot/)
  integer, parameter, dimension(3) :: n  = (/imax,jmax,ktot/)
  real(8), parameter, dimension(3) :: l   = (/lx,ly,lz/)
  real(8), parameter, dimension(3) :: dl  = (/dx,dy,dz/)
  real(8), parameter, dimension(3) :: dli = (/dxi,dyi,dzi/)
  !real(8), dimension(0:imax+1,0:jmax+1,0:ktot+1) :: u,v,w,p,up,vp,wp,pp
  real(8), dimension(:,:,:),allocatable :: u,v,w,p,up,vp,wp,pp
#ifdef USE_CUDA
  attributes(managed):: u,v,w,p,up,vp,wp,pp
#endif
  !real(8), dimension(imax,jmax,ktot)    :: dudtrko,dvdtrko,dwdtrko
  real(8), dimension(:,:,:),allocatable, target  :: dudtrk_A,dvdtrk_A,dwdtrk_A
  real(8), dimension(:,:,:),allocatable, target  :: dudtrk_B,dvdtrk_B,dwdtrk_B
  real(8), dimension(:,:,:), pointer :: dudtrko,dvdtrko,dwdtrko
  real(8), dimension(:,:,:), pointer :: dudtrk,dvdtrk,dwdtrk, rk_tmp
  real(8), dimension(3) :: tauxo,tauyo,tauzo
  real(8), dimension(3) :: f
  type(C_PTR), dimension(2,2) :: arrplanp
  !real(8), dimension(imax,jmax) :: lambdaxyp
  real(8), dimension(:,:),allocatable :: lambdaxyp
  !real(8), dimension(ktot) :: ap,bp,cp
  real(8), dimension(:),allocatable :: ap,bp,cp
  real(8) :: normfftp
  type rhs_bound
    !real(8), dimension(n(2),n(3),0:1) :: x
    !real(8), dimension(n(1),n(3),0:1) :: y
    !real(8), dimension(n(1),n(2),0:1) :: z
#ifdef USE_CUDA
    real(8), dimension(:,:,:),allocatable,managed :: x,y,z
#else
    real(8), dimension(:,:,:),allocatable :: x,y,z
#endif
  end type rhs_bound 
  integer :: i,j,k,im,ip,jm,jp,km,kp
#ifdef IMPDIFF
  type(C_PTR), dimension(2,2) :: arrplanu,arrplanv,arrplanw
  real(8), dimension(:,:),allocatable :: lambdaxyu,lambdaxyv,lambdaxyw
  real(8), dimension(:)  ,allocatable :: au,av,aw,bu,bv,bw,bb,cu,cv,cw
  real(8) :: normfftu,normfftv,normfftw
  real(8) :: alpha,alphai
  type(rhs_bound), allocatable :: rhsbu,rhsbv,rhsbw
  real(8), dimension(:,:,:),allocatable    :: dudtrkd,dvdtrkd,dwdtrkd
  #ifdef USE_CUDA
  attributes(managed):: dudtrkd,dvdtrkd,dwdtrkd,lambdaxyu,lambdaxyv,lambdaxyw
  attributes(managed):: au,av,aw,bu,bv,bw,bb,cu,cv,cw
  attributes(managed):: rhsbu,rhsbv,rhsbw
  #endif
#endif
  type(rhs_bound), allocatable :: rhsbp
  real(8) :: ristep
  real(8) :: dt,dti,dtmax,time,dtrk,dtrki,divtot,divmax
  integer :: irk,istep
  !real(8), dimension(0:ktot+1) :: dzc,dzf,zc,zf,dzci,dzfi
  real(8), dimension(:),allocatable :: dzc,dzf,zc,zf,dzci,dzfi,dzflzi,dzclzi,zclzi
#ifdef USE_CUDA
  integer :: istat
  integer(kind=cuda_count_kind):: freeMem,totMem
  integer(8):: totEle
  attributes(managed):: dzc,dzf,zc,zf,dzci,dzfi,dzflzi,dzclzi,zclzi,dudtrk_A,dvdtrk_A,dwdtrk_A,lambdaxyp,ap,bp,cp,rhsbp
  attributes(managed):: dudtrko,dvdtrko,dwdtrko,dudtrk,dvdtrk,dwdtrk,dudtrk_B,dvdtrk_B,dwdtrk_B
#endif
  real(8) :: meanvel
  real(8), dimension(3) :: dpdl
  !real(8), allocatable, dimension(:) :: var
  real(8), dimension(10) :: var
#ifdef TIMING
  real(8) :: dt12,dt12av,dt12min,dt12max
#endif
  real(8):: f1d,f2d,f3d
  character(len=7) :: fldnum
  integer :: lenr,kk
  logical :: kill
  !
  !$call omp_set_num_threads(nthreadsmax)
  call initmpi(ng,cbcpre)

!
! Allocate memory.
!
#ifdef USE_CUDA
  istat = cudaMemAdvise( u, size(u), cudaMemAdviseSetPreferredLocation, mydev )
  istat = cudaMemAdvise( v, size(v), cudaMemAdviseSetPreferredLocation, mydev )
  istat = cudaMemAdvise( w, size(w), cudaMemAdviseSetPreferredLocation, mydev )
  istat = cudaMemAdvise( p, size(p), cudaMemAdviseSetPreferredLocation, mydev )
  istat = cudaMemAdvise( up, size(up), cudaMemAdviseSetPreferredLocation, mydev )
  istat = cudaMemAdvise( vp, size(vp), cudaMemAdviseSetPreferredLocation, mydev )
  istat = cudaMemAdvise( wp, size(wp), cudaMemAdviseSetPreferredLocation, mydev )
  istat = cudaMemAdvise( pp, size(pp), cudaMemAdviseSetPreferredLocation, mydev )
  istat = cudaMemAdvise( dudtrk_A, size(dudtrk_A), cudaMemAdviseSetPreferredLocation, mydev )
  istat = cudaMemAdvise( dvdtrk_A, size(dvdtrk_A), cudaMemAdviseSetPreferredLocation, mydev )
  istat = cudaMemAdvise( dwdtrk_A, size(dwdtrk_A), cudaMemAdviseSetPreferredLocation, mydev )
  istat = cudaMemAdvise( dudtrk_B, size(dudtrk_B), cudaMemAdviseSetPreferredLocation, mydev )
  istat = cudaMemAdvise( dvdtrk_B, size(dvdtrk_B), cudaMemAdviseSetPreferredLocation, mydev )
  istat = cudaMemAdvise( dwdtrk_B, size(dwdtrk_B), cudaMemAdviseSetPreferredLocation, mydev )
#endif

  allocate( u(0:imax+1,0:jmax+1,0:ktot+1)) 
  allocate( v(0:imax+1,0:jmax+1,0:ktot+1))
  allocate( w(0:imax+1,0:jmax+1,0:ktot+1)) 
  allocate( p(0:imax+1,0:jmax+1,0:ktot+1)) 
  allocate(up(0:imax+1,0:jmax+1,0:ktot+1)) 
  allocate(vp(0:imax+1,0:jmax+1,0:ktot+1)) 
  allocate(wp(0:imax+1,0:jmax+1,0:ktot+1)) 
  allocate(pp(0:imax+1,0:jmax+1,0:ktot+1)) 
  allocate(dudtrk_A(imax,jmax,ktot))
  allocate(dvdtrk_A(imax,jmax,ktot))
  allocate(dwdtrk_A(imax,jmax,ktot))
  allocate(dudtrk_B(imax,jmax,ktot))
  allocate(dvdtrk_B(imax,jmax,ktot))
  allocate(dwdtrk_B(imax,jmax,ktot))
#ifdef USE_CUDA
  istat = cudaMemPrefetchAsync( u, size(u), cudaCpuDeviceId, 0 )
  istat = cudaMemPrefetchAsync( v, size(v), cudaCpuDeviceId, 0 )
  istat = cudaMemPrefetchAsync( w, size(w), cudaCpuDeviceId, 0 )
  istat = cudaMemPrefetchAsync( p, size(p), cudaCpuDeviceId, 0 )
  istat = cudaMemPrefetchAsync( up, size(up), mydev, 0 )
  istat = cudaMemPrefetchAsync( vp, size(vp), mydev, 0 )
  istat = cudaMemPrefetchAsync( wp, size(wp), mydev, 0 )
  istat = cudaMemPrefetchAsync( pp, size(pp), mydev, 0 )
  istat = cudaMemPrefetchAsync( dudtrk_A, size(dudtrk_A), mydev, 0 )
  istat = cudaMemPrefetchAsync( dvdtrk_A, size(dvdtrk_A), mydev, 0 )
  istat = cudaMemPrefetchAsync( dwdtrk_A, size(dwdtrk_A), mydev, 0 )
  istat = cudaMemPrefetchAsync( dudtrk_B, size(dudtrk_B), mydev, 0 )
  istat = cudaMemPrefetchAsync( dvdtrk_B, size(dvdtrk_B), mydev, 0 )
  istat = cudaMemPrefetchAsync( dwdtrk_B, size(dwdtrk_B), mydev, 0 )
#endif

  dudtrko => dudtrk_A
  dvdtrko => dvdtrk_A
  dwdtrko => dwdtrk_A
  dudtrk => dudtrk_B
  dvdtrk => dudtrk_B
  dwdtrk => dudtrk_B

#ifdef IMPDIFF
  allocate(dudtrkd(imax,jmax,ktot))    
  allocate(dvdtrkd(imax,jmax,ktot))   
  allocate(dwdtrkd(imax,jmax,ktot))    

  allocate(rhsbu,rhsbv,rhsbw)
  allocate(rhsbu%x(n(2),n(3),0:1))
  allocate(rhsbu%y(n(1),n(3),0:1))
  allocate(rhsbu%z(n(1),n(2),0:1))
  allocate(rhsbv%x(n(2),n(3),0:1))
  allocate(rhsbv%y(n(1),n(3),0:1))
  allocate(rhsbv%z(n(1),n(2),0:1))
  allocate(rhsbw%x(n(2),n(3),0:1))
  allocate(rhsbw%y(n(1),n(3),0:1))
  allocate(rhsbw%z(n(1),n(2),0:1))
#ifdef USE_CUDA
  allocate(lambdaxyu(itot/dims(2),jtot/dims(1)))
  allocate(lambdaxyv(itot/dims(2),jtot/dims(1)))
  allocate(lambdaxyw(itot/dims(2),jtot/dims(1)))
#else
  allocate(lambdaxyu(imax,jmax))
  allocate(lambdaxyv(imax,jmax))
  allocate(lambdaxyw(imax,jmax))
#endif
  allocate(au(ktot)) 
  allocate(bu(ktot))
  allocate(cu(ktot)) 
  allocate(av(ktot)) 
  allocate(bv(ktot))
  allocate(cv(ktot)) 
  allocate(aw(ktot)) 
  allocate(bw(ktot))
  allocate(cw(ktot)) 
  allocate(bb(ktot))
#endif

  allocate(rhsbp)
  allocate(rhsbp%x(n(2),n(3),0:1))
  allocate(rhsbp%y(n(1),n(3),0:1))
  allocate(rhsbp%z(n(1),n(2),0:1))
#ifdef USE_CUDA
  allocate(lambdaxyp(itot/dims(2),jtot/dims(1)))
#else
  allocate(lambdaxyp(imax,jmax))
#endif
  allocate(ap(ktot))
  allocate(bp(ktot))
  allocate(cp(ktot))

  allocate( dzc(0:ktot+1))
  allocate( dzf(0:ktot+1))
  allocate ( zc(0:ktot+1))
  allocate(  zf(0:ktot+1))
  allocate(dzci(0:ktot+1))
  allocate(dzfi(0:ktot+1))
  allocate(dzflzi(0:ktot+1))
  allocate(dzclzi(0:ktot+1))
  allocate(zclzi(0:ktot+1))

#ifdef USE_CUDA
  istat = cudaMemAdvise( rhsbp%x, size(rhsbp%x), cudaMemAdviseSetPreferredLocation, mydev )
  istat = cudaMemAdvise( rhsbp%y, size(rhsbp%y), cudaMemAdviseSetPreferredLocation, mydev )
  istat = cudaMemAdvise( rhsbp%z, size(rhsbp%z), cudaMemAdviseSetPreferredLocation, mydev )
  istat = cudaMemPrefetchAsync( rhsbp%x, size(rhsbp%x), cudaCpuDeviceId, 0 )
  istat = cudaMemPrefetchAsync( rhsbp%y, size(rhsbp%y), cudaCpuDeviceId, 0 )
  istat = cudaMemPrefetchAsync( rhsbp%z, size(rhsbp%z), cudaCpuDeviceId, 0 )

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

!!!!!!!!!
  if(myid.eq.0) print*, '******************************'
  if(myid.eq.0) print*, '*** Beginning of simulation ***'
  if(myid.eq.0) print*, '******************************'
  if(myid.eq.0) print*, ''

#ifdef USE_CUDA
  if(myid.eq.0) then 
      print*, ' GPU accelerated version, grid size:',n(1)*dims(1),n(2)*dims(2),n(3)
      #ifndef IMPDIFF
      totEle=(8*int((imax+2)*(jmax+2)*(ktot+2),8) + 6*int(n(1)*n(2)*n(3),8) + 6*(n(1)*n(2)+n(2)*n(3)+n(1)*n(3)) )
      #else
      totEle=(8*int((imax+2)*(jmax+2)*(ktot+2),8) + 9*int(n(1)*n(2)*n(3),8) +24*(n(1)*n(2)+n(2)*n(3)+n(1)*n(3)) )
      #endif
      totEle = totEle+ int((n(1)*dims(1)*n(2)*dims(2)/dims(1)*n(3)/dims(2)),8) + &
                       int((n(1)*n(2)*dims(2)*n(3)/dims(2)),8) + int(n(1)*n(2)*ng(3),8) + &
                       int((n(1)*dims(1)+2)* n(2)*dims(2)/dims(1)* n(3)/dims(2),8) +  &
                       int((n(2)*dims(2)+2)* n(1)* n(3)/dims(2) , 8 ) + &
                       int(n(2)*dims(2)* n(1)* n(3)/dims(2), 8 ) 
      print*, ' Estimated Memory usage (GB) per MPI task: ',totEle*8./(1024.**3)
      ! Missing fft plan
      istat=cudaMemGetInfo(freeMem,totMem)
      print *," Memory on GPU (GB)                      : ",totMem/(1024.**3)
      print *," "
  endif
#endif 

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
    inquire (iolength=lenr) dzc(1)
    open(99,file=trim(datadir)//'grid.bin',access='direct',recl=4*n(3)*lenr)
    write(99,rec=1) dzc(1:n(3)),dzf(1:n(3)),zc(1:n(3)),zf(1:n(3))
    close(99)
    open(99,file=trim(datadir)//'grid.out')
    do kk=0,ktot+1
      write(99,'(5E15.7)') 0.d0,zf(kk),zc(kk),dzf(kk),dzc(kk)
    enddo
    close(99)
  endif
  !
  ! test input files before proceeding with the calculation
  !
#ifdef USE_NVTX
  call nvtxStartRange("sanity", 3)
#endif
  call test_sanity(ng,n,dims,cbcvel,cbcpre,bcvel,bcpre,is_outflow,is_forced, &
                   dli,dzci,dzfi)
#ifdef USE_NVTX
  call nvtxEndRange
#endif
  !
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
     dzci(k)=1.d0/dzc(k)
     dzfi(k)=1.d0/dzf(k)
     dzflzi(k)=dzf(k)/lz
     dzclzi(k)=dzc(k)/lz
     zclzi(k)=zc(k)/lz
  enddo
  !@cuf istat=cudaDeviceSynchronize()

#ifdef USE_CUDA
  istat = cudaMemPrefetchAsync( dzci, size(dzci), cudaCpuDeviceId, 0 )
  istat = cudaMemPrefetchAsync( dzfi, size(dzfi), cudaCpuDeviceId, 0 )
  istat = cudaMemPrefetchAsync( dzflzi, size(dzflzi), cudaCpuDeviceId, 0 )
  istat = cudaMemPrefetchAsync( dzclzi, size(dzclzi), cudaCpuDeviceId, 0 )
  istat = cudaMemPrefetchAsync( zclzi, size(zclzi), cudaCpuDeviceId, 0 )
#endif

  if(.not.restart) then
    istep = 0
    time = 0.d0
#ifdef USE_NVTX
      call nvtxStartRange("initflow", 4)
#endif
    call initflow(inivel,n,zclzi,dzclzi,dzflzi,visc,uref,u,v,w,p)
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
  call bounduvw(cbcvel,n,bcvel,is_outflow,dl,dzc,dzf,u,v,w)
  call boundp(cbcpre,n,bcpre,dl,dzc,dzf,p)
#ifdef USE_NVTX
      call nvtxStartRange("chkdt", 1)
#endif
  call chkdt(n,dl,dzci,dzfi,visc,u,v,w,dtmax)
#ifdef USE_NVTX
      call nvtxEndRange
#endif
  dt = cfl*dtmax
  if(myid.eq.0) print*, 'dtmax = ', dtmax, 'dt = ',dt
  dti = 1.d0/dt
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
  do while(istep.lt.nstep)

 #ifdef USE_NVTX
      call nvtxStartRange("timestep", istep)
 #endif

#ifdef TIMING
    dt12 = MPI_WTIME()
#endif
    istep = istep + 1
    time = time + dt
    if(myid.eq.0) print*, 'Timestep #', istep, 'Time = ', time
    dpdl(:)  = 0.d0
    tauxo(:) = 0.d0
    tauyo(:) = 0.d0
    tauzo(:) = 0.d0
    do irk=1,3
      dtrk = sum(rkcoeff(:,irk))*dt
      dtrki = 1.d0/dtrk
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
        end do
        end do
       end do
  !@cuf istat=cudaDeviceSynchronize()
      end if

      if(is_forced(2)) then
      f2d=f(2)
      !$cuf kernel do(3) <<<*,*>>> 
       do k=1,n(3)
        do j=1,n(2)
         do i=1,n(1)
           vp(i,j,k) = vp(i,j,k) + f2d
        end do
        end do
       end do
  !@cuf istat=cudaDeviceSynchronize()
      end if

      if(is_forced(3)) then
      !$cuf kernel do(3) <<<*,*>>> 
      f3d=f(3)
       do k=1,n(3)
        do j=1,n(2)
         do i=1,n(1)
           wp(i,j,k) = wp(i,j,k) + f3d
        end do
        end do
       end do
  !@cuf istat=cudaDeviceSynchronize()
      end if
      #ifdef USE_NVTX
      call nvtxEndRange
      #endif

#ifdef IMPDIFF
      alpha = -1.d0/(.5d0*visc*dtrk)
      #ifdef USE_CUDA
      !$cuf kernel do(3) <<<*,*>>> 
       do k=1,n(3)
        do j=1,n(2)
         do i=1,n(1)
           up(i,j,k) = up(i,j,k) *alpha
        end do
        end do
       end do
  !@cuf istat=cudaDeviceSynchronize()
      #else
      !$OMP WORKSHARE
      up(1:n(1),1:n(2),1:n(3)) = up(1:n(1),1:n(2),1:n(3))*alpha
      !$OMP END WORKSHARE
      #endif
      
      !$cuf kernel do(1) <<<*,*>>> 
      do k=1,ktot
       bb(k) = bu(k) + alpha
      end do 
  !@cuf istat=cudaDeviceSynchronize()

      #ifdef USE_NVTX
      call nvtxStartRange("solver_u", irk+3)
      #endif
      call updt_rhs_b((/'f','c','c'/),cbcvel(:,:,1),n,rhsbu%x,rhsbu%y,rhsbu%z,up)
      call solver(n,arrplanu,normfftu,lambdaxyu,au,bb,cu,cbcvel(:,3,1),(/'f','c','c'/),up)
      #ifdef USE_NVTX
      call nvtxEndRange
      #endif
      #ifdef USE_CUDA
      !$cuf kernel do(3) <<<*,*>>> 
       do k=1,n(3)
        do j=1,n(2)
         do i=1,n(1)
           vp(i,j,k) = vp(i,j,k) *alpha
        end do
        end do
       end do
  !@cuf istat=cudaDeviceSynchronize()
      #else
      !$OMP WORKSHARE
      vp(1:n(1),1:n(2),1:n(3)) = vp(1:n(1),1:n(2),1:n(3))*alpha
      !$OMP END WORKSHARE
      #endif
      !$cuf kernel do(1) <<<*,*>>> 
      do k=1,ktot
       bb(k) = bv(k) + alpha
      end do 
  !@cuf istat=cudaDeviceSynchronize()

      #ifdef USE_NVTX
      call nvtxStartRange("solver_v", irk+4)
      #endif
      call updt_rhs_b((/'c','f','c'/),cbcvel(:,:,2),n,rhsbv%x,rhsbv%y,rhsbv%z,vp)
      call solver(n,arrplanv,normfftv,lambdaxyv,av,bb,cv,cbcvel(:,3,2),(/'c','f','c'/),vp)
      #ifdef USE_NVTX
      call nvtxEndRange
      #endif
      #ifdef USE_CUDA
      !$cuf kernel do(3) <<<*,*>>> 
       do k=1,n(3)
        do j=1,n(2)
         do i=1,n(1)
           wp(i,j,k) = wp(i,j,k) *alpha
        end do
        end do
       end do
  !@cuf istat=cudaDeviceSynchronize()
      #else
      !$OMP WORKSHARE
      wp(1:n(1),1:n(2),1:n(3)) = wp(1:n(1),1:n(2),1:n(3))*alpha
      !$OMP END WORKSHARE
      #endif
      !$cuf kernel do(1) <<<*,*>>> 
      do k=1,ktot
       bb(k) = bw(k) + alpha
      end do 
  !@cuf istat=cudaDeviceSynchronize()

      #ifdef USE_NVTX
      call nvtxStartRange("solver_w", irk+5)
      #endif
      call updt_rhs_b((/'c','c','f'/),cbcvel(:,:,3),n,rhsbw%x,rhsbw%y,rhsbw%z,wp)
      call solver(n,arrplanw,normfftw,lambdaxyw,aw,bb,cw,cbcvel(:,3,3),(/'c','c','f'/),wp)
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
      endif
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
      call bounduvw(cbcvel,n,bcvel,no_outflow,dl,dzc,dzf,up,vp,wp) ! outflow BC only at final velocity
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
      !call solver(n,arrplanp,normfftp,lambdaxyp,ap,bp,cp,cbcpre(:,3),(/'c','c','c'/),pp(1:n(1),1:n(2),1:n(3)))
      call solver(n,arrplanp,normfftp,lambdaxyp,ap,bp,cp,cbcpre(:,3),(/'c','c','c'/),pp)
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
      call bounduvw(cbcvel,n,bcvel,is_outflow,dl,dzc,dzf,u,v,w)
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
      !$OMP SHARED(p,pp,dzfi,dzci,alphai)
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
                        (pp(ip,j,k)-2.d0*pp(i,j,k)+pp(im,j,k))*(dxi**2) + &
                        (pp(i,jp,k)-2.d0*pp(i,j,k)+pp(i,jm,k))*(dyi**2) + &
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
        end do
        end do
       end do
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
    if(mod(istep,icheck).eq.0) then
      if(myid.eq.0) print*, 'Checking stability and divergence...'
      call chkdt(n,dl,dzci,dzfi,visc,u,v,w,dtmax)
      dt  = cfl*dtmax
      if(myid.eq.0) print*, 'dtmax = ', dtmax, 'dt = ',dt
      if(dtmax.lt.small) then
        if(myid.eq.0) print*, 'ERROR: timestep is too small.'
        if(myid.eq.0) print*, 'Aborting ...'
        istep = nstep + 1 ! i.e. exit main loop
        kill = .true.
      endif
      dti = 1.d0/dt
      call chkdiv(n,dli,dzfi,u,v,w,divtot,divmax)
      if(divmax.gt.small.or.divtot.ne.divtot) then
        if(myid.eq.0) print*, 'ERROR: maximum divergence is too large.'
        if(myid.eq.0) print*, 'Aborting ...'
        istep = nstep + 1 ! i.e. exit main loop
        kill = .true.
      endif
    endif
    !
    ! output routines below
    !
    if(mod(istep,iout0d).eq.0) then
      !allocate(var(4))
      var(1) = 1.d0*istep
      var(2) = dt
      var(3) = time
      call out0d(trim(datadir)//'time.out',3,var)
      if(any(is_forced(:))) then
        var(1)   = time
        var(2:4) = dpdl(1:3)
        call out0d(trim(datadir)//'forcing.out',4,var)
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
      include 'out3d.h90'
    endif
    if(mod(istep,isave ).eq.0) then
      ristep = 1.d0*istep
      call load('w',trim(datadir)//'fld.bin',n,u(1:n(1),1:n(2),1:n(3)), &
                                               v(1:n(1),1:n(2),1:n(3)), &
                                               w(1:n(1),1:n(2),1:n(3)), &
                                               p(1:n(1),1:n(2),1:n(3)), &
                                               time,ristep)
    if(myid.eq.0) print*, '*** Checkpoint saved at time = ', time, 'time step = ', istep, '. ***'
    endif
#ifdef TIMING
      dt12 = MPI_WTIME()-dt12
      call MPI_ALLREDUCE(dt12,dt12av ,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(dt12,dt12min,1,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(dt12,dt12max,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
      if(myid.eq.0) print*, 'Avrg, min & max elapsed time: '
      if(myid.eq.0) print*, dt12av/(1.d0*product(dims)),dt12min,dt12max
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
! Deallocate memory.
!
  deallocate( u,v,w,p,up,vp,wp,pp)
  deallocate(dudtrk_A,dvdtrk_A,dwdtrk_A)
  deallocate(dudtrk_B,dvdtrk_B,dwdtrk_B)
#ifdef IMPDIFF
  deallocate(dudtrkd,dvdtrkd,dwdtrkd)
  deallocate(rhsbu%x,rhsbu%y,rhsbu%z)
  deallocate(rhsbv%x,rhsbv%y,rhsbv%z)
  deallocate(rhsbw%x,rhsbw%y,rhsbw%z)
  deallocate(lambdaxyu,lambdaxyv,lambdaxyw)
  deallocate(au,bu,cu)
  deallocate(av,bv,cv)
  deallocate(aw,bw,cw)
  deallocate(bb)
#endif
  deallocate(rhsbp%x,rhsbp%y,rhsbp%z)
  deallocate(lambdaxyp,ap,bp,cp)
  deallocate(dzc,dzf,zc,zf,dzci,dzfi,dzflzi,dzclzi)

  if(myid.eq.0.and.(.not.kill)) print*, '*** Fim ***'
  call decomp_2d_finalize
  call MPI_FINALIZE(ierr)
  call exit
end program cans
