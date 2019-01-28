module mod_rk
  use mpi
  use mod_common_mpi, only: ierr
  use mod_param, only: is_forced,velf,dims
  use mod_debug, only: chkmean
  use mod_mom  , only: momxad,momyad,momzad,momxp,momyp,momzp,momp,momxyzad
  use mod_momd , only: momxpd,momypd,momzpd,momxa,momya,momza
  use mod_moms , only: momsad
#ifdef USE_NVTX
  use nvtx
#endif 
  implicit none
  private
  public rk,rk_id
  contains
  subroutine rk(rkpar,n,dli,dzci,dzfi,dzflzi,dzclzi,visc,dt,l,u,v,w,p,dudtrk,dvdtrk,dwdtrk,dudtrko,dvdtrko,dwdtrko,tauxo,tauyo,tauzo,up,vp,wp,f)
    !
    ! low-storage 3rd-order Runge-Kutta scheme 
    ! for time integration of the momentum equations.
    !
    !@cuf use cudafor
    implicit none
    real(8), intent(in), dimension(2) :: rkpar
    integer, intent(in), dimension(3) :: n
    real(8), intent(in) :: visc,dt
    real(8), intent(in   ), dimension(3) :: dli,l 
    real(8), intent(in   ), dimension(0:) :: dzci,dzfi,dzflzi,dzclzi
    real(8), intent(in   ), dimension(0:,0:,0:) :: u ,v ,w,p
    real(8), intent(inout), dimension(:,:,:) :: dudtrko,dvdtrko,dwdtrko
    real(8), intent(inout), dimension(3) :: tauxo,tauyo,tauzo
    real(8), intent(out), dimension(0:,0:,0:) :: up,vp,wp
    real(8), intent(out), dimension(3) :: f
    !real(8),              dimension(n(1),n(2),n(3)) ::          dudtrk, dvdtrk, dwdtrk
    real(8),              dimension(:,:,:) ::          dudtrk, dvdtrk, dwdtrk
    real(8) :: factor1,factor2,factor12
    real(8), dimension(3) :: taux,tauy,tauz
#ifdef USE_CUDA
   attributes(managed):: dzci,dzfi,dzflzi,dzclzi,u ,v ,w,p, dudtrk, dvdtrk, dwdtrk, up,vp,wp
   attributes(managed):: dudtrko,dvdtrko,dwdtrko
   integer:: istat
#endif
    real(8) :: dxi, dyi
    integer :: nx,ny,nz,im,ip,jm,jp,km,kp,i,j,k,nxg,nyg,nzg
    real(8) :: uuip,uuim,uvjp,uvjm,uwkp,uwkm
    real(8) :: dudxp,dudxm,dudyp,dudym,dudzp,dudzm
    real(8) :: taux2d,taux3d
    real(8) :: uvip,uvim,vvjp,vvjm,wvkp,wvkm
    real(8) :: dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm
    real(8) :: tauy1d,tauy3d
    real(8) :: uwip,uwim,vwjp,vwjm,wwkp,wwkm
    real(8) :: dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm
    real(8) :: tauz1d,tauz2d,dudtrk_temp,dvdtrk_temp,dwdtrk_temp
    real(8) :: mean
    !
    factor1 = rkpar(1)*dt
    factor2 = rkpar(2)*dt
    factor12 = factor1 + factor2
    !
    dxi = dli(1)
    dyi = dli(2)
    nx = n(1)
    ny = n(2)
    nz = n(3)
    !
#ifdef USE_NVTX
    call nvtxStartRange("momxyzad",5)
#endif
    !
#ifdef USE_CUDA
    !$cuf kernel do(3) <<<*,(8,8,8)>>>
#else
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k,im,jm,km,ip,jp,kp) &
    !$OMP PRIVATE(uuip,uuim,uvjp,uvjm,uwkp,uwkm) &
    !$OMP PRIVATE(dudxp,dudxm,dudyp,dudym,dudzp,dudzm) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,dzi,visc,u,v,w,dudt,dzci,dzfi)
#endif
    do k=1,nz
      kp = k + 1
      km = k - 1
      do j=1,ny
        jp = j + 1
        jm = j - 1
        do i=1,nx
          ip = i + 1
          im = i - 1
          ! touch u
          dudxp = (u(ip,j,k)-u(i,j,k))*dxi
          dudxm = (u(i,j,k)-u(im,j,k))*dxi
          dudyp = (u(i,jp,k)-u(i,j,k))*dyi
          dudym = (u(i,j,k)-u(i,jm,k))*dyi
          dudzp = (u(i,j,kp)-u(i,j,k))*dzci(k)
          dudzm = (u(i,j,k)-u(i,j,km))*dzci(km)
          uuip  = 0.25d0*( u(ip,j,k)+u(i,j,k) )*( u(ip,j ,k )+u(i,j ,k ) )
          uuim  = 0.25d0*( u(im,j,k)+u(i,j,k) )*( u(im,j ,k )+u(i,j ,k ) )

          dudtrk_temp = dxi*(     -uuip + uuim ) + (dudxp-dudxm)*visc*dxi + &
                                               + (dudyp-dudym)*visc*dyi + &
                                               + (dudzp-dudzm)*visc*dzfi(k)


          uvjp  = 0.25d0*( u(i,jp,k)+u(i,j,k) ) !*( v(ip,j ,k )+v(i,j ,k ) )
          uvjm  = 0.25d0*( u(i,jm,k)+u(i,j,k) ) !*( v(ip,jm,k )+v(i,jm,k ) )
          uwkp  = 0.25d0*( u(i,j,kp)+u(i,j,k) ) !*( w(ip,j ,k )+w(i,j ,k ) )
          uwkm  = 0.25d0*( u(i,j,km)+u(i,j,k) ) !*( w(ip,j ,km)+w(i,j ,km) )
          uvip  = 0.25d0*( u(i ,j,k)+u(i ,jp,k) ) !*( v(i,j,k )+v(ip,j ,k) )
          uvim  = 0.25d0*( u(im,j,k)+u(im,jp,k) ) !*( v(i,j,k )+v(im,j ,k) )
          uwip  = 0.25d0*( u(i ,j ,k)+u(i ,j ,kp) ) !*( w(i,j,k)+w(ip,j,k) )
          uwim  = 0.25d0*( u(im,j ,k)+u(im,j ,kp) ) !*( w(i,j,k)+w(im,j,k) )

          ! touch v
          vvjp  = 0.25d0*( v(i,j,k )+v(i,jp,k)  )*( v(i,j,k )+v(i ,jp,k) )
          vvjm  = 0.25d0*( v(i,j,k )+v(i,jm,k)  )*( v(i,j,k )+v(i ,jm,k) )
          dvdxp = (v(ip,j,k)-v(i,j,k))*dxi
          dvdxm = (v(i,j,k)-v(im,j,k))*dxi
          dvdyp = (v(i,jp,k)-v(i,j,k))*dyi
          dvdym = (v(i,j,k)-v(i,jm,k))*dyi
          dvdzp = (v(i,j,kp)-v(i,j,k))*dzci(k)
          dvdzm = (v(i,j,k)-v(i,j,km))*dzci(km)
          uvip  = uvip*( v(i,j,k )+v(ip,j ,k) )
          uvim  = uvim*( v(i,j,k )+v(im,j ,k) )


          dvdtrk_temp =   dxi*(     -uvip + uvim ) + (dvdxp-dvdxm)*visc*dxi+ &
                        dyi*(     -vvjp + vvjm ) + (dvdyp-dvdym)*visc*dyi+ &
                                                   (dvdzp-dvdzm)*visc*dzfi(k)

          uvjp  = uvjp*( v(ip,j ,k )+v(i,j ,k ) )
          uvjm  = uvjm*( v(ip,jm,k )+v(i,jm,k ) )
          dudtrk_temp = dudtrk_temp + dyi*(     -uvjp + uvjm )

          wvkp  = 0.25d0*( v(i ,j ,kp)+v(i ,j ,k ) ) !*( w(i,j,k )+w(i,jp,k ) )
          wvkm  = 0.25d0*( v(i ,j ,km)+v(i ,j ,k ) ) !*( w(i,j,km)+w(i,jp,km) )
          vwjp  = 0.25d0*( v(i ,j ,k )+v(i ,j ,kp) ) !*( w(i,j,k )+w(i,jp,k ) )
          vwjm  = 0.25d0*( v(i ,jm,k )+v(i ,jm,kp) ) !*( w(i,j,k )+w(i,jm,k ) )

          !touch w
          wwkp  = 0.25d0*( w(i,j,k)+w(i,j,kp) )*( w(i ,j ,k)+w(i ,j ,kp) )
          wwkm  = 0.25d0*( w(i,j,k)+w(i,j,km) )*( w(i ,j ,k)+w(i ,j ,km) )
          dwdxp = (w(ip,j,k)-w(i,j,k))*dxi
          dwdxm = (w(i,j,k)-w(im,j,k))*dxi
          dwdyp = (w(i,jp,k)-w(i,j,k))*dyi
          dwdym = (w(i,j,k)-w(i,jm,k))*dyi
          dwdzp = (w(i,j,kp)-w(i,j,k))*dzfi(kp)
          dwdzm = (w(i,j,k)-w(i,j,km))*dzfi(k)
          uwip  = uwip*( w(i,j,k)+w(ip,j,k) )
          uwim  = uwim*( w(i,j,k)+w(im,j,k) )
          vwjp  = vwjp*( w(i,j,k )+w(i,jp,k ) )
          vwjm  = vwjm*( w(i,j,k )+w(i,jm,k ) )


          dwdtrk_temp =   dxi*(     -uwip + uwim ) + (dwdxp-dwdxm)*visc*dxi+ &
                          dyi*(     -vwjp + vwjm ) + (dwdyp-dwdym)*visc*dyi+ &
                          dzci(k)*( -wwkp + wwkm ) + (dwdzp-dwdzm)*visc*dzci(k)

          uwkp  = uwkp*( w(ip,j ,k )+w(i,j ,k ) )
          uwkm  = uwkm*( w(ip,j ,km)+w(i,j ,km) )
          dudtrk_temp = dudtrk_temp + dzfi(k)*( -uwkp + uwkm )

          wvkp  = wvkp*( w(i,j,k )+w(i,jp,k ) )
          wvkm  = wvkm*( w(i,j,km)+w(i,jp,km) )
          dvdtrk_temp = dvdtrk_temp + dzfi(k)*( -wvkp + wvkm )

          !
          ! Momentum balance
          !
          dudtrk(i,j,k) = dudtrk_temp
          dvdtrk(i,j,k) = dvdtrk_temp
          dwdtrk(i,j,k) = dwdtrk_temp

        enddo
      enddo
    enddo
#ifndef USE_CUDA
    !$OMP END PARALLEL DO
#endif
    !@cuf istat=cudaDeviceSynchronize()
    taux(:) = 0.
    taux2d=0.d0
    taux3d=0.d0
#ifdef USE_CUDA
    !$cuf kernel do(2) <<<*,*>>>
#endif
    do k=1,nz
      do i=1,nx
        dudyp = (u(i,1 ,k)-u(i,0   ,k))*dyi*visc*dzflzi(k)
        dudym = (u(i,ny,k)-u(i,ny+1,k))*dyi*visc*dzflzi(k)
        taux2d = taux2d + (dudyp+dudym)
      enddo
    enddo
    !@cuf istat=cudaDeviceSynchronize()
    taux(2)=taux2d
#ifdef USE_CUDA
    !$cuf kernel do(2) <<<*,*>>>
#endif
    do j=1,ny
      do i=1,nx
        dudzp = (u(i,j,1 )-u(i,j,0   ))*dzci(0)*visc
        dudzm = (u(i,j,nz)-u(i,j,nz+1))*dzci(nz)*visc
        taux3d = taux3d + (dudzp+dudzm)
      enddo
    enddo
    !@cuf istat=cudaDeviceSynchronize()
    taux(3)=taux3d

    call mpi_allreduce(MPI_IN_PLACE,taux(1),3,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    nxg = nx*dims(1)
    nyg = ny*dims(2)
    nzg = nz
    taux(1) = taux(1)/(1.d0*nyg)
    taux(2) = taux(2)/(1.d0*nxg)
    taux(3) = taux(3)/(1.d0*nxg*nyg)

    tauy(:) = 0.
    tauy1d=0.d0
    tauy3d=0.d0
#ifdef USE_CUDA
    !$cuf kernel do(2) <<<*,*>>>
#endif
    do k=1,nz
      do j=1,ny
        dvdxp = (v(1 ,j,k)-v(0   ,j,k))*dxi*visc*dzflzi(k)
        dvdxm = (v(nx,j,k)-v(nx+1,j,k))*dxi*visc*dzflzi(k)
        tauy1d = tauy1d + (dvdxp+dvdxm)
      enddo
    enddo
    !@cuf istat=cudaDeviceSynchronize()
    tauy(1)=tauy1d
#ifdef USE_CUDA
    !$cuf kernel do(2) <<<*,*>>>
#endif
    do j=1,ny
      do i=1,nx
        dvdzp = (v(i,j,1 )-v(i,j,0   ))*dzci(0)*visc
        dvdzm = (v(i,j,nz)-v(i,j,nz+1))*dzci(nz)*visc
        tauy3d = tauy3d + (dvdzp+dvdzm)
      enddo
    enddo
    !@cuf istat=cudaDeviceSynchronize()
    tauy(3)=tauy3d
    call mpi_allreduce(MPI_IN_PLACE,tauy(1),3,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    nxg = nx*dims(1)
    nyg = ny*dims(2)
    nzg = nz
    tauy(1) = tauy(1)/(1.d0*nyg)
    tauy(2) = tauy(2)/(1.d0*nxg)
    tauy(3) = tauy(3)/(1.d0*nxg*nyg)

    tauz(:) = 0.
    tauz1d = 0.d0
    tauz2d = 0.d0
#ifdef USE_CUDA
    !$cuf kernel do(2) <<<*,*>>>
#endif
    do k=1,nz
      do j=1,ny
        dwdxp = (w(1 ,j,k)-w(0   ,j,k))*dxi*visc*dzflzi(k)
        dwdxm = (w(nx,j,k)-w(nx+1,j,k))*dxi*visc*dzflzi(k)
        tauz1d = tauz1d + (dwdxp+dwdxm)
      enddo
    enddo
    !@cuf istat=cudaDeviceSynchronize()
    tauz(1)=tauz1d
#ifdef USE_CUDA
    !$cuf kernel do(2) <<<*,*>>>
#endif
    do k=1,nz
      do i=1,nx
        dwdyp = (w(i,1,k )-w(i,0   ,k))*dyi*visc*dzflzi(k)
        dwdym = (w(i,ny,k)-w(i,ny+1,k))*dyi*visc*dzflzi(k)
        tauz2d = tauz2d + (dwdyp+dwdym)
      enddo
    enddo
    !@cuf istat=cudaDeviceSynchronize()
    tauz(2)=tauz2d
    call mpi_allreduce(MPI_IN_PLACE,tauz(1),3,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    nxg = nx*dims(1)
    nyg = ny*dims(2)
    nzg = nz
    tauz(1) = tauz(1)/(1.d0*nyg)
    tauz(2) = tauz(2)/(1.d0*nxg)
    tauz(3) = tauz(3)/(1.d0*nxg*nyg)

#ifdef USE_NVTX
      call nvtxEndRange
#endif

    f(1) = (factor1*sum(taux(:)/l(:)) + factor2*sum(tauxo(:)/l(:)))
    f(2) = (factor1*sum(tauy(:)/l(:)) + factor2*sum(tauyo(:)/l(:)))
    f(3) = (factor1*sum(tauz(:)/l(:)) + factor2*sum(tauzo(:)/l(:)))
    tauxo(:) = taux(:)
    tauyo(:) = tauy(:)
    tauzo(:) = tauz(:)
#ifdef USE_NVTX
    call nvtxStartRange("update_cuf2",1)
#endif
#ifdef USE_CUDA
    !$cuf kernel do(3) <<<*,*>>>
#else
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(n,factor1,factor2,up,vp,wp,dudtrk,dvdtrk,dwdtrk,dudtrko,dvdtrko,dwdtrko)
#endif
    do k=1,n(3)
      kp = k + 1
      do j=1,n(2)
        jp = j + 1
        do i=1,n(1)
          ip = i + 1
          ! could be split in two loops, because factor2=0 for istep=1, but like this reads nicer
          !up(i,j,k) = up(i,j,k) + factor1*dudtrk(i,j,k) + factor2*dudtrko(i,j,k)
          !vp(i,j,k) = vp(i,j,k) + factor1*dvdtrk(i,j,k) + factor2*dvdtrko(i,j,k)
          !wp(i,j,k) = wp(i,j,k) + factor1*dwdtrk(i,j,k) + factor2*dwdtrko(i,j,k)

          ! momp and update_cuf1 in place:
          up(i,j,k) = u(i,j,k) + factor12*( - dxi*    ( p(ip,j,k)-p(i,j,k) )  ) + factor1*dudtrk(i,j,k) + factor2*dudtrko(i,j,k)
          vp(i,j,k) = v(i,j,k) + factor12*( - dyi*    ( p(i,jp,k)-p(i,j,k) )  ) + factor1*dvdtrk(i,j,k) + factor2*dvdtrko(i,j,k)
          wp(i,j,k) = w(i,j,k) + factor12*( - dzci(k)*( p(i,j,kp)-p(i,j,k) )  ) + factor1*dwdtrk(i,j,k) + factor2*dwdtrko(i,j,k)

          ! d*dtrk buffers swapped with pointers in main
          !dudtrko(i,j,k) = dudtrk(i,j,k)
          !dvdtrko(i,j,k) = dvdtrk(i,j,k)
          !dwdtrko(i,j,k) = dwdtrk(i,j,k)
        enddo
      enddo
    enddo
#ifndef USE_CUDA
    !$OMP END PARALLEL DO
#endif

 #ifdef USE_NVTX
      call nvtxEndRange
      call nvtxStartRange("chkmean",1)
 #endif
    !
    ! bulk velocity forcing
    !
    f(:) = 0.d0
    if(is_forced(1)) then
      call chkmean(n,dzflzi,up,mean)
      f(1) = velf(1) - mean
    endif
    if(is_forced(2)) then
      call chkmean(n,dzflzi,vp,mean)
      f(2) = velf(2) - mean
    endif
    if(is_forced(3)) then
      call chkmean(n,dzclzi,wp,mean)
      f(3) = velf(3) - mean
    endif
 #ifdef USE_NVTX
      call nvtxEndRange
 #endif
    return
  end subroutine rk
  subroutine rk_id(rkpar,n,dli,dzci,dzfi,dzflzi,dzclzi,visc,dt,l,u,v,w,p,& 
                   dudtrk,dvdtrk,dwdtrk,dudtrkd,dvdtrkd,dwdtrkd, &
                   dudtrko,dvdtrko,dwdtrko,tauxo,tauyo,tauzo,up,vp,wp,f)
    !
    ! low-storage 3rd-order Runge-Kutta scheme 
    ! for time integration of the momentum equations with implicit diffusion.
    !
    implicit none
    real(8), intent(in), dimension(2) :: rkpar
    integer, intent(in), dimension(3) :: n
    real(8), intent(in) :: visc,dt
    real(8), intent(in   ), dimension(3) :: dli,l 
    real(8), intent(in   ), dimension(0:) :: dzci,dzfi,dzflzi,dzclzi
    real(8), intent(in   ), dimension(0:,0:,0:) :: u ,v ,w,p
    real(8), intent(inout), dimension(:,:,:) :: dudtrko,dvdtrko,dwdtrko
    real(8), intent(inout), dimension(3) :: tauxo,tauyo,tauzo
    real(8), intent(out), dimension(0:,0:,0:) :: up,vp,wp
    real(8), intent(out), dimension(3) :: f
    real(8),              dimension(:,:,:) ::          dudtrk, dvdtrk, dwdtrk
    real(8),              dimension(:,:,:) ::          dudtrkd, dvdtrkd, dwdtrkd
    real(8) :: factor1,factor2,factor12
    real(8), dimension(3) :: taux,tauy,tauz
#ifdef USE_CUDA
   attributes(managed):: dzci,dzfi,dzflzi,dzclzi,u ,v ,w,p, dudtrk, dvdtrk, dwdtrk, up,vp,wp
   attributes(managed):: dudtrko,dvdtrko,dwdtrko,dudtrkd, dvdtrkd, dwdtrkd
   integer:: istat
#endif
    integer :: i,j,k
    real(8) :: mean
    !
    factor1 = rkpar(1)*dt
    factor2 = rkpar(2)*dt
    factor12 = factor1 + factor2
    !
#ifdef USE_NVTX
      call nvtxStartRange("momxpd",1)
#endif
    call momxpd(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,dzflzi,visc,p,u,dudtrk,dudtrkd,taux)
#ifdef USE_NVTX
      call nvtxEndRange
      call nvtxStartRange("momypd",2)
#endif
    call momypd(n(1),n(2),n(3),dli(2),dli(2),dzci,dzfi,dzflzi,visc,p,v,dvdtrk,dvdtrkd,tauy)
#ifdef USE_NVTX
      call nvtxEndRange
      call nvtxStartRange("momyzd",2)
#endif
    call momzpd(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,dzflzi,visc,p,w,dwdtrk,dwdtrkd,tauz)
    f(1) = factor12*sum(taux(:)/l(:))
    f(2) = factor12*sum(tauy(:)/l(:))
    f(3) = factor12*sum(tauz(:)/l(:))
#ifdef USE_NVTX
      call nvtxEndRange
#endif
    ! alternatively, calculate force from the mean velocity directly
#ifdef USE_CUDA
    !$cuf kernel do(3) <<<*,*>>>
#else
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(n,factor12,u,v,w,up,vp,wp,dudtrk,dvdtrk,dwdtrk)
#endif
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          up(i,j,k) = u(i,j,k) + factor12*dudtrk(i,j,k)
          vp(i,j,k) = v(i,j,k) + factor12*dvdtrk(i,j,k)
          wp(i,j,k) = w(i,j,k) + factor12*dwdtrk(i,j,k)
        enddo
      enddo
    enddo
#ifndef USE_CUDA
    !$OMP END PARALLEL DO
#endif
    !@cuf istat=cudaDeviceSynchronize()

#ifdef USE_NVTX
      call nvtxStartRange("momxa",1)
#endif
    call momxa(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,u,v,w,dudtrk)
#ifdef USE_NVTX
      call nvtxEndRange
      call nvtxStartRange("momya",2)
#endif
    call momya(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,u,v,w,dvdtrk)
#ifdef USE_NVTX
      call nvtxEndRange
      call nvtxStartRange("momza",2)
#endif
    call momza(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,u,v,w,dwdtrk)
#ifdef USE_NVTX
      call nvtxEndRange
#endif

#ifdef USE_CUDA
    !$cuf kernel do(3) <<<*,*>>>
#else
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(n,factor1,factor2,up,vp,wp,dudtrk,dvdtrk,dwdtrk,dudtrko,dvdtrko,dwdtrko)
#endif
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          ! could be split in two loops, because factor2=0 for istep=1, but like this reads nicer
          up(i,j,k) = up(i,j,k) + factor1*dudtrk(i,j,k) + factor2*dudtrko(i,j,k)
          vp(i,j,k) = vp(i,j,k) + factor1*dvdtrk(i,j,k) + factor2*dvdtrko(i,j,k)
          wp(i,j,k) = wp(i,j,k) + factor1*dwdtrk(i,j,k) + factor2*dwdtrko(i,j,k)
          dudtrko(i,j,k) = dudtrk(i,j,k)
          dvdtrko(i,j,k) = dvdtrk(i,j,k)
          dwdtrko(i,j,k) = dwdtrk(i,j,k)
        enddo
      enddo
    enddo
#ifndef USE_CUDA
    !$OMP END PARALLEL DO
#endif
    !@cuf istat=cudaDeviceSynchronize()
    !
    ! bulk velocity forcing
    !
    f(:) = 0.d0
    if(is_forced(1)) then
      call chkmean(n,dzflzi,up,mean)
      f(1) = velf(1) - mean
    endif
    if(is_forced(2)) then
      call chkmean(n,dzflzi,vp,mean)
      f(2) = velf(2) - mean
    endif
    if(is_forced(3)) then
      call chkmean(n,dzclzi,wp,mean)
      f(3) = velf(3) - mean
    endif
    !
    ! compute rhs of helmholtz equation
    !
#ifdef USE_CUDA
    !$cuf kernel do(3) <<<*,*>>>
#else
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(n,factor12,factor2,visc,up,vp,wp,dudtrkd,dvdtrkd,dwdtrkd)
#endif
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          up(i,j,k) = up(i,j,k) - .5d0*factor12*dudtrkd(i,j,k)
          vp(i,j,k) = vp(i,j,k) - .5d0*factor12*dvdtrkd(i,j,k)
          wp(i,j,k) = wp(i,j,k) - .5d0*factor12*dwdtrkd(i,j,k)
        enddo
      enddo
    enddo
#ifndef USE_CUDA
    !$OMP END PARALLEL DO
#endif
    !@cuf istat=cudaDeviceSynchronize()
    return
  end subroutine rk_id
  subroutine rk_scal(rkpar,n,dli,dzci,dzfi,dzflzi,dzclzi,visc,dt,u,v,w,dsdtrko,s)
    !
    ! low-storage 3rd-order Runge-Kutta scheme 
    ! for time integration of the scalar field.
    !
    implicit none
    real(8), intent(in   ), dimension(2) :: rkpar
    integer, intent(in   ), dimension(3) :: n
    real(8), intent(in   ), dimension(3) :: dli
    real(8), intent(in   ), dimension(0:) :: dzci,dzfi,dzflzi,dzclzi
    real(8), intent(in   ) :: visc,dt
    real(8), intent(in   ), dimension(0:,0:,0:) :: u,v,w
    real(8), intent(inout), dimension(:,:,:) :: dsdtrko
    real(8), intent(inout), dimension(0:,0:,0:) :: s
    real(8),              dimension(n(1),n(2),n(3)) :: dsdtrk
    real(8) :: factor1,factor2,factor12
    integer :: i,j,k
    !
    factor1 = rkpar(1)*dt
    factor2 = rkpar(2)*dt
    factor12 = factor1 + factor2
    call momsad(n(1),n(2),n(3),dli(1),dli(2),dli(3),dzci,dzfi,dzflzi,visc,u,v,w,s,dsdtrk)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(n,factor1,factor2,s,dsdtrk,dsdtrko)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          s(i,j,k) = s(i,j,k) + factor1*dsdtrk(i,j,k) + factor2*dsdtrko(i,j,k)
          dsdtrko(i,j,k) = dsdtrk(i,j,k)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    return
  end subroutine rk_scal
end module mod_rk
