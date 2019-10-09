module mod_momd
  use mpi
  use mod_param     , only: dims
  use mod_common_mpi, only: ierr
  use mod_types
  !@cuf use cudafor
  implicit none
  private
  public momxa,momya,momza,momxpd,momypd,momzpd
  contains
  subroutine momxa(nx,ny,nz,dxi,dyi,dzci,dzfi,u,v,w,dudt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    real(rp), dimension(0:,0:,0:), intent(in) :: u,v,w
    real(rp), dimension(:,:,:), intent(out) :: dudt
    integer :: im,ip,jm,jp,km,kp,i,j,k
    real(rp) :: uuip,uuim,uvjp,uvjm,uwkp,uwkm
#ifdef USE_CUDA
    attributes(managed) :: u,v,w,dudt,dzci,dzfi
    integer :: istat
#endif
    !
#ifdef USE_CUDA
    !$cuf kernel do(3) <<<*,*>>>
#else
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k,im,jm,km,ip,jp,kp) &
    !$OMP PRIVATE(uuip,uuim,uvjp,uvjm,uwkp,uwkm) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,u,v,w,dudt,dzci,dzfi)
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
          uuip  = 0.25*( u(ip,j,k)+u(i,j,k) )*( u(ip,j ,k )+u(i,j ,k ) )
          uuim  = 0.25*( u(im,j,k)+u(i,j,k) )*( u(im,j ,k )+u(i,j ,k ) )
          uvjp  = 0.25*( u(i,jp,k)+u(i,j,k) )*( v(ip,j ,k )+v(i,j ,k ) )
          uvjm  = 0.25*( u(i,jm,k)+u(i,j,k) )*( v(ip,jm,k )+v(i,jm,k ) )
          uwkp  = 0.25*( u(i,j,kp)+u(i,j,k) )*( w(ip,j ,k )+w(i,j ,k ) )
          uwkm  = 0.25*( u(i,j,km)+u(i,j,k) )*( w(ip,j ,km)+w(i,j ,km) )
          !
          ! momentum balance
          !
          dudt(i,j,k) = dxi*(     -uuip + uuim ) + &
                        dyi*(     -uvjp + uvjm ) + &
                        dzfi(k)*( -uwkp + uwkm )
        enddo
      enddo
    enddo
#ifndef USE_CUDA
    !$OMP END PARALLEL DO
#endif
    !@cuf istat=cudaDeviceSynchronize()
    return
  end subroutine momxa
  !
  subroutine momya(nx,ny,nz,dxi,dyi,dzci,dzfi,u,v,w,dvdt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    real(rp), dimension(0:,0:,0:), intent(in) :: u,v,w
    real(rp), dimension(:,:,:), intent(out) :: dvdt
    integer :: im,ip,jm,jp,km,kp,i,j,k
    real(rp) :: uvip,uvim,vvjp,vvjm,wvkp,wvkm
#ifdef USE_CUDA
    attributes(managed) :: u,v,w,dvdt,dzci,dzfi
    integer :: istat
#endif
    !
#ifdef USE_CUDA
    !$cuf kernel do(3) <<<*,*>>>
#else
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k,im,jm,km,ip,jp,kp) &
    !$OMP PRIVATE(uvip,uvim,vvjp,vvjm,wvkp,wvkm) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,dzci,dzfi,u,v,w,dvdt)
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
          uvip  = 0.25*( u(i ,j,k)+u(i ,jp,k) )*( v(i,j,k )+v(ip,j ,k) )
          uvim  = 0.25*( u(im,j,k)+u(im,jp,k) )*( v(i,j,k )+v(im,j ,k) )
          vvjp  = 0.25*( v(i,j,k )+v(i,jp,k)  )*( v(i,j,k )+v(i ,jp,k) )
          vvjm  = 0.25*( v(i,j,k )+v(i,jm,k)  )*( v(i,j,k )+v(i ,jm,k) )
          wvkp  = 0.25*( w(i,j,k )+w(i,jp,k)  )*( v(i,j,kp)+v(i ,j ,k) )
          wvkm  = 0.25*( w(i,j,km)+w(i,jp,km) )*( v(i,j,km)+v(i ,j ,k) )
          !
          ! momentum balance
          !
          dvdt(i,j,k) = dxi*(     -uvip + uvim ) + &
                        dyi*(     -vvjp + vvjm ) + &
                        dzfi(k)*( -wvkp + wvkm )
        enddo
      enddo
    enddo
#ifndef USE_CUDA
    !$OMP END PARALLEL DO
#endif
    !@cuf istat=cudaDeviceSynchronize()
    return
  end subroutine momya
  !
  subroutine momza(nx,ny,nz,dxi,dyi,dzci,dzfi,u,v,w,dwdt)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    real(rp), dimension(0:,0:,0:), intent(in) :: u,v,w
    real(rp), dimension(:,:,:), intent(out) :: dwdt
    integer :: im,ip,jm,jp,km,kp,i,j,k
    real(rp) :: uwip,uwim,vwjp,vwjm,wwkp,wwkm
#ifdef USE_CUDA
    attributes(managed) :: u,v,w,dwdt,dzci,dzfi
    integer :: istat
#endif
    !
#ifdef USE_CUDA
    !$cuf kernel do(3) <<<*,*>>>
#else
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k,im,jm,km,ip,jp,kp) &
    !$OMP PRIVATE(uwip,uwim,vwjp,vwjm,wwkp,wwkm) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,dzci,dzfi,u,v,w,dwdt)
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
          uwip  = 0.25*( w(i,j,k)+w(ip,j,k) )*( u(i ,j ,k)+u(i ,j ,kp) )
          uwim  = 0.25*( w(i,j,k)+w(im,j,k) )*( u(im,j ,k)+u(im,j ,kp) )
          vwjp  = 0.25*( w(i,j,k)+w(i,jp,k) )*( v(i ,j ,k)+v(i ,j ,kp) )
          vwjm  = 0.25*( w(i,j,k)+w(i,jm,k) )*( v(i ,jm,k)+v(i ,jm,kp) )
          wwkp  = 0.25*( w(i,j,k)+w(i,j,kp) )*( w(i ,j ,k)+w(i ,j ,kp) )
          wwkm  = 0.25*( w(i,j,k)+w(i,j,km) )*( w(i ,j ,k)+w(i ,j ,km) )
          !
          ! momentum balance
          !
          dwdt(i,j,k) = dxi*(     -uwip + uwim ) + &
                        dyi*(     -vwjp + vwjm ) + &
                        dzci(k)*( -wwkp + wwkm )
        enddo
      enddo
    enddo
#ifndef USE_CUDA
    !$OMP END PARALLEL DO
#endif
    !@cuf istat=cudaDeviceSynchronize()
    return
  end subroutine momza
  !
  subroutine momxpd(nx,ny,nz,dxi,dyi,dzci,dzfi,dzflzi,visc,p,u,dudt,dudtd,taux)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi,visc
    real(rp), intent(in), dimension(0:) :: dzci,dzfi,dzflzi
    real(rp), dimension(0:,0:,0:), intent(in) :: p,u
    real(rp), dimension(:,:,:), intent(out) :: dudt,dudtd
    real(rp), dimension(3)  , intent(out) :: taux
    real(rp) :: dudxp,dudxm,dudyp,dudym,dudzp,dudzm
    integer  :: im,ip,jm,jp,km,kp,i,j,k
    integer  :: nxg,nyg,nzg
    real(rp) :: taux2d,taux3d
#ifdef USE_CUDA
    attributes(managed) :: u,p,dudt,dudtd,dzci,dzfi,dzflzi
    integer :: istat
#endif
 
#ifdef USE_CUDA
    !$cuf kernel do(3) <<<*,*>>>
#else
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k,im,jm,km,ip,jp,kp) &
    !$OMP PRIVATE(dudxp,dudxm,dudyp,dudym,dudzp,dudzm) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,dzci,dzfi,u,p,dudt,dudtd,visc)
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
          dudxp = (u(ip,j,k)-u(i,j,k))*dxi
          dudxm = (u(i,j,k)-u(im,j,k))*dxi
          dudyp = (u(i,jp,k)-u(i,j,k))*dyi
          dudym = (u(i,j,k)-u(i,jm,k))*dyi
          dudzp = (u(i,j,kp)-u(i,j,k))*dzci(k)
          dudzm = (u(i,j,k)-u(i,j,km))*dzci(km)
          dudtd(i,j,k) = (dudxp-dudxm)*visc*dxi + &
                         (dudyp-dudym)*visc*dyi + &
                         (dudzp-dudzm)*visc*dzfi(k)
          dudt(i,j,k) = - dxi*( p(ip,j,k)-p(i,j,k) ) + dudtd(i,j,k)
        enddo
      enddo
    enddo
#ifndef USE_CUDA
    !$OMP END PARALLEL DO
#endif
    !@cuf istat=cudaDeviceSynchronize()
    taux(:) = 0.
    taux2d  = 0.
    taux3d  = 0.
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
    taux(2) = taux2d
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
    taux(3) = taux3d
    !@cuf istat=cudaDeviceSynchronize()
    call mpi_allreduce(MPI_IN_PLACE,taux(1),3,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    nxg = nx*dims(1)
    nyg = ny*dims(2)
    nzg = nz
    taux(1) = taux(1)/(1.*nyg)
    taux(2) = taux(2)/(1.*nxg)
    taux(3) = taux(3)/(1.*nxg*nyg)
    return
  end subroutine momxpd
  !
  subroutine momypd(nx,ny,nz,dxi,dyi,dzci,dzfi,dzflzi,visc,p,v,dvdt,dvdtd,tauy)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi,visc
    real(rp), intent(in), dimension(0:) :: dzci,dzfi,dzflzi
    real(rp), dimension(0:,0:,0:), intent(in) :: p,v
    real(rp), dimension(:,:,:), intent(out) :: dvdt,dvdtd
    real(rp), dimension(3)  , intent(out) :: tauy
    real(rp) :: dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm
    integer  :: im,ip,jm,jp,km,kp,i,j,k
    integer  :: nxg,nyg,nzg
    real(rp) :: tauy1d,tauy3d
#ifdef USE_CUDA
    attributes(managed) :: v,p,dvdt,dvdtd,dzci,dzfi,dzflzi
    integer :: istat
#endif
    !
#ifdef USE_CUDA
    !$cuf kernel do(3) <<<*,*>>>
#else
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k,im,jm,km,ip,jp,kp) &
    !$OMP PRIVATE(dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,dzci,dzfi,v,p,dvdt,dvdtd,visc)
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
          dvdxp = (v(ip,j,k)-v(i,j,k))*dxi
          dvdxm = (v(i,j,k)-v(im,j,k))*dxi
          dvdyp = (v(i,jp,k)-v(i,j,k))*dyi
          dvdym = (v(i,j,k)-v(i,jm,k))*dyi
          dvdzp = (v(i,j,kp)-v(i,j,k))*dzci(k)
          dvdzm = (v(i,j,k)-v(i,j,km))*dzci(km)
          dvdtd(i,j,k) = (dvdxp-dvdxm)*visc*dxi+ &
                         (dvdyp-dvdym)*visc*dyi+ &
                         (dvdzp-dvdzm)*visc*dzfi(k)
          dvdt(i,j,k) = - dyi*( p(i,jp,k)-p(i,j,k) ) + dvdtd(i,j,k) 
        enddo
      enddo
    enddo
#ifndef USE_CUDA
    !$OMP END PARALLEL DO
#endif
    !@cuf istat=cudaDeviceSynchronize()
    tauy(:) = 0.
    tauy1d  = 0.
    tauy3d  = 0.
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
    tauy(1) = tauy1d
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
    tauy(3) = tauy3d
    !@cuf istat=cudaDeviceSynchronize()
    call mpi_allreduce(MPI_IN_PLACE,tauy(1),3,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    nxg = nx*dims(1)
    nyg = ny*dims(2)
    nzg = nz
    tauy(1) = tauy(1)/(1.*nyg)
    tauy(2) = tauy(2)/(1.*nxg)
    tauy(3) = tauy(3)/(1.*nxg*nyg)
    return
  end subroutine momypd
  !
  subroutine momzpd(nx,ny,nz,dxi,dyi,dzci,dzfi,dzflzi,visc,p,w,dwdt,dwdtd,tauz)
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi,visc
    real(rp), intent(in), dimension(0:) :: dzci,dzfi,dzflzi
    real(rp), dimension(0:,0:,0:), intent(in) :: p,w
    real(rp), dimension(:,:,:), intent(out) :: dwdt,dwdtd
    real(rp), dimension(3)  , intent(out) :: tauz
    real(rp) :: dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm
    integer  :: im,ip,jm,jp,km,kp,i,j,k
    integer  :: nxg,nyg,nzg
    real(rp) :: tauz1d,tauz2d
#ifdef USE_CUDA
    attributes(managed) :: w,p,dwdt,dwdtd,dzci,dzfi,dzflzi
    integer :: istat
#endif
    !
#ifdef USE_CUDA
    !$cuf kernel do(3) <<<*,*>>>
#else
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k,im,jm,km,ip,jp,kp) &
    !$OMP PRIVATE(dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,dzci,dzfi,w,p,dwdt,dwdtd,visc)
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
          dwdxp = (w(ip,j,k)-w(i,j,k))*dxi
          dwdxm = (w(i,j,k)-w(im,j,k))*dxi
          dwdyp = (w(i,jp,k)-w(i,j,k))*dyi
          dwdym = (w(i,j,k)-w(i,jm,k))*dyi
          dwdzp = (w(i,j,kp)-w(i,j,k))*dzfi(kp)
          dwdzm = (w(i,j,k)-w(i,j,km))*dzfi(k)
          dwdtd(i,j,k) = (dwdxp-dwdxm)*visc*dxi+ &
                         (dwdyp-dwdym)*visc*dyi+ &
                         (dwdzp-dwdzm)*visc*dzci(k)
          dwdt(i,j,k) = - dzci(k)*( p(i,j,kp)-p(i,j,k) ) + dwdtd(i,j,k)
        enddo
      enddo
    enddo
#ifndef USE_CUDA
    !$OMP END PARALLEL DO
#endif
    !@cuf istat=cudaDeviceSynchronize()
    tauz(:) = 0.
    tauz1d  = 0.
    tauz2d  = 0.
#ifdef USE_CUDA
    !$cuf kernel do(2) <<<*,*>>>
#endif
    do k=1,nz
      do j=1,ny
        dwdxp = (w(1 ,j,k)-w(0   ,j,k))*dxi*visc*dzflzi(k)
        dwdxm = (w(nx,j,k)-w(nx+1,j,k))*dxi*visc*dzflzi(k)
        tauz1d = tauz1d+ (dwdxp+dwdxm)
      enddo
    enddo
    tauz(1) = tauz1d
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
    tauz(2) = tauz2d
    !@cuf istat=cudaDeviceSynchronize()
    call mpi_allreduce(MPI_IN_PLACE,tauz(1),3,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    nxg = nx*dims(1)
    nyg = ny*dims(2)
    nzg = nz
    tauz(1) = tauz(1)/(1.*nyg)
    tauz(2) = tauz(2)/(1.*nxg)
    tauz(3) = tauz(3)/(1.*nxg*nyg)
    return
  end subroutine momzpd
end module mod_momd
