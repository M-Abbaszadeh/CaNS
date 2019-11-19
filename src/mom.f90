module mod_mom
  use mpi
  use mod_param     , only: dims, bforce
  use mod_common_mpi, only: ierr
  use mod_types
  implicit none
  private
  public momxad,momyad,momzad,momxp,momyp,momzp,momp,momxyzad
  contains
  subroutine momxyzad(nx,ny,nz,dxi,dyi,dzi,dzci,dzfi,dzflzi,visc,u,v,w,dudt,dvdt,dwdt,taux,tauy,tauz)
    !@cuf use cudafor
    implicit none
    integer, intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi,dzi,visc
    real(rp), intent(in), dimension(0:) :: dzci,dzfi,dzflzi
    real(rp), dimension(0:,0:,0:), intent(in) :: u,v,w
    real(rp), dimension(:,:,:), intent(out) :: dudt,dvdt,dwdt
    real(rp), dimension(3)  , intent(out) :: taux,tauy,tauz
    integer :: im,ip,jm,jp,km,kp,i,j,k
    real(rp) :: uuip,uuim,uvjp,uvjm,uwkp,uwkm
    real(rp) :: dudxp,dudxm,dudyp,dudym,dudzp,dudzm
    real(rp) :: taux2d,taux3d
    real(rp) :: uvip,uvim,vvjp,vvjm,wvkp,wvkm
    real(rp) :: dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm
    real(rp) :: bforcex,bforcey,bforcez
    real(rp) :: tauy1d,tauy3d
    real(rp) :: uwip,uwim,vwjp,vwjm,wwkp,wwkm
    real(rp) :: dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm
    real(rp) :: tauz1d,tauz2d,dudt_temp,dvdt_temp,dwdt_temp
#ifdef USE_CUDA
    attributes(managed) :: u,v,w,dudt,dvdt,dwdt,dzci,dzfi,dzflzi
    integer :: istat
#endif
    integer :: nxg,nyg,nzg
    !
    bforcex = bforce(1)
    bforcey = bforce(2)
    bforcez = bforce(3)
#ifdef USE_CUDA
    !$cuf kernel do(3) <<<*,(8,8,8)>>>
#else
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k,im,jm,km,ip,jp,kp) &
    !$OMP PRIVATE(uuip,uuim,uvjp,uvjm,uwkp,uwkm) &
    !$OMP PRIVATE(dudxp,dudxm,dudyp,dudym,dudzp,dudzm) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,dzi,visc,u,v,w,dudt,dzci,dzfi,bforcex,bforcey,bforcez)
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
          !
          ! touch u
          !
          dudxp = (u(ip,j,k)-u(i,j,k))*dxi
          dudxm = (u(i,j,k)-u(im,j,k))*dxi
          dudyp = (u(i,jp,k)-u(i,j,k))*dyi
          dudym = (u(i,j,k)-u(i,jm,k))*dyi
          dudzp = (u(i,j,kp)-u(i,j,k))*dzci(k)
          dudzm = (u(i,j,k)-u(i,j,km))*dzci(km)
          uuip  = 0.25*( u(ip,j,k)+u(i,j,k) )*( u(ip,j ,k )+u(i,j ,k ) )
          uuim  = 0.25*( u(im,j,k)+u(i,j,k) )*( u(im,j ,k )+u(i,j ,k ) )
          !
          dudt_temp = dxi*(     -uuip + uuim ) + (dudxp-dudxm)*visc*dxi + &
                                               + (dudyp-dudym)*visc*dyi + &
                                               + (dudzp-dudzm)*visc*dzfi(k)
          !
          uvjp  = 0.25*( u(i,jp,k)+u(i,j,k) ) !*( v(ip,j ,k )+v(i,j ,k ) )
          uvjm  = 0.25*( u(i,jm,k)+u(i,j,k) ) !*( v(ip,jm,k )+v(i,jm,k ) )
          uwkp  = 0.25*( u(i,j,kp)+u(i,j,k) ) !*( w(ip,j ,k )+w(i,j ,k ) )
          uwkm  = 0.25*( u(i,j,km)+u(i,j,k) ) !*( w(ip,j ,km)+w(i,j ,km) )
          uvip  = 0.25*( u(i ,j,k)+u(i ,jp,k) ) !*( v(i,j,k )+v(ip,j ,k) )
          uvim  = 0.25*( u(im,j,k)+u(im,jp,k) ) !*( v(i,j,k )+v(im,j ,k) )
          uwip  = 0.25*( u(i ,j ,k)+u(i ,j ,kp) ) !*( w(i,j,k)+w(ip,j,k) )
          uwim  = 0.25*( u(im,j ,k)+u(im,j ,kp) ) !*( w(i,j,k)+w(im,j,k) )
          !
          ! touch v
          !
          vvjp  = 0.25*( v(i,j,k )+v(i,jp,k)  )*( v(i,j,k )+v(i ,jp,k) )
          vvjm  = 0.25*( v(i,j,k )+v(i,jm,k)  )*( v(i,j,k )+v(i ,jm,k) )
          dvdxp = (v(ip,j,k)-v(i,j,k))*dxi
          dvdxm = (v(i,j,k)-v(im,j,k))*dxi
          dvdyp = (v(i,jp,k)-v(i,j,k))*dyi
          dvdym = (v(i,j,k)-v(i,jm,k))*dyi
          dvdzp = (v(i,j,kp)-v(i,j,k))*dzci(k)
          dvdzm = (v(i,j,k)-v(i,j,km))*dzci(km)
          uvip  = uvip*( v(i,j,k )+v(ip,j ,k) )
          uvim  = uvim*( v(i,j,k )+v(im,j ,k) )
          !
          dvdt_temp =   dxi*(     -uvip + uvim ) + (dvdxp-dvdxm)*visc*dxi+ &
                        dyi*(     -vvjp + vvjm ) + (dvdyp-dvdym)*visc*dyi+ &
                                                   (dvdzp-dvdzm)*visc*dzfi(k)
          uvjp  = uvjp*( v(ip,j ,k )+v(i,j ,k ) )
          uvjm  = uvjm*( v(ip,jm,k )+v(i,jm,k ) )
          dudt_temp = dudt_temp + dyi*(     -uvjp + uvjm )
          !
          wvkp  = 0.25*( v(i ,j ,kp)+v(i ,j ,k ) ) !*( w(i,j,k )+w(i,jp,k ) )
          wvkm  = 0.25*( v(i ,j ,km)+v(i ,j ,k ) ) !*( w(i,j,km)+w(i,jp,km) )
          vwjp  = 0.25*( v(i ,j ,k )+v(i ,j ,kp) ) !*( w(i,j,k )+w(i,jp,k ) )
          vwjm  = 0.25*( v(i ,jm,k )+v(i ,jm,kp) ) !*( w(i,j,k )+w(i,jm,k ) )
          !
          ! touch w
          !
          wwkp  = 0.25*( w(i,j,k)+w(i,j,kp) )*( w(i ,j ,k)+w(i ,j ,kp) )
          wwkm  = 0.25*( w(i,j,k)+w(i,j,km) )*( w(i ,j ,k)+w(i ,j ,km) )
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
          !
          ! momentum balance
          !
          dwdt_temp =   dxi*(     -uwip + uwim ) + (dwdxp-dwdxm)*visc*dxi+ &
                        dyi*(     -vwjp + vwjm ) + (dwdyp-dwdym)*visc*dyi+ &
                        dzci(k)*( -wwkp + wwkm ) + (dwdzp-dwdzm)*visc*dzci(k)
          !
          uwkp  = uwkp*( w(ip,j ,k )+w(i,j ,k ) )
          uwkm  = uwkm*( w(ip,j ,km)+w(i,j ,km) )
          dudt_temp = dudt_temp + dzfi(k)*( -uwkp + uwkm )
          !
          wvkp  = wvkp*( w(i,j,k )+w(i,jp,k ) )
          wvkm  = wvkm*( w(i,j,km)+w(i,jp,km) )
          dvdt_temp = dvdt_temp + dzfi(k)*( -wvkp + wvkm )
          !
          dudt(i,j,k) = dudt_temp + bforcex
          dvdt(i,j,k) = dvdt_temp + bforcey
          dwdt(i,j,k) = dwdt_temp + bforcez
        enddo
      enddo
    enddo
#ifndef USE_CUDA
    !$OMP END PARALLEL DO
#endif
    !@cuf istat=cudaDeviceSynchronize()
    !
    taux(:) = 0.
    taux2d = 0.
    taux3d = 0.
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
    !@cuf istat=cudaDeviceSynchronize()
    taux(3) = taux3d
    call mpi_allreduce(MPI_IN_PLACE,taux(1),3,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    nxg = nx*dims(1)
    nyg = ny*dims(2)
    nzg = nz
    taux(1) = taux(1)/(1.*nyg)
    taux(2) = taux(2)/(1.*nxg)
    taux(3) = taux(3)/(1.*nxg*nyg)
    tauy(:) = 0.
    tauy1d = 0.
    tauy3d = 0.
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
    !@cuf istat=cudaDeviceSynchronize()
    tauy(3) = tauy3d
    call mpi_allreduce(MPI_IN_PLACE,tauy(1),3,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    nxg = nx*dims(1)
    nyg = ny*dims(2)
    nzg = nz
    tauy(1) = tauy(1)/(1.*nyg)
    tauy(2) = tauy(2)/(1.*nxg)
    tauy(3) = tauy(3)/(1.*nxg*nyg)
    !
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
        tauz1d = tauz1d + (dwdxp+dwdxm)
      enddo
    enddo
    !@cuf istat=cudaDeviceSynchronize()
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
    !@cuf istat=cudaDeviceSynchronize()
    tauz(2)=tauz2d
    call mpi_allreduce(MPI_IN_PLACE,tauz(1),3,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    nxg = nx*dims(1)
    nyg = ny*dims(2)
    nzg = nz
    tauz(1) = tauz(1)/(1.*nyg)
    tauz(2) = tauz(2)/(1.*nxg)
    tauz(3) = tauz(3)/(1.*nxg*nyg)
    return
  end subroutine momxyzad
  !
  subroutine momxad(nx,ny,nz,dxi,dyi,dzi,dzci,dzfi,dzflzi,visc,u,v,w,dudt,taux)
    !@cuf use cudafor
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi,dzi,visc
    real(rp), intent(in), dimension(0:) :: dzci,dzfi,dzflzi
    real(rp), dimension(0:,0:,0:), intent(in) :: u,v,w
    real(rp), dimension(:,:,:), intent(out) :: dudt
    real(rp), dimension(3)  , intent(out) :: taux
    integer :: im,ip,jm,jp,km,kp,i,j,k
    real(rp) :: uuip,uuim,uvjp,uvjm,uwkp,uwkm
    real(rp) :: dudxp,dudxm,dudyp,dudym,dudzp,dudzm
    real(rp) :: bforcex
    real(rp) :: taux2d,taux3d
#ifdef USE_CUDA
    attributes(managed) :: u,v,w,dudt,dzci,dzfi,dzflzi
    integer :: istat
#endif
    integer :: nxg,nyg,nzg
    !
    bforcex = bforce(1)
#ifdef USE_CUDA
    !$cuf kernel do(3) <<<*,*>>>
#else
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k,im,jm,km,ip,jp,kp) &
    !$OMP PRIVATE(uuip,uuim,uvjp,uvjm,uwkp,uwkm) &
    !$OMP PRIVATE(dudxp,dudxm,dudyp,dudym,dudzp,dudzm) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,dzi,visc,u,v,w,dudt,dzci,dzfi)
    !$OMP SHARED(nx,ny,nz,dxi,dyi,dzi,visc,u,v,w,dudt,dzci,dzfi,bforcex)
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
          dudxp = (u(ip,j,k)-u(i,j,k))*dxi
          dudxm = (u(i,j,k)-u(im,j,k))*dxi
          dudyp = (u(i,jp,k)-u(i,j,k))*dyi
          dudym = (u(i,j,k)-u(i,jm,k))*dyi
          dudzp = (u(i,j,kp)-u(i,j,k))*dzci(k)
          dudzm = (u(i,j,k)-u(i,j,km))*dzci(km)
          !
          ! Momentum balance
          !
          dudt(i,j,k) = dxi*(     -uuip + uuim ) + (dudxp-dudxm)*visc*dxi + &
                        dyi*(     -uvjp + uvjm ) + (dudyp-dudym)*visc*dyi + &
                        dzfi(k)*( -uwkp + uwkm ) + (dudzp-dudzm)*visc*dzfi(k) + &
                        bforcex
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
    !@cuf istat=cudaDeviceSynchronize()
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
    !@cuf istat=cudaDeviceSynchronize()
    taux(3) = taux3d
    call mpi_allreduce(MPI_IN_PLACE,taux(1),3,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,taux(1),3,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    nxg = nx*dims(1)
    nyg = ny*dims(2)
    nzg = nz
    taux(1) = taux(1)/(1.*nyg)
    taux(2) = taux(2)/(1.*nxg)
    taux(3) = taux(3)/(1.*nxg*nyg)
    return
  end subroutine momxad
  !
  subroutine momyad(nx,ny,nz,dxi,dyi,dzi,dzci,dzfi,dzflzi,visc,u,v,w,dvdt,tauy)
    !@cuf use cudafor
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi,dzi,visc
    real(rp), intent(in), dimension(0:) :: dzci,dzfi,dzflzi
    real(rp), dimension(0:,0:,0:), intent(in) :: u,v,w
    real(rp), dimension(:,:,:), intent(out) :: dvdt
    real(rp), dimension(3), intent(out) :: tauy
    integer :: im,ip,jm,jp,km,kp,i,j,k
    real(rp) :: uvip,uvim,vvjp,vvjm,wvkp,wvkm
    real(rp) :: dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm
    real(rp) :: bforcey
    real(rp) :: tauy1d,tauy3d
#ifdef USE_CUDA
    attributes(managed) :: u,v,w,dvdt,dzci,dzfi,dzflzi
    integer :: istat
#endif
    integer :: nxg,nyg,nzg
    !
    bforcey = bforce(2)
#ifdef USE_CUDA
    !$cuf kernel do(3) <<<*,*>>>
#else
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k,im,jm,km,ip,jp,kp) &
    !$OMP PRIVATE(uvip,uvim,vvjp,vvjm,wvkp,wvkm) &
    !$OMP PRIVATE(dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,dzi,visc,u,v,w,dvdt,dzci,dzfi)
    !$OMP SHARED(nx,ny,nz,dxi,dyi,dzi,visc,u,v,w,dvdt,dzci,dzfi,bforcey)
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
          dvdxp = (v(ip,j,k)-v(i,j,k))*dxi
          dvdxm = (v(i,j,k)-v(im,j,k))*dxi
          dvdyp = (v(i,jp,k)-v(i,j,k))*dyi
          dvdym = (v(i,j,k)-v(i,jm,k))*dyi
          dvdzp = (v(i,j,kp)-v(i,j,k))*dzci(k)
          dvdzm = (v(i,j,k)-v(i,j,km))*dzci(km)
          !
          ! Momentum balance
          !
          dvdt(i,j,k) = dxi*(     -uvip + uvim ) + (dvdxp-dvdxm)*visc*dxi+ &
                        dyi*(     -vvjp + vvjm ) + (dvdyp-dvdym)*visc*dyi+ &
                        dzfi(k)*( -wvkp + wvkm ) + (dvdzp-dvdzm)*visc*dzfi(k)+ &
                        bforcey
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
    !@cuf istat=cudaDeviceSynchronize()
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
    !@cuf istat=cudaDeviceSynchronize()
    tauy(3) = tauy3d
    call mpi_allreduce(MPI_IN_PLACE,tauy(1),3,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    nxg = nx*dims(1)
    nyg = ny*dims(2)
    nzg = nz
    tauy(1) = tauy(1)/(1.*nyg)
    tauy(2) = tauy(2)/(1.*nxg)
    tauy(3) = tauy(3)/(1.*nxg*nyg)
    return
  end subroutine momyad
  !
  subroutine momzad(nx,ny,nz,dxi,dyi,dzi,dzci,dzfi,dzflzi,visc,u,v,w,dwdt,tauz)
    !@cuf use cudafor
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi,dzi,visc
    real(rp), intent(in), dimension(0:) :: dzci,dzfi,dzflzi
    real(rp), dimension(0:,0:,0:), intent(in) :: u,v,w
    real(rp), dimension(:,:,:), intent(out) :: dwdt
    real(rp), dimension(3), intent(out) :: tauz
    integer :: im,ip,jm,jp,km,kp,i,j,k
    real(rp) :: uwip,uwim,vwjp,vwjm,wwkp,wwkm
    real(rp) :: dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm
    real(rp) :: bforcez
    real(rp) :: tauz1d,tauz2d
#ifdef USE_CUDA
    attributes(managed) :: u,v,w,dwdt,dzci,dzfi,dzflzi
    integer :: istat
#endif
    integer :: nxg,nyg,nzg
    !
    bforcez = bforce(3)
#ifdef USE_CUDA
    !$cuf kernel do(3) <<<*,*>>>
#else
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k,im,jm,km,ip,jp,kp) &
    !$OMP PRIVATE(uwip,uwim,vwjp,vwjm,wwkp,wwkm) &
    !$OMP PRIVATE(dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm) &
    !$OMP SHARED(nx,ny,nz,dxi,dyi,dzi,visc,u,v,w,dwdt,dzci,dzfi,bforcez)
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
          dwdxp = (w(ip,j,k)-w(i,j,k))*dxi
          dwdxm = (w(i,j,k)-w(im,j,k))*dxi
          dwdyp = (w(i,jp,k)-w(i,j,k))*dyi
          dwdym = (w(i,j,k)-w(i,jm,k))*dyi
          dwdzp = (w(i,j,kp)-w(i,j,k))*dzfi(kp)
          dwdzm = (w(i,j,k)-w(i,j,km))*dzfi(k)
          !
          ! Momentum balance
          !
          dwdt(i,j,k) = dxi*(     -uwip + uwim ) + (dwdxp-dwdxm)*visc*dxi+ &
                        dyi*(     -vwjp + vwjm ) + (dwdyp-dwdym)*visc*dyi+ &
                        dzci(k)*( -wwkp + wwkm ) + (dwdzp-dwdzm)*visc*dzci(k)+ &
                        bforcez
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
    call mpi_allreduce(MPI_IN_PLACE,tauz(1),3,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    nxg = nx*dims(1)
    nyg = ny*dims(2)
    nzg = nz
    tauz(1) = tauz(1)/(1.*nyg)
    tauz(2) = tauz(2)/(1.*nxg)
    tauz(3) = tauz(3)/(1.*nxg*nyg)
    return
  end subroutine momzad
  !
  subroutine momxp(nx,ny,nz,dxi,p,dudt)
    !@cuf use cudafor
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi 
    real(rp), dimension(0:,0:,0:), intent(in) :: p
    real(rp), dimension(:,:,:), intent(out) :: dudt
#ifdef USE_CUDA
    attributes(managed) :: p,dudt
    integer :: istat
#endif
    integer :: i,j,k
    integer :: ip
#ifdef USE_CUDA
    !$cuf kernel do(3) <<<*,*>>>
#else
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k,ip) &
    !$OMP SHARED(nx,ny,nz,dxi,p,dudt)
#endif
    do k=1,nz
      do j=1,ny
        do i=1,nx
          ip = i + 1
          dudt(i,j,k) = - dxi*( p(ip,j,k)-p(i,j,k) )
        enddo
      enddo
    enddo
#ifndef USE_CUDA
    !$OMP END PARALLEL DO
#endif
    !@cuf istat=cudaDeviceSynchronize()
    return
  end subroutine momxp
  !
  subroutine momyp(nx,ny,nz,dyi,p,dvdt)
    !@cuf use cudafor
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dyi 
    real(rp), dimension(0:,0:,0:), intent(in) :: p
    real(rp), dimension(:,:,:), intent(out) :: dvdt
#ifdef USE_CUDA
    attributes(managed) :: p,dvdt
    integer :: istat
#endif
    integer :: i,j,k
    integer :: jp
#ifdef USE_CUDA
    !$cuf kernel do(3) <<<*,*>>>
#else
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k,jp) &
    !$OMP SHARED(nx,ny,nz,dyi,p,dvdt)
#endif
    do k=1,nz
      do j=1,ny
        jp = j + 1
        do i=1,nx
          dvdt(i,j,k) = - dyi*( p(i,jp,k)-p(i,j,k) )
        enddo
      enddo
    enddo
#ifndef USE_CUDA
    !$OMP END PARALLEL DO
#endif
    !@cuf istat=cudaDeviceSynchronize()
    return
  end subroutine momyp
  !
  subroutine momzp(nx,ny,nz,dzci,p,dwdt)
    !@cuf use cudafor
    implicit none
    integer , intent(in) :: nx,ny,nz
    real(rp), intent(in), dimension(0:) :: dzci
    real(rp), dimension(0:,0:,0:), intent(in) :: p
    real(rp), dimension(:,:,:), intent(out) :: dwdt
#ifdef USE_CUDA
    attributes(managed) :: p,dwdt,dzci
    integer :: istat
#endif
    integer :: kp
    integer :: i,j,k
    !
#ifdef USE_CUDA
    !$cuf kernel do(3) <<<*,*>>>
#else
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k,kp) &
    !$OMP SHARED(nx,ny,nz,p,dwdt,dzci)
#endif
    do k=1,nz
      kp = k + 1
      do j=1,ny
        do i=1,nx
          dwdt(i,j,k) = - dzci(k)*( p(i,j,kp)-p(i,j,k) )
        enddo
      enddo
    enddo
#ifndef USE_CUDA
    !$OMP END PARALLEL DO
#endif
    !@cuf istat=cudaDeviceSynchronize()
    return
  end subroutine momzp
  subroutine momp(nx,ny,nz,dxi,dyi,dzci,p,dudt,dvdt,dwdt)
    !@cuf use cudafor
    !
    ! combined momxp,momyp,momzp
    !
    implicit none
    integer, intent(in) :: nx,ny,nz
    real(rp), intent(in) :: dxi,dyi 
    real(rp), intent(in), dimension(0:) :: dzci
    real(rp), dimension(0:,0:,0:), intent(in) :: p
    real(rp), dimension(:,:,:), intent(out) :: dudt,dvdt,dwdt
#ifdef USE_CUDA
    attributes(managed) :: p,dudt,dvdt,dwdt,dzci
    integer :: istat
#endif
    integer :: ip,jp,kp
    integer :: i,j,k
    !
#ifdef USE_CUDA
    !$cuf kernel do(3) <<<*,*>>>
#else
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k,ip,jp,kp) &
    !$OMP SHARED(nx,ny,nz,p,dudt,dvdt,dwdt,dzci,dxi,dyi)
#endif
    do k=1,nz
      kp = k + 1
      do j=1,ny
        jp = j + 1
        do i=1,nx
          ip = i + 1
          dudt(i,j,k) = - dxi*    ( p(ip,j,k)-p(i,j,k) )
          dvdt(i,j,k) = - dyi*    ( p(i,jp,k)-p(i,j,k) )
          dwdt(i,j,k) = - dzci(k)*( p(i,j,kp)-p(i,j,k) )
        enddo
      enddo
    enddo
#ifndef USE_CUDA
    !$OMP END PARALLEL DO
#endif
    !@cuf istat=cudaDeviceSynchronize()
    return
  end subroutine momp
end module mod_mom
