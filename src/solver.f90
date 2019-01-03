module mod_solver
  use iso_c_binding, only: C_PTR
  use decomp_2d
  use mod_fft   , only: fftd,ffti
  use mod_param , only: dims
#ifdef USE_CUDA
  use cudafor
  use mod_fftw_param
#endif
  implicit none
  private
  public solver
  contains
  subroutine solver(n,arrplan,normfft,lambdaxy,a,b,c,bcz,c_or_f,pz_pad)
    implicit none
    integer, intent(in), dimension(3) :: n
    type(C_PTR), intent(in), dimension(2,2) :: arrplan
    real(8), intent(in) :: normfft
    real(8), intent(in), dimension(n(1),n(2)) :: lambdaxy
    real(8), intent(in), dimension(n(3)) :: a,b,c
    character(len=1), dimension(0:1), intent(in) :: bcz
    character(len=1), intent(in), dimension(3) :: c_or_f
    real(8), intent(inout), dimension(0:,0:,0:) :: pz_pad
    !real(8), dimension(n(1)*dims(1),n(2)*dims(2)/dims(1),n(3)/dims(2)) :: px
    !real(8), dimension(n(1)*dims(1)/dims(1),n(2)*dims(2),n(3)/dims(2)) :: py
    real(8), allocatable, dimension(:,:,:) :: px,py,pz
    integer :: i,j,k
#ifdef USE_CUDA
    attributes(managed) :: px,py,pz,pz_pad,lambdaxy,a,b,c
    integer :: istat,ii
    real(8), allocatable, dimension(:,:,:), device :: py_t
    complex(8), allocatable, dimension(:,:,:), device :: pxc, pyc, pyc_t
#endif
    integer :: dotrans
    integer, dimension(3) :: ng
    integer :: q

    dotrans = 0
    if( dims(1)*dims(2) .ne. 1) dotrans = 1

    ng(:) = n(:)
    ng(1:2) = ng(1:2)*dims(1:2)
    allocate(px(ng(1),ng(2)/dims(1),ng(3)/dims(2)))
    allocate(py(ng(1)/dims(1),ng(2),ng(3)/dims(2)))
    allocate(pz(ng(1)/dims(1),ng(2)/dims(2),ng(3)))

#ifdef USE_CUDA
    allocate( pxc( ng(1)/2 + 1, ng(2)/dims(1), ng(3)/dims(2) ) )
    allocate( pyc_t( ng(2)/2 + 1, ng(1)/dims(1), ng(3)/dims(2) ) )
    allocate( py_t( ng(2), ng(1)/dims(1), ng(3)/dims(2) ) )

    istat = cudaMemAdvise( px, size(px), cudaMemAdviseSetPreferredLocation, 0 )
    istat = cudaMemAdvise( py, size(py), cudaMemAdviseSetPreferredLocation, 0 )
    istat = cudaMemAdvise( pz, size(pz), cudaMemAdviseSetPreferredLocation, 0 )

    istat = cudaMemPrefetchAsync( pz_pad, size(pz_pad), 0, 0)
    istat = cudaMemPrefetchAsync( pz, size(pz), 0, 0)

    if(dotrans) then
      istat = cudaMemPrefetchAsync( px, size(px), cudaCpuDeviceId, 0)
      istat = cudaMemPrefetchAsync( py, size(py), cudaCpuDeviceId, 0)
    else
      istat = cudaMemPrefetchAsync( px, size(px), 0, 0)
      istat = cudaMemPrefetchAsync( py, size(py), 0, 0)
    endif

    !$cuf kernel do(3) <<<*,*>>>
#endif
    do k=1,ng(3)
    do j=1,ng(2)/dims(2)
    do i=1,ng(1)/dims(1)
      pz(i,j,k) = pz_pad(i,j,k)
    enddo
    enddo
    enddo
    !
    if(dotrans) then
#ifdef USE_CUDA
      istat = cudaMemPrefetchAsync(pz,size(pz),cudaCpuDeviceId,0)
      !@cuf istat = cudaDeviceSynchronize()
#endif
      !call transpose_z_to_x(pz,px)
      call transpose_z_to_y(pz,py)
      call transpose_y_to_x(py,px)
    else
#ifdef USE_CUDA
      !$cuf kernel do(3) <<<*,*>>>
#endif
      do k=1,ng(3)
      do j=1,ng(2)/dims(2)
      do i=1,ng(1)/dims(1)
        px(i,j,k) = pz(i,j,k)
      enddo
      enddo
      enddo
      !istat = cudaMemcpy( px, pz, size(px), cudaMemcpyDeviceToDevice )
    endif

#ifdef USE_CUDA
    if( dotrans ) istat = cudaMemPrefetchAsync(px,size(px),0,0)
    istat = cufftExecD2Z(cufft_plan_fwd_x, px, pxc)
!
!$cuf kernel do(2) <<<*,*>>>
    do k=1,ng(3)/dims(2)
    do j=1,ng(2)/dims(1)
      do i=1,(ng(1)/2)+1
        px(i,j,k) = REAL(pxc(i,j,k))
      end do
      ii=2
      do i=ng(1),(ng(1)/2)+2,-1
        px(i,j,k) = AIMAG(pxc(ii,j,k))
        ii = ii + 1
      end do
    end do
    end do

    if( dotrans ) istat = cudaMemPrefetchAsync(px,size(px),cudaCpuDeviceId,0)
    !@cuf istat = cudaDeviceSynchronize()
#else
    call fftd(arrplan(1,1),px) ! fwd transform in x
#endif
    !
    if(dotrans) then
      call transpose_x_to_y(px,py)
    else
#ifdef USE_CUDA
      !$cuf kernel do(3) <<<*,*>>>
#endif
      do k=1,ng(3)/dims(2)
      do j=1,ng(2)
      do i=1,ng(1)/dims(1)
        py(i,j,k) = px(i,j,k)
      enddo
      enddo
      enddo
      !istat = cudaMemcpy( py, px, size(px), cudaMemcpyDeviceToDevice )
    endif

#ifdef USE_CUDA
    if( dotrans ) istat = cudaMemPrefetchAsync(py,size(py),0,0)

    !$cuf kernel do(3) <<<*,*>>>
    do k=1,ng(3)/dims(2)
    do j=1,ng(1)/dims(1)
    do i=1,ng(2)
      py_t(i,j,k) = py(j,i,k)
    enddo
    enddo
    enddo

    istat = cufftExecD2Z(cufft_plan_fwd_y, py_t, pyc_t)

!$cuf kernel do(2) <<<*,*>>>
    do k=1,ng(3)/dims(2)
    do j=1,ng(1)/dims(1)
      do i=1,(ng(2)/2)+1
        py(j,i,k) = REAL(pyc_t(i,j,k))
      end do
      ii=2
      do i=ng(2),(ng(2)/2)+2,-1
        py(j,i,k) = AIMAG(pyc_t(ii,j,k))
        ii = ii + 1
      end do
    end do
    end do

    if( dotrans ) istat = cudaMemPrefetchAsync(py,size(py),cudaCpuDeviceId,0)
    !@cuf istat = cudaDeviceSynchronize()
#else
    call fftd(arrplan(1,2),py) ! fwd transform in y
#endif
    !
    if(dotrans) then
      call transpose_y_to_z(py,pz)
    else
#ifdef USE_CUDA
      !$cuf kernel do(3) <<<*,*>>>
#endif
      do k=1,ng(3)
      do j=1,ng(2)/dims(2)
      do i=1,ng(1)/dims(1)
        pz(i,j,k) = py(i,j,k)
      enddo
      enddo
      enddo
      !istat = cudaMemcpy( pz, py, size(pz), cudaMemcpyDeviceToDevice )
    endif


    q = 0
    if(c_or_f(3).eq.'f'.and.bcz(1).eq.'D') q = 1
    if(bcz(0)//bcz(1).eq.'PP') then
      call gaussel_periodic(n(1),n(2),n(3)-q,a,b,c,lambdaxy,pz)
    else
#ifdef USE_CUDA
      istat = cudaMemPrefetchAsync(lambdaxy,size(lambdaxy),0,0)
      if( dotrans ) istat = cudaMemPrefetchAsync(pz,size(pz),0,0)
      call gaussel_gpu(         n(1),n(2),n(3)-q,a,b,c,lambdaxy,pz)
      if( dotrans ) istat = cudaMemPrefetchAsync(pz,size(pz),cudaCpuDeviceId,0)
      !@cuf istat = cudaDeviceSynchronize()
#else
      call gaussel(         n(1),n(2),n(3)-q,a,b,c,lambdaxy,pz)
#endif
    endif
    !
    if(dotrans) then
      call transpose_z_to_y(pz,py)
    else
#ifdef USE_CUDA
      !$cuf kernel do(3) <<<*,*>>>
#endif
      do k=1,ng(3)/dims(2)
      do j=1,ng(2)
      do i=1,ng(1)/dims(1)
        py(i,j,k) = pz(i,j,k)
      enddo
      enddo
      enddo
      !istat = cudaMemcpy( py, pz, size(py), cudaMemcpyDeviceToDevice )
    endif

#ifdef USE_CUDA
    if( dotrans ) istat = cudaMemPrefetchAsync(py,size(py),0,0)

    pyc_t = (0.d0,0.d0)
!$cuf kernel do(2) <<<*,*>>>
    do k=1,ng(3)/dims(2)
    do j=1,ng(1)/dims(1)
      do i=1,(ng(2)/2)+1
        pyc_t(i,j,k)%re = py(j,i,k)
      end do
      ii=2
      do i=ng(2),(ng(2)/2)+2,-1
        pyc_t(ii,j,k)%im = py(j,i,k)
        ii = ii + 1
      end do
    end do
    end do

    istat = cufftExecZ2D(cufft_plan_bwd_y, pyc_t, py_t)

    !$cuf kernel do(3) <<<*,*>>>
    do k=1,ng(3)/dims(2)
    do j=1,ng(1)/dims(1)
    do i=1,ng(2)
      py(j,i,k) = py_t(i,j,k)
    enddo
    enddo
    enddo

    if( dotrans ) istat = cudaMemPrefetchAsync(py,size(py),cudaCpuDeviceId,0)
    !@cuf istat = cudaDeviceSynchronize()
#else
    call ffti(arrplan(2,2),py) ! bwd transform in y
#endif
    !
    if(dotrans) then
      call transpose_y_to_x(py,px)
    else
#ifdef USE_CUDA
      !$cuf kernel do(3) <<<*,*>>>
#endif
      do k=1,ng(3)/dims(2)
      do j=1,ng(2)/dims(1)
      do i=1,ng(1)
        px(i,j,k) = py(i,j,k)
      enddo
      enddo
      enddo
      !istat = cudaMemcpy( px, py, size(px), cudaMemcpyDeviceToDevice )
    endif

#ifdef USE_CUDA
    if( dotrans ) istat = cudaMemPrefetchAsync(px,size(px),0,0)

!$cuf kernel do(2) <<<*,*>>>
    do k=1,ng(3)/dims(2)
    do j=1,ng(2)/dims(1)
      do i=1,(ng(1)/2)+1
        pxc(i,j,k)%re = px(i,j,k)
      end do
      ii=2
      do i=ng(1),(ng(1)/2)+2,-1
        pxc(ii,j,k)%im = px(i,j,k)
        ii = ii + 1
      end do
    end do
    end do

    istat = cufftExecZ2D(cufft_plan_bwd_x, pxc, px)

    if( dotrans ) istat = cudaMemPrefetchAsync(px,size(px),cudaCpuDeviceId,0)
    !@cuf istat = cudaDeviceSynchronize()
#else
    call ffti(arrplan(2,1),px) ! bwd transform in x
#endif
    !
    if(dotrans) then
      !call transpose_x_to_z(px,pz)
      call transpose_x_to_y(px,py)
      call transpose_y_to_z(py,pz)
    else
#ifdef USE_CUDA
      !$cuf kernel do(3) <<<*,*>>>
#endif
      do k=1,ng(3)
      do j=1,ng(2)/dims(2)
      do i=1,ng(1)/dims(1)
        pz(i,j,k) = px(i,j,k)
      enddo
      enddo
      enddo
      !istat = cudaMemcpy( pz, px, size(pz), cudaMemcpyDeviceToDevice )
    endif

#ifdef USE_CUDA
    if( dotrans ) istat = cudaMemPrefetchAsync(pz,size(pz),0,0)

!$cuf kernel do(3) <<<*,*>>>
#endif
    do k=lbound(pz,3),ubound(pz,3)
    do j=lbound(pz,2),ubound(pz,2)
    do i=lbound(pz,1),ubound(pz,1)
      pz_pad(i,j,k) = pz(i,j,k)*normfft
    enddo
    enddo
    enddo

    deallocate(px,py,pz)
#ifdef USE_CUDA
    deallocate(pxc,pyc_t,py_t)
#endif

    return
  end subroutine solver
  !
#ifdef USE_CUDA
  subroutine gaussel_gpu(nx,ny,n,a,b,c,lambdaxy,p)
    implicit none
    integer, intent(in) :: nx,ny,n
    real(8), intent(in), dimension(:), managed :: a,b,c
    real(8), intent(in), dimension(nx,ny), managed :: lambdaxy
    real(8), intent(inout), dimension(:,:,:), managed :: p
    real(8), dimension(nx,ny,n), device :: bb
    real(8), dimension(nx,ny,n), device :: d
    real(8) :: z
    integer :: i,j,k
    !
    !solve tridiagonal system
    !

    !$cuf kernel do(2) <<<*,*>>>
    do j=1,ny
      do i=1,nx
        do k=1,n
          bb(i,j,k) = b(k) + lambdaxy(i,j)
        enddo
        z = 1.d0/bb(i,j,1)
        d(i,j,1) = c(1)*z
        p(i,j,1) = p(i,j,1)*z
        do k=2,n-1
          z = 1.d0/(bb(i,j,k)-a(k)*d(i,j,k-1))
          d(i,j,k) = c(k)*z
          p(i,j,k) = (p(i,j,k)-a(k)*p(i,j,k-1))*z
        enddo
        z = bb(i,j,n)-a(n)*d(i,j,n-1)
        if(z.ne.0.d0) then
          p(i,j,n) = (p(i,j,n)-a(n)*p(i,j,n-1))/z
        else
          p(i,j,n) = 0.d0
        endif
        !
        ! backward substitution
        !
        do k=n-1,1,-1
          p(i,j,k) = p(i,j,k) - d(i,j,k)*p(i,j,k+1)
        enddo

      enddo
    enddo

    return
  end subroutine gaussel_gpu
#endif

  subroutine gaussel(nx,ny,n,a,b,c,lambdaxy,p)
    implicit none
    integer, intent(in) :: nx,ny,n
    real(8), intent(in), dimension(:) :: a,b,c
    real(8), intent(in), dimension(nx,ny) :: lambdaxy
    real(8), intent(inout), dimension(:,:,:) :: p
    real(8), dimension(n) :: bb
    integer :: i,j
    !
    !solve tridiagonal system
    !
    !$OMP PARALLEL DEFAULT(none) &
    !$OMP PRIVATE(i,j,bb) &
    !$OMP SHARED(nx,ny,n,a,b,c,lambdaxy,p)
    !$OMP DO COLLAPSE(2)
    do j=1,ny
      do i=1,nx
        bb(:) = b(1:n) + lambdaxy(i,j)
        call dgtsv_homebrewed(n,a,bb,c,p(i,j,1:n))
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    return
  end subroutine gaussel
  !
  subroutine gaussel_periodic(nx,ny,n,a,b,c,lambdaxy,p)
    implicit none
    integer, intent(in) :: nx,ny,n
    real(8), intent(in), dimension(:) :: a,b,c
    real(8), intent(in), dimension(nx,ny) :: lambdaxy
    real(8), intent(inout), dimension(:,:,:) :: p
    real(8), dimension(n) :: bb,p1,p2
    integer :: i,j,info
    !
    !solve tridiagonal system
    !
    !$OMP PARALLEL DEFAULT(none) &
    !$OMP PRIVATE(i,j,bb,p1,p2) &
    !$OMP SHARED(nx,ny,n,a,b,c,lambdaxy,p)
    !$OMP DO COLLAPSE(2)
    do j=1,ny
      do i=1,nx
        bb(:)  = b(:) + lambdaxy(i,j)
        p1(1:n-1) = p(i,j,1:n-1)
        call dgtsv_homebrewed(n-1,a(1:n-1),bb(1:n-1),c(1:n-2),p1(1:n-1))
        p2(:) = 0.d0
        p2(1  ) = -a(1  )
        p2(n-1) = -c(n-1)
        call dgtsv_homebrewed(n-1,a(2:n-1),bb(1:n-1),c(1:n-2),p2(1:n-1))
        p(i,j,n) = (p(i,j,n) - c(n)*p1(1) - a(n)*p1(n-1)) / &
                   (bb(   n) + c(n)*p2(1) + a(n)*p2(n-1))
        p(i,j,1:n-1) = p1(1:n-1) + p2(1:n-1)*p(i,j,n)
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    return
  end subroutine gaussel_periodic
  subroutine dgtsv_homebrewed(n,a,b,c,p)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in   ), dimension(:) :: a,b,c
    real(8), intent(inout), dimension(:) :: p
    real(8), dimension(n) :: d
    real(8) :: z
    integer :: l
    !
    ! Gauss elimination
    !
    z = 1.d0/b(1)
    d(1) = c(1)*z
    p(1) = p(1)*z
    do l=2,n-1
      z    = 1.d0/(b(l)-a(l)*d(l-1))
      d(l) = c(l)*z
      p(l) = (p(l)-a(l)*p(l-1))*z
    enddo
    z = b(n)-a(n)*d(n-1)
    if(z.ne.0.d0) then
      p(n) = (p(n)-a(n)*p(n-1))/z
    else
      p(n) = 0.d0
    endif
    !
    ! backward substitution
    !
    do l=n-1,1,-1
      p(l) = p(l) - d(l)*p(l+1)
    enddo
    return
  end subroutine dgtsv_homebrewed
end module mod_solver
