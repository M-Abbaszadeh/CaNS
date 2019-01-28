module mod_solver
  use iso_c_binding, only: C_PTR
  use decomp_2d
  use mod_fft   , only: fftd,ffti
  use mod_param , only: dims
#ifdef USE_CUDA
  use cudafor
  use mod_fftw_param
  use mod_common_mpi, only: mydev
#endif
  implicit none
  real(8), allocatable, dimension(:,:,:) :: px,py,pz
#ifdef USE_CUDA
  real(8), allocatable, dimension(:,:,:) :: pw
  attributes(managed) :: px,py,pz,pw
  real(8), allocatable, dimension(:,:,:), device :: py_t
#ifdef EPHC
  real(8), allocatable, dimension(:,:,:), device :: pxc, pyc_t
#else
  complex(8), allocatable, dimension(:,:,:), device :: pxc, pyc_t
#endif
#endif
  private
  public solver
  contains
  subroutine solver(n,arrplan,normfft,lambdaxy,a,b,c,bcz,c_or_f,pz_pad)
    implicit none
    integer, intent(in), dimension(3) :: n
    type(C_PTR), intent(in), dimension(2,2) :: arrplan
    real(8), intent(in) :: normfft
#ifdef USE_CUDA
    real(8), intent(in), dimension(n(1)*dims(1)/dims(2),n(2)*dims(2)/dims(1)) :: lambdaxy
#else
    real(8), intent(in), dimension(n(1),n(2)) :: lambdaxy
#endif
    real(8), intent(in), dimension(n(3)) :: a,b,c
    character(len=1), dimension(0:1), intent(in) :: bcz
    character(len=1), intent(in), dimension(3) :: c_or_f
    real(8), intent(inout), dimension(0:,0:,0:) :: pz_pad
    !real(8), dimension(n(1)*dims(1),n(2)*dims(2)/dims(1),n(3)/dims(2)) :: px
    !real(8), dimension(n(1)*dims(1)/dims(1),n(2)*dims(2),n(3)/dims(2)) :: py
    !real(8), allocatable, dimension(:,:,:) :: px,py,pz
    integer :: i,j,k
#ifdef USE_CUDA
    !attributes(managed) :: px,py,pz,pz_pad,lambdaxy,a,b,c
    attributes(managed) :: pz_pad,lambdaxy,a,b,c
    integer :: istat,ii,ng1,ng2
    !real(8), allocatable, dimension(:,:,:), device :: py_t
    !complex(8), allocatable, dimension(:,:,:), device :: pxc, pyc, pyc_t
#endif
    integer, dimension(3) :: ng
    integer :: q

    ng(:) = n(:)
    ng(1:2) = ng(1:2)*dims(1:2)
    if ( .not. allocated(px) ) allocate(px(ng(1),ng(2)/dims(1),ng(3)/dims(2)))
    if ( .not. allocated(py) ) allocate(py(ng(1)/dims(1),ng(2),ng(3)/dims(2)))
    if ( .not. allocated(pz) ) allocate(pz(ng(1)/dims(1),ng(2)/dims(2),ng(3)))
#ifdef USE_CUDA
    if ( .not. allocated(pw) ) allocate(pw(ng(1)/dims(2),ng(2)/dims(1),ng(3)))
#ifdef EPHC
    if ( .not. allocated(pxc  )) allocate( pxc( 2*(ng(1)/2 + 1), ng(2)/dims(1), ng(3)/dims(2) ) )
    if ( .not. allocated(pyc_t)) allocate( pyc_t( 2*(ng(2)/2 + 1), ng(1)/dims(1), ng(3)/dims(2) ) )
#else
    if ( .not. allocated(pxc  )) allocate( pxc( ng(1)/2 + 1, ng(2)/dims(1), ng(3)/dims(2) ) )
    if ( .not. allocated(pyc_t)) allocate( pyc_t( ng(2)/2 + 1, ng(1)/dims(1), ng(3)/dims(2) ) )
#endif
    if ( .not. allocated(py_t )) then
      allocate( py_t( ng(2), ng(1)/dims(1), ng(3)/dims(2) ) )
      istat = cudaMemAdvise( px, size(px), cudaMemAdviseSetPreferredLocation, mydev )
      istat = cudaMemAdvise( py, size(py), cudaMemAdviseSetPreferredLocation, mydev )
      istat = cudaMemAdvise( pz, size(pz), cudaMemAdviseSetPreferredLocation, mydev )
      istat = cudaMemAdvise( pw, size(pw), cudaMemAdviseSetPreferredLocation, mydev )
      istat = cudaMemPrefetchAsync( px, size(px), mydev, 0)
      istat = cudaMemPrefetchAsync( py, size(py), mydev, 0)
      istat = cudaMemPrefetchAsync( pz, size(pz), mydev, 0)
      istat = cudaMemPrefetchAsync( pw, size(pw), mydev, 0)
      istat = cudaMemPrefetchAsync( pz_pad, size(pz_pad), mydev, 0)
      istat = cudaMemPrefetchAsync(lambdaxy,size(lambdaxy),mydev,0)
      istat = cudaMemPrefetchAsync( a, size(a), mydev, 0)
      istat = cudaMemPrefetchAsync( b, size(b), mydev, 0)
      istat = cudaMemPrefetchAsync( c, size(c), mydev, 0)
    endif
#endif

#ifdef EPHC
    if(dims(1) * dims(2) .eq. 1) then

#define Y_B4_X
#ifdef Y_B4_X
    !$cuf kernel do(3) <<<*,(8,8,8)>>>
    do k=1,ng(3)/dims(2)
    do j=1,ng(1)/dims(1)
    do i=1,ng(2)
      py_t(i,j,k) = pz_pad(j,i,k)
    enddo
    enddo
    enddo

    istat = cufftExecD2Z(cufft_plan_fwd_y, py_t, pyc_t)

    ng2 = ng(2)
    !$cuf kernel do(3) <<<*,(8,8,8)>>>
    do k=1,ng(3)/dims(2)
    do i=1,ng(2)
      do j=1,ng(1)/dims(1)
        if( i .eq. 1 ) then
          px(j,i,k) = pyc_t(i,j,k)
        else
          px(j,i,k) = pyc_t(i+1,j,k)
        endif
      end do
    end do
    end do
    !
    !call transpose_y_to_x(py,px)
    !
    istat = cufftExecD2Z(cufft_plan_fwd_x, px, pxc)

    ng1 = ng(1)
    !$cuf kernel do(3) <<<*,*>>>
    do k=1,ng(3)/dims(2)
    do j=1,ng(2)/dims(1)
      do i=1,ng(1)
        if( i .eq. 1)  then
          pw(i,j,k) = pxc(i,j,k)
        else
          pw(i,j,k) = pxc(i+1,j,k)
        endif
      end do
    end do
    end do
    !
    !call transpose_x_to_z(px,pw)
    !
    q = 0
    if(c_or_f(3).eq.'f'.and.bcz(1).eq.'D') q = 1
    if(bcz(0)//bcz(1).eq.'PP') then
      call gaussel_periodic_gpu(ng(1)/dims(2),ng(2)/dims(1),n(3)-q,a,b,c,lambdaxy,pw,px,py,pxc,pyc_t)
    else
      call gaussel_gpu(         ng(1)/dims(2),ng(2)/dims(1),n(3)-q,a,b,c,lambdaxy,pw,px,py)
    endif
    !
    !call transpose_z_to_x(pw,px)
    !
    ng1 = ng(1)
    !$cuf kernel do(3) <<<*,*>>>
    do k=1,ng(3)/dims(2)
    do j=1,ng(2)/dims(1)
      do i=1,ng(1)
        if( i .eq. 1)  then
          pxc(i,j,k) = pw(i,j,k)
        else
          pxc(i+1,j,k) = pw(i,j,k)
        endif
      end do
    end do
    end do

    istat = cufftExecZ2D(cufft_plan_bwd_x, pxc, px)
    !
    !call transpose_x_to_y(px,py)
    !
    ng2 = ng(2)
    !$cuf kernel do(3) <<<*,*>>>
    do k=1,ng(3)/dims(2)
    do j=1,ng(1)/dims(1)
      do i=1,ng(2)
        if( i .eq. 1 ) then
          pyc_t(i,j,k) = px(j,i,k)
        else
          pyc_t(i+1,j,k) = px(j,i,k)
        endif
      end do
    end do
    end do

    istat = cufftExecZ2D(cufft_plan_bwd_y, pyc_t, py_t)
    !
    !$cuf kernel do(3) <<<*,(8,8,8)>>>
    do k=1,ng(3)
    do j=1,ng(2)/dims(2)
    do i=1,ng(1)/dims(1)
      pz_pad(i,j,k) = py_t(j,i,k)*normfft
    enddo
    enddo
    enddo
#else
    ! X first then Y (not as fast as Y first)
    !$cuf kernel do(3) <<<*,(8,8,8)>>>
    do k=1,ng(3)/dims(2)
    do j=1,ng(2)/dims(1)
    do i=1,ng(1)
      px(i,j,k) = pz_pad(i,j,k)
    enddo
    enddo
    enddo

    istat = cufftExecD2Z(cufft_plan_fwd_x, px, pxc)

    !$cuf kernel do(3) <<<*,(8,8,8)>>>
    do k=1,ng(3)/dims(2)
    do j=1,ng(1)/dims(1)
    do i=1,ng(2)
      if(j .eq. 1) then
        py_t(i,j,k) = pxc(j,i,k)
      else
        py_t(i,j,k) = pxc(j+1,i,k)
      endif
    enddo
    enddo
    enddo

    istat = cufftExecD2Z(cufft_plan_fwd_y, py_t, pyc_t)

    ng2 = ng(2)
    !$cuf kernel do(3) <<<*,(8,8,8)>>>
    do k=1,ng(3)/dims(2)
    do i=1,ng(2)
      do j=1,ng(1)/dims(1)
        if( i .eq. 1 ) then
          pw(j,i,k) = pyc_t(i,j,k)
        else
          pw(j,i,k) = pyc_t(i+1,j,k)
        endif
      end do
    end do
    end do
    !
    q = 0
    if(c_or_f(3).eq.'f'.and.bcz(1).eq.'D') q = 1
    if(bcz(0)//bcz(1).eq.'PP') then
      call gaussel_periodic_gpu(ng(1)/dims(2),ng(2)/dims(1),n(3)-q,a,b,c,lambdaxy,pw,px,py,pxc,pyc_t)
    else
      call gaussel_gpu(         ng(1)/dims(2),ng(2)/dims(1),n(3)-q,a,b,c,lambdaxy,pw,px,py)
    endif
    !
    ng2 = ng(2)
    !$cuf kernel do(3) <<<*,(8,8,8)>>>
    do k=1,ng(3)/dims(2)
    do j=1,ng(1)/dims(1)
      do i=1,ng(2)
        if( i .eq. 1 ) then
          pyc_t(i,j,k) = pw(j,i,k)
        else
          pyc_t(i+1,j,k) = pw(j,i,k)
        endif
      end do
    end do
    end do

    istat = cufftExecZ2D(cufft_plan_bwd_y, pyc_t, py_t)

    ng1 = ng(1)
    !$cuf kernel do(3) <<<*,(8,8,8)>>>
    do k=1,ng(3)/dims(2)
    do j=1,ng(2)/dims(1)
      do i=1,ng(1)
        if( i .eq. 1)  then
          pxc(i,j,k) = py_t(j,i,k)
        else
          pxc(i+1,j,k) = py_t(j,i,k)
        endif
      end do
    end do
    end do

    istat = cufftExecZ2D(cufft_plan_bwd_x, pxc, px)

    !$cuf kernel do(3) <<<*,*>>>
    do k=lbound(pz,3),ubound(pz,3)
    do j=lbound(pz,2),ubound(pz,2)
    do i=lbound(pz,1),ubound(pz,1)
      pz_pad(i,j,k) = px(i,j,k)*normfft
    enddo
    enddo
    enddo

#endif

    else ! if single rank

#endif

#ifdef USE_CUDA
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
    call transpose_z_to_y(pz,py)
    !
#ifdef USE_CUDA
    !$cuf kernel do(3) <<<*,*>>>
    do k=1,ng(3)/dims(2)
    do j=1,ng(1)/dims(1)
    do i=1,ng(2)
      py_t(i,j,k) = py(j,i,k)
    enddo
    enddo
    enddo

    istat = cufftExecD2Z(cufft_plan_fwd_y, py_t, pyc_t)

#ifdef EPHC
    ng2 = ng(2)
    !$cuf kernel do(3) <<<*,(8,8,8)>>>
    do k=1,ng(3)/dims(2)
    do i=1,ng(2)
      do j=1,ng(1)/dims(1)
        if( i .eq. 1 ) then
          py(j,i,k) = pyc_t(i,j,k)
        else
          py(j,i,k) = pyc_t(i+1,j,k)
        endif
      end do
    end do
    end do
#else
    ng2 = ng(2)
    !$cuf kernel do(3) <<<*,*>>>
    do k=1,ng(3)/dims(2)
    do j=1,ng(1)/dims(1)
      do i=1,ng(2)
        if( i .le. (ng2/2)+1 )  then
          py(j,i,k) = REAL( pyc_t(i,j,k) )
        else
          py( j, ng2 - (i - (ng2/2 + 2)),k) = AIMAG( pyc_t(i-(ng2/2),j,k) )
        endif
      end do
    end do
    end do
#endif

#else
    call fftd(arrplan(1,2),py) ! fwd transform in y
#endif
    !
    call transpose_y_to_x(py,px)
    !
#ifdef USE_CUDA
    istat = cufftExecD2Z(cufft_plan_fwd_x, px, pxc)

#ifdef EPHC
    ng1 = ng(1)
    !$cuf kernel do(3) <<<*,*>>>
    do k=1,ng(3)/dims(2)
    do j=1,ng(2)/dims(1)
      do i=1,ng(1)
        if( i .eq. 1)  then
          px(i,j,k) = pxc(i,j,k)
        else
          px(i,j,k) = pxc(i+1,j,k)
        endif
      end do
    end do
    end do
#else
    ng1 = ng(1)
    !$cuf kernel do(3) <<<*,*>>>
    do k=1,ng(3)/dims(2)
    do j=1,ng(2)/dims(1)
      do i=1,ng(1)
        if( i .le. (ng1/2)+1 )  then
          px(i,j,k) = REAL(pxc(i,j,k))
        else
          px( ng1 - (i - (ng1/2 + 2)),j,k) = AIMAG(pxc(i-(ng1/2),j,k))
        endif
      end do
    end do
    end do
#endif

#else
    call fftd(arrplan(1,1),px) ! fwd transform in x
#endif
    !
#ifdef USE_CUDA
    call transpose_x_to_z(px,pw)
#else
    call transpose_x_to_y(px,py)
    call transpose_y_to_z(py,pz)
#endif
    !
    q = 0
    if(c_or_f(3).eq.'f'.and.bcz(1).eq.'D') q = 1
    if(bcz(0)//bcz(1).eq.'PP') then
#ifdef USE_CUDA
      call gaussel_periodic_gpu(ng(1)/dims(2),ng(2)/dims(1),n(3)-q,a,b,c,lambdaxy,pw,px,py,pxc,pyc_t)
#else
      call gaussel_periodic(n(1),n(2),n(3)-q,a,b,c,lambdaxy,pz)
#endif
    else
#ifdef USE_CUDA
      call gaussel_gpu(         ng(1)/dims(2),ng(2)/dims(1),n(3)-q,a,b,c,lambdaxy,pw,px,py)
#else
      call gaussel(         n(1),n(2),n(3)-q,a,b,c,lambdaxy,pz)
#endif
    endif
    !
#ifdef USE_CUDA
    call transpose_z_to_x(pw,px)
#else
    call transpose_z_to_y(pz,py)
    call transpose_y_to_x(py,px)
#endif
    !
#ifdef USE_CUDA

#ifdef EPHC
    ng1 = ng(1)
    !$cuf kernel do(3) <<<*,*>>>
    do k=1,ng(3)/dims(2)
    do j=1,ng(2)/dims(1)
      do i=1,ng(1)
        if( i .eq. 1)  then
          pxc(i,j,k) = px(i,j,k)
        else
          pxc(i+1,j,k) = px(i,j,k)
        endif
      end do
    end do
    end do
#else
    ng1 = ng(1)
    !$cuf kernel do(3) <<<*,*>>>
    do k=1,ng(3)/dims(2)
    do j=1,ng(2)/dims(1)
    do i=1,(ng1/2)+1
       pxc(i,j,k)%re = px(i,j,k)
    end do
    end do
    end do

    !$cuf kernel do(3) <<<*,*>>>
    do k=1,ng(3)/dims(2)
    do j=1,ng(2)/dims(1)
    do i=(ng1/2)+2,ng1
       pxc(i-(ng1/2),j,k)%im = px( ng1 - (i - (ng1/2 + 2)),j,k)
    end do
    end do
    end do
#endif
    istat = cufftExecZ2D(cufft_plan_bwd_x, pxc, px)

#else
    call ffti(arrplan(2,1),px) ! bwd transform in x
#endif
    !
    call transpose_x_to_y(px,py)
    !
#ifdef USE_CUDA

#ifdef EPHC
    !pyc_t = 0.d0

    ng2 = ng(2)
    !$cuf kernel do(3) <<<*,*>>>
    do k=1,ng(3)/dims(2)
    do j=1,ng(1)/dims(1)
      do i=1,ng(2)
        if( i .eq. 1 ) then
          pyc_t(i,j,k) = py(j,i,k)
        else
          pyc_t(i+1,j,k) = py(j,i,k)
        endif
      end do
    end do
    end do
#else
    !pyc_t = (0.d0,0.d0)

    ng2 = ng(2)
    !$cuf kernel do(3) <<<*,*>>>
    do k=1,ng(3)/dims(2)
    do j=1,ng(1)/dims(1)
    do i=1,(ng2/2)+1
       pyc_t(i,j,k)%re = py(j,i,k)
    end do
    end do
    end do

    !$cuf kernel do(3) <<<*,*>>>
    do k=1,ng(3)/dims(2)
    do j=1,ng(1)/dims(1)
    do i=(ng2/2)+2,ng2
       pyc_t(i-(ng2/2),j,k)%im = py( j, ng2 - (i - (ng2/2 + 2)),k)
    end do
    end do
    end do
#endif
    istat = cufftExecZ2D(cufft_plan_bwd_y, pyc_t, py_t)

    !$cuf kernel do(3) <<<*,*>>>
    do k=1,ng(3)/dims(2)
    do i=1,ng(2)
    do j=1,ng(1)/dims(1)
      py(j,i,k) = py_t(i,j,k)
    enddo
    enddo
    enddo
#else
    call ffti(arrplan(2,2),py) ! bwd transform in y
#endif
    !
    call transpose_y_to_z(py,pz)
    !
#ifdef USE_CUDA
    !$cuf kernel do(3) <<<*,*>>>
#endif
    do k=lbound(pz,3),ubound(pz,3)
    do j=lbound(pz,2),ubound(pz,2)
    do i=lbound(pz,1),ubound(pz,1)
      pz_pad(i,j,k) = pz(i,j,k)*normfft
    enddo
    enddo
    enddo

#ifdef EPHC
    endif
#endif

    !deallocate(px,py,pz,pw)
#ifdef USE_CUDA
    !deallocate(pxc,pyc_t,py_t)
#endif

    return
  end subroutine solver
  !
#ifdef USE_CUDA
  subroutine gaussel_gpu(nx,ny,n,a,b,c,lambdaxy,p,bb,d)
    implicit none
    integer, intent(in) :: nx,ny,n
    real(8), intent(in), dimension(:), managed :: a,b,c
    real(8), intent(in), dimension(nx,ny), managed :: lambdaxy
    real(8), intent(inout), dimension(:,:,:), managed :: p
    real(8), intent(inout), dimension(nx,ny,n), device :: bb
    real(8), intent(inout), dimension(nx,ny,n), device :: d
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

  subroutine gaussel_periodic_gpu(nx,ny,n,a,b,c,lambdaxy,p,bb,d,p1,p2)
    implicit none
    integer, intent(in) :: nx,ny,n
    real(8), intent(in), dimension(:), managed :: a,b,c
    real(8), intent(in), dimension(nx,ny), managed :: lambdaxy
    real(8), intent(inout), dimension(:,:,:), managed :: p
    !DIR$ IGNORE_TKR p1,p2
    real(8), dimension(nx,ny,n), device :: bb,p1,p2,d
    real(8) :: z
    integer :: i,j,k
    !
    !solve tridiagonal system
    !

    !$cuf kernel do(2) <<<*,*>>>
    do j=1,ny
      do i=1,nx
        do k=1,n
          bb(i,j,k)  = b(k) + lambdaxy(i,j)
        enddo
        do k=1,n-1
          p1(i,j,k) = p(i,j,k)
        enddo

        !call dgtsv_homebrewed(n-1,a(1:n-1),bb(1:n-1),c(1:n-2),p1(1:n-1))
        z = 1.d0/bb(i,j,1)
        d(i,j,1) = c(1)*z
        p1(i,j,1) = p1(i,j,1)*z
        do k=2,n-2
          z = 1.d0/(bb(i,j,k)-a(k)*d(i,j,k-1))
          d(i,j,k) = c(k)*z
          p1(i,j,k) = (p1(i,j,k)-a(k)*p1(i,j,k-1))*z
        enddo
        z = bb(i,j,n-1)-a(n-1)*d(i,j,n-2)
        if(z.ne.0.d0) then
          p1(i,j,n-1) = (p1(i,j,n-1)-a(n-1)*p1(i,j,n-2))/z
        else
          p1(i,j,n-1) = 0.d0
        endif
        !
        ! backward substitution
        !
        do k=n-2,1,-1
          p1(i,j,k) = p1(i,j,k) - d(i,j,k)*p1(i,j,k+1)
        enddo


        do k=1,n
          p2(i,j,k) = 0.d0
        enddo

        p2(i,j,1  ) = -a(1  )
        p2(i,j,n-1) = -c(n-1)

        !call dgtsv_homebrewed(n-1,a(2:n-1),bb(1:n-1),c(1:n-2),p2(1:n-1))
        z = 1.d0/bb(i,j,1)
        d(i,j,1) = c(1)*z
        p2(i,j,1) = p2(i,j,1)*z
        do k=2,n-2
          z = 1.d0/(bb(i,j,k)-a(k+1)*d(i,j,k-1))
          d(i,j,k) = c(k)*z
          p2(i,j,k) = (p2(i,j,k)-a(k+1)*p2(i,j,k-1))*z
        enddo
        z = bb(i,j,n-1)-a(n)*d(i,j,n-2)
        if(z.ne.0.d0) then
          p2(i,j,n-1) = (p2(i,j,n-1)-a(n)*p2(i,j,n-2))/z
        else
          p2(i,j,n-1) = 0.d0
        endif
        !
        ! backward substitution
        !
        do k=n-2,1,-1
          p2(i,j,k) = p2(i,j,k) - d(i,j,k)*p2(i,j,k+1)
        enddo

        p(i,j,n) = (p(i,j,n) - c(n)*p1(i,j,1) - a(n)*p1(i,j,n-1)) / &
                   (bb(i,j,n) + c(n)*p2(i,j,1) + a(n)*p2(i,j,n-1))
        do k=1,n-1
          p(i,j,k) = p1(i,j,k) + p2(i,j,k)*p(i,j,n)
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    return
  end subroutine gaussel_periodic_gpu

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
