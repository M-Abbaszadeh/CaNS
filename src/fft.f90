module mod_fft
  use iso_c_binding , only: C_INT
  use mod_common_mpi, only: ierr
  use mod_fftw_param
  use mod_param     , only: pi,dims
  use mod_types
  !$ use omp_lib
  private
  public fftini,fftend,fft
#ifdef USE_CUDA
  public fftf_gpu,fftb_gpu,signal_processing
#endif
  contains
  subroutine fftini(nx,ny,nz,bcxy,c_or_f,arrplan,normfft)
    implicit none
    integer, intent(in) :: nx,ny,nz
    character(len=1), intent(in), dimension(0:1,2) :: bcxy
    character(len=1), intent(in), dimension(2) :: c_or_f
    type(C_PTR), intent(out), dimension(2,2) :: arrplan
    real(rp), intent(out) :: normfft
    real(rp), dimension(nx,ny/dims(1),nz/dims(2))  :: arrx
    real(rp), dimension(nx/dims(1),ny,nz/dims(2))  :: arry
    type(C_PTR) :: plan_fwd_x,plan_bwd_x, &
                   plan_fwd_y,plan_bwd_y
    type(fftw_iodim), dimension(1) :: iodim
    type(fftw_iodim), dimension(2) :: iodim_howmany
    integer :: kind_fwd,kind_bwd
    real(rp), dimension(2) :: norm
    integer(C_INT) :: nx_x,ny_x,nz_x, &
                      nx_y,ny_y,nz_y
    integer :: ix,iy
#ifdef USE_CUDA
    integer :: istat
    integer(int_ptr_kind()) :: worksize, max_worksize
    integer, pointer :: null_fptr
    call c_f_pointer( c_null_ptr, null_fptr )
    max_worksize = 0
#endif
#ifdef SINGLE_PRECISION
    !$ call sfftw_init_threads(ierr)
    !$ call sfftw_plan_with_nthreads(omp_get_max_threads())
#else
    !$ call dfftw_init_threads(ierr)
    !$ call dfftw_plan_with_nthreads(omp_get_max_threads())
#endif
    !
    ! fft in x
    !
    ! prepare plans with guru interface
    !
    nx_x = nx
    ny_x = ny/dims(1)
    nz_x = nz/dims(2)
    nx_y = nx/dims(1)
    ny_y = ny
    nz_y = nz/dims(2)
    !
    normfft = 1.
    ix = 0 
    ! size of transform reduced by 1 point with Dirichlet BC in face
    if(bcxy(0,1)//bcxy(1,1).eq.'DD'.and.c_or_f(1).eq.'f') ix = 1
    iodim(1)%n  = nx_x-ix
    iodim(1)%is = 1
    iodim(1)%os = 1
    iodim_howmany(1)%n  = ny_x
    iodim_howmany(1)%is = nx_x
    iodim_howmany(1)%os = nx_x
    iodim_howmany(2)%n  = nz_x
    iodim_howmany(2)%is = nx_x*ny_x
    iodim_howmany(2)%os = nx_x*ny_x
    call find_fft(bcxy(:,1),c_or_f(1),kind_fwd,kind_bwd,norm)
    plan_fwd_x=fftw_plan_guru_r2r(1,iodim,2,iodim_howmany,arrx,arrx,kind_fwd,FFTW_ESTIMATE)
    plan_bwd_x=fftw_plan_guru_r2r(1,iodim,2,iodim_howmany,arrx,arrx,kind_bwd,FFTW_ESTIMATE)
    normfft = normfft*norm(1)*(nx_x+norm(2)-ix)
#ifdef USE_CUDA
    if( .not. allocated( cufft_workspace ) ) then
      batch = ny_x*nz_x
      !istat = cufftPlan1D(cufft_plan_fwd_x,nx_x,CUFFT_FWD_TYPE,batch)
      istat = cufftCreate( cufft_plan_fwd_x )
      istat = cufftSetAutoAllocation( cufft_plan_fwd_x, 0 )
      istat = cufftMakePlanMany(cufft_plan_fwd_x,1,nx_x,null_fptr,1,nx_x,null_fptr,1,nx_x,CUFFT_FWD_TYPE,batch,worksize)
      max_worksize = max(worksize,max_worksize)
      !
      !istat = cufftPlan1D(cufft_plan_bwd_x,nx_x,CUFFT_BWD_TYPE,batch);
      istat = cufftCreate( cufft_plan_bwd_x )
      istat = cufftSetAutoAllocation( cufft_plan_bwd_x, 0 )
      istat = cufftMakePlanMany(cufft_plan_bwd_x,1,nx_x,null_fptr,1,nx_x,null_fptr,1,nx_x,CUFFT_BWD_TYPE,batch,worksize)
      max_worksize = max(worksize,max_worksize)
    endif
#endif
    !
    ! fft in y
    !
    ! prepare plans with guru interface
    !
    iy = 0
    ! size of transform reduced by 1 point with Dirichlet BC in face
    if(bcxy(0,2)//bcxy(1,2).eq.'DD'.and.c_or_f(2).eq.'f') iy = 1
    iodim(1)%n  = ny_y-iy
    iodim(1)%is = nx_y
    iodim(1)%os = nx_y
    iodim_howmany(1)%n  = nx_y
    iodim_howmany(1)%is = 1
    iodim_howmany(1)%os = 1
    iodim_howmany(2)%n  = nz_y
    iodim_howmany(2)%is = nx_y*ny_y
    iodim_howmany(2)%os = nx_y*ny_y
    call find_fft(bcxy(:,2),c_or_f(2),kind_fwd,kind_bwd,norm)
    plan_fwd_y=fftw_plan_guru_r2r(1,iodim,2,iodim_howmany,arry,arry,kind_fwd,FFTW_ESTIMATE)
    plan_bwd_y=fftw_plan_guru_r2r(1,iodim,2,iodim_howmany,arry,arry,kind_bwd,FFTW_ESTIMATE)
    normfft = normfft*norm(1)*(ny_y+norm(2)-iy)
#ifdef USE_CUDA
    if( .not. allocated( cufft_workspace ) ) then
      batch = nx_y*nz_y
      !istat = cufftPlan1D(cufft_plan_fwd_y,ny_y,CUFFT_FWD_TYPE,batch)
      istat = cufftCreate( cufft_plan_fwd_y )
      istat = cufftSetAutoAllocation( cufft_plan_fwd_y, 0 )
      istat = cufftMakePlanMany(cufft_plan_fwd_y,1,ny_y,null_fptr,1,ny_y,null_fptr,1,ny_y,CUFFT_FWD_TYPE,batch,worksize)
      max_worksize = max(worksize,max_worksize)
      !
      !istat = cufftPlan1D(cufft_plan_bwd_y,ny_y,CUFFT_BWD_TYPE,batch)
      istat = cufftCreate( cufft_plan_bwd_y )
      istat = cufftSetAutoAllocation( cufft_plan_bwd_y, 0 )
      istat = cufftMakePlanMany(cufft_plan_bwd_y,1,ny_y,null_fptr,1,ny_y,null_fptr,1,ny_y,CUFFT_BWD_TYPE,batch,worksize)
      max_worksize = max(worksize,max_worksize)
      !
      allocate(cufft_workspace(max_worksize/(2*sizeof(1._rp))))
      !
      istat = cufftSetWorkArea( cufft_plan_fwd_x, cufft_workspace )
      istat = cufftSetWorkArea( cufft_plan_bwd_x, cufft_workspace )
      istat = cufftSetWorkArea( cufft_plan_fwd_y, cufft_workspace )
      istat = cufftSetWorkArea( cufft_plan_bwd_y, cufft_workspace )
    endif
#endif
    !
    normfft = normfft**(-1)
    arrplan(1,1) = plan_fwd_x
    arrplan(2,1) = plan_bwd_x
    arrplan(1,2) = plan_fwd_y
    arrplan(2,2) = plan_bwd_y
    return
  end subroutine fftini
  !
  subroutine fftend(arrplan)
    implicit none
    type(C_PTR), intent(in), dimension(2,2) :: arrplan
#ifdef SINGLE_PRECISION
    call sfftw_destroy_plan(arrplan(1,1))
    call sfftw_destroy_plan(arrplan(1,2))
    call sfftw_destroy_plan(arrplan(2,1))
    call sfftw_destroy_plan(arrplan(2,2))
    !$call sfftw_cleanup_threads(ierr)
#else
    call dfftw_destroy_plan(arrplan(1,1))
    call dfftw_destroy_plan(arrplan(1,2))
    call dfftw_destroy_plan(arrplan(2,1))
    call dfftw_destroy_plan(arrplan(2,2))
    !$call dfftw_cleanup_threads(ierr)
#endif
    return
  end subroutine fftend
  !
  subroutine fft(plan,arr)
    implicit none
    type(C_PTR), intent(in) :: plan 
    real(rp), intent(inout), dimension(:,:,:) :: arr
#ifdef SINGLE_PRECISION
    call sfftw_execute_r2r(plan,arr,arr)
#else
    call dfftw_execute_r2r(plan,arr,arr)
#endif
    return
  end subroutine fft
  !
#ifdef USE_CUDA
  subroutine fftf_gpu(plan,arrin,arrout)
    implicit none
    integer , intent(in) :: plan 
    real(rp), intent(in ), dimension(:,:,:), device :: arrin
    real(rp), intent(out), dimension(:,:,:), device :: arrout
    integer :: istat
#ifdef SINGLE_PRECISION
    istat = cufftExecR2C(plan,arrin,arrout)
#else
    istat = cufftExecD2Z(plan,arrin,arrout)
#endif
    return
  end subroutine fftf_gpu
  subroutine fftb_gpu(plan,arrin,arrout)
    implicit none
    integer , intent(in) :: plan 
    real(rp), intent(in ), dimension(:,:,:), device :: arrin
    real(rp), intent(out), dimension(:,:,:), device :: arrout
    integer :: istat
#ifdef SINGLE_PRECISION
    istat = cufftExecC2R(plan,arrin,arrout)
#else
    istat = cufftExecZ2D(plan,arrin,arrout)
#endif
    return
  end subroutine fftb_gpu
#endif
  !
  subroutine find_fft(bc,c_or_f,kind_fwd,kind_bwd,norm)
  implicit none
  character(len=1), intent(in), dimension(0:1) :: bc
  character(len=1), intent(in) :: c_or_f
  integer , intent(out) :: kind_fwd,kind_bwd
  real(rp), intent(out), dimension(2) :: norm
  if(c_or_f.eq.'c') then
    select case(bc(0)//bc(1))
    case('PP')
      kind_fwd = FFTW_R2HC
      kind_bwd = FFTW_HC2R
      norm = (/1.,0./)
    case('NN')
      kind_fwd = FFTW_REDFT10
      kind_bwd = FFTW_REDFT01
      norm = (/2.,0./)
    case('DD')
      kind_fwd = FFTW_RODFT10
      kind_bwd = FFTW_RODFT01
      norm = (/2.,0./)
    case('ND')
      kind_fwd = FFTW_REDFT11
      kind_bwd = FFTW_REDFT11
      norm = (/2.,0./)
    case('DN')
      kind_fwd = FFTW_RODFT11
      kind_bwd = FFTW_RODFT11
      norm = (/2.,0./)
    end select
  elseif(c_or_f.eq.'f') then
    select case(bc(0)//bc(1))
    case('PP')
      kind_fwd = FFTW_R2HC
      kind_bwd = FFTW_HC2R
      norm = (/1.,0./)
    case('NN')
      kind_fwd = FFTW_REDFT00
      kind_bwd = FFTW_REDFT00
      norm = (/2.,-1./)
    case('DD')
      kind_fwd = FFTW_RODFT00
      kind_bwd = FFTW_RODFT00
      norm = (/2.,1./)
    case('ND')
      kind_fwd = FFTW_REDFT10
      kind_bwd = FFTW_REDFT01
      norm = (/2.,0./)
    case('DN')
      kind_fwd = FFTW_RODFT01
      kind_bwd = FFTW_RODFT10
      norm = (/2.,0./)
    end select
  endif
  return
  end subroutine find_fft
#ifdef USE_CUDA
  subroutine posp_fftf(n,idir,arr)
    !
    ! post-processing of a signal following a forward FFT
    ! to order the data as follows:
    ! (r[0],r[n],r[1],i[1],...,r[n-1],i[n-1])
    !
    implicit none
    integer , intent(in   ), dimension(3) :: n       ! dimensions of input/output array
    integer , intent(in   ) :: idir                  ! direction where the transform is taken
    real(rp), intent(inout), dimension(:,:,:) :: arr ! input/output array
    integer :: j,k,nn
    attributes(device) :: arr
    !
    nn = n(idir)-2
    select case(idir)
    case(1)
      !
      !$cuf kernel do(2) <<<*,*>>>
      do k=1,n(3)
        do j=1,n(2)
          arr(2,j,k) = arr(nn+1,j,k)
        enddo
      enddo
    end select
    return
  end subroutine posp_fftf
  subroutine prep_fftb(n,idir,arr)
    !
    ! pre-processing of a signal preciding a backward FFT
    ! to order the data as follows:
    ! (r[0],i[0],r[1],i[1],...,r[n-1],i[n-1],r[n],i[n])
    ! note that i[0] = i[n] = 0
    !
    implicit none
    integer , intent(in   ), dimension(3) :: n       ! dimensions of input/output array
    integer , intent(in   ) :: idir                  ! direction where the transform is taken
    real(rp), intent(inout), dimension(:,:,:) :: arr ! input/output array
    integer :: j,k,nn
    attributes(device) :: arr
    !
    nn = n(idir)-2
    select case(idir)
    case(1)
      !$cuf kernel do(2) <<<*,*>>>
      do k=1,n(3)
        do j=1,n(2)
          arr(nn+1,j,k) = arr(2,j,k)
          arr(2   ,j,k) = 0.
        enddo
      enddo
    end select
    return
  end subroutine prep_fftb
  subroutine prep_dctiif(n,idir,arr)
    !
    ! pre-processing of a signal to perform a fast forward cosine transform
    ! with FFTs (see Makhoul 1980)
    ! 
    ! the input signal x(n) is pre-processed into a signal v(n)
    ! as follows:
    !
    ! v(n) = x(2n       ),              0 <= n <= floor((N-1)/2)
    !      = x(2N -2n -1), floor((N+1)/2) <= n <= N-1
    ! with n = 0,...,N-1 and N is the total number of elements 
    ! of the signal.
    !
    implicit none
    integer , intent(in   ), dimension(3) :: n       ! dimensions of input/output array
    integer , intent(in   ) :: idir                  ! array direction where the transform is taken
    real(rp), intent(inout), dimension(:,:,:) :: arr ! input/output array
    real(rp), allocatable, dimension(:,:,:) :: arr_tmp   ! needs to be alocatable actually!
    integer :: i,j,k,nn,ii
    attributes(device) :: arr
    attributes(device) :: arr_tmp
    !
    nn = n(idir)
    select case(idir)
    case(1)
      if(.not.allocated(arr_tmp)) allocate(arr_tmp(0:n(idir)-1,n(2),n(3)))
      !$cuf kernel do(2) <<<*,*>>>
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            ii = i-1
            if(    ii.le.(nn-1)/2) then
              arr_tmp(ii,j,k) = arr(2*ii+1       ,j,k)
            else
              arr_tmp(ii,j,k) = arr(2*(nn-ii)-1+1,j,k)
            endif
          enddo
          do i=1,n(1)
            ii = i-1
            arr(i,j,k) = arr_tmp(ii,j,k)
          enddo
        enddo
      enddo
      !deallocate(arr_tmp)
    end select
    return
  end subroutine prep_dctiif
  subroutine posp_dctiif(n,idir,arr)
    ! 
    ! post-processing of a signal to perform a fast cosine transform
    ! with FFTs (see Makhoul 1980)
    !
    implicit none
    integer , intent(in   ), dimension(3) :: n 
    integer , intent(in   ) :: idir
    real(rp), intent(inout), dimension(:,:,:) :: arr
    real(rp), allocatable, dimension(:,:,:) :: arr_tmp
    integer :: i,j,k,ii,nn
    attributes(device) :: arr
    attributes(device) :: arr_tmp
    real(rp) :: arg
    !
    nn = n(idir)-2
    select case(idir)
    case(1)
      if(.not.allocated(arr_tmp)) allocate(arr_tmp(0:n(idir)-1+1,n(2),n(3)))
      !$cuf kernel do(2) <<<*,*>>>
      do k=1,n(3)
        do j=1,n(2)
          do i=1,nn+2,2
            ii = (i-1)/2
            !arr_tmp(ii   ,j,k) =    real( &
            !                         2.*exp(-ri_unit*pi*ii/(2.*nn))*cmplx(arr(i,j,k),arr(i+1,j,k),rp) &
            !                        )
            !arr_tmp(nn-ii,j,k) = - aimag( &
            !                         2.*exp(-ri_unit*pi*ii/(2.*nn))*cmplx(arr(i,j,k),arr(i+1,j,k),rp) &
            !                        ) ! = 0 for ii=0
            arg = -pi*ii/(2.*nn)
            arr_tmp(ii   ,j,k) =  2.*(cos(arg)*arr(i,j,k) - sin(arg)*arr(i+1,j,k))
            arr_tmp(nn-ii,j,k) = -2.*(sin(arg)*arr(i,j,k) + cos(arg)*arr(i+1,j,k))
          enddo
          do i=1,nn
            ii = i-1
            arr(i,j,k) = arr_tmp(ii,j,k)
          enddo
        enddo
      enddo
      !deallocate(arr_tmp)
    end select
    return
  end subroutine posp_dctiif
  subroutine prep_dctiib(n,idir,arr)
    !
    ! pre-processing of a signal to perform a fast backward cosine transform
    ! with FFTs (see Makhoul 1980)
    !
    implicit none
    integer , intent(in   ), dimension(3) :: n 
    integer , intent(in   ) :: idir
    real(rp), intent(inout), dimension(:,:,:) :: arr
    real(rp), allocatable, dimension(:,:,:) :: arr_tmp
    integer :: i,j,k,nn,ii
    attributes(device) :: arr
    attributes(device) :: arr_tmp
    real(rp) :: arg
    !
    nn = n(idir)-2
    select case(idir)
    case(1)
      if(.not.allocated(arr_tmp)) allocate(arr_tmp(0:n(idir)-1+1,n(2),n(3)))
      !$cuf kernel do(2) <<<*,*>>>
      do k=1,n(3)
        do j=1,n(2)
          arr(nn+1,j,k) = 0.
          arr(nn+2,j,k) = 0.
          do i=1,nn+2,2
            ii = (i-1)/2
            !arr_tmp(2*ii  ,j,k)  = real( 1.*exp(ri_unit*pi*ii/(2.*nn))*(arr(ii+1,j,k)-ri_unit*arr(nn-ii+1,j,k)))
            !arr_tmp(2*ii+1,j,k)  = aimag(1.*exp(ri_unit*pi*ii/(2.*nn))*(arr(ii+1,j,k)-ri_unit*arr(nn-ii+1,j,k)))
            arg = pi*ii/(2.*nn)
            arr_tmp(2*ii  ,j,k) = 1.*(cos(arg)*arr(ii+1,j,k) + sin(arg)*arr(nn-ii+1,j,k))
            arr_tmp(2*ii+1,j,k) = 1.*(sin(arg)*arr(ii+1,j,k) - cos(arg)*arr(nn-ii+1,j,k))
          enddo
          do i=1,nn+2
            ii = i-1
            arr(i,j,k) = arr_tmp(ii,j,k)
          enddo
        enddo
      enddo
      !deallocate(arr_tmp)
    end select
    return
  end subroutine prep_dctiib
  subroutine posp_dctiib(n,idir,arr)
    !
    ! post-processing of a signal to perform a fast forward cosine transform
    ! with FFTs (see Makhoul 1980)
    ! 
    ! the input signal v(n) is post-processed into a signal x(n)
    ! as follows:
    !
    ! v(n) = x(2n       ),              0 <= n <= floor((N-1)/2)
    !      = x(2N -2n -1), floor((N+1)/2) <= n <= N-1
    ! with n = 0,...,N-1 and N is the total number of elements 
    ! of the signal.
    !
    implicit none
    integer , intent(in   ), dimension(3) :: n       ! dimensions of input/output array
    integer , intent(in   ) :: idir                  ! array direction where the transform is taken
    real(rp), intent(inout), dimension(:,:,:) :: arr ! input/output array
    real(rp), allocatable, dimension(:,:,:) :: arr_tmp ! needs to be alocatable actually!
    integer :: i,j,k,nn,ii
    attributes(device) :: arr
    attributes(device) :: arr_tmp
    !
    nn = n(idir)
    select case(idir)
    case(1)
      if(.not.allocated(arr_tmp)) allocate(arr_tmp(0:n(idir)-1,n(2),n(3)))
      !$cuf kernel do(2) <<<*,*>>>
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            ii = i-1
            arr_tmp(ii,j,k) = arr(i,j,k)
          enddo
          do i=1,n(1)
            ii = i-1
            if(    ii.le.(nn-1)/2) then
              arr(2*ii+1       ,j,k) = arr_tmp(ii,j,k)
            else
              arr(2*(nn-ii)-1+1,j,k) = arr_tmp(ii,j,k)
            endif
          enddo
        enddo
      enddo
      !deallocate(arr_tmp)
    end select
    return
  end subroutine posp_dctiib
  subroutine signal_processing(pre_or_pos,f_or_b,cbc,c_or_f,n,idir,arr)
    implicit none
    !
    ! wrapper subroutine for signal processing to compute FFT-based transforms
    ! (can also be done with pointers to a subroutine like in initgrid.f90)
    !
    integer,          intent(in) :: pre_or_pos ! prior (0) or after (1) fft
    character(len=1), intent(in) :: f_or_b     ! forward or backward transform
    character(len=2), intent(in) :: cbc        ! type of boundary condition
    character(len=1), intent(in) :: c_or_f     ! cell- or face-centred BC?
    integer, intent(in), dimension(3)         :: n
    integer, intent(in)                       :: idir
    real(rp), intent(inout), dimension(:,:,:) :: arr
    attributes(device) :: arr
    integer :: istat
    select case(cbc)
    case('PP')
        select case(f_or_b)
        case('F')
          if(    pre_or_pos.eq.0) then
            return
          elseif(pre_or_pos.eq.1) then
            call posp_fftf(n,idir,arr)
          else
          endif
        case('B')
          if(    pre_or_pos.eq.0) then
            call prep_fftb(n,idir,arr)
          elseif(pre_or_pos.eq.1) then
            return
          else
          endif
        end select
      return
    case('NN')
      if(c_or_f.eq.'c') then
        select case(f_or_b)
        case('F')
          if(    pre_or_pos.eq.0) then
            call prep_dctiif(n,idir,arr)
          elseif(pre_or_pos.eq.1) then
            call posp_dctiif(n,idir,arr)
          else
          endif
        case('B')
          if(    pre_or_pos.eq.0) then
            call prep_dctiib(n,idir,arr)
          elseif(pre_or_pos.eq.1) then
            call posp_dctiib(n,idir,arr)
          else
          endif
        case default
        end select
      endif
    case default
      ! ERROR  (trap this in sanity check in sanity.f90)
    end select
    return
  end subroutine signal_processing
#endif
end module mod_fft
