module mod_bound
  use mpi
  use mod_common_mpi, only: ierr,status,comm_cart,left,right,front,back,xhalo,yhalo
  use mod_common_mpi, only: xsl_buf, xrl_buf, ysl_buf, yrl_buf, xsr_buf, xrr_buf, ysr_buf, yrr_buf
  use mod_param     , only: dims
#ifdef USE_CUDA
  use cudafor
  use nvtx
#endif
  use mod_types
  implicit none
  private
  public boundp,bounduvw,updt_rhs_b
  contains
  subroutine bounduvw(cbc,n,bc,is_correc,dl,dzc,dzf,u,v,w)
    !
    ! imposes velocity boundary conditions
    !
    implicit none
    character(len=1), intent(in), dimension(0:1,3,3) :: cbc
    integer , intent(in), dimension(3) :: n 
    real(rp), intent(in), dimension(0:1,3,3) :: bc
    logical , intent(in)                   :: is_correc
    real(rp), intent(in), dimension(3) :: dl
    real(rp), intent(in), dimension(0:) :: dzc,dzf
    real(rp), intent(inout), dimension(0:,0:,0:) :: u,v,w
#ifdef USE_CUDA
    attributes(managed) :: u,v,w
#endif
    logical :: impose_norm_bc
    !
    call updthalo((/n(1),n(2)/),1,u)
    call updthalo((/n(1),n(2)/),2,u)
    call updthalo((/n(1),n(2)/),1,v)
    call updthalo((/n(1),n(2)/),2,v)
    call updthalo((/n(1),n(2)/),1,w)
    call updthalo((/n(1),n(2)/),2,w)
    !
    impose_norm_bc = (.not.is_correc).or.(cbc(0,1,1)//cbc(1,1,1).eq.'PP')
    if(left .eq.MPI_PROC_NULL) then
      if(impose_norm_bc) call set_bc(cbc(0,1,1),0,n(1),1,.false.,bc(0,1,1),dl(1),u)
                         call set_bc(cbc(0,1,2),0,n(1),1,.true. ,bc(0,1,2),dl(1),v)
                         call set_bc(cbc(0,1,3),0,n(1),1,.true. ,bc(0,1,3),dl(1),w)
    endif
    if(right.eq.MPI_PROC_NULL) then
      if(impose_norm_bc) call set_bc(cbc(1,1,1),1,n(1),1,.false.,bc(1,1,1),dl(1),u)
                         call set_bc(cbc(1,1,2),1,n(1),1,.true. ,bc(1,1,2),dl(1),v)
                         call set_bc(cbc(1,1,3),1,n(1),1,.true. ,bc(1,1,3),dl(1),w)
    endif
    impose_norm_bc = (.not.is_correc).or.(cbc(0,2,2)//cbc(1,2,2).eq.'PP')
    if(front.eq.MPI_PROC_NULL) then
                         call set_bc(cbc(0,2,1),0,n(2),2,.true. ,bc(0,2,1),dl(2),u)
      if(impose_norm_bc) call set_bc(cbc(0,2,2),0,n(2),2,.false.,bc(0,2,2),dl(2),v)
                         call set_bc(cbc(0,2,3),0,n(2),2,.true. ,bc(0,2,3),dl(2),w)
     endif
    if(back .eq.MPI_PROC_NULL) then
                         call set_bc(cbc(1,2,1),1,n(2),2,.true. ,bc(1,2,1),dl(2),u)
      if(impose_norm_bc) call set_bc(cbc(1,2,2),1,n(2),2,.false.,bc(1,2,2),dl(2),v)
                         call set_bc(cbc(1,2,3),1,n(2),2,.true. ,bc(1,2,3),dl(2),w)
    endif
    impose_norm_bc = (.not.is_correc).or.(cbc(0,3,3)//cbc(1,3,3).eq.'PP')
                       call set_bc(cbc(0,3,1),0,n(3),3,.true. ,bc(0,3,1),dzc(0)   ,u)
                       call set_bc(cbc(0,3,2),0,n(3),3,.true. ,bc(0,3,2),dzc(0)   ,v)
    if(impose_norm_bc) call set_bc(cbc(0,3,3),0,n(3),3,.false.,bc(0,3,3),dzf(0)   ,w)
                       call set_bc(cbc(1,3,1),1,n(3),3,.true. ,bc(1,3,1),dzc(n(3)),u)
                       call set_bc(cbc(1,3,2),1,n(3),3,.true. ,bc(1,3,2),dzc(n(3)),v)
    if(impose_norm_bc) call set_bc(cbc(1,3,3),1,n(3),3,.false.,bc(1,3,3),dzf(n(3)),w)
    return
  end subroutine bounduvw
  !
  subroutine boundp(cbc,n,bc,dl,dzc,dzf,p)
    !
    ! imposes pressure boundary conditions
    !
    implicit none
    character(len=1), intent(in), dimension(0:1,3) :: cbc
    integer , intent(in), dimension(3) :: n 
    real(rp)         , intent(in), dimension(0:1,3) :: bc
    real(rp), intent(in), dimension(3) :: dl
    real(rp), intent(in), dimension(0:) :: dzc,dzf
    real(rp), intent(inout), dimension(0:,0:,0:) :: p
#ifdef USE_CUDA
    attributes(managed) :: p
#endif
    !
    call updthalo((/n(1),n(2)/),1,p)
    call updthalo((/n(1),n(2)/),2,p)
    !
    if(left .eq.MPI_PROC_NULL) then
      call set_bc(cbc(0,1),0,n(1),1,.true.,bc(0,1),dl(1),p)
    endif
    if(right.eq.MPI_PROC_NULL) then
      call set_bc(cbc(1,1),1,n(1),1,.true.,bc(1,1),dl(1),p)
    endif
    if(front.eq.MPI_PROC_NULL) then
      call set_bc(cbc(0,2),0,n(2),2,.true.,bc(0,2),dl(2),p)
     endif
    if(back .eq.MPI_PROC_NULL) then
      call set_bc(cbc(1,2),1,n(2),2,.true.,bc(1,2),dl(2),p)
    endif
    call set_bc(cbc(0,3),0,n(3),3,.true.,bc(0,3),dzc(0)   ,p)
    call set_bc(cbc(1,3),1,n(3),3,.true.,bc(1,3),dzc(n(3)),p)
    return
  end subroutine boundp
  !
  subroutine set_bc(ctype,ibound,n,idir,centered,rvalue,dr,p)
    implicit none
    character(len=1), intent(in) :: ctype
    integer , intent(in) :: ibound,n,idir
    logical , intent(in) :: centered
    real(rp), intent(in) :: rvalue,dr
    real(rp), intent(inout), dimension(0:,0:,0:) :: p
    real(rp) :: factor,sgn
#ifdef USE_CUDA
    attributes(managed) :: p
    integer ::  i,j,k,istat
#endif
    !
    factor = rvalue
    if(ctype.eq.'D'.and.centered) then
      factor = 2.*factor
      sgn    = -1.
    endif
    if(ctype.eq.'N'.and.centered) then
      if(    ibound.eq.0) then
        factor = -dr*factor
      elseif(ibound.eq.1) then
        factor =  dr*factor
      endif
      sgn    = 1.
    endif
    !
    select case(ctype)
    case('P')
      select case(idir)
      case(1)
        !p(0  ,:,:) = p(n,:,:)
        !p(n+1,:,:) = p(1,:,:)
      case(2)
        !p(:,0  ,:) = p(:,n,:)
        !p(:,n+1,:) = p(:,1,:)
      case(3)
        #ifdef USE_CUDA
        !$cuf kernel do(2) <<<*,*>>>
        do j=lbound(p,2),ubound(p,2)
          do i=lbound(p,1),ubound(p,1)
            p(i,j,0  ) = p(i,j,n)
            p(i,j,n+1) = p(i,j,1)
          enddo
        enddo
        #else
        !$OMP WORKSHARE
        p(:,:,0  ) = p(:,:,n)
        p(:,:,n+1) = p(:,:,1)
        !$OMP END WORKSHARE
        #endif
      end select
    case('D','N')
      if(centered) then
        select case(idir)
        case(1)
          if    (ibound.eq.0) then
            #ifdef USE_CUDA
            !$cuf kernel do(2) <<<*,*>>>
            do k=lbound(p,3),ubound(p,3)
              do j=lbound(p,2),ubound(p,2)
                p(0  ,j,k) = factor+sgn*p(1,j,k)
              enddo
            enddo
            #else
            !$OMP WORKSHARE
             p(0  ,:,:) = factor+sgn*p(1,:,:)
            !$OMP END WORKSHARE
            #endif
          elseif(ibound.eq.1) then
            #ifdef USE_CUDA
            !$cuf kernel do(2) <<<*,*>>>
              do k=lbound(p,3),ubound(p,3)
                do j=lbound(p,2),ubound(p,2)
                  p(n+1,j,k) = factor+sgn*p(n,j,k)
                enddo
              enddo
            #else
            !$OMP WORKSHARE
            p(n+1,:,:) = factor+sgn*p(n,:,:)
            !$OMP END WORKSHARE
            #endif
          endif
        case(2)
          if    (ibound.eq.0) then
            #ifdef USE_CUDA
            !$cuf kernel do(2) <<<*,*>>>
            do k=lbound(p,3),ubound(p,3)
              do i=lbound(p,1),ubound(p,1)
                p(i,0  ,k) = factor+sgn*p(i,1,k)
              enddo
            enddo
            #else
            !$OMP WORKSHARE
            p(:,0  ,:) = factor+sgn*p(:,1,:)
            !$OMP END WORKSHARE
            #endif
          elseif(ibound.eq.1) then
            #ifdef USE_CUDA
            !$cuf kernel do(2) <<<*,*>>>
            do k=lbound(p,3),ubound(p,3)
              do i=lbound(p,1),ubound(p,1)
                p(i,n+1,k) = factor+sgn*p(i,n,k)
              enddo
            enddo
            #else
            !$OMP WORKSHARE
            p(:,n+1,:) = factor+sgn*p(:,n,:)
            !$OMP END WORKSHARE
            #endif
          endif
        case(3)
          if    (ibound.eq.0) then
            #ifdef USE_CUDA
            !$cuf kernel do(2) <<<*,*>>>
            do j=lbound(p,2),ubound(p,2)
              do i=lbound(p,1),ubound(p,1)
                p(i,j,0  ) = factor+sgn*p(i,j,1)
              enddo
            enddo
            #else
            !$OMP WORKSHARE
            p(:,:,0  ) = factor+sgn*p(:,:,1)
            !$OMP END WORKSHARE
            #endif
          elseif(ibound.eq.1) then
            #ifdef USE_CUDA
            !$cuf kernel do(2) <<<*,*>>>
            do j=lbound(p,2),ubound(p,2)
              do i=lbound(p,1),ubound(p,1)
                p(i,j,n+1) = factor+sgn*p(i,j,n)
              enddo
            enddo
            #else
            !$OMP WORKSHARE
            p(:,:,n+1) = factor+sgn*p(:,:,n)
            !$OMP END WORKSHARE
            #endif
          endif
        end select
      elseif(.not.centered.and.ctype.eq.'D') then
        select case(idir)
        case(1)
          if    (ibound.eq.0) then
            #ifdef USE_CUDA
            !$cuf kernel do(2) <<<*,*>>>
            do k=lbound(p,3),ubound(p,3)
              do j=lbound(p,2),ubound(p,2)
                p(0,j,k) = factor 
              enddo
            enddo
            #else
            !$OMP WORKSHARE
            p(0,:,:) = factor 
            !$OMP END WORKSHARE
            #endif
          elseif(ibound.eq.1) then
            #ifdef USE_CUDA
            !$cuf kernel do(2) <<<*,*>>>
            do k=lbound(p,3),ubound(p,3)
              do j=lbound(p,2),ubound(p,2)
                p(n,j,k) = factor 
                p(n+1,j,k) = p(n-1,j,k)
              enddo
            enddo
            #else
            !$OMP WORKSHARE
            p(n  ,:,:) = factor
            p(n+1,:,:) = p(n-1,:,:)
            !$OMP END WORKSHARE
            #endif
          endif
        case(2)
          if    (ibound.eq.0) then
            #ifdef USE_CUDA
            !$cuf kernel do(2) <<<*,*>>>
            do k=lbound(p,3),ubound(p,3)
              do i=lbound(p,1),ubound(p,1)
                p(i,0  ,k) = factor
              enddo
            enddo
            #else
            !$OMP WORKSHARE
            p(:,0,:) = factor 
            !$OMP END WORKSHARE
            #endif
          elseif(ibound.eq.1) then
            #ifdef USE_CUDA
            !$cuf kernel do(2) <<<*,*>>>
            do k=lbound(p,3),ubound(p,3)
              do i=lbound(p,1),ubound(p,1)
                p(i,n  ,k) = factor
                p(i,n+1,k) = p(i,n-1,k)
              enddo
            enddo
            #else
            !$OMP WORKSHARE
            p(:,n  ,:) = factor
            p(:,n+1,:) = p(:,n-1,:)
            !$OMP END WORKSHARE
            #endif
          endif
        case(3)
          if    (ibound.eq.0) then
            #ifdef USE_CUDA
            !$cuf kernel do(2) <<<*,*>>>
            do j=lbound(p,2),ubound(p,2)
              do i=lbound(p,1),ubound(p,1)
                p(i,j,0) = factor
              enddo
            enddo
            #else
            !$OMP WORKSHARE
            p(:,:,0) = factor 
            !$OMP END WORKSHARE
            #endif
          elseif(ibound.eq.1) then
            #ifdef USE_CUDA
            !$cuf kernel do(2) <<<*,*>>>
            do j=lbound(p,2),ubound(p,2)
              do i=lbound(p,1),ubound(p,1)
                p(i,j,n) = factor
                p(i,j,n+1) = p(i,j,n-1)
              enddo
            enddo
            #else
            !$OMP WORKSHARE
            p(:,:,n  ) = factor
            p(:,:,n+1) = p(:,:,n-1)
            !$OMP END WORKSHARE
            #endif
          endif
        end select
      elseif(.not.centered.and.ctype.eq.'N') then
        select case(idir)
        case(1)
          if    (ibound.eq.0) then
            #ifdef USE_CUDA
            !$cuf kernel do(2) <<<*,*>>>
            do k=lbound(p,3),ubound(p,3)
              do j=lbound(p,2),ubound(p,2)
                p(0,j,k) = 1.*factor + p(1  ,j,k)
              enddo
            enddo
            #else
            !$OMP WORKSHARE
            !p(0,:,:) = 1./3.*(-2.*factor+4.*p(1  ,:,:)-p(2  ,:,:))
            p(0,:,:) = 1.*factor + p(1  ,:,:)
            !$OMP END WORKSHARE
            #endif
          elseif(ibound.eq.1) then
            #ifdef USE_CUDA
            !$cuf kernel do(2) <<<*,*>>>
            do k=lbound(p,3),ubound(p,3)
              do j=lbound(p,2),ubound(p,2)
                p(n,j,k) = 1.*factor + p(n-1,j,k)
                p(n+1,j,k) = p(n,j,k) ! not needed
              enddo
            enddo
            #else
            !$OMP WORKSHARE
            !p(n,:,:) = 1./3.*(-2.*factor+4.*p(n-1,:,:)-p(n-2,:,:))
            p(n,:,:) = 1.*factor + p(n-1,:,:)
            p(n+1,:,:) = p(n,:,:) ! not needed
            !$OMP END WORKSHARE
            #endif
          endif
        case(2)
          if    (ibound.eq.0) then
            #ifdef USE_CUDA
            !$cuf kernel do(2) <<<*,*>>>
            do k=lbound(p,3),ubound(p,3)
              do i=lbound(p,1),ubound(p,1)
                p(i,0  ,k) = 1.*factor + p(i,1  ,k) 
              enddo
            enddo
            #else
            !$OMP WORKSHARE
            !p(:,0  ,:) = 1./3.*(-2.*factor+4.*p(:,1,:)-p(:,2  ,:))
            p(:,0,:) = 1.*factor + p(:,1  ,:)
            !$OMP END WORKSHARE
            #endif
          elseif(ibound.eq.1) then
            #ifdef USE_CUDA
            !$cuf kernel do(2) <<<*,*>>>
            do k=lbound(p,3),ubound(p,3)
              do i=lbound(p,1),ubound(p,1)
                p(i,n,k) = 1.*factor + p(i,n-1,k)
                p(i,n+1,k) = p(i,n,k) ! not needed
              enddo
            enddo
            #else
            !$OMP WORKSHARE
            !p(:,n,:) = 1./3.*(-2.*factor+4.*p(:,n-1,:)-p(:,n-2,:))
            p(:,n,:) = 1.*factor + p(:,n-1,:)
            p(:,n+1,:) = p(:,n,:) ! not needed
            !$OMP END WORKSHARE
            #endif
          endif
        case(3)
          if    (ibound.eq.0) then
            #ifdef USE_CUDA
            !$cuf kernel do(2) <<<*,*>>>
            do j=lbound(p,2),ubound(p,2)
              do i=lbound(p,1),ubound(p,1)
                p(i,j,0  ) = 1.*factor + p(i,j,1  )
              enddo
            enddo
            #else
            !$OMP WORKSHARE
            !p(:,:,0) = 1./3.*(-2.*factor+4.*p(:,:,1  )-p(:,:,2  ))
            p(:,:,0) = 1.*factor + p(:,:,1  )
            !$OMP END WORKSHARE
            #endif
          elseif(ibound.eq.1) then
            #ifdef USE_CUDA
            !$cuf kernel do(2) <<<*,*>>>
            do j=lbound(p,2),ubound(p,2)
              do i=lbound(p,1),ubound(p,1)
                p(i,j,n) = 1.*factor + p(i,j,n-1)
                p(i,j,n+1) = p(i,j,n) ! not needed
              enddo
            enddo
            #else
            !$OMP WORKSHARE
            !p(:,:,n) = 1./3.*(-2.*factor+4.*p(:,:,n-1)-p(:,:,n-2))
            p(:,:,n) = 1.*factor + p(:,:,n-1)
            p(:,:,n+1) = p(:,:,n) ! not needed
            !$OMP END WORKSHARE
            #endif
          endif
        end select
      endif
    end select
    !
    ! n.b. need to add this sync for pre-Pascal cards
    !      since a managed variable will be touched 
    !      on the host memory before synchronization
    !
    !@cuf istat=cudaDeviceSynchronize()
    return
  end subroutine set_bc
  !
  subroutine inflow(n,idir,dl,dzf,vel2d,u,v,w)
    implicit none
    integer, intent(in), dimension(3) :: n
    integer, intent(in) :: idir
    real(rp), intent(in), dimension(3) :: dl
    real(rp), intent(in), dimension(0:) :: dzf
    real(rp), dimension(0:,0:), intent(in) :: vel2d
    real(rp), dimension(0:,0:,0:), intent(inout) :: u,v,w
    real(rp) :: dx,dy,dxi,dyi
    real(rp), dimension(0:n(3)+1) :: dzfi
    integer :: i,j,k
    !
    dx   = dl(1)
    dxi  = dl(1)**(-1)
    dy   = dl(2)
    dyi  = dl(2)**(-1)
    dzfi = dzf**(-1)
    select case(idir)
      case(1) ! x direction
        if(left.eq.MPI_PROC_NULL) then
          i = 0
          do k=1,n(3)
            do j=1,n(2)
              u(i,j,k) = vel2d(j,k)
            enddo
          enddo 
        endif
      case(2) ! y direction
        j = 0
        if(front.eq.MPI_PROC_NULL) then
          do k=1,n(3)
            do i=1,n(1)
              v(i,j,k) = vel2d(i,k)
            enddo
          enddo 
        endif
      case(3) ! z direction
        k = 0
        do j=1,n(2)
          do i=1,n(1)
            w(i,j,k) = vel2d(i,j)
          enddo
        enddo 
    end select
    return
  end subroutine inflow
  !
  subroutine updt_rhs_b(c_or_f,cbc,n,rhsbx,rhsby,rhsbz,p)
    implicit none
    character, intent(in), dimension(3) :: c_or_f
    character(len=1), intent(in), dimension(0:1,3) :: cbc
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(:,:,0:) :: rhsbx,rhsby,rhsbz
    real(rp), intent(inout), dimension(0:,0:,0:) :: p
    integer , dimension(3) :: q
    integer :: idir
#ifdef USE_CUDA
    attributes(managed) :: p,rhsbx,rhsby,rhsbz
    integer:: i,j,k,ind
#endif
    q(:) = 0
    do idir = 1,3
      if(c_or_f(idir).eq.'f'.and.cbc(1,idir).eq.'D') q(idir) = 1
    enddo
    if(left.eq.MPI_PROC_NULL) then
      #ifdef USE_CUDA
      !$cuf kernel do(2) <<<*,*>>>
      do k=1,n(3)
        do j=1,n(2)
          p(1  ,j,k) = p(1  ,j,k) + rhsbx(j,k,0)
        enddo
      enddo
      #else
      !$OMP WORKSHARE
      p(1        ,1:n(2),1:n(3)) = p(1        ,1:n(2),1:n(3)) + rhsbx(:,:,0)
      !$OMP END WORKSHARE
      #endif
    endif  
    if(right.eq.MPI_PROC_NULL) then
      #ifdef USE_CUDA
      ind=n(1)-q(1)
      !$cuf kernel do(2) <<<*,*>>>
      do k=1,n(3)
        do j=1,n(2)
          p(ind,j,k) = p(ind,j,k) + rhsbx(j,k,1)
        enddo
      enddo
      #else
      !$OMP WORKSHARE
      p(n(1)-q(1),1:n(2),1:n(3)) = p(n(1)-q(1),1:n(2),1:n(3)) + rhsbx(:,:,1)
      !$OMP END WORKSHARE
      #endif
    endif
    if(front.eq.MPI_PROC_NULL) then
      #ifdef USE_CUDA
      !$cuf kernel do(2) <<<*,*>>>
      do k=1,n(3)
        do i=1,n(1)
          p(i,1  ,k) = p(i,1  ,k) + rhsby(i,k,0)
        enddo
      enddo
      #else
      !$OMP WORKSHARE
      p(1:n(1),1        ,1:n(3)) = p(1:n(1),1        ,1:n(3)) + rhsby(:,:,0)
      !$OMP END WORKSHARE
      #endif
    endif
    if(back.eq.MPI_PROC_NULL) then
      #ifdef USE_CUDA
      ind=n(2)-q(2)
      !$cuf kernel do(2) <<<*,*>>>
      do k=1,n(3)
        do i=1,n(1)
          p(i,ind,k) = p(i,ind,k) + rhsby(i,k,1)
        enddo
      enddo
      #else
      !$OMP WORKSHARE
      p(1:n(1),n(2)-q(2),1:n(3)) = p(1:n(1),n(2)-q(2),1:n(3)) + rhsby(:,:,1)
      !$OMP END WORKSHARE
      #endif
    endif
    #ifdef USE_CUDA
    ind=n(3)-q(3)
    !$cuf kernel do(2) <<<*,*>>>
    do j=1,n(2)
      do i=1,n(1)
        p(i,j,1  ) = p(i,j,1  ) + rhsbz(i,j,0)
        p(i,j,ind) = p(i,j,ind) + rhsbz(i,j,1)
      enddo
    enddo
    #else
    !$OMP WORKSHARE
    p(1:n(1),1:n(2),1        ) = p(1:n(1),1:n(2),1        ) + rhsbz(:,:,0)
    p(1:n(1),1:n(2),n(3)-q(3)) = p(1:n(1),1:n(2),n(3)-q(3)) + rhsbz(:,:,1)
    !$OMP END WORKSHARE
    #endif
    return
  end subroutine updt_rhs_b
  !
  subroutine updthalo(n,idir,p)
    implicit none
    integer , dimension(2), intent(in) :: n
    integer , intent(in) :: idir
    real(rp), dimension(0:,0:,0:), intent(inout) :: p
#ifdef USE_CUDA
    attributes(managed) :: p
    integer :: istat
#endif
    integer :: i,j,k,n_1,n_2
    !integer :: requests(4), statuses(MPI_STATUS_SIZE,4)
    !
    !  this subroutine updates the halos that store info
    !  from the neighboring computational sub-domain
    !
#ifdef USE_NVTX
      call nvtxStartRange("updthalo", 1)
#endif
    select case(idir)
    case(1) ! x direction
      !if( .false. ) then
      if(dims(1) .eq.  1) then
#ifdef USE_CUDA
        n_1=n(1)
        !$cuf kernel do(2) <<<*,*>>>
        do k=lbound(p,3),ubound(p,3)
          do j=lbound(p,2),ubound(p,2)
            p(n_1+1 ,j,k) = p(  1,j,k)
            p(0     ,j,k) = p(n_1,j,k)
          enddo
        enddo
#else
        !$OMP WORKSHARE
        p(n(1)+1,:,:) = p(   1,:,:)
        p(0     ,:,:) = p(n(1),:,:)
        !$OMP END WORKSHARE
#endif
      else
        n_1=n(1)
#ifdef USE_CUDA
        !$cuf kernel do(2) <<<*,*>>>
#endif
        do k=lbound(p,3),ubound(p,3)
          do j=lbound(p,2),ubound(p,2)
            xsl_buf(j,k) = p(  1,j,k)
            xsr_buf(j,k) = p(n_1,j,k)
          enddo
        enddo
        !@cuf istat = cudaDeviceSynchronize()
        !
        call MPI_SENDRECV(xsl_buf(0,0), size( xsl_buf ),MPI_REAL_RP,left ,0, &
                          xrr_buf(0,0), size( xrr_buf ),MPI_REAL_RP,right,0, &
                          comm_cart,status,ierr)
        call MPI_SENDRECV(xsr_buf(0,0), size( xsr_buf ),MPI_REAL_RP,right,0, &
                          xrl_buf(0,0), size( xrl_buf ),MPI_REAL_RP,left ,0, &
                          comm_cart,status,ierr)
#ifdef USE_CUDA
        !$cuf kernel do(2) <<<*,*>>>
#endif
        do k=lbound(p,3),ubound(p,3)
          do j=lbound(p,2),ubound(p,2)
            p(n_1+1,j,k) = xrr_buf(j,k)
            p(0    ,j,k) = xrl_buf(j,k)
          enddo
        enddo
        !@cuf istat = cudaDeviceSynchronize()
        !
        !call MPI_SENDRECV(p(1     ,0,0),1,xhalo,left ,0, &
        !                  p(n(1)+1,0,0),1,xhalo,right,0, &
        !                  comm_cart,status,ierr)
        !call MPI_SENDRECV(p(n(1),0,0),1,xhalo,right,0, &
        !                  p(0   ,0,0),1,xhalo,left ,0, &
        !                  comm_cart,status,ierr)
        !!call MPI_IRECV(p(0     ,0,0),1,xhalo,left ,1, &
        !!               comm_cart,requests(2),error)
        !!call MPI_IRECV(p(n(1)+1,0,0),1,xhalo,right,0, &
        !!               comm_cart,requests(1),error)
        !!call MPI_ISSEND(p(n(1),0,0),1,xhalo,right,1, &
        !!               comm_cart,requests(4),error)
        !!call MPI_ISSEND(p(1   ,0,0),1,xhalo,left ,0, &
        !!               comm_cart,requests(3),error)
        !!call MPI_WAITALL(4, requests, statuses, error)
      endif
    case(2) ! y direction
      !if( .false. ) then
      if( dims(2) .eq.  1 ) then
#ifdef USE_CUDA
        n_2=n(2)
        !$cuf kernel do(2) <<<*,*>>>
        do k=lbound(p,3),ubound(p,3)
          do i=lbound(p,1),ubound(p,1)
            p(i,n_2+1,k) = p(i,  1,k)
            p(i,    0,k) = p(i,n_2,k)
          enddo
        enddo
#else
        !$OMP WORKSHARE
        p(:,n(2)+1,:)  = p(:, 1  ,:)
        p(:, 0     ,:) = p(:,n(2),:)
        !$OMP END WORKSHARE
#endif
      else
        n_2=n(2)
#ifdef USE_CUDA
        !$cuf kernel do(2) <<<*,*>>>
#endif
        do k=lbound(p,3),ubound(p,3)
          do i=lbound(p,1),ubound(p,1)
            ysl_buf(i,k) = p(i,  1,k)
            ysr_buf(i,k) = p(i,n_2,k)
          enddo
        enddo
        !@cuf istat = cudaDeviceSynchronize()
        !
        call MPI_SENDRECV(ysl_buf(0,0), size( ysl_buf ),MPI_REAL_RP,front,0, &
                          yrr_buf(0,0), size( yrr_buf ),MPI_REAL_RP,back ,0, &
                          comm_cart,status,ierr)
        call MPI_SENDRECV(ysr_buf(0,0), size( ysr_buf ),MPI_REAL_RP,back ,0, &
                          yrl_buf(0,0), size( yrl_buf ),MPI_REAL_RP,front,0, &
                          comm_cart,status,ierr)
#ifdef USE_CUDA
        !$cuf kernel do(2) <<<*,*>>>
#endif
        do k=lbound(p,3),ubound(p,3)
          do i=lbound(p,1),ubound(p,1)
            p(i,n_2+1,k) = yrr_buf(i,k)
            p(i,    0,k) = yrl_buf(i,k)
          enddo
        enddo
        !@cuf istat = cudaDeviceSynchronize()
        !call MPI_SENDRECV(p(0,1     ,0),1,yhalo,front,0, &
        !                  p(0,n(2)+1,0),1,yhalo,back ,0, &
        !                  comm_cart,status,ierr)
        !call MPI_SENDRECV(p(0,n(2),0),1,yhalo,back ,0, &
        !                  p(0,0   ,0),1,yhalo,front,0, &
        !                  comm_cart,status,ierr)
        !!call MPI_IRECV(p(0,n(2)+1,0),1,yhalo,back ,0, &
        !!               comm_cart,requests(1),error)
        !!call MPI_IRECV(p(0,0     ,0),1,yhalo,front,1, &
        !!               comm_cart,requests(2),error)
        !!call MPI_ISSEND(p(0,1   ,0),1,yhalo,front,0, &
        !!               comm_cart,requests(3),error)
        !!call MPI_ISSEND(p(0,n(2),0),1,yhalo,back ,1, &
        !!               comm_cart,requests(4),error)
        !!call MPI_WAITALL(4, requests, statuses, error)
      endif
    end select
#ifdef USE_NVTX
      call nvtxEndRange
#endif
    return
  end subroutine updthalo
end module mod_bound
