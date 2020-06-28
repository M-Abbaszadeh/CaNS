module mod_correc
  use mod_types
  implicit none
  private
  public correc
  contains
  subroutine correc(n,dli,dzci,dt,p,up,vp,wp,u,v,w)
    !@cuf use cudafor
    !
    ! corrects the velocity so that it is divergence free
    !
    implicit none
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(3) :: dli
    real(rp), intent(in), dimension(0:) :: dzci
    real(rp), intent(in) :: dt
    real(rp), intent(in) , dimension(0:,0:,0:) :: p,up,vp,wp
    real(rp), intent(out), dimension(0:,0:,0:) :: u,v,w
    real(rp) :: factori,factorj
    !real(rp), dimension(0:n(3)+1) :: factork
#ifdef USE_CUDA
    attributes(managed) :: p,up,vp,wp,u,v,w,dzci
    integer:: istat
#endif
    integer :: i,j,k,ip,jp,kp
    !
    !factor = rkcoeffab(rkiter)*dt
    !
    factori = dt*dli(1)
    factorj = dt*dli(2)
#ifdef USE_CUDA
    !$cuf kernel do(3) <<<*,*>>>
#else
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(n,factori,u,up,p) &
    !$OMP PRIVATE(i,j,k,ip)
#endif
    do k=0,n(3)+1
      do j=0,n(2)+1
        do i=0,n(1)
          ip = i+1
          u(i,j,k) = up(i,j,k) - factori*(   p(ip,j,k)-p(i,j,k))
        enddo
      enddo
    enddo
#ifndef USE_CUDA
    !$OMP END PARALLEL DO
#endif
#ifdef USE_CUDA
    !$cuf kernel do(3) <<<*,*>>>
#else
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(n,factorj,v,vp,p) &
    !$OMP PRIVATE(i,j,k,jp)
#endif
    do k=0,n(3)+1
      do j=0,n(2)
        jp = j+1
        do i=0,n(1)+1
          v(i,j,k) = vp(i,j,k) - factorj*(   p(i,jp,k)-p(i,j,k))
        enddo
      enddo
    enddo
#ifndef USE_CUDA
    !$OMP END PARALLEL DO
#endif
#ifdef USE_CUDA
    !$cuf kernel do(3) <<<*,*>>>
#else
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(n,dt,dzci,w,wp,p) &
    !$OMP PRIVATE(i,j,k,kp)
#endif
    do k=0,n(3)
      kp = k+1
      do j=0,n(2)+1
        do i=0,n(1)+1
          w(i,j,k) = wp(i,j,k) - dt*dzci(k)*(p(i,j,kp)-p(i,j,k))
        enddo
      enddo
    enddo
#ifndef USE_CUDA
    !$OMP END PARALLEL DO
#endif
    !@cuf istat=cudaDeviceSynchronize()
    return
  end subroutine correc
end module mod_correc
