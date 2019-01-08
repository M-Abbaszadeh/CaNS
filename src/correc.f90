module mod_correc
  implicit none
  private
  public correc
  contains
  subroutine correc(n,dli,dzci,dt,p,up,vp,wp,u,v,w)
    !
    ! corrects the velocity so that it is divergence free
    !
    !@cuf use cudafor
    implicit none
    integer, intent(in), dimension(3) :: n
    real(8), intent(in), dimension(3) :: dli
    real(8), intent(in), dimension(0:) :: dzci
    real(8), intent(in) :: dt
    real(8), intent(in) , dimension(0:,0:,0:) :: p,up,vp,wp
    real(8), intent(out), dimension(0:,0:,0:) :: u,v,w
    real(8) :: factori,factorj
#ifdef USE_CUDA
    attributes(managed):: p,up,vp,wp,u,v,w,dzci
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
    !$OMP SHARED(n,factori,factorj,dzci,dt,u,v,w,up,vp,wp,p) &
    !$OMP PRIVATE(i,j,k,ip,jp,kp)
#endif
    do k=1,n(3)
      kp = k+1
      do j=1,n(2)
        jp = j+1
        do i=1,n(1)
          ip = i+1
          u(i,j,k) = up(i,j,k) - factori*(   p(ip,j,k)-p(i,j,k))
          v(i,j,k) = vp(i,j,k) - factorj*(   p(i,jp,k)-p(i,j,k))
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
