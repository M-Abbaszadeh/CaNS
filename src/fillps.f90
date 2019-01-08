module mod_fillps
  implicit none
  private
  public fillps
  contains
  subroutine fillps(n,dli,dzfi,dti,up,vp,wp,p)
    !
    !  fill the right-hand side of the Poisson equation for the correction pressure.
    !
    !  the discrete divergence is:
    !
    !  w(i,j,k)-w(i,j,k-1)   v(i,j,k)-v(i,j-1,k)   u(i,j,k)-u(i-1,j,k)
    !  ------------------- + ------------------- + -------------------  = div
    !          dz                    dy                    dx
    !
    !@cuf use cudafor
    implicit none
    integer, intent(in), dimension(3) :: n
    real(8), intent(in), dimension(3) :: dli
    real(8), intent(in), dimension(0:) :: dzfi
    real(8), intent(in) :: dti
    real(8), intent(in ), dimension(0:,0:,0:) :: up,vp,wp
    real(8), intent(out), dimension(0:,0:,0:) :: p
    real(8) :: dtidxi,dtidyi!,dtidzi
#ifdef USE_CUDA
    attributes(managed):: up,vp,wp,p,dtidzfi,dzfi
    integer:: istat
#endif
    integer :: i,j,k,im,jm,km
    !
    dtidxi = dti*dli(1)
    dtidyi = dti*dli(2)
#ifdef USE_CUDA
   !$cuf kernel do(3) <<<*,*>>>
#else
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(n,p,up,vp,wp,dti,dzfi,dtidyi,dtidxi) &
    !$OMP PRIVATE(i,j,k,im,jm,km)
#endif
    do k=1,n(3)
      km = k-1
      do j=1,n(2)
        jm = j-1
        do i=1,n(1)
          im = i-1
          p(i,j,k) = ( &
                      (wp(i,j,k)-wp(i,j,km))*dti*dzfi(k)+ &
                      (vp(i,j,k)-vp(i,jm,k))*dtidyi    + &
                      (up(i,j,k)-up(im,j,k))*dtidxi    )
        enddo
      enddo
    enddo
#ifndef USE_CUDA
    !$OMP END PARALLEL DO
#endif
!@cuf istat=cudaDeviceSynchronize()
    !
    return
  end subroutine fillps
end module mod_fillps
