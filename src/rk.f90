module mod_rk
  use mod_param, only: is_forced,velf
  use mod_debug, only: chkmean
  use mod_mom  , only: momxad,momyad,momzad,momxp,momyp,momzp
  use mod_momd , only: momxpd,momypd,momzpd,momxa,momya,momza
  use mod_moms , only: momsad
#ifdef USE_NVTX
  use nvtx
#endif 
  implicit none
  private
  public rk,rk_id
  contains
  subroutine rk(rkpar,n,dli,dzci,dzfi,dzflzi,dzclzi,visc,dt,l,u,v,w,p,dudtrko,dvdtrko,dwdtrko,tauxo,tauyo,tauzo,up,vp,wp,f)
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
    real(8),              dimension(n(1),n(2),n(3)) ::          dudtrk, dvdtrk, dwdtrk
    real(8) :: factor1,factor2,factor12
    real(8), dimension(3) :: taux,tauy,tauz
#ifdef USE_CUDA
   attributes(managed):: dzci,dzfi,dzflzi,dzclzi,u ,v ,w,p, dudtrk, dvdtrk, dwdtrk, up,vp,wp
   attributes(managed):: dudtrko,dvdtrko,dwdtrko
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
      call nvtxStartRange("momxp",1)
 #endif
    call momxp(n(1),n(2),n(3),dli(1),p,dudtrk)
 #ifdef USE_NVTX
      call nvtxEndRange
      call nvtxStartRange("momyp",2)
 #endif
    call momyp(n(1),n(2),n(3),dli(2),p,dvdtrk)
 #ifdef USE_NVTX
      call nvtxEndRange
      call nvtxStartRange("momzp",3)
 #endif
    call momzp(n(1),n(2),n(3),dzci  ,p,dwdtrk)
 #ifdef USE_NVTX
      call nvtxEndRange
 #endif

!@cuf istat=cudaDeviceSynchronize()

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
      call nvtxStartRange("momxad",4)
 #endif
    call momxad(n(1),n(2),n(3),dli(1),dli(2),dli(3),dzci,dzfi,dzflzi,visc,u,v,w,dudtrk,taux)
 #ifdef USE_NVTX
      call nvtxEndRange
      call nvtxStartRange("momyad",5)
 #endif
    call momyad(n(1),n(2),n(3),dli(1),dli(2),dli(3),dzci,dzfi,dzflzi,visc,u,v,w,dvdtrk,tauy)
 #ifdef USE_NVTX
      call nvtxEndRange
      call nvtxStartRange("momzad",6)
 #endif
    call momzad(n(1),n(2),n(3),dli(1),dli(2),dli(3),dzci,dzfi,dzflzi,visc,u,v,w,dwdtrk,tauz)
 #ifdef USE_NVTX
      call nvtxEndRange
 #endif

    f(1) = (factor1*sum(taux(:)/l(:)) + factor2*sum(tauxo(:)/l(:)))
    f(2) = (factor1*sum(tauy(:)/l(:)) + factor2*sum(tauyo(:)/l(:)))
    f(3) = (factor1*sum(tauz(:)/l(:)) + factor2*sum(tauzo(:)/l(:)))
    tauxo(:) = taux(:)
    tauyo(:) = tauy(:)
    tauzo(:) = tauz(:)
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
    return
  end subroutine rk
  subroutine rk_id(rkpar,n,dli,dzci,dzfi,dzflzi,dzclzi,visc,dt,l,u,v,w,p,dudtrko,dvdtrko,dwdtrko,tauxo,tauyo,tauzo,up,vp,wp,f)
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
    real(8),              dimension(n(1),n(2),n(3)) ::          dudtrk, dvdtrk, dwdtrk
    real(8),              dimension(n(1),n(2),n(3)) ::          dudtrkd, dvdtrkd, dwdtrkd
    real(8) :: factor1,factor2,factor12
    real(8), dimension(3) :: taux,tauy,tauz
#ifdef USE_CUDA
   attributes(managed):: dzci,dzfi,dzflzi,dzclzi,u ,v ,w,p, dudtrk, dvdtrk, dwdtrk, up,vp,wp
   attributes(managed):: dudtrko,dvdtrko,dwdtrko
   integer:: istat
#endif
    integer :: i,j,k
    real(8) :: mean
    !
    factor1 = rkpar(1)*dt
    factor2 = rkpar(2)*dt
    factor12 = factor1 + factor2
    !
    call momxpd(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,dzflzi,visc,p,u,dudtrk,dudtrkd,taux)
    call momypd(n(1),n(2),n(3),dli(2),dli(2),dzci,dzfi,dzflzi,visc,p,v,dvdtrk,dvdtrkd,tauy)
    call momzpd(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,dzflzi,visc,p,w,dwdtrk,dwdtrkd,tauz)
    f(1) = factor12*sum(taux(:)/l(:))
    f(2) = factor12*sum(tauy(:)/l(:))
    f(3) = factor12*sum(tauz(:)/l(:))
    ! alternatively, calculate force from the mean velocity directly
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(n,factor12,u,v,w,up,vp,wp,dudtrk,dvdtrk,dwdtrk)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          up(i,j,k) = u(i,j,k) + factor12*dudtrk(i,j,k)
          vp(i,j,k) = v(i,j,k) + factor12*dvdtrk(i,j,k)
          wp(i,j,k) = w(i,j,k) + factor12*dwdtrk(i,j,k)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    call momxa(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,u,v,w,dudtrk)
    call momya(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,u,v,w,dvdtrk)
    call momza(n(1),n(2),n(3),dli(1),dli(2),dzci,dzfi,u,v,w,dwdtrk)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(n,factor1,factor2,up,vp,wp,dudtrk,dvdtrk,dwdtrk,dudtrko,dvdtrko,dwdtrko)
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
    !$OMP END PARALLEL DO
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
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(n,factor12,factor2,visc,up,vp,wp,dudtrkd,dvdtrkd,dwdtrkd)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          up(i,j,k) = up(i,j,k) - .5d0*factor12*dudtrkd(i,j,k)
          vp(i,j,k) = vp(i,j,k) - .5d0*factor12*dvdtrkd(i,j,k)
          wp(i,j,k) = wp(i,j,k) - .5d0*factor12*dwdtrkd(i,j,k)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
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
