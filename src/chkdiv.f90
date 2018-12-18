module mod_chkdiv
  use mpi
  use mod_common_mpi, only: myid,coord,ierr
  implicit none
  private
  public chkdiv
  contains
  subroutine chkdiv(n,dli,dzfi,u,v,w,divtot,divmax)
    !
    ! checks the divergence of the velocity field
    !
     !@cuf use cudafor
    implicit none
    integer, intent(in), dimension(3) :: n
    real(8), intent(in), dimension(3) :: dli
    real(8), intent(in), dimension(0:) :: dzfi
    real(8), intent(in), dimension(0:,0:,0:) :: u,v,w
    real(8), intent(out) :: divtot,divmax
    real(8) :: dxi,dyi,div!,dzi,div
#ifdef USE_CUDA
    attributes(managed):: u,v,w,dzfi
    integer:: istat
#endif
    integer :: i,j,k,im,jm,km
    !integer :: ii,jj
    !
    dxi = dli(1)
    dyi = dli(2)
    !dzi = dli(3)
    divtot = 0.d0
    divmax = 0.d0
#ifdef USE_CUDA
    !$cuf kernel do(3) <<<*,*>>>
    do k=1,n(3)
       do j=1,n(2)
          do i=1,n(1)
           km = k-1
           jm = j-1
           im = i-1
#else
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP  SHARED(n,u,v,w,dxi,dyi,dzfi) &
    !$OMP  PRIVATE(i,j,k,im,jm,km,div) &
    !$OMP  REDUCTION(+:divtot) &
    !$OMP  REDUCTION(max:divmax)
    do k=1,n(3)
       km = k-1
       do j=1,n(2)
          jm = j-1
          do i=1,n(1)
             im = i-1
#endif
             div = (w(i,j,k)-w(i,j,km))*dzfi(k) + &
                   (v(i,j,k)-v(i,jm,k))*dyi     + &
                   (u(i,j,k)-u(im,j,k))*dxi
             divmax = max(divmax,abs(div))
             divtot = divtot + div
             !ii = coord(1)*n(1)+i
             !jj = coord(2)*n(2)+j
             !if(abs(div).ge.1.e-12) print*,div,'Large divergence at grid cell: ',ii,jj,k,div
          enddo
       enddo
    enddo
#ifndef USE_CUDA
    !$OMP END PARALLEL DO
#endif
!@cuf istat=cudaDeviceSynchronize()
    call mpi_allreduce(MPI_IN_PLACE,divtot,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,divmax,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
    if(myid.eq.0) print*, 'Total divergence = ', divtot, '| Maximum divergence = ', divmax
    return
  end subroutine chkdiv
end module mod_chkdiv
