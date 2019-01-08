module mod_common_mpi
  use mpi
  use mod_param, only: dims
  implicit none
  integer :: myid
  integer :: left,right,front,back
  integer, dimension(2) :: coord
  integer :: comm_cart,ierr
  integer :: xhalo,yhalo
  integer :: status(MPI_STATUS_SIZE)
  real(8), allocatable, dimension(:,:) :: xsl_buf, xrl_buf, ysl_buf, yrl_buf, xsr_buf, xrr_buf, ysr_buf, yrr_buf
#ifdef USE_CUDA
  integer :: mydev
  attributes( managed ) :: xsl_buf, xrl_buf, ysl_buf, yrl_buf, xsr_buf, xrr_buf, ysr_buf, yrr_buf
#endif
end module mod_common_mpi
