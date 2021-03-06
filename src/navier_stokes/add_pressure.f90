!-------------------------------------------------------------------------------
! add pressure gradient to right hand side of navier-stokes
! INPUT:
!   nlk   pre-computed source terms (penalization, non-linear) FOURIER space
! OUTPUT:
!   nlk   now divergence-free right hand side
!-------------------------------------------------------------------------------
subroutine add_pressure (nlk)
  use vars
  implicit none
  real(kind=pr),dimension(0:nx-1,0:ny-1,1:2), intent (inout) :: nlk

  real(kind=pr),dimension(:,:),allocatable :: work1, work2, work3
  integer :: iy

  ! allocate work memory
  allocate ( work1(0:nx-1,0:ny-1), work2(0:nx-1,0:ny-1), work3(0:nx-1,0:ny-1) )

  !---------------------------------------------------------------------------
  ! classic pressure projects all source terms
  !---------------------------------------------------------------------------
  ! divergence of non-linear terms
  call cofdx ( nlk(:,:,1), work1 )
  call cofdy ( nlk(:,:,2), work2 )

  ! solve poisson eqn to get the pressure (in work3)
  call poisson ( work1+work2, work3)

  ! compute gradient of pressure
  call cofdx (work3, work1)
  call cofdy (work3, work2)

  ! add gradient to RHS
  !$omp parallel do private(iy)
  do iy=0,ny-1
    nlk(:,iy,1) = nlk(:,iy,1) + work1(:,iy)
    nlk(:,iy,2) = nlk(:,iy,2) + work2(:,iy)
  enddo
  !$omp end parallel do

  ! release work memory
  deallocate(work1, work2, work3)
end subroutine add_pressure
