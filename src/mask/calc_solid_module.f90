module calc_solid_module
  implicit none
  contains

!===============================================================================
!solid_initialisation
!===============================================================================
subroutine solid_initialisation (solid)
  use vars
  implicit none
  type(solid_data_struct), intent(inout) :: solid

  !  fill the structure for the solid parameters with ini data
    solid%position(1)     = position_x
    solid%position(2)     = position_y
    solid%velocity(1)     = velocity_x
    solid%velocity(2)     = velocity_y
    solid%acceleration(1) = acceleration_x
    solid%acceleration(2) = acceleration_x
    solid%aeroForce(1)    = 0.d0
    solid%aeroForce(2)    = 0.d0

end subroutine solid_initialisation

!===============================================================================
! The subroutine for the RHS of the solid.
!===============================================================================
subroutine RK2_rhs_solid (solid, solid_tmp, dt, step)
!-------------------------------------------------------------------------------
! This embedded routine performs two steps of RK2
! "step" can be "1" or "2"
! step 1 with "solid" is for the euler step
! step 2 with "solid_tmp for the second RK2 step (RHS evaluation with the
! argument defined above)
!-------------------------------------------------------------------------------
  use vars
  implicit none
  type(solid_data_struct), intent(inout) :: solid
  type(solid_data_struct), intent(inout) :: solid_tmp
  real(kind=pr),intent(in) :: dt
  integer,intent(in) :: step

  select case (step)
    ! the euler step -------------------------------------------------------------
    case (1)
      solid%acceleration(1) = 0.d0
      solid%acceleration(2) = g

      solid_tmp%position(1:2) = solid%position(1:2) + dt*solid%velocity(1:2)
      solid_tmp%velocity(1:2) = solid%velocity(1:2) + dt*solid%acceleration(1:2)

      ! acceleration is constant in time. force the solid_tmp%acceleration
      ! for the next RK step. Forcing not necessary.
      solid_tmp%acceleration(1) = 0.d0
      solid_tmp%acceleration(2) = g

      ! Output: solid_tmp
    ! the euler step END ---------------------------------------------------------

    ! the second RK2 step---------------------------------------------------------
    case (2)
      solid%position(1:2) = solid%position(1:2) + 0.5d0*dt*(solid%velocity(1:2)     + solid_tmp%velocity(1:2)     )
      solid%velocity(1:2) = solid%velocity(1:2) + 0.5d0*dt*(solid%acceleration(1:2) + solid_tmp%acceleration(1:2) )

      ! Output: solid at the time t+1 with RK2
    ! the second RK2 step END-----------------------------------------------------

    ! unknown step : error--------------------------------------------------------
    case default
      write (*,*) step
      write (*,*) '??? ERROR: Invalid initial condition'
      stop
    ! unknown step : error END ---------------------------------------------------
  end select

end subroutine RK2_rhs_solid

!===============================================================================
! This subroutine calculates the forces
!===============================================================================
subroutine calc_forces (mask, u, us, solid)
  use vars
  implicit none
  type(solid_data_struct), intent(inout) :: solid
  real(kind=pr),dimension(0:nx-1,0:ny-1),intent(in) :: mask
  real(kind=pr),dimension(0:nx-1,0:ny-1,1:2),intent(in) :: u, us
  ! not optimal to use the whole field...(?)
  real(kind=pr),dimension(0:nx-1, 0:ny-1,1:2) :: u_diff
  integer :: ix,iy

  u_diff = 0.d0 ! correkt?

  ! have to think about it for the best performens

  !$omp parallel do private(ix,iy)
  do ix=0,nx-1
    do iy=0,ny-1

      if (mask(ix,iy) > 0.d0) then
        ! calc the velocity difference
        u_diff(ix,iy,1) = u(ix,iy,1) - solid%velocity(1)
        u_diff(ix,iy,2) = u(ix,iy,2) - solid%velocity(2)
      endif

    enddo
  enddo
  !$omp end parallel do

  solid%aeroForce(1) = dx * dy * sum(sum(u_diff(:,:,1), DIM = 1), DIM =1 ) / eps
  solid%aeroForce(2) = dx * dy * sum(sum(u_diff(:,:,2), DIM = 1), dim =1 ) / eps

end subroutine calc_forces

end module calc_solid_module
