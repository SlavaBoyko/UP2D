module RK2_module
  implicit none
  contains

subroutine RK2 (time, dt,it, u, uk, p, vort, nlk, mask, us, mask_sponge, solid)
  use vars
  use rhs
  use calc_solid_module
  implicit none
  type(solid_data_struct), intent(inout) :: solid
  type(solid_data_struct) :: solid_tmp ! for the RK time step
  integer,intent(in) :: it
  real(kind=pr),intent(out) :: dt
  real(kind=pr),intent(in) :: time
  real(kind=pr),dimension(0:nx-1,0:ny-1), intent(inout) :: mask, mask_sponge
  real(kind=pr),dimension(0:nx-1,0:ny-1), intent (inout) :: vort, p
  real(kind=pr),dimension(0:nx-1,0:ny-1,1:2), intent (inout) :: u, uk, nlk
  real(kind=pr),dimension(0:nx-1,0:ny-1,1:2), intent(inout) :: us
  real(kind=pr),dimension(:,:), allocatable :: workvis
  real(kind=pr),dimension(:,:,:), allocatable :: nlk2, uk_tmp, u_tmp

  integer :: iy
  real(kind=pr) :: adjust_dt

  ! allocate memory for work arrays. note we should not use them from stack, but
  ! allocateable instead, to avoid issues with stack limitations.
  allocate(nlk2(0:nx-1, 0:ny-1,1:2), uk_tmp(0:nx-1, 0:ny-1,1:2), u_tmp(0:nx-1, 0:ny-1,1:2))
  allocate(workvis(0:nx-1, 0:ny-1))

  !-- determine time step
  dt = adjust_dt (time,it,u)

  !-- compute integrating factor
  call cal_vis (dt, workvis)
  !-- mask and us
  call create_mask (time, mask, us, solid)
  !-- RHS and pressure
  call cal_nlk (time, u, uk, vort, nlk, mask, us, mask_sponge)
  call add_pressure (nlk)
  !-- Calculate forces
  call calc_forces (solid,   mask, u, us)
  !-- first RK2 step for solid
  call RK2_rhs_solid(solid, solid_tmp, dt, 1) !<- this one need forces. they are in solid

  !-- do the euler ste
  !$omp parallel do private (iy)
  do iy=0,ny-1
    uk_tmp(:,iy,1) = (uk(:,iy,1)+dt*nlk(:,iy,1))*dealiase(:,iy)*workvis(:,iy)
    uk_tmp(:,iy,2) = (uk(:,iy,2)+dt*nlk(:,iy,2))*dealiase(:,iy)*workvis(:,iy)
  enddo
  !$omp end parallel do

  !-- mean flow forcing
  call mean_flow (uk_tmp,time)

  !-- velocity in phys. space
  call ifft (uk_tmp(:,:,1), u_tmp(:,:,1))
  call ifft (uk_tmp(:,:,2), u_tmp(:,:,2))

  !---------------------------------------------------------------------------------
  ! do second RK2 step (RHS evaluation with the argument defined above)
  !---------------------------------------------------------------------------------
  !-- mask and us
  call create_mask (time+dt, mask, us, solid_tmp)
  !-- RHS and pressure
  call cal_nlk (time+dt, u_tmp, uk_tmp, vort, nlk2, mask, us, mask_sponge)
  call add_pressure (nlk2)
  !-- Calculate forces
  call calc_forces (solid,   mask, u_tmp, us) !<- the u_tmp is different
  !-- second RK2 step for solid
  call RK2_rhs_solid(solid, solid_tmp, dt, 2) !<- this one need forces. they are in solid

  !-- sum up all the terms (final step)
  !$omp parallel do private (iy)
  do iy=0,ny-1
    uk(:,iy,1) = ( uk(:,iy,1)*workvis(:,iy) + 0.5d0*dt*( nlk(:,iy,1)*workvis(:,iy) + nlk2(:,iy,1) ) )*dealiase(:,iy)
    uk(:,iy,2) = ( uk(:,iy,2)*workvis(:,iy) + 0.5d0*dt*( nlk(:,iy,2)*workvis(:,iy) + nlk2(:,iy,2) ) )*dealiase(:,iy)
  enddo
  !$omp end parallel do

  !-- mean flow forcing
  call mean_flow (uk,time + dt)

  !-- velocity in phys. space
  !-- note: we now advanced in fourier space (uk) and transform back to phys. space (u)
  !-- to ensure that bith fields are synchronous. this is useful, for example, for drag
  !-- computation
  call ifft (uk(:,:,1), u(:,:,1))
  call ifft (uk(:,:,2), u(:,:,2))

  deallocate(nlk2, uk_tmp, u_tmp,workvis)

!  at the end of the time step, we consistently return u and uk
end subroutine RK2

end module
