!===============================================================================
subroutine create_mask (time, mask, us, u, solid)
  use vars
  implicit none
  type(solid_data_struct), intent(inout) :: solid
  real(kind=pr), intent (in) :: time
  real(kind=pr),dimension(0:nx-1,0:ny-1),intent(inout) :: mask
  real(kind=pr),dimension(0:nx-1,0:ny-1,1:2),intent(inout) :: us, u
  real(kind=pr) :: R
  integer :: ix, iy

  mask = 0.d0
  us   = 0.d0

  select case (iMask)
    case('cylinder')
      call cylinder(mask, us)
    case('ellipse')
      call ellipse(mask, us)
    case('free_ellipse')
      call free_ellipse(time, mask, us, u, solid)
    case('moving_cylinder')
      call moving_cylinder(time, mask, us)
    case('none')
      mask = 0.d0
    case default
      write (*,*) "mask not defnd", iMask
      stop
  end select

  !$omp parallel do private(iy)
  do iy=0,ny-1
    mask(:,iy) = mask(:,iy) / eps
  enddo
  !$omp end parallel do

end subroutine create_mask

!===============================================================================

subroutine cylinder(mask, us)
  use vars
  implicit none
  real(kind=pr),dimension(0:nx-1,0:ny-1),intent(inout) :: mask
  real(kind=pr),dimension(0:nx-1,0:ny-1,1:2),intent(inout) :: us
  integer :: ix,iy
  real(kind=pr)::R,R0,smooth

  x0 = xl/2.d0
  y0 = yl/2.d0
  R0 = 1.0d0
  smooth = 2.d0*max(dx,dy)

  !$omp parallel do private(ix,iy,R)
  do ix=0,nx-1
    do iy=0,ny-1
      R = dsqrt( (dble(ix)*dx-x0)**2 +(dble(iy)*dy-y0)**2 )
      call SmoothStep(mask(ix,iy), R, R0, smooth)
    enddo
  enddo
  !$omp end parallel do
end subroutine cylinder

!===============================================================================

subroutine ellipse(mask, us)
  use vars
  implicit none
  real(kind=pr),dimension(0:nx-1,0:ny-1),intent(inout) :: mask
  real(kind=pr),dimension(0:nx-1,0:ny-1,1:2),intent(inout) :: us
  integer :: ix,iy
  real(kind=pr)::R,R0,x,y

  x0 = xl/2.d0
  y0 = yl/2.d0

  !$omp parallel do private(ix,iy,R,x,y)
  do ix=0,nx-1
    do iy=0,ny-1
      x = dble(ix)*dx-x0
      y = dble(iy)*dy-y0

      R = (x/0.5d0)**2  +  (y/0.1d0)**2
      if (R<= 1.d0) then
        mask(ix,iy) = 1.d0
        us(ix,iy,2) = -1.d0
      endif

    enddo
  enddo
  !$omp end parallel do
end subroutine ellipse

!===============================================================================
! NO SMOOTHING FUNKTION !
subroutine free_ellipse(time, mask, us, u, solid)
  use vars
  use calc_solid_module
  implicit none
  type(solid_data_struct), intent(inout) :: solid
  real(kind=pr),intent(in) :: time
  real(kind=pr),dimension(0:nx-1,0:ny-1),intent(inout) :: mask
  real(kind=pr),dimension(0:nx-1,0:ny-1,1:2),intent(inout) :: us, u
  integer :: ix,iy
  real(kind=pr)::R,R0,x,y
  real(kind=pr)::x_tmp,y_tmp ! used for rotation
  real(kind=pr) :: Fx, Fy, cross_p, u_diff_x, u_diff_y ! used for the Forces

  mask = 0.d0
  cross_p = 0.d0
  Fx = 0.d0
  Fy = 0.d0

  call periodize_solid_cog (solid)

  x0 = solid%position(1)
  y0 = solid%position(2)

  ! Bounding container speed up
  ix_start = ceiling( (x0 - (a + 0.1d0) )/dx )
  ix_end   = ceiling( (x0 + (a + 0.1d0) )/dx )
  iy_start = ceiling( (y0 - (a + 0.1d0) )/dy )
  iy_end   = ceiling( (y0 + (a + 0.1d0) )/dy )

  if (iy_end >= ny-1) then; iy_end = ny-1; iy_start = 0;  endif;
  if (iy_start <= 0 ) then; iy_end = ny-1; iy_start = 0;  endif;
  if (ix_start <= 0  ) then; ix_start = 0; ix_end = nx-1; endif;
  if (ix_end >= nx-1 ) then; ix_start = 0; ix_end = nx-1; endif;

  !$omp parallel do private(ix,iy, R, x,y, x_tmp,y_tmp, u_diff_x, u_diff_y) &
  !$omp& reduction(+:Fx, Fy, cross_p)
!  do ix=0,nx-1
!    do iy=0,ny-1
  do ix=ix_start,ix_end
    do iy=iy_start,iy_end
      x = dble(ix)*dx-x0
      y = dble(iy)*dy-y0

      call periodize_solid_coordinate (x,y)

      x_tmp = cos(solid%ang_position)*x - sin(solid%ang_position)*y
      y_tmp = sin(solid%ang_position)*x + cos(solid%ang_position)*y

      R = (x_tmp/a)**2  +  (y_tmp/b)**2
      if (R<= 1.d0) then
        mask(ix,iy) = 1.d0
        ! update the solid velocity
        ! this ist the solid velocity
        us(ix,iy,1) = solid%velocity(1) - solid%ang_velocity * y
        us(ix,iy,2) = solid%velocity(2) + solid%ang_velocity * x
        ! we can also compute the Forces here (!)
        ! they look like this in seperate routine
        ! This is handy becouse we are here in the mask region
        u_diff_x = u(ix,iy,1)- us(ix,iy,1)
        u_diff_y = u(ix,iy,2)- us(ix,iy,2)

          !Fx = Fx + mask(ix,iy)*( u_diff_x )
          !Fy = Fy + mask(ix,iy)*( u_diff_y )

          !Cross product: r x u
          !cross_p = cross_p + mask(ix,iy) * (x * u_diff_y - y * u_diff_x)

        Fx = Fx + u_diff_x
        Fy = Fy + u_diff_y

        cross_p = cross_p + (x * u_diff_y - y * u_diff_x)
      endif
    enddo
  enddo
  !$omp end parallel do

  solid%aeroForce(1) = Fx * dx * dy /eps
  solid%aeroForce(2) = Fy * dx * dy /eps

  solid%momentum = cross_p * dx * dy /eps

end subroutine free_ellipse

!===============================================================================

subroutine moving_cylinder(time,mask, us)
  use vars
  implicit none
  real(kind=pr),intent(in) :: time
  real(kind=pr),dimension(0:nx-1,0:ny-1),intent(inout) :: mask
  real(kind=pr),dimension(0:nx-1,0:ny-1,1:2),intent(inout) :: us

  integer :: ix,iy
  real(kind=pr) :: R, R0, smooth

  R0 = 0.5d0
  smooth = 2.d0*max(dx,dy)
  x0 = xl/2.d0 + sin(2.d0*pi*time)
  y0 = yl/2.d0

  !$omp parallel do private(ix,iy,R)
  do ix=0,nx-1
    do iy=0,ny-1
      R = dsqrt( (dble(ix)*dx-x0)**2 +(dble(iy)*dy-y0)**2 )
      call SmoothStep(mask(ix,iy), R, R0, smooth)
      if (mask(ix,iy) > 0.d0) then
        us(ix,iy,1) = 2.d0*pi*cos(2.d0*pi*time)
        us(ix,iy,2) = 0.d0
      endif
    enddo
  enddo
  !$omp end parallel do
end subroutine moving_cylinder


!===============================================================================

subroutine SmoothStep (f,x,t,h)
  use vars
  implicit none
  real(kind=pr), intent (out) :: f
  real(kind=pr), intent (in)  :: x,t,h

  if (x<=t-h) then
    f = 1.d0
  elseif (((t-h)<x).and.(x<(t+h))) then
    f = 0.5d0*(1.0d0+cos((x-t+h)*pi/(2.d0*h)) )
  else
    f = 0.0d0
  endif
end subroutine SmoothStep
