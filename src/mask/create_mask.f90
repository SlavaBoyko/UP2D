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
    case('free_hat')
      call free_hat(time, mask, us, u, solid)
    case('free_triangle')
      call free_triangle(time, mask, us, u, solid)
    case('moving_cylinder')
      call moving_cylinder(time, mask, us)
    case('none')
      mask = 0.d0
    case default
      write (*,*) "mask not defnd: ", iMask
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
  integer :: ix,iy, lb, rb, bb, tb
  real(kind=pr)::R,R0,x,y
  real(kind=pr)::x_tmp,y_tmp
  real(kind=pr) :: Fx, Fy, cross_p, u_diff_x, u_diff_y ! used for the Forces


  mask = 0.d0
  cross_p = 0.d0
  Fx = 0.d0
  Fy = 0.d0

  call periodize_solid_cog (solid)

  x0 = solid%position(1)
  y0 = solid%position(2)

  !Bounding container speed up
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
  do ix=ix_start, ix_end
    do iy=iy_start, iy_end
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

subroutine free_hat(time, mask, us, u, solid) ! <- do we need to pass the mask into the routine ?
  use vars
  use calc_solid_module
  use bound_container_module
  implicit none
  type(solid_data_struct), intent(inout) :: solid
  real(kind=pr),intent(in) :: time
  real(kind=pr),dimension(0:nx-1,0:ny-1),intent(inout) :: mask
  real(kind=pr),dimension(0:nx-1,0:ny-1,1:2),intent(inout) :: us, u
  integer :: ix,iy, lb, rb, bb, tb
  real(kind=pr)::R,R0,x,y
  real(kind=pr)::x_tmp,y_tmp,theta ! used for rotation . theta = position angle in this routine
  real(kind=pr) :: Fx, Fy, cross_p, u_diff_x, u_diff_y ! used for the Forces
  real(kind=pr),dimension(1:2,1) :: CS_leg_1,    & ! coordinate system of the leg_1. If the hat points up, leg_1 is the left one.
                                    CS_leg_2,    & ! leg_2
                                    CS_hat         ! the global coordinate system

  us = 0.d0   ! <- imprtant ?
  mask = 0.d0
  cross_p = 0.d0
  Fx = 0.d0
  Fy = 0.d0

  call periodize_solid_cog (solid)

  x0 = solid%position(1)
  y0 = solid%position(2)

  call get_bounding_container_index (rb,lb,bb,tb,solid)
  ! rb = right bounding; lb = left bounding ; bb = bottom bounding; tb = top bounding

  !$omp parallel do private(ix,iy, R, u_diff_x, u_diff_y, CS_leg_1,CS_leg_2,CS_hat) &
  !$omp& reduction(+:Fx, Fy, cross_p)
   do ix=lb, rb
     do iy=bb ,tb
       ! |-------CS_hat--------------|----shift to the legs CS--|
        CS_hat(1,1) = dble(ix)*dx - x0             ! x coordinate
        CS_hat(2,1) = dble(iy)*dy - y0 - cg_shift  ! y coordinate

        call periodize_solid_coordinate (CS_hat(1,1),CS_hat(2,1))

        CS_leg_1 = matmul( rotate_leg(:,:,1), CS_hat) ! <- rotate the CS of the leg_1
        CS_leg_2 = matmul( rotate_leg(:,:,2), CS_hat) ! <- rotate the CS of the leg_2

        call rotate_hat_leg_cog ( CS_leg_1, solid%ang_position, 1 ) ! <- index for the leg
        call rotate_hat_leg_cog ( CS_leg_2, solid%ang_position, 2 )

        call build_smooth_hat (CS_leg_1, CS_leg_2, mask, ix, iy) !Output: mask


        us(ix,iy,1) = solid%velocity(1) - solid%ang_velocity * ( CS_hat(2,1) + cg_shift )
        us(ix,iy,2) = solid%velocity(2) + solid%ang_velocity *  CS_hat(1,1)

        !u_diff_x = ( u(ix,iy,1)- us(ix,iy,1) ) * mask(ix,iy)
        !u_diff_y = ( u(ix,iy,2)- us(ix,iy,2) ) * mask(ix,iy)

        !Fx = Fx + u_diff_x
        !Fy = Fy + u_diff_y


        !cross_p = cross_p + (  CS_hat(1,1) * u_diff_y   -  (CS_hat(2,1)+cg_shift) * u_diff_x  )
    enddo
  enddo
  !$omp end parallel do

        !write(*,*) cross_p
  !solid%aeroForce(1) = Fx * dx * dy /eps
  !solid%aeroForce(2) = Fy * dx * dy /eps
  !write(*,*) solid%aeroForce(2)
  !solid%momentum = cross_p * dx * dy /eps

end subroutine free_hat

!===============================================================================

subroutine free_triangle(time, mask, us, u, solid) ! <- do we need to pass the mask into the routine ?
  use vars
  use calc_solid_module
  use bound_container_module
  implicit none
  type(solid_data_struct), intent(inout) :: solid
  real(kind=pr),intent(in) :: time
  real(kind=pr),dimension(0:nx-1,0:ny-1),intent(inout) :: mask
  real(kind=pr),dimension(0:nx-1,0:ny-1,1:2),intent(inout) :: us, u
  integer :: ix,iy, lb, rb, bb, tb
  real(kind=pr)::R,R0,x,y
  real(kind=pr)::x_tmp,y_tmp,theta ! used for rotation . theta = position angle in this routine
  real(kind=pr) :: Fx, Fy, cross_p, u_diff_x, u_diff_y ! used for the Forces
  real(kind=pr),dimension(1:2,1) :: CS_leg_1,    & ! coordinate system of the leg_1. If the hat points up, leg_1 is the left one.
                                    CS_leg_2,    & ! leg_2
                                    CS_hat,      & ! the global coordinate system
                                    CS_r
  us = 0.d0   ! <- imprtant ?
  mask = 0.d0
  cross_p = 0.d0
  Fx = 0.d0
  Fy = 0.d0

  call periodize_solid_cog (solid)

  x0 = solid%position(1)
  y0 = solid%position(2)

  call get_bounding_container_index (rb,lb,bb,tb,solid)
  ! rb = right bounding; lb = left bounding ; bb = bottom bounding; tb = top bounding

  !$omp parallel do private(ix,iy, R, u_diff_x, u_diff_y, CS_leg_1,CS_leg_2,CS_hat,CS_r) &
  !$omp& reduction(+:Fx, Fy, cross_p)
   do ix=lb, rb
     do iy=bb ,tb
       ! |-------CS_hat--------------|----shift to the legs CS--|
        CS_hat(1,1) = dble(ix)*dx - x0             ! x coordinate
        CS_hat(2,1) = dble(iy)*dy - y0             ! y coordinate

        call periodize_solid_coordinate (CS_hat(1,1),CS_hat(2,1))

        ! rotate around the cog
        CS_r(1,1) =  cos(solid%ang_position)*CS_hat(1,1) + sin(solid%ang_position)*CS_hat(2,1);
        CS_r(2,1) = -sin(solid%ang_position)*CS_hat(1,1) + cos(solid%ang_position)*CS_hat(2,1);

        if (   1.d0/3.d0 * leg_l * cos(alpha) + CS_r(2,1)                       >= 0.d0 .and. &
             + CS_r(2,1)  - 1.d0 / tan(alpha)*CS_r(1,1) - 2.d0/3.d0 * leg_l * cos(alpha) <= 0.d0 .and. &
             + CS_r(2,1)  + 1.d0 / tan(alpha)*CS_r(1,1) - 2.d0/3.d0 * leg_l * cos(alpha) <= 0.d0      &
             ) then
             mask(ix,iy) = 1.d0
        endif

        us(ix,iy,1) = ( solid%velocity(1) - solid%ang_velocity *  CS_r(2,1) ) * mask(ix,iy)
        us(ix,iy,2) = ( solid%velocity(2) + solid%ang_velocity *  CS_r(1,1) )* mask(ix,iy)

        us(ix,iy,1) = ( solid%velocity(1) - solid%ang_velocity * (CS_hat(2,1)+cg_shift)  ) * mask(ix,iy)
        us(ix,iy,2) = ( solid%velocity(2) + solid%ang_velocity *  CS_hat(1,1)            ) * mask(ix,iy)

        u_diff_x = ( u(ix,iy,1)- us(ix,iy,1) ) * mask(ix,iy)
        u_diff_y = ( u(ix,iy,2)- us(ix,iy,2) ) * mask(ix,iy)

        Fx = Fx + u_diff_x
        Fy = Fy + u_diff_y

        cross_p = cross_p + (  CS_r(1,1) * u_diff_y   -  (CS_r(2,1) ) * u_diff_x  )
        !cross_p = cross_p + (  CS_hat(1,1) * u_diff_y   -  (CS_hat(2,1)+cg_shift) * u_diff_x  )

    enddo
  enddo
  !$omp end parallel do
  solid%aeroForce(1) = Fx * dx * dy /eps
  solid%aeroForce(2) = Fy * dx * dy /eps

  solid%momentum = cross_p * dx * dy /eps

end subroutine free_triangle


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

!===============================================================================

subroutine build_smooth_hat (CS_1, CS_2, mask, ix, iy)
  use vars
  implicit none
  real(kind=pr)::R, R_ll, R_rl
  integer,intent(in) :: ix,iy
  real(kind=pr),dimension(0:nx-1,0:ny-1),intent(inout) :: mask
  real(kind=pr),dimension(1:2,1) :: CS_1,    & ! coordinate system of the leg_1. If the hat points up, leg_1 is the left one.
                                    CS_2       ! leg_2

  R = sqrt( CS_1(1,1)**2 + CS_1(2,1)**2 )
  if (R <= leg_h / 2.d0) then
    mask(ix,iy) = 1.d0
  endif

  if (R >= leg_h / 2.d0 .and. R <= leg_h / 2.d0 + smooth_length) then
    mask(ix,iy) = 0.5*(1 + cos((sqrt(CS_1(1,1)**2 +CS_1(2,1)**2) - leg_h/2.d0)*pi/smooth_length) );
  endif

  ! round the left leg
  R_ll = sqrt( ( CS_1(1,1)-leg_l )**2 + CS_1(2,1)**2 )
  if (R_ll <= leg_h / 2.d0) then
    mask(ix,iy) = 1.d0
  endif

  if (R_ll >= leg_h / 2.d0 .and. R_ll <= leg_h / 2.d0 + smooth_length) then
    mask(ix,iy) = 0.5*(1 + cos((sqrt( ( CS_1(1,1) - leg_l)**2 +CS_1(2,1)**2) - leg_h/2.d0)*pi/smooth_length) );
  endif

  ! round the right leg
  R_rl = sqrt( ( CS_2(1,1)-leg_l )**2 + CS_2(2,1)**2 )
  if (R_rl <= leg_h / 2.d0) then
    mask(ix,iy) = 1.d0
  endif

  if (R_rl >= leg_h / 2.d0 .and. R_rl <= leg_h / 2.d0 + smooth_length) then
    mask(ix,iy) = 0.5*(1 + cos((sqrt( ( CS_2(1,1) - leg_l)**2 + CS_2(2,1)**2) - leg_h/2.d0)*pi/smooth_length) );
  endif

  !smooth inside left
  if ( CS_1(1,1) <=  leg_l                        .and.  &
       CS_1(1,1) >=  0.d0                         .and.  & ! in x
       CS_1(2,1) <=  leg_h / 2.d0 + smooth_length .and.  & ! in y
       CS_1(2,1) >=  leg_h / 2.d0                        &
      ) then
        mask(ix,iy) = 0.5d0 * cos((CS_1(2,1) - leg_h/2.d0)*pi / smooth_length) + 0.5;
  endif

  !smooth inside right
  if ( CS_2(1,1) <=  leg_l                        .and.  &
       CS_2(1,1) >=  leg_h / 2.d0                 .and.  & ! in x
       CS_2(2,1) >= -leg_h / 2.d0 - smooth_length .and.  & ! in y
       CS_2(2,1) <= -leg_h / 2.d0                        &
      ) then
        mask(ix,iy) = 0.5d0 * cos((CS_2(2,1) + leg_h/2.d0)*pi / smooth_length) + 0.5;
  endif

  if ( CS_1(1,1) <=  leg_l        .and.  &
       CS_1(1,1) >=  0.d0         .and.  &
       CS_1(2,1) <=  leg_h / 2.d0 .and.  &
       CS_1(2,1) >= -leg_h / 2.d0        &
      ) then
        mask(ix,iy) = 1.d0
  endif

  if ( CS_2(1,1) <=  leg_l        .and.  &
       CS_2(1,1) >=  0.d0         .and.  &
       CS_2(2,1) <=  leg_h / 2.d0 .and.  &
       CS_2(2,1) >= -leg_h / 2.d0        &
      ) then
        mask(ix,iy) = 1.d0
  endif

  !smooth outside left
  if ( CS_1(1,1) <=  leg_l                        .and.  &
       CS_1(1,1) >=  0                            .and.  & ! in x
       CS_1(2,1) >= -leg_h / 2.d0 - smooth_length .and.  & ! in y
       CS_1(2,1) <= -leg_h / 2.d0                        &
      ) then
        mask(ix,iy) = 0.5d0 * cos((CS_1(2,1) + leg_h/2.d0)*pi / smooth_length) + 0.5;
  endif

  !smooth outside right
  if ( CS_2(1,1) <=  leg_l                        .and.  &
       CS_2(1,1) >=  0                            .and.  & ! in x
       CS_2(2,1) <=  leg_h / 2.d0 + smooth_length .and.  & ! in y
       CS_2(2,1) >=  leg_h / 2.d0                        &
      ) then
        mask(ix,iy) = 0.5d0 * cos((CS_2(2,1) - leg_h/2.d0)*pi / smooth_length) + 0.5;
  endif

end subroutine build_smooth_hat
