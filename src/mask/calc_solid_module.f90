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
    !character, intent(in) :: iMask
    !real(kind=pr), intent(out) :: J
    !real(kind=pr), intent(in)  :: a, b, Mass

    solid%aeroForce(1)    = 0.d0
    solid%aeroForce(2)    = 0.d0

    ! if none J
    J = 1.d0;
    select case (iMask)
      !-free_ellipse------------------------------------------------------------
      case('free_ellipse')
        a = 0.5d0
        b = 0.1d0
        J = Mass * (a**2 + b**2) / 5
      !-free_ellipse END--------------------------------------------------------

      !-free_hat------------------------------------------------------------
      case('free_hat')
        write(*,*) 'Preprocessing the solid: hat...'
        ! Preprocessing is done here. WE ARE IN THE hat BODY SYSTEM NOW.
        ! The hat have two "legs". Every leg have there own coordinate system (CS).
        ! CS(leg_1) and CS(leg_2) collapse in point where the legs meet (definition).
        ! If the hat points up, leg_1 is the left one.
        ! let us shift CS's up so, that the CS of the whole hat is in the
        ! center of gravity (cg)
        cg_shift = leg_l/2.d0 * cos(alpha); ! <- used in create_mask

        !The size (length) of the smoothing is controlled with n_cell_smooth
        smooth_length = n_cell_smooth * dx;

        ! Alpha is the 1/2 opening angle. To build the hat we need to transform this
        ! angle in CS(leg_1) and CS(leg_2).
        alpha_leg_1 = 3.d0*pi/2.d0 - alpha;
        alpha_leg_2 = 3.d0*pi/2.d0 + alpha;

        ! Let us compute the rotation matrix for each leg
        rotate_leg(:,:,1) = reshape((/cos(alpha_leg_1),-sin(alpha_leg_1), &
                                      sin(alpha_leg_1), cos(alpha_leg_1)/),(/2,2/)); ! <- you see NOT the structure of the matrix, LOOK UP "reshape" for more info

        rotate_leg(:,:,2) = reshape((/cos(alpha_leg_2),-sin(alpha_leg_2), &
                                      sin(alpha_leg_2), cos(alpha_leg_2)/),(/2,2/)); ! <- you see NOT the structure of the matrix, LOOK UP "reshape" for more info

        ! The rotation around the cg is made with the "solid%position_angle"
        ! To make the rotation we need first translate the center of the CS(leg_1)
        ! and the CS(leg_2) to the cg of the whole hat. Make the rotation with
        ! "solid%position_angle" and then translate the CS(leg_1) and CS(leg_2)
        ! back. For that operation we need a constant translation distance.
        ! The translation distance for the CS(leg_1) and the CS(_2) is the same,
        ! due to the definition that they collapse in a point, where the legs meet.
        ! The modulus is computed as follow:
        ! in x direction in CS(leg)
        cg_rot_dist(1) = leg_l / 2.d0 * cos(alpha)*cos(alpha);
        ! in y direction in CS(leg)
        cg_rot_dist(2) = leg_l / 2.d0 * cos(alpha)*sin(alpha);

        ! compute the Moment of inertia for iMask
        !J = Mass * leg_l**2 * ( 1 - 0.75d0 * cos(alpha)**2 ) / 3.d0
        J = Mass * ( 1 - 0.75d0 * cos(alpha)**2 ) / 3.d0

        write(*,*) 'Preprocessing the solid: hat...DONE'
      !-free_hat END--------------------------------------------------------

      !-free_triangle------------------------------------------------------------
    case('free_triangle')
        write(*,*) 'Preprocessing the solid: triangle...'

        cg_shift = 2.d0 * leg_l/3.d0 * cos(alpha); ! <- used in create_mask

        !The size (length) of the smoothing is controlled with n_cell_smooth
        smooth_length = n_cell_smooth * dx;

        ! Alpha is the 1/2 opening angle. To build the hat we need to transform this
        ! angle in CS(leg_1) and CS(leg_2).
        alpha_leg_1 = 3.d0*pi/2.d0 - alpha;
        alpha_leg_2 = 3.d0*pi/2.d0 + alpha;

        ! Let us compute the rotation matrix for each leg
        rotate_leg(:,:,1) = reshape((/cos(alpha_leg_1),-sin(alpha_leg_1), &
                                      sin(alpha_leg_1), cos(alpha_leg_1)/),(/2,2/)); ! <- you see NOT the structure of the matrix, LOOK UP "reshape" for more info

        rotate_leg(:,:,2) = reshape((/cos(alpha_leg_2),-sin(alpha_leg_2), &
                                      sin(alpha_leg_2), cos(alpha_leg_2)/),(/2,2/)); ! <- you see NOT the structure of the matrix, LOOK UP "reshape" for more info


        cg_rot_dist(1) = leg_l / 2.d0 * cos(alpha)*cos(alpha);
        ! in y direction in CS(leg)
        cg_rot_dist(2) = leg_l / 2.d0 * cos(alpha)*sin(alpha);

        ! compute the Moment of inertia for iMask
        !J = Mass * leg_l**2 * ( 1 - 0.75d0 * cos(alpha)**2 ) / 3.d0
        J = Mass * ( leg_l**2.d0 + ( leg_l* cos(alpha) )**2.d0 ) / 36.d0

        write(*,*) 'Preprocessing the solid: triangle...DONE'
      !-free_triangle END--------------------------------------------------------

      !-cylinder------------------------------------------------------------
      case('cylinder') ! <- need to be set here for the test unit

      !-cylinder END--------------------------------------------------------

      ! unknown step : error--------------------------------------------------------
      case default
        write (*,*) iMask
        write (*,*) '??? ERROR: this initial condition is not defined in '
        stop
      ! unknown step : error END ---------------------------------------------------

    end select

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
        ! momentum balance
        solid%acceleration(1) = solid%aeroForce(1) / Mass
        solid%acceleration(2) = solid%aeroForce(2) / Mass + g

        solid_tmp%position(1:2) = solid%position(1:2) + dt*solid%velocity(1:2)
        solid_tmp%velocity(1:2) = solid%velocity(1:2) + dt*solid%acceleration(1:2)
        solid_tmp%acceleration(1:2) = solid%acceleration(1:2)
        ! angular momentum balance
        solid%ang_acceleration = solid%momentum / J

        solid_tmp%ang_position = solid%ang_position + dt*solid%ang_velocity
        solid_tmp%ang_velocity = solid%ang_velocity + dt*solid%ang_acceleration
        solid_tmp%ang_acceleration = solid%ang_acceleration
        ! Output: solid_tmp
      ! the euler step END ---------------------------------------------------------

      ! the second RK2 step---------------------------------------------------------
      case (2)
        ! momentum balance
        ! acceleration is constant in time. force the solid_tmp%acceleration
        ! for the next RK step. Forcing not necessary.
        solid%acceleration(1) = solid%aeroForce(1) / Mass
        solid%acceleration(2) = solid%aeroForce(2) / Mass + g
        !write(*,*) solid_tmp%acceleration(2)
        solid%position(1:2) = solid%position(1:2) + 0.5d0*dt*(solid%velocity(1:2)     + solid_tmp%velocity(1:2)     )
        solid%velocity(1:2) = solid%velocity(1:2) + 0.5d0*dt*(solid%acceleration(1:2) + solid_tmp%acceleration(1:2) )

        ! angular momentum balance
        solid%ang_acceleration = solid%momentum / J

        solid%ang_position = solid%ang_position + 0.5d0*dt*(solid%ang_velocity     + solid_tmp%ang_velocity     )
        solid%ang_velocity = solid%ang_velocity + 0.5d0*dt*(solid%ang_acceleration + solid_tmp%ang_acceleration )

        !write (*,*) 'This is ang velocity'
        !write (*,*) solid%ang_velocity
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
  subroutine calc_forces (solid,   mask, u, us)
    use vars
    implicit none
    type(solid_data_struct), intent(inout) :: solid
    real(kind=pr),dimension(0:nx-1,0:ny-1),intent(in) :: mask
    real(kind=pr),dimension(0:nx-1,0:ny-1,1:2),intent(in) :: u, us
    real(kind=pr) :: Fx, Fy, cross_p, u_diff_x, u_diff_y, x, y
    integer :: ix,iy

    cross_p = 0.d0
    Fx = 0.d0
    Fy = 0.d0
    !speed up with Bounding box ?

    !!!$omp parallel do private(ix,iy,x,y,u_diff_x,u_diff_y) &
    !!!$omp& reduction(+:Fx, Fy, cross_p)
    do ix=0,nx-1
      do iy=0,ny-1
        x = dble(ix)*dx-x0
        y = dble(iy)*dy-y0

        !call periodize_solid_coordinate (x,y)

        u_diff_x = u(ix,iy,1)-us(ix,iy,1)
        u_diff_y = u(ix,iy,2)-us(ix,iy,2)

        Fx = Fx + mask(ix,iy)*( u_diff_x )
        Fy = Fy + mask(ix,iy)*( u_diff_y )

        !Cross product: r x u
        cross_p = cross_p + mask(ix,iy) * (x * u_diff_y - y * u_diff_x)
      enddo
    enddo
    !!!$omp end parallel do

    solid%aeroForce(1) = Fx * dx * dy
    solid%aeroForce(2) = Fy * dx * dy

    solid%momentum = cross_p * dx * dy
  end subroutine calc_forces

!===============================================================================
! Keep the solid in the domain
!===============================================================================
  ! This routine checks if the center of gravity is in the fluid domain and puts
  ! it in if not so. This is the periodic condition for the solid
  subroutine periodize_solid_cog (solid)
    use vars
    implicit none
    type(solid_data_struct), intent(inout) :: solid

    if ( y0 < 0.d0 ) then;  solid%position(2) = y0 + yl;  endif
    if ( y0 > yl   ) then;  solid%position(2) = y0 - yl;  endif
    if ( x0 < 0.d0 ) then;  solid%position(1) = x0 + xl;  endif
    if ( x0 > xl   ) then;  solid%position(1) = x0 - xl;  endif

  end subroutine periodize_solid_cog

!===============================================================================
!Keep the solid continuously (intact) at the domain boundary
!===============================================================================
  subroutine periodize_solid_coordinate (x,y)
    use vars
    implicit none
    real(kind=pr),intent(inout) :: x, y

    if ( x > xl/2.d0  ) then;  x = x - xl; endif
    if ( y > yl/2.d0  ) then;  y = y - yl; endif
    if ( x < -xl/2.d0 ) then;  x = x + xl; endif
    if ( y < -yl/2.d0 ) then;  y = y + yl; endif

  end subroutine periodize_solid_coordinate

!===============================================================================
! Rotate the legs of the hat around the center of gravity  (by rotating the CS)
!===============================================================================
  subroutine rotate_hat_leg_cog (CS,theta,leg_number)
    use vars
    implicit none
    real(kind=pr),intent(in) :: theta
    integer, intent(in) :: leg_number
    real(kind=pr),dimension(1:2,1),intent(inout) :: CS   ! <- coordinate system
    real(kind=pr),dimension(1:2,1)               :: CS_r ! temorary variable of the rotation
    select case (leg_number)
      !-rotating the leg_1 around the cog---------------------------------------
      case(1)
        ! translate to the cog in CS(leg_1)
        CS(1,1) = CS(1,1) - cg_rot_dist(1);
        CS(2,1) = CS(2,1) - cg_rot_dist(2);

        ! rotate around the cog
        CS_r(1,1) =  cos(theta)*CS(1,1) + sin(theta)*CS(2,1);
        CS_r(2,1) = -sin(theta)*CS(1,1) + cos(theta)*CS(2,1);

        CS(1,1) = CS_r(1,1) + cg_rot_dist(1);
        CS(2,1) = CS_r(2,1) + cg_rot_dist(2);
        !Output : CS

      !-rotating the leg_1 around the cog END-----------------------------------

      !-rotating the leg_2 around the cog---------------------------------------
      case(2)
        ! translate to the cog in CS(leg_2)
        CS(1,1) = CS(1,1) - cg_rot_dist(1);
        CS(2,1) = CS(2,1) + cg_rot_dist(2);

        ! rotate around the cog
        CS_r(1,1) =  cos(theta)*CS(1,1) + sin(theta)*CS(2,1);
        CS_r(2,1) = -sin(theta)*CS(1,1) + cos(theta)*CS(2,1);

        CS(1,1) = CS_r(1,1) + cg_rot_dist(1);
        CS(2,1) = CS_r(2,1) - cg_rot_dist(2);
        !Output : CS
      !-rotating the leg_2 around the cog END-----------------------------------

      ! unknown leg_number : error----------------------------------------------
      case default
        write (*,*) leg_number
        write (*,*) '??? ERROR: no such leg_number in "rotate_hat_leg_cog" '
        stop
      ! unknown leg_number : error END -----------------------------------------

    end select

  end subroutine rotate_hat_leg_cog

!===============================================================================
! Rotate gurney of the hat around the center of gravity  (by rotating the CS)
!===============================================================================
  subroutine rotate_hat_gurney (CS,theta,leg_number)
    use vars
    implicit none
    real(kind=pr),intent(in) :: theta
    integer, intent(in) :: leg_number
    real(kind=pr),dimension(1:2,1),intent(inout) :: CS   ! <- coordinate system
    real(kind=pr),dimension(1:2,1)               :: CS_r ! temorary variable of the rotation
    select case (leg_number)
      !-rotating the leg_1 around the cog---------------------------------------
      case(1)
        ! translate to the cog in CS(leg_1)
        CS(1,1) = CS(1,1) - ( - leg_l * 0.5d0 * cos(alpha) * sin(alpha) )
        CS(2,1) = CS(2,1) - ( - leg_l + 0.5d0 * leg_l * cos(alpha) * cos(alpha) )

        ! rotate around the cog
        CS_r(1,1) =  cos(theta)*CS(1,1) + sin(theta)*CS(2,1);
        CS_r(2,1) = -sin(theta)*CS(1,1) + cos(theta)*CS(2,1);

        CS(1,1) = CS_r(1,1) + ( - leg_l * 0.5d0 * cos(alpha) * sin(alpha) )
        CS(2,1) = CS_r(2,1) + ( - leg_l + 0.5d0 * leg_l * cos(alpha) * cos(alpha))
        !Output : CS

      !-rotating the leg_1 around the cog END-----------------------------------

      !-rotating the leg_2 around the cog---------------------------------------
      case(2)
        ! translate to the cog in CS(leg_2)
        CS(1,1) = CS(1,1) - ( - leg_l * 0.5d0 * cos(alpha) * sin(alpha) );
        CS(2,1) = CS(2,1) - ( leg_l - 0.5d0 * leg_l * cos(alpha) * cos(alpha) );

        ! rotate around the cog
        CS_r(1,1) =  cos(theta)*CS(1,1) + sin(theta)*CS(2,1);
        CS_r(2,1) = -sin(theta)*CS(1,1) + cos(theta)*CS(2,1);

        CS(1,1) = CS_r(1,1) + ( - leg_l * 0.5d0 * cos(alpha) * sin(alpha) )
        CS(2,1) = CS_r(2,1) + ( leg_l - 0.5d0 * leg_l * cos(alpha) * cos(alpha) )
        !Output : CS
      !-rotating the leg_2 around the cog END-----------------------------------

      ! unknown leg_number : error----------------------------------------------
      case default
        write (*,*) leg_number
        write (*,*) '??? ERROR: no such leg_number in "rotate_hat_leg_cog" '
        stop
      ! unknown leg_number : error END -----------------------------------------

    end select

  end subroutine rotate_hat_gurney

end module calc_solid_module
