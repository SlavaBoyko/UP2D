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

    ! compute the Moment of inertia for iMask
    select case (iMask)
      !-free_ellipse------------------------------------------------------------
      case('free_ellipse')
        a = 0.5d0
        b = 0.1d0
        J = Mass * (a**2 + b**2) / 5
      !-free_ellipse END--------------------------------------------------------

      !-free_hut------------------------------------------------------------
      case('free_hut')

      !-free_hut END--------------------------------------------------------

      !-free_hut------------------------------------------------------------
      case('cylinder') ! <- need to be set here for the test unit

      !-free_hut END--------------------------------------------------------

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

        ! angular momentum balance
        solid%ang_acceleration = solid%momentum / J

        solid_tmp%ang_position = solid%ang_position + dt*solid%ang_velocity
        solid_tmp%ang_velocity = solid%ang_velocity + dt*solid%ang_acceleration
        ! Output: solid_tmp
      ! the euler step END ---------------------------------------------------------

      ! the second RK2 step---------------------------------------------------------
      case (2)
        ! momentum balance
        ! acceleration is constant in time. force the solid_tmp%acceleration
        ! for the next RK step. Forcing not necessary.
        solid_tmp%acceleration(1) = solid%aeroForce(1) / Mass
        solid_tmp%acceleration(2) = solid%aeroForce(2) / Mass + g

        solid%position(1:2) = solid%position(1:2) + 0.5d0*dt*(solid%velocity(1:2)     + solid_tmp%velocity(1:2)     )
        solid%velocity(1:2) = solid%velocity(1:2) + 0.5d0*dt*(solid%acceleration(1:2) + solid_tmp%acceleration(1:2) )

        ! angular momentum balance
        solid_tmp%ang_acceleration = solid%momentum / J

        solid%ang_position = solid%ang_position + 0.5d0*dt*(solid%ang_velocity     + solid_tmp%ang_velocity     )
        solid%ang_velocity = solid%ang_velocity + 0.5d0*dt*(solid%ang_acceleration + solid_tmp%ang_acceleration )

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

        call periodize_solid_coordinate (x,y)

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

end module calc_solid_module
