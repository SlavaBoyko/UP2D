module bound_container_module
  implicit none
  contains

!===============================================================================
! main call routine
!===============================================================================
  subroutine get_bounding_container_index (rb,lb,bb,tb,solid)
    use vars
    implicit none
    real(kind=pr), dimension(1:3,1:2)   :: bc_nodes
    type(solid_data_struct), intent(in) :: solid
    integer,intent(out) :: rb,  &  ! right  bounding
                           lb,  &  ! left   bounding
                           bb,  &  ! bottom bounding
                           tb      ! top    bounding

    select case (BC)
      !-no BC------------------------------------------------------------
      case('no')
        lb = 0
        rb = nx-1
        bb = 0
        tb = ny-1
      !-no BC END--------------------------------------------------------

      !- BC------------------------------------------------------------
      case('yes')
        call get_3_bounding_points (bc_nodes,solid)
        call get_array_index (rb,lb,bb,tb,bc_nodes,solid)
      !- BC END--------------------------------------------------------

      ! unknown BC case : error--------------------------------------------------------
      case default
        write (*,*) BC
        write (*,*) '??? ERROR: no such BC case. bound_container_modul '
        stop
      ! unknown BC case : error END ---------------------------------------------------

    end select

  end subroutine get_bounding_container_index

!===============================================================================
! find the 3 bounding points of the hat
!===============================================================================
  subroutine get_3_bounding_points (bc_nodes,solid)
    use vars
    implicit none
    type(solid_data_struct),           intent(in)    :: solid
    real(kind=pr), dimension(1:3,1:2), intent(out)   :: bc_nodes
    real(kind=pr), dimension(1:2,1:2)                :: points, points_r
    ! Point one is the top one. we are in the domain coordinae system
    ! we add the center of gravity of the hat later
    ! alpha = 1/2 opening angle

    !                      |------------|  <- buffer,to get out of the solid
    bc_nodes(1,1) = -(leg_l + 2.d0*leg_h) / 2.d0*cos(alpha)*sin(solid%ang_position)
    bc_nodes(1,2) =  (leg_l + 2.d0*leg_h) / 2.d0*cos(alpha)*cos(solid%ang_position)
    !                |------------------------------------| <-  center of gravity

    ! The sidepoints
    points(1:2,1) = leg_l * 1.d0
    points(1:2,2) = 0.d0

    ! rotate (positiv) the basis vektor in respekt to the leg_1 coordinate system
    ! and shift it up by the cg distance
    points_r(1,1) = cos(alpha_leg_1) * points(1,1)                         ! first point, x coordinate
    points_r(1,2) = sin(alpha_leg_1) * points(1,1) + leg_l/2.d0*cos(alpha) ! first point, y coordinate

    ! do the same for the right leg point
    points_r(2,1) = cos(alpha_leg_2) * points(2,1)                         ! second point, x coordinate
    points_r(2,2) = sin(alpha_leg_2) * points(2,1) + leg_l/2.d0*cos(alpha) ! second point, y coordinate

    ! perform the last rotation in respect to the cg
    ! first point
    !                                    |------------| <- x coordinate of the leg_1
    points(1,1) = cos(solid%ang_position)*points_r(1,1) - sin(solid%ang_position)*points_r(1,2);
    points(1,2) = sin(solid%ang_position)*points_r(1,1) + cos(solid%ang_position)*points_r(1,2);
    !                                                                           |---------------| <- y coordinate of the leg_1
    ! second point
    !                                    |------------| <- x coordinate of the leg_2
    points(2,1) = cos(solid%ang_position)*points_r(2,1) - sin(solid%ang_position)*points_r(2,2);
    points(2,2) = sin(solid%ang_position)*points_r(2,1) + cos(solid%ang_position)*points_r(2,2);
    !                                                                           |---------------| <- y coordinate of the leg_2

    ! set the poits to write out
    bc_nodes(2,1) = points(1,1) ! the leg_1 point "x"
    bc_nodes(2,2) = points(1,2) ! the leg_1 point "y"

    bc_nodes(3,1) = points(2,1) ! the leg_2 point "x"
    bc_nodes(3,2) = points(2,2) ! the leg_2 point "y"

  end subroutine get_3_bounding_points

!===============================================================================
! find the array index of the 3 bounding points of the hat
!===============================================================================
  subroutine get_array_index (rb,lb,bb,tb,bc_nodes,solid)
    use vars
    implicit none
    type(solid_data_struct),           intent(in)   :: solid
    real(kind=pr), dimension(1:3,1:2), intent(in)   :: bc_nodes
    real(kind=pr)                                   :: rb_r,lb_r,bb_r,tb_r ! <- index r for "real number"
    integer, intent(out)                            :: rb,lb,bb,tb

    ! we have the bounding box. Now we set it arround the hat
    rb_r = solid%position(1) + abs(maxval(bc_nodes(:,1)))  + buffer + smooth_length;
    lb_r = solid%position(1) - abs(minval(bc_nodes(:,1)))  - buffer - smooth_length;
    tb_r = solid%position(2) + abs(maxval(bc_nodes(:,2)))  + buffer + smooth_length;
    bb_r = solid%position(2) - abs(minval(bc_nodes(:,2)))  - buffer - smooth_length;

    ! convert to integer
    lb = ceiling( lb_r / dx )
    rb = ceiling( rb_r / dx )
    bb = ceiling( bb_r / dy )
    tb = ceiling( tb_r / dy )

    ! if we are at the boundary of the domain turn of the bounding_container
    if (tb >= ny-1 ) then; tb = ny-1; bb = 0   ;  endif;
    if (bb <= 0    ) then; tb = ny-1; bb = 0   ;  endif;
    if (lb <= 0    ) then; lb = 0   ; rb = nx-1;  endif;
    if (rb >= nx-1 ) then; lb = 0   ; rb = nx-1;  endif;

    !Output: rb,lb,bb,tb
  end subroutine get_array_index

end module bound_container_module
