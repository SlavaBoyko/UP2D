module statistic_module
  implicit none
  contains

!==============================================================================
! write out ekin
!==============================================================================
  subroutine write_out_ekin (time, u)
    use vars
    implicit none
    real(kind=pr),intent(in) :: time
    real(kind=pr),dimension(0:nx-1,0:ny-1,1:2), intent(inout) :: u

    open  (14, file = 'ekin.t', status = 'unknown', access = 'append')
      write (14,'(2(es15.8,1x))') time, sum( (u(:,:,1)**2 + u(:,:,2)**2)/2.d0 )*dx*dy
    close (14)

  end subroutine write_out_ekin

!==============================================================================
! write out the forces and momentum
!==============================================================================
  subroutine write_out_solid_forces (time, solid)
    use vars
    implicit none
    real(kind=pr), intent(in) :: time
    type(solid_data_struct), intent(in) :: solid
    logical :: there

    !open (14, file = 'hat_data.txt', status = 'new' )
    !  write(14,*) 'hi'
    !close (14)
    inquire( file = 'hat_data.txt' , exist= there )

    if (there) then
      open (14, file = 'hat_data.txt', status = 'unknown', access = 'append')
        write (14,'(6(es15.8,1x))') time, solid%aeroForce(1), &
                                          solid%aeroForce(2), &
                                          solid%momentum    , &
                                          solid%position(1) , &
                                          solid%position(2)
      close (14)
    else
      open (14, file = 'hat_data.txt', status = 'new', access = 'append')
        write(14,'(6(1x,a))') '# Time', 'Fx', 'Fy', 'Momentum', 'Postition_x', 'Position_y'
      close(14)
    endif
    ! open (14, file = 'momentum.txt', status = 'unknown', access = 'append')
    !   write (14,'(2(es15.8,1x))') time, solid%momentum
    ! close (14)
    !
    ! open (14, file = 'y_position.txt', status = 'unknown', access = 'append')
    !   write (14,'(2(es15.8,1x))') time, solid%position(2)
    ! close (14)
    !
    ! open (14, file = 'x_position.txt', status = 'unknown', access = 'append')
    !   write (14,'(2(es15.8,1x))') time, solid%position(1)
    ! close (14)

  end subroutine write_out_solid_forces

!==============================================================================
! Output how much time remains in the simulation.
!==============================================================================
  subroutine are_we_there_yet(time, wtime_tstart, dt1, it)
    use vars
    use timing_module
    implicit none

    real(kind=pr),intent(inout) :: time,wtime_tstart,dt1
    integer,intent(inout) :: it
    real(kind=pr):: time_left, t2, time_nt

    ! elapsed time since time stepping started
    t2 = MPI_wtime() - wtime_tstart
    ! estimate remaining time until we reach tmax
    time_left = (tmax-time) * (t2/(time-tstart))
    ! estimate remaining time until we real nt time steps
    time_nt = 9e9*dble(nt-it) * (t2/dble(it))
    ! remaining time is minimum of both
    time_left = min( time_left, time_nt )

    write(*,'("time left: ",i2,"d ",i2,"h ",i2,"m ",i2,"s dt=",es10.2,"s t=",g10.2)') &
    floor(time_left/(24.d0*3600.d0))   ,&
    floor(mod(time_left,24.d0*3600.d0)/3600.d0),&
    floor(mod(time_left,3600.d0)/60.d0),&
    floor(mod(mod(time_left,3600.d0),60.d0)),&
    dt1,time
  end subroutine are_we_there_yet

end module statistic_module
