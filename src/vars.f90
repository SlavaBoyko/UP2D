module vars
  implicit none
  integer, parameter :: pr = kind (0.d0)
  integer, parameter :: strlen=80

  integer, save :: nx, ny
  integer, save :: nt
  integer, save :: iDealias

  real(kind=pr), save :: xl,yl,dx,dy,x0,y0
  real(kind=pr), save :: Tmax, CFL, tsave, tdrag, dt_fixed=0.d0, dt_max=0.d0, tsave_first=0.d0
  real(kind=pr), save :: nu, eps, pi, tstart, g, Mass
  real(kind=pr), save :: ux_mean, uy_mean
  integer, save :: itsave
  character(len=strlen),save :: intelligent_dt = "yes"
  character(len=strlen),save :: inicond, iMask, iMeanFlow, iMethod

  integer, save :: iSaveVelocity, iSaveVorticity, iSaveMask, iSavePressure

  ! deliberately reduce code to second order FD?
  logical, save :: FD_2nd = .false.

  ! sponge term
  character (len=strlen), save :: iSpongeType
  real(kind=pr), save :: eps_sponge
  integer, save :: use_sponge = 0

  ! memory
  real(kind=pr), dimension(:,:), allocatable, save :: dealiase

  ! solid ======================================================================
  ! rigid solid with one center of gravity (cg). Parameters of cg

  ! define the structure for the solid parameters
  type solid_data_struct
    ! the first index ("1") is the "x" and the ("2") "y"
    real(kind=pr), dimension(1:2) :: position,         &
                                     velocity,         &
                                     acceleration,     &
                                     aeroForce

    real(kind=pr)                 :: ang_position,     &
                                     ang_velocity,     &
                                     ang_acceleration, &
                                     momentum

  end type solid_data_struct

  !Moment of inertia
  real(kind=pr),save :: J

  !this params only for the free_ellipse
  real(kind=pr),save :: a, b

  !Bounding box
  integer, save:: ix_start, ix_end, iy_start, iy_end
!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!

!-----------------------------------------------------------------------------
! Condition for output conditions.
! return true after tfrequ time units or itfrequ time steps or if we're at the
! and of the simulation
!-----------------------------------------------------------------------------
logical function time_for_output( time, dt, it, tfrequ, ifrequ, tmax, tfirst )
  implicit none
  real(kind=pr), intent(in) :: time, dt, tfrequ, tfirst
  real(kind=pr), intent(in) :: tmax ! final time (if we save at the end of run)
  integer, intent(in) :: it ! time step counter
  integer, intent(in) :: ifrequ ! save every ifrequ time steps

  real(kind=pr) :: tnext1, tnext2

  time_for_output = .false.

  ! we never save before tfirst
  if (time<tfirst) return

  if (intelligent_dt=="yes") then
    ! with intelligent time stepping activated, the time step is adjusted not
    ! to pass by tsave,tintegral,tmax,tslice
    ! this is the next instant we want to save
    tnext1 = dble(ceiling(time/tfrequ))*tfrequ
    tnext2 = dble(floor  (time/tfrequ))*tfrequ
    ! please note that the time actually is very close to the next instant we
    ! want to save. however, it may be slightly less or larger. therefore, we
    ! cannot just check (time-tnext), since tnext may be wrong
    if ((abs(time-tnext1)<=1.0d-6).or.(abs(time-tnext2)<=1.0d-6).or.&
        (modulo(it,ifrequ)==0).or.(abs(time-tmax)<=1.0d-6)) then
      time_for_output = .true.
    endif
  else
    ! without intelligent time stepping, we save output when we're close enough
    if ( (modulo(time,tfrequ)<dt).or.(modulo(it,ifrequ)==0).or.(time==tmax) ) then
      time_for_output = .true.
    endif
  endif
end function

end module vars
