subroutine mean_flow (uk,time)
  use vars
  implicit none
  real(kind=pr), dimension(0:nx-1,0:ny-1,1:2), intent(inout) :: uk
  real(kind=pr), intent(in) :: time

  select case (iMeanFlow)
    case ('none')
      ! don't do anything
    case ('constant')
      uk(0:1,0:1,1) = ux_mean
      uk(0:1,0:1,2) = uy_mean
    case ('oscillating')
      uk(0:1,0:1,1) = ux_mean
      uk(0:1,0:1,2) = uy_mean * pi * sin(2.d0*pi*time)
  end select

end subroutine mean_flow
