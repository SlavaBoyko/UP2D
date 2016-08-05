!-------------------------------------------------------------------------------
! set up the mask for the sponge term
!-------------------------------------------------------------------------------
subroutine sponge_mask(time, mask_sponge, solid)
  use vars
  implicit none
  type(solid_data_struct), intent(in) :: solid
  real(kind=pr), intent (in) :: time
  real(kind=pr),dimension(0:nx-1,0:ny-1), intent(inout) :: mask_sponge
  real(kind=pr) :: R, x0_s, y0_s
  integer :: ix, iy

  if (use_sponge == 0) return

  mask_sponge = 0.d0

  ! read out data for the moving sponge_mask
  x0_s = solid%position(1)
  y0_s = solid%position(2)

  select case (iSpongeType)
    case ('everywhere')
      mask_sponge = 1.d0
    case('none')
      mask_sponge = 0.d0
    case('allEdges')
      call sponge_all_four_edges(mask_sponge, x0_s, y0_s)
    case default
      write (*,*) "mask not defnd", iSpongeType
      stop
  end select

end subroutine sponge_mask

!-------------------------------------------------------------------------------
! Add the sponge term to the non-linear terms in fourier space
!-------------------------------------------------------------------------------
subroutine add_sponge_term(time, nlk, vor, mask_sponge, work1, work2)
  use vars
  implicit none
  real(kind=pr),intent(in) :: time
  real(kind=pr),dimension(0:nx-1,0:ny-1,1:2),intent(inout) :: nlk
  real(kind=pr),dimension(0:nx-1,0:ny-1),intent(inout) :: vor, mask_sponge, work1, work2

  real(kind=pr),dimension(:,:,:), allocatable :: sp_tmp
  integer :: iy

  if (use_sponge == 1) then
    ! allocate temporary array
    allocate( sp_tmp(0:nx-1, 0:ny-1,1:2) )

    ! apply sponge penalization to vorticity
    !$omp parallel do private(iy)
    do iy=0,ny-1
      work1(:,iy) = -mask_sponge(:,iy)*vor(:,iy)/eps_sponge
    enddo
    !$omp end parallel do

    call fft(work1,work2)

    ! obtain the velocity
    call vorticity2velocity( work2, sp_tmp )
    ! add sponge term to NL terms (in F-space)
    nlk = nlk + sp_tmp

    deallocate(sp_tmp)
  endif

end subroutine

!-------------------------------------------------------------------------------
! Add the sponge term to the non-linear terms in fourier space
!-------------------------------------------------------------------------------

subroutine sponge_all_four_edges(mask_sponge, x0_s, y0_s)
  use vars
  use calc_solid_module ! we use one functtion of this here
  implicit none
  real(kind=pr), intent(in) :: x0_s, y0_s
  real(kind=pr),dimension(0:nx-1,0:ny-1), intent(inout) :: mask_sponge
  integer :: ix,iy
  real(kind=pr)::x,y

  !$omp parallel do private(ix,iy,x,y)
  do ix=0,nx-1
    do iy=0,ny-1
      x = dble(ix)*dx - x0_s
      y = dble(iy)*dy - y0_s

       ! here we can use the function from the calc_solid_module
       ! to periodize the sponge_mask in case it is moving
       call periodize_solid_coordinate (x,y)

       if ( x >= xl/2.d0 - Sp_thickness * dx .and. x <= xl/2.d0  .or. &
            x >= -xl/2.d0 .and. x <= - xl/2.d0 + Sp_thickness * dx  .or. &
            y >= yl/2.d0 - Sp_thickness * dy .and. y <= yl/2.d0 .or. &
            y >= -yl/2.d0 .and. y <= - yl/2.d0 + Sp_thickness * dx ) then

         mask_sponge(ix,iy) = 1.d0
       endif

    enddo
  enddo
  !$omp end parallel do

end subroutine sponge_all_four_edges
