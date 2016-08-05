! Wrapper to read parameters from an ini file for fsi.  Also reads
! parameters which are set in the vars module.
subroutine get_params(paramsfile,solid)
  use vars
  use ini_files_parser
  type(inifile) :: PARAMS
  character(len=*) :: paramsfile
  type(solid_data_struct) :: solid

  ! Read the paramsfile and put the length i and the text in PARAMS
  call read_ini_file(PARAMS, paramsfile, .true.)

  ! Resolution section
  call read_param(PARAMS,"Resolution","nx",nx, 4)
  call read_param(PARAMS,"Resolution","ny",ny, 4)
  call read_param(PARAMS,"Time","nt",nt, 9999999)
  call read_param(PARAMS,"Time","Tmax",Tmax,1.d9)
  call read_param(PARAMS,"Time","CFL",cfl,0.1d0)
  call read_param(PARAMS,"Time","iMethod",iMethod,"RK2")
  call read_param(PARAMS,"Time","dt_fixed",dt_fixed,0.0d0)
  call read_param(PARAMS,"Time","dt_max",dt_max,0.0d0)
  call read_param(PARAMS,"ReynoldsNumber","nu",nu,1.d-2)

  ! Initial conditions section
  call read_param(PARAMS,"InitialCondition","inicond",inicond, "none")

  ! Dealasing section
  call read_param(PARAMS,"Dealiasing","iDealias",iDealias, 1)

  ! Penalization section
  call read_param(PARAMS,"Penalization","iMask",iMask, "none")
  call read_param(PARAMS,"Penalization","position_x",solid%position(1), 2.d0)
  call read_param(PARAMS,"Penalization","position_y",solid%position(2), 2.d0)
  call read_param(PARAMS,"Penalization","velocity_x",solid%velocity(1), 0.d0)
  call read_param(PARAMS,"Penalization","velocity_y",solid%velocity(2), 0.d0)
  call read_param(PARAMS,"Penalization","acceleration_x",solid%acceleration(1), 0.d0)
  call read_param(PARAMS,"Penalization","acceleration_y",solid%acceleration(2), 0.d0)
  call read_param(PARAMS,"Penalization","position_angle",solid%ang_position, 0.d0)
  call read_param(PARAMS,"Penalization","angular_velocity",solid%ang_velocity, 0.d0)
  call read_param(PARAMS,"Penalization","angular_acceleration",solid%ang_acceleration, 0.d0)
  call read_param(PARAMS,"Penalization","g",g, 0.d0)
  call read_param(PARAMS,"Penalization","Mass",Mass, 1.d0)
  call read_param(PARAMS,"Penalization","alpha",alpha, 1.57d0) ! pi/2 standard
  call read_param(PARAMS,"Penalization","leg_l",leg_l, 1.d0)
  call read_param(PARAMS,"Penalization","leg_h",leg_h, 0.2d0)
  call read_param(PARAMS,"Penalization","n_cell_smooth",n_cell_smooth, 3)
  call read_param(PARAMS,"Penalization","bounding_container",BC, "no")
  call read_param(PARAMS,"Penalization","buffer",buffer, 0.1d0)

  !Penalization parameter
  call read_param(PARAMS,"PenalizationParameter","C_k",C_k, 0.d0)
  call read_param(PARAMS,"PenalizationParameter","eps",eps, 1.d-2)

  ! Geometry section
  call read_param(PARAMS,"Geometry","xl",xl, 1.d0)
  call read_param(PARAMS,"Geometry","yl",yl, 1.d0)

  ! saving section
  call read_param(PARAMS,"Saving","tsave",tsave, 1.d0)
  call read_param(PARAMS,"Saving","itsave",itsave, 999999)
  call read_param(PARAMS,"Saving","tdrag",tdrag, 1.d0)
  call read_param(PARAMS,"Saving","iSaveVorticity",iSaveVorticity, 1)
  call read_param(PARAMS,"Saving","iSaveVelocity",iSaveVelocity, 1)
  call read_param(PARAMS,"Saving","iSavePressure",iSavePressure, 1)
  call read_param(PARAMS,"Saving","iSaveMask",iSaveMask, 1)
  call read_param(PARAMS,"Saving","iSaveSponge",iSaveSponge, 1)

  ! sponge
  call read_param(PARAMS,"Sponge","iSpongeType",iSpongeType, "none")
  call read_param(PARAMS,"Sponge","use_sponge",use_sponge, 0)
  call read_param(PARAMS,"Sponge","moving_sponge",moving_sponge, 0)
  call read_param(PARAMS,"Sponge","eps_sponge",eps_sponge, 1.0d0)
  call read_param(PARAMS,"Sponge","Sp_thickness",Sp_thickness, 10)

  ! mean flow
  call read_param(PARAMS,"MeanFlow","iMeanFlow",iMeanFlow, "none")
  call read_param(PARAMS,"MeanFlow","ux_mean",ux_mean, 0.d0)
  call read_param(PARAMS,"MeanFlow","uy_mean",uy_mean, 0.d0)

  ! Set other parameters (all procs)
  pi = 4.d0 * datan(1.d0)
  ! lattice spacing is global
  dx = xl/dble(nx)
  dy = yl/dble(ny)

  ! clean ini file
  call clean_ini_file(PARAMS)

  ! PenalizationParameter case
  if (C_k > 1.d-10)then
    eps = (C_k * dx)**2 / nu
    write(*,*) '*** INFO: C_k is set, the new eps is:', eps
  endif
end subroutine get_params
