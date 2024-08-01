#define DEBUG_OUT
!------------------------------------------------------------------------------------------------------------------
! MODULE: module_fr_sfire_fields
!   Store global fire variables in alphabetical order
!   
!   Each variable MUST be defined with short description and unit MUST be specified in scientific notation.
!   The short description MUST start with "!<" string to be compliant with fortran language server.
!------------------------------------------------------------------------------------------------------------------

module module_fr_sfire_fields

  implicit none
  
  ! fire params type
  type fire_params
    real, pointer,  dimension(:, :)     :: vx,vy                  !< wind velocity (m/s)
    real, pointer,  dimension(:, :)     :: zsf                    !< terrain height (m)
    real, pointer,  dimension(:, :)     :: dzdxf,dzdyf            !< terrain grad (1)
    real, pointer,  dimension(:, :)     :: bbb,phisc,phiwc,r_0    !< spread formula coefficients
    real, pointer,  dimension(:, :)     :: fgip                   !< init mass of surface fuel (kg/m^2)
    real, pointer,  dimension(:, :)     :: ischap                 !< 1 if chapparal
    real, pointer,  dimension(:, :)     :: fuel_time              !< time to burn to 1/e (s)
    real, pointer,  dimension(:, :)     :: fmc_g                  !< fuel moisture contents, ground (1)
    real, pointer,  dimension(:, :)     :: nfuel_cat              !< fuel category (integer values)
    real, pointer,  dimension(:, :)     :: tempf                  !< air temperature (K) (Balbi model only)
    real, pointer,  dimension(:, :)     :: rhof                   !< air pressure (Pa) (Balbi model only)
    real, pointer,  dimension(:, :)     :: ffwidth                !< fire front width (m) (Balbi model only)
    real, pointer,  dimension(:, :, :)  :: fmc_gc_f               !< fuel moisture contents, ground by class (1)
  end type fire_params
  
  
  ! 2d fire grid
  real, pointer,  public, dimension(:, :) :: can_bottom         !< canopy bottom height [m]
  real, pointer,  public, dimension(:, :) :: can_bulk_density   !< canopy bulk density [kg m-3]
  real, pointer,  public, dimension(:, :) :: can_fuel_burnt     !< canopy fuel burnt durint time step [kg m-2]
  real, pointer,  public, dimension(:, :) :: can_top            !< canopy top height [m]
  real, pointer,  public, dimension(:, :) :: cfb                !< crown fraction burnt [-]
  real, pointer,  public, dimension(:, :) :: mlhfx_eff_time     !< effective time after ignition for MLHFX scheme [-]
  real, pointer,  public, dimension(:, :) :: surf_moisture_frac !< surface moisture fraction for each layer for heat fluxes computation [-]
  
  ! 2d fire grid + layers
  real, pointer,  public, dimension(:, :, :)  :: mlhfx_fuel_mass_frac !< surface biomass mass fraction for each layer for MLHFX scheme [-]
  
  ! 2d atm grid
  real, pointer,  public, dimension(:, :) :: can_fx_z0          !< height of injection for canopy heat fluxes [m]  
  
  ! Multilayers Heat Fluxes scheme (MLHFX) data for 4 layers (1h, 10h, 100h, 1000h fuel size)
  real, pointer,  protected, dimension(:, :, :) :: mlhfx_consum_rate !< fuel fraction consumption rate data (4 layers, length of data given by mlhfx_dim, 2 for time and consumption) [s, s-1]
  integer,        protected, dimension(4)       :: mlhfx_len_data    !< length of time series for 1h, 10h, 100h, 1000h fuels

contains


!------------------------------------------------------------------------------------------------------------------
! SUBROUTINE: set_sfire_fields_from_grid
!>   Define pointers of SFIRE field module
!------------------------------------------------------------------------------------------------------------------
subroutine set_sfire_fields_from_grid(grid)

  use module_domain, only: domain

  implicit none

  type(domain),   intent(in)  :: grid

  ! pointer to grid variable
  cfb => grid%cfb
  can_top => grid%can_top
  can_bottom => grid%can_bottom
  can_bulk_density => grid%can_bulk_density
  can_fuel_burnt => grid%can_fuel_burnt
  can_fx_z0 => grid%can_fx_z0
  mlhfx_fuel_mass_frac => grid%mlhfx_fuel_mass_frac
  mlhfx_eff_time => grid%mlhfx_eff_time
  surf_moisture_frac => grid%surf_moisture_frac
  mlhfx_consum_rate => grid%mlhfx_consum_rate

end subroutine set_sfire_fields_from_grid


!------------------------------------------------------------------------------------------------------------------
! SUBROUTINE: initialize_sfire_fields
!>   Initialize fields defined in sfire_fields module
!!
!!  Set a default value for each field in the module  
!------------------------------------------------------------------------------------------------------------------
subroutine initialize_sfire_fields(restart, i_fire_can_fuel_read)
  
  implicit none

  ! 0. Declarations
  !    ============

  logical,  intent(in   ) :: restart              !< flag for restart
  integer,  intent(in   ) :: i_fire_can_fuel_read !< source of data for canopy data

  ! 1. Initialize value on array
  !    =========================

  if (.not. restart) then
    cfb(:, :)                 = 0.
    if (i_fire_can_fuel_read /= 1) then
      ! set default value only if variable not in input file
      can_top(:, :)             = 0.
      can_bottom(:, :)          = 0.
      can_bulk_density(:, :)    = 0.
    endif
    can_fuel_burnt(:, :)      = 0.
    surf_moisture_frac(:, :)  = 1.
    mlhfx_fuel_mass_frac(:, :, :) = 1.
    mlhfx_eff_time(:, :) = 0.
  endif

  can_fx_z0(:, :)           = 0.

end subroutine initialize_sfire_fields

!------------------------------------------------------------------------------------------------------------------
! SUBROUTINE: read_mlhf_data
!>   Read content of the file SFIRE_MLHF_BioSmoke.TBL and put content in mlhf_biosmoke_data structure
!! 
!! Read the necessary data for MultiLayers Heat Flux scheme (MLHF)
!! Load the fuel fraction consumption rate for the size classes (1h, 10, 100h, 1000h corresponding 
!! to 0.25, 1, 3, 8 in) 
!------------------------------------------------------------------------------------------------------------------
subroutine read_mlhf_data()
  
  use module_fr_sfire_util, only: crash, message, mlhfx_dim

  implicit none
  
  ! 0. Declarations
  !    ------------
  
  ! 0.1 Arguments

  ! 0.2 Locals
  character(len=*), parameter  :: filename = "SFIRE_MLHF_BioSmoke.TBL" !< name of the file containing the data
  integer, parameter :: io = 23 !< iounit
  integer :: istat, k, i_layer
  character(len=256) :: line, ierr, msg
  character(len=17) :: model
  
  ! 1. Initialize arrays
  !    -----------------
  
  mlhfx_consum_rate(:, :, :) = 0.
  
  ! 2. Open .TBL file and read content
  !    -------------------------------

  open(unit=io, file=filename, status="old", action="read", iostat=istat, iomsg=ierr)
  if (istat /= 0) then
    write(msg, *) "read_mlhf_data: error reading data file: ", filename
    call message(msg, level=0)
    call crash(ierr)
  endif
  
  ! read headers
  read(io, *) line
  read(io, *) line
  ! read data name
  read(io, '(A17)') model
  ! check if name is correct
  if (model == "BioSmoke 4 layers") then
    do i_layer = 1, 4
      ! dummy line
      read(io, *) line
      ! read number of class
      read(io, *) mlhfx_len_data(i_layer)
      ! write(*,*) "len data:", mlhfx_len_data(i_layer)
      ! check if size if ok
      if (mlhfx_len_data(i_layer) > mlhfx_dim) then
        call crash("read_mlhf_data: mlhfx_dim is too small for 1h data")
      endif
      ! dummy line
      read(io, *) line
      read(io, *) line
      ! read data
      do k = 1, mlhfx_len_data(i_layer)
        read(io, '(F10.1, 3X, E15.7)') mlhfx_consum_rate(i_layer, k, 1), mlhfx_consum_rate(i_layer, k, 2)
        ! write(*, *) "data ", k, mlhfx_consum_rate(i_layer, k, 1), mlhfx_consum_rate(i_layer, k, 2)
      enddo
    enddo
  endif

  close(unit=io)
  
end subroutine read_mlhf_data

end module module_fr_sfire_fields
