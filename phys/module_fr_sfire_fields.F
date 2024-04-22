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

  ! 2d fire grid
  real, pointer,  public, dimension(:, :) :: can_bottom         !< canopy bottom height [m]
  real, pointer,  public, dimension(:, :) :: can_bulk_density   !< canopy bulk density [kg m-3]
  real, pointer,  public, dimension(:, :) :: can_fuel_burnt     !< canopy fuel burnt durint time step [kg m-2]
  real, pointer,  public, dimension(:, :) :: can_top            !< canopy top height [m]
  real, pointer,  public, dimension(:, :) :: cfb                !< crown fraction burnt [-]
  real, pointer,  public, dimension(:, :) :: surf_moisture_frac !< surface moisture fraction for each layer for heat fluxes computation [-]
  
  ! 2d fire grid + layers
  real, pointer,  public, dimension(:, :, :)  :: surf_mass_frac_b !< surface biomass mass fraction for each layer for heat fluxes computation [-]
  real, pointer,  public, dimension(:, :, :)  :: surf_mass_frac_c !< surface char mass fraction for each layer for heat fluxes computation [-]
  
  ! 2d atm grid
  real, pointer,  public, dimension(:, :) :: can_fx_z0          !< height of injection for canopy heat fluxes [m]
  
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

end module module_fr_sfire_fields
