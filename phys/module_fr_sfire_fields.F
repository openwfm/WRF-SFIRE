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
  
  ! 2d fire grid + layers
  real, pointer,  public, dimension(:, :, :)  :: surf_mass_frac !< surface mass fraction for each layer for heat fluxes computation [-]
  
  ! 2d atm grid
  real, pointer,  public, dimension(:, :) :: can_fx_z0          !< height of injection for canopy heat fluxes [m]
            
end module module_fr_sfire_fields
