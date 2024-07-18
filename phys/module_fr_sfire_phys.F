!*** Jan Mandel August-October 2007 email: jmandel@ucar.edu or Jan.Mandel@gmail.com
!
! This file contains parts copied and/or adapted from earlier codes by 
! Terry Clark, Janice Coen, Don Latham, and Net Patton.

module module_fr_sfire_phys

use module_model_constants, only: cp,xlv
use module_fr_sfire_util
use module_fr_sfire_fields, only: fire_params
#ifdef USE_STARFIRE
use module_fr_starfire_phys, only: fire_ros_wirc
#endif

implicit none
PRIVATE

! verbosity and diagnostic switches, can be activated when this module is USEd
logical, public:: compare_rothermel=.false.       ! run comparison in fire_ros
logical, public:: call_write_fuels_m=.true.       ! create file fuels_m for graphics postprocessing 
integer, public:: balbi_msglevel=2                ! prints in fire_ros_balbi, 0=nothing 1=errors 2=warnings 3=diagnostics
integer, public:: balbi_nc_max=10                 ! max number of balbi no convergence cases
integer, public:: balbi_nc_count=0                ! counter of balbi no convergence cases
integer, public:: behave_msglevel=2               ! prints in fire_ros_behave, 0=nothing 1=errors 2=warnings 3=diagnostics

! subroutines and functions 
PUBLIC:: init_fuel_cats,fire_ros,heat_fluxes,set_nfuel_cat,set_fire_params, &
write_fuels_m,fire_risk,fire_intensity,fuel_moisture,advance_moisture,fuel_name,&
fire_rate_of_spread, set_canopy_data, initialize_sfire_fields

! types
public:: fire_params

! variables
PUBLIC:: fire_wind_height,fcz0,fcwh,have_fuel_cats,nfuelcats, &
  no_fuel_cat,no_fuel_cat2,windrf,moisture_classes,itree,fueldepthm,ibeh
PUBLIC:: mfuelcats
! NOTE: fcwh and fcz0 are called fwh and fz0 in read/write statements


! moisture behavior, see Mandel et al EGU 2012

!! To add moisture classes:
!  1. change parameter max_moisture_classes below
!  2. change the default of nfmc to the same value in Registry/registry.fire
!  3. add the appropriate lines real::fmc_gw<number>= <default values>
!  4. add dfault 

!*** dimensions
   INTEGER, PARAMETER :: mfuelcats = 100     ! max number of fuel categories 
  integer, parameter::max_moisture_classes=6
!***

  integer, parameter::zm=max_moisture_classes - 3
  integer:: moisture_classes=3
  real, dimension(max_moisture_classes):: drying_lag,wetting_lag,saturation_moisture,saturation_rain, &
         rain_threshold,rec_drying_lag_sec,rec_wetting_lag_sec
  integer, dimension(max_moisture_classes):: drying_model,wetting_model,fmc_gc_initialization,is_dead
   ! relative weights of moisture class for each fuel category
   integer::itmp
   CHARACTER (len=80), DIMENSION(max_moisture_classes), save :: moisture_class_name
   real, dimension(mfuelcats):: &  ! should sum up to one
     fmc_gw01=(/ (1.0, itmp=1,mfuelcats) /), & 
     fmc_gw02=(/ (0.0, itmp=1,mfuelcats) /), & 
     fmc_gw03=(/ (0.0, itmp=1,mfuelcats) /), & 
     fmc_gw04=(/ (0.0, itmp=1,mfuelcats) /), & 
     fmc_gw05=(/ (0.0, itmp=1,mfuelcats) /), &
     fmc_gw06=(/ (0.0, itmp=1,mfuelcats) /) 
   real, dimension(mfuelcats):: savr_gc01, savr_gc02, savr_gc03, savr_gc04, savr_gc05, savr_gc06
   real, dimension(mfuelcats):: fuelload_gc01, fuelload_gc02, fuelload_gc03, fuelload_gc04, fuelload_gc05, fuelload_gc06
    
  data moisture_class_name /'1-hour fuel','10-hour fuel','100-hour fuel',zm*'NOT USED'/
  data drying_lag                 /1., 10., 100. , zm*0./  ! time lag (h) approaching equilibrium moisture 
  data wetting_lag                   /14, 140., 1400., zm*0./  ! time lag (h) for approaching saturation in rain
  data saturation_moisture        /2.5, 2.5, 2.5 , zm*0./  ! saturation moisture contents (1) in rain
  data saturation_rain            /8.0, 8.0, 8.0 , zm*0./  ! stronger rain matters only in duration (mm/h) 
  data rain_threshold             /0.05, 0.05, 0.05, zm*0 /! rain intensity this small is same as nothing
  data drying_model               /1,   1,   1,    zm*1 /
  data wetting_model              /1,   1,   1,    zm*1 /
  data fmc_gc_initialization      /2,   2,   2,    zm*2 /  ! initialization 0=input, 1=from fuelmc_g, 2=from equilibrium
  data is_dead                    /1,   1,   1,    zm*0 /  ! identify moisture classes representing dead fuel
  real, dimension(7)::eq_p
  data eq_p/  1.035e-09, &   !(3.893e-10, 1.681e-09) ! coefficients of the equilibrium fuel moisture polynomial
              -2.62e-07, &   !(-4.593e-07, -6.473e-08) ! fitted from the graph in Schroeder and Buck (1970)
              2.507e-05, &   !(2.194e-06, 4.795e-05)
              -0.001107, &   !(-0.002353, 0.000139)
                0.02245, &   !(-0.009188, 0.05409)
               -0.05901, &   !(-0.3721, 0.254)
                  3.043/   !(2.17, 3.915)
   
! =========================================================================

! Following table copied from module_fr_cawfe_fuel by Ned Patton with minor changes. 
! Based on:           Clark, T. L., J. L. Coen and D. Latham: 2004, 
!                      "Description of a coupled atmosphere-fire model",
!                      International Journal of Wildland Fire, 13, 49-63.
!
! edited by Jan Mandel   jmandel@ucar.edu  September 2007
!
! - moved all fuel related constants and the initialization subroutine here
! - copied descriptions for fuel categories from fire_sfc.m4 in the original CAWFE code 
! This file had to be copied under a new name because packages in wrf physics
! layer are not allowed to call each other.

!D in col 2 means quantity derived from the others
!
!  Scalar constants (data same for all fuel categories):
!       HFGL           SURFACE FIRE HEAT FLUX THRESHOLD TO IGNITE CANOPY (W/m^2)
!       CMBCNST        JOULES PER KG OF DRY FUEL
!       FUELHEAT       FUEL PARTICLE LOW HEAT CONTENT, BTU/LB
!       FUELMC_G       FUEL PARTICLE (SURFACE) MOISTURE CONTENT
!D      BMST           RATIO OF LATENT TO SENSIBLE HEAT FROM SFC BURN:
!                        % of total fuel mass that is water (not quite
!                        = % fuel moisture).    BMST= (H20)/(H20+DRY)
!                        so BMST = FUELMC_G / (1 + FUELMC_G)  where
!                        FUELMC_G = ground fuel moisture
!
!  Data arrays indexed by fuel category:
!       FGI            INITIAL TOTAL MASS OF SURFACE FUEL (KG/M**2)
!       FUELDEPTHM     FUEL DEPTH, IN M  (CONVERTED TO FT)              
!       SAVR           FUEL PARTICLE SURFACE-AREA-TO-VOLUME RATIO, 1/FT
!                         GRASS: 3500., 10 hr fuel: 109.,  100 hr fuel: 30.
!       FUELMCE        MOISTURE CONTENT OF EXTINCTION; 0.30 FOR MANY DEAD FUELS; 0.15 FOR GRASS
!       FUELDENS       OVENDRY PARTICLE DENSITY, LB/FT^3
!       ST             FUEL PARTICLE TOTAL MINERAL CONTENT
!       SE             FUEL PARTICLE EFFECTIVE MINERAL CONTENT
!       WEIGHT         WEIGHTING PARAMETER THAT DETERMINES THE SLOPE OF THE MASS LOSS CURVE
!                        RANGES FROM ~5 (FAST BURNUP) TO 1000 ( ~40% DECR OVER 10 MIN).
!       FCI_D          INITIAL DRY MASS OF CANOPY FUEL
!       FCT            BURN OUT TIME FOR CANOPY FUEL, AFTER DRY (S)
!       ichap          1 if chaparral, 0 if not
!       itree          Tree type index, manually set, which corresposponds to Massman et al. 2017 
!                      tree profiles where 1 is Aspen, (2) Scots Pine, (3) Jack Pine, 
!                      (4) Loblolly Pine, (5) Hardwood, (6) Spruce, (7)Corn and (8) Rice. Fuel 
!                      categories not associated with trees (i.e tall grass) are set to zero.
!                      Added by DVM 6/30/2019. 
!       idynamic       1 if dynamic fuel, 0 if not
!D      FCI            INITIAL TOTAL MASS OF CANOPY FUEL
!D      FCBR           FUEL CANOPY BURN RATE (KG/M**2/S) 

!

! scalar fuel coefficients
   REAL, SAVE:: cmbcnst,hfgl,fuelmc_g,fuelmc_c, fire_wind_height,ffwidth_thresh
! computed values
   REAL, SAVE:: fuelheat

! defaults, may be changed in init_fuel_cats
   DATA cmbcnst  / 17.433e+06/             ! J/kg dry fuel
   DATA hfgl     / 17.e4 /                ! W/m^2
   DATA fuelmc_g / 0.08  /                ! set = 0 for dry ground fuel
   DATA fuelmc_c / 1.00  /                ! set = 0 for dry canopy
   DATA fire_wind_height/6.096/           ! m, 6.096m Behave, 0 to use fwh in each category
   DATA ffwidth_thresh /0.0/              ! Balbi ROS threshold for convective cooling (W_00)
!  REAL, PARAMETER :: bmst     = fuelmc_g/(1+fuelmc_g)
!  REAL, PARAMETER :: fuelheat = cmbcnst * 4.30e-04     ! convert J/kg to BTU/lb
!  real, parameter :: xlv      = 2.5e6                  ! to make it selfcontained
!  real, parameter :: cp      =  7.*287./2              ! to make it selfcontained


! fuel categorytables
   INTEGER, PARAMETER :: nf=14              ! fuel cats in data stmts, for fillers only`
   INTEGER, SAVE      :: nfuelcats = 13     ! number of fuel categories, can be reset from namelist.fire
   INTEGER, PARAMETER :: zf = mfuelcats-nf  ! number of zero fillers in data stmt 
   INTEGER, SAVE      :: no_fuel_cat = 14   ! special no fuel category outside of 1:nfuelcats
   INTEGER, SAVE      :: no_fuel_cat2 = 99999999  ! all categories between no_fuel_cat and no_fuel_cat2 are no fuel
   INTEGER, SAVE      :: ibeh=1             ! spread formula 0 = MacArthur, 1 = Rothermel, 2 = Balbi
   CHARACTER (len=80), DIMENSION(mfuelcats ), save :: fuel_name
   INTEGER, DIMENSION( mfuelcats ), save :: ichap
   INTEGER, DIMENSION( mfuelcats ), save :: idynamic
   INTEGER, DIMENSION( mfuelcats ), save :: itree !< tree type
   REAL   , DIMENSION( mfuelcats ), save :: windrf,weight,fgi,fci,fci_d,fct,fcbr, &
   fueldepthm,fueldens,fuelmce,   &
   fcwh,fcz0, ffw, &
   savr,st,se,adjr0,adjrw,adjrs, &
   fmc_gl_stdev,fmc_gl_ndwi_0,fmc_gl_ndwi_rate,fmc_gl_ndwi_stdev, k_tc
   REAL   , DIMENSION( mfuelcats , max_moisture_classes), save :: fmc_gw, savr_gc, fuelload_gc

   PUBLIC :: n_tree_type
   INTEGER, PARAMETER                       :: n_tree_type = 8          !< number of tree type categories
   REAL   , DIMENSION(n_tree_type), SAVE    :: tree_z_max               !< canopy top height (AGL) as function of tree type.
   REAL   , DIMENSION(n_tree_type), SAVE    :: tree_z_min               !< canopy bottom height (AGL) as function of tree type.
   REAL   , DIMENSION(n_tree_type), SAVE    :: tree_burn_time           !< canopy characteristic burn time

   ! LLNL fuel load data for 1h, 10h, 100h, live fuel and associated burn time
  integer,  parameter :: n_layer_hf_surf = 4        !< number of layers used for surface heat fluxes computation (for now fixed to 4)
  integer,  parameter :: n_fuel_cat_A13 = 13        !< number of fuel categories in Anderson fuel classification
  ! fuel_load_layers_hf is from: https://www.fs.usda.gov/Internet/FSE_DOCUMENTS/stelprdb5111790.pdf + conversion coef from [ton acre-1] to [kg m-2]
  real, dimension(n_fuel_cat_A13, n_layer_hf_surf), save  :: fuel_load_layers_hf = reshape(& 
    (/0.1658858, 0.       , 0.       , 0.       , &
      0.44834  , 0.22417  , 0.112085 , 0.112085 , &
      0.6747517, 0.       , 0.       , 0.       , &
      1.1230917, 0.8989217, 0.44834  , 1.1230917, &
      0.22417  , 0.112085 , 0.       , 0.44834  , &
      0.336255 , 0.560425 , 0.44834  , 0.       , &
      0.2533121, 0.4191979, 0.336255 , 0.0829429, &
      0.336255 , 0.22417  , 0.560425 , 0.       , &
      0.6545764, 0.0919097, 0.0336255, 0.       , &
      0.6747517, 0.44834  , 1.1230917, 0.44834  , &
      0.336255 , 1.0110067, 1.2351767, 0.       , &
      0.8989217, 3.1451051, 3.7055301, 0.       , &
      1.5714317, 5.1648768, 6.2879685, 0.        /), &
      shape(fuel_load_layers_hf), order=(/2, 1/)) !< fuel load for each fuel class and for each size class [kg m-2]
  real, dimension(n_layer_hf_surf), save  :: layer_burn_time_b  !< biomass burn time for 1h, 10h, 100h, and live fuel [s]
  real, dimension(n_layer_hf_surf), save  :: layer_burn_time_c  !< char burn time for 1h, 10h, 100h, and live fuel [s]
  real, dimension(n_layer_hf_surf), save  :: layer_char_fraction  !< char fraction for 1h, 10h, 100h, and live fuel [-] 
! =============================================================================
! Standard 13 fire behavior fuel models (for surface fires), along with some
!          estimated canopy properties (for crown fire).
! =============================================================================
   DATA fuel_name /  &
     'FUEL MODEL 1: Short grass (1 ft)', &
     'FUEL MODEL 2: Timber (grass and understory)', &
     'FUEL MODEL 3: Tall grass (2.5 ft)', &
     'FUEL MODEL 4: Chaparral (6 ft)', &
     'FUEL MODEL 5: Brush (2 ft) ', &
     'FUEL MODEL 6: Dormant brush, hardwood slash', &
     'FUEL MODEL 7: Southern rough', &
     'FUEL MODEL 8: Closed timber litter', &
     'FUEL MODEL 9: Hardwood litter', &
     'FUEL MODEL 10: Timber (litter + understory)', &
     'FUEL MODEL 11: Light logging slash', &
     'FUEL MODEL 12: Medium logging slash', &
     'FUEL MODEL 13: Heavy logging slash', &
     'FUEL MODEL 14: no fuel', &
     zf*' '/
   DATA windrf /0.36, 0.36, 0.44,  0.55,  0.42,  0.44,  0.44, &     
                0.36, 0.36, 0.36,  0.36,  0.43,  0.46,  1e-7, zf*0 / ! added jmandel October 2010
   DATA fgi / 0.166, 0.897, 0.675, 2.468, 0.785, 1.345, 1.092, &
              1.121, 0.780, 2.694, 2.582, 7.749, 13.024, 1.e-7, zf*0.  /
   DATA fueldepthm /0.305,  0.305,  0.762, 1.829, 0.61,  0.762,0.762, &
                    0.0610, 0.0610, 0.305, 0.305, 0.701, 0.914, 0.305,zf*0. /
   DATA savr / 3500., 2784., 1500., 1739., 1683., 1564., 1562.,  &
               1889., 2484., 1764., 1182., 1145., 1159., 3500., zf*0. /
   DATA savr_gc01 / 3500., 3000., 1500., 2000., 2000., 1750., 1750.,  &
               2000., 2500., 2000., 1500., 1500., 1500., 2000., zf*0. /
   DATA savr_gc02 / nf * 109., zf*9999. /
   DATA savr_gc03 / nf * 30., zf*9999. /
   DATA savr_gc04 / nf * 9999., zf*9999. / ! 1000h not use yet
   DATA savr_gc05 / nf * 9999., zf*9999. /
   DATA savr_gc06 / 9999., 1500., 9999., 1500., 1500., 9999., 1500.,  &
               9999., 9999., 1500., 9999., 9999., 9999., 9999., zf*0. /
   DATA fuelload_gc01 / 0.034, 0.092, 0.138, 0.23, 0.046, 0.069, 0.052, &
               0.069, 0.134, 0.138, 0.069, 0.184, 0.322, 0.001, zf*0. /
   DATA fuelload_gc02 / 0., 0.046, 0., 0.184, 0.023, 0.115, 0.086, &
               0.046, 0.019, 0.092, 0.207, 0.644, 1.058, 0., zf*0. /
   DATA fuelload_gc03 / 0., 0.023, 0., 0.092, 0., 0.092, 0.069, &
               0.115, 0.007, 0.23, 0.253, 0.759, 1.288, 0., zf*0. /
   DATA fuelload_gc04 / nf * 0., zf*0. /
   DATA fuelload_gc05 / nf * 0., zf*0. /
   DATA fuelload_gc06 / 0., 0.023, 0., 0.23, 0.092, 0., 0.017, &
               0., 0., 0.092, 0., 0., 0., 0., zf*0. /
   DATA fuelmce / 0.12, 0.15, 0.25, 0.20, 0.20, 0.25, 0.40,  &
                  0.30, 0.25, 0.25, 0.15, 0.20, 0.25, 0.12 , zf*0. / 
   DATA fueldens / nf * 32., zf*0. /   ! 32 if solid, 19 if rotten.
   DATA st / nf* 0.0555 , zf*0./
   DATA se / nf* 0.010 , zf*0./
! ----- Notes on weight: (4) - best fit of Latham data;
!                 (5)-(7) could be 60-120; (8)-(10) could be 300-1600;
!                 (11)-(13) could be 300-1600
   DATA weight / 7.,  7.,  7., 180., 100., 100., 100.,  &
              900., 900., 900., 900., 900., 900., 7. , zf*0./ 
! ----- 1.12083 is 5 tons/acre.  5-50 tons/acre orig., 100-300 after blowdown
   DATA fci_d / 0., 0., 0., 1.123, 0., 0., 0.,  &
            1.121, 1.121, 1.121, 1.121, 1.121, 1.121, 0., zf*0./
   DATA fct / 60., 60., 60., 60., 60., 60., 60.,  &
            60., 120., 180., 180., 180., 180. , 60. , zf*0.   /
   DATA ichap / 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 , zf*0/
   DATA itree / 0, 5, 0, 0, 0, 5, 5, 5, 5, 5, 0, 0, 0, 0 , zf*0/  
   DATA idynamic / 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 , zf*0/  
   DATA fcwh  / 6.096, 6.096, 6.096, 6.096, 6.096, 6.096, 6.096, &
         6.096, 6.096, 6.096, 6.096, 6.096, 6.096, 6.096, zf*0. / ! consistent with BEHAVE
! roughness length 0.13*fueldepthm except cat 3 fz0=0.1 for consistency with landuse
 ! fz0 =   0.0396,0.0396,0.0991,0.2378,0.0793,0.0991,0.0991,
   DATA fcz0 /  0.0396,0.0396,0.1000,0.2378,0.0793,0.0991,0.0991, &
         0.0079,0.0079,0.0396,0.0396,0.0911,0.1188,0.0396, zf * 0. /
   !DATA fcwh /0.6 , 0.6,  1.5, 36,  1.2,  1.5, 1.5, &
   !           0.12, 0.12, 0.6, 0.6, 1.38, 1.8, 1.8, zf*0 / ! fuel wind height, added jm 2/23/11
   !DATA fcz0 /0.3, 0.3, 0.75, 18.,  0.6,  0.75, 0.75, &
   !           0.06, 0.06, 0.3, 0.3, 0.69, 0.9, 0.9, zf*0 /  ! fuel roughness height, added jm 2/23/11
   DATA ffw /nf* 0.9, zf*0/
   DATA fmc_gl_ndwi_0 /nf*0.1, zf*0./
   DATA fmc_gl_ndwi_rate /nf*0.6, zf*0./
   DATA fmc_gl_ndwi_stdev /nf*0.2, zf*0./
   DATA fmc_gl_stdev /nf*0.2, zf*0./
   DATA adjr0 /mfuelcats*1./
   DATA adjrw /mfuelcats*1./
   DATA adjrs /mfuelcats*1./
   DATA k_tc / 4.025e-02, 1.621e-02, 3.752e-02, 3.000e-02, 1.812e-02, 3.621e-02, 3.778e-02, &
    8.770e-04, 1.781e-03, 8.216e-03, 7.405e-03, 1.772e-02, 2.169e-02, 3.000e-02, zf*0 / ! fuel thermal conductivity

   DATA tree_z_max / 10., 20., 15., 24., 22.5, 10., 2.2, 0.73 /         ! canopy top height as function of tree type [m]
   DATA tree_z_min / 3. , 6. , 5. , 12., 10. , 5. , 1.,  0.73 /         ! canopy bottom height as function of tree type [m]
   DATA tree_burn_time / 10. , 10. , 10. , 10., 10. , 10. , 10., 10. /  ! canopy burning time function of tree type [m]
   
   DATA layer_burn_time_b / 7., 500., 3600., 900. /   ! burn time for biomass [s]
   DATA layer_burn_time_c / 7., 500., 3600., 900. /   ! burn time for char [s]
   DATA layer_char_fraction / 0.2, 0.3, 0.4, 0.5 / ! char fraction [-]
   ! =========================================================================

  logical, save :: have_fuel_cats=.false.

contains

subroutine fuel_moisture(                &
    id,                                  & ! for prints and maybe file names
    nfmc,                                &
    ids,ide, jds,jde,                    & ! atm grid dimensions
    ims,ime, jms,jme,                    &
    ips,ipe,jps,jpe,                     &
    its,ite,jts,jte,                     &
    ifds, ifde, jfds, jfde,              & ! fire grid dimensions
    ifms, ifme, jfms, jfme,              &
    ifts,ifte,jfts,jfte,                 &
    ir,jr,                               & ! atm/fire grid ratio
    nfuel_cat,                           & ! fuel data
    fndwi,                                & ! satellite sensing on fire grid
    fmc_gc,                              & ! moisture contents by class on atmospheric grid
    fmc_g,                               & ! weighted fuel moisture contents on fire grid
    fmc_gc_f                              & ! moisture contents by class on fire grid
    )

implicit none

!**** arguments
integer, intent(in)::                    &
    id,nfmc,                             &
    ids,ide, jds,jde,                    & ! atm grid dimensions
    ims,ime, jms,jme,                    &
    ips,ipe,jps,jpe,                     &
    its,ite,jts,jte,                     &
    ifds, ifde, jfds, jfde,              & ! fire grid dimensions
    ifms, ifme, jfms, jfme,              &
    ifts,ifte,jfts,jfte,                 &
    ir,jr                                  ! atm/fire grid ratio


real,intent(in),dimension(ifms:ifme,jfms:jfme):: nfuel_cat, & ! fuel data
                                       fndwi                  ! satellite sensing interpolated to fire grid
real,intent(inout),dimension(ims:ime,nfmc,jms:jme):: fmc_gc
real,intent(out),dimension(ifms:ifme,jfms:jfme):: fmc_g ! fuel data
real,intent(out),dimension(ifms:ifme,nfmc,jfms:jfme):: fmc_gc_f ! fuel data

!**** local
real, dimension(its-1:ite+1,jts-1:jte+1):: fmc_k  ! copy of fmc_gc(:,k,:) 
real, dimension(ifts:ifte,jfts:jfte):: fmc_f,  &   ! interpolation of fmc_gc(:,k,:) to the fire grid 
                                       nwdi_f      ! inerpolation of nwdi to the fire grid
integer::i,j,k,n
integer::ibs,ibe,jbs,jbe
real::f1,w1,w2,f2,fa,fc

character(len=128)::msg

call check_mesh_2dim(ifts,ifte,jfts,jfte,ifds,ifde,jfds,jfde) ! check if fire tile fits into domain
call check_mesh_2dim(ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme) ! check if fire tile fits into memory

do j=jfts,jfte
    do i=ifts,ifte
        fmc_g(i,j)=0.               ! initialize sum over classes
    enddo
enddo

! one beyond the tile but not beyond the domain boundary
ibs=max(ids,its-1)
ibe=min(ide,ite+1)
jbs=max(jds,jts-1)
jbe=min(jde,jte+1)

call check_mesh_2dim(ibs,ibe,jbs,jbe,ims,ime,jms,jme) ! check if tile with halo fits into memory

do k=1,moisture_classes

    ! copy halo beyond the tile but not beyond the domain boundary
    do j=jbs,jbe
        do i=ibs,ibe
            fmc_k(i,j)=fmc_gc(i,k,j)      ! copy slice to 2d array
        enddo
    enddo

    call print_2d_stats(ibs,ibe,jbs,jbe,its-1,ite+1,jts-1,jte+1,fmc_k,'fuel_moisture: fmc_k')

    ! interpolate moisture contents in the class k to the fire mesh
    call interpolate_z2fire(id,0,        & ! for debug output, <= 0 no output
    ids,ide,jds,jde,                     & ! atm grid dimensions
    its-1,ite+1,jts-1,jte+1,             & ! memory dimensions
    ips,ipe,jps,jpe,                     &
    its,ite,jts,jte,                     &
    ifds, ifde, jfds, jfde,              & ! fire grid dimensions
    ifts, ifte, jfts, jfte,              &
    ifts,ifte,  jfts,jfte,               &
    ir,jr,                               & ! atm/fire grid ratio
    fmc_k,                               & ! atm grid arrays in
    fmc_f)                                 ! fire grid arrays out

    call print_2d_stats(ifts,ifte,jfts,jfte,ifts,ifte,jfts,jfte,fmc_f,'fuel_moisture: fmc_f')
    
    do j=jfts,jfte
        do i=ifts,ifte
            fmc_gc_f(i,k,j)=fmc_f(i,j)
        enddo
    enddo

    if(k .eq. kfmc_ndwi)then  ! if live moisture, assimilate ndwi
        call print_2d_stats(ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme,fndwi,'fuel_moisture: fndwi')
        write(msg,'(a,i4)')'Assimilating NDWI in fuel moisture class ',k
        call message(msg)
    endif

    ! add moisture contents for class k to the fuel moisture
    do j=jfts,jfte
        do i=ifts,ifte
            n = nfuel_cat(i,j)
            if(n > 0)then
                if(k .ne. kfmc_ndwi)then 
                   fmc_g(i,j)=fmc_g(i,j)+fmc_gw(n,k)*fmc_f(i,j)      ! add to sum over classes
                else                                                    ! if live moisture, assimilate ndwi
                   f1=fmc_f(i,j)
                   w1 = fmc_gl_stdev(n)                            
                   w1 = 1./(w1*w1)                                      ! weight of forecast
                   w2 = fmc_gl_ndwi_stdev(n)
                   w2 = 1./(w2*w2)                                      ! weight of update
                   f2 = fmc_gl_ndwi_0(n) + fmc_gl_ndwi_rate(n) * fndwi(i,j)   ! from regression
                   fa = (w1*f1 + w2*f2) / (w1 + w2)                     ! updated value
                   fc = fmc_gw(n,k)*fa                                  ! times proportion of live fuel
                   fmc_g(i,j)=fmc_g(i,j)+fc      ! add to sum over classes
                   ! write(*,*)'NDWI:',i,j,f1,f2,w1,w2,f1,fa,fmc_gw(n,k),fc,fmc_g(i,j)
                endif
            endif
        enddo
    enddo

enddo

call print_2d_stats(ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme,fmc_g,'fuel_moisture: fmc_g')
call print_3d_stats_by_slice(ifts,ifte,1,moisture_classes,jfts,jfte,ifms,ifme,1,nfmc,jfms,jfme,fmc_gc_f,'fmc_gc_f')

end subroutine fuel_moisture

subroutine advance_moisture(    &
    initialize,                 & ! initialize timestepping. true on the first call at time 0, then false
    ims,ime,  jms,jme,          & ! memory dimensions
    its,ite,  jts,jte,          & ! tile dimensions
    nfmc,                       & ! dimension of moisture fields
    moisture_dt,                & ! timestep = time step time elapsed from the last call
    fmep_decay_tlag,            & ! moisture extended model assimilated diffs. decay time lag
    rainc, rainnc,              & ! accumulated rain 
    t2, q2, psfc,               & ! temperature (K), vapor contents (kg/kg), pressure (Pa) at the surface
    rain_old,                   & ! previous value of accumulated rain
    t2_old, q2_old, psfc_old,   & ! previous values of the atmospheric state at surface
    rh_fire,                    & ! relative humidity at surface, for diagnostic only
    fmc_gc,                     & ! fuel moisture by class, updated
    fmep,                       & ! fuel moisture extended model parameters
    fmc_equi,                   & ! fuel moisture equilibrium by class, for diagnostics only
    fmc_lag                     & ! fuel moisture tendency by classe, for diagnostics only
    )

implicit none

!*** arguments
logical, intent(in):: initialize
integer, intent(in)::           &
    ims,ime,  jms,jme,          & ! memory dimensions
    its,ite,  jts,jte,          & ! tile dimensions
    nfmc                          ! number of moisture fields
real, intent(in):: moisture_dt, fmep_decay_tlag
real, intent(in), dimension(ims:ime,jms:jme):: t2, q2, psfc, rainc, rainnc
real, intent(inout), dimension(ims:ime,jms:jme):: t2_old, q2_old, psfc_old, rain_old 
real, intent(inout), dimension(ims:ime,nfmc,jms:jme):: fmc_gc
real, intent(inout), dimension(ims:ime,2,jms:jme):: fmep
real, intent(out), dimension(ims:ime,nfmc,jms:jme):: fmc_equi, fmc_lag
real, intent(out), dimension(ims:ime,jms:jme)::rh_fire 

!*** global
! fuel properties moisture set by init_fuel_cats

!*** local
integer:: i,j,k
real::rain_int, T, P, Q, QRS, ES, RH, tend, EMC_d, EMC_w, EMC, R, rain_diff, fmc, rlag, equi, &
    d, w, rhmax, rhmin, change, rainmax,rainmin, fmc_old, H, deltaS, deltaE
real, parameter::tol=1e-2 ! relative change larger than that will switch to exponential ode solver 
character(len=256)::msg
logical::bad_wrf
integer::msglevel=2
logical, parameter::check_rh=.false.
integer::check_data=2  ! 0=nothing, 1=replace quietly, 2=warning (also printed if msglevel>2), 3=crash
real::epsilon,Pws,Pw,t2_min,q2_min,psfc_min
real::t2_floor=200.    ! minimum allowed temperature (K)
real::q2_floor=1e-8    ! miniumu allowed moisture contents (kg/kg)
real::psfc_floor=1000. ! mimimum allowed surface pressiure (Pa)


!*** executable

! check arguments

if(msglevel>1)then
!$OMP CRITICAL(SFIRE_PHYS_CRIT)
    write(msg,'(a,f10.2,a,i4,a,i4)')'advance_moisture dt=',moisture_dt,'s using ',moisture_classes,' classes from possible ',nfmc
!$OMP END CRITICAL(SFIRE_PHYS_CRIT)
    call message(msg,level=2)
endif

if(moisture_classes > nfmc .or. moisture_classes > max_moisture_classes)then
!$OMP CRITICAL(SFIRE_PHYS_CRIT)
    write(msg,*)'advance_moisture: moisture_classes=',moisture_classes, &
       ' > nfmc=',nfmc,' or >  max_moisture_classes=',max_moisture_classes
!$OMP END CRITICAL(SFIRE_PHYS_CRIT)
    call crash(msg)
endif

call print_2d_stats(its,ite,jts,jte,ims,ime,jms,jme,t2,'T2')
call print_2d_stats(its,ite,jts,jte,ims,ime,jms,jme,q2,'Q2')
call print_2d_stats(its,ite,jts,jte,ims,ime,jms,jme,psfc,'PSFC')

if(initialize) then 
    call message('advance_moisture: initializing, copying surface variables to old')
    call copy2old
else
    call print_3d_stats_by_slice(its,ite,1,moisture_classes,jts,jte,ims,ime,1,nfmc,jms,jme,fmc_gc,'before advance fmc_gc')
endif

if(check_data.ge.2 .or. msglevel.ge.2)then
    t2_min = huge(t2_min)
    q2_min = huge(q2_min)
    psfc_min = huge(psfc_min)
    do j=jts,jte
        do i=its,ite
            t2_min=min(t2(i,j),t2_min)
            q2_min=min(q2(i,j),q2_min)
            psfc_min=min(psfc(i,j),psfc_min)
        enddo
    enddo
    bad_wrf = ( t2_min<t2_floor .or. psfc_min<psfc_floor .or. q2_min<q2_floor)
    if(bad_wrf .or. msglevel.ge.2)then
!$OMP CRITICAL(SFIRE_PHYS_CRIT)
91      format(a,3(2x,a,e11.3))
        write(msg,91)'minimal  ','t2',t2_min,'q2',q2_min,'psfc',psfc_min
        call message(msg,level=0) 
        write(msg,91)'floor    ','t2',t2_floor,'q2',q2_floor,'psfc',psfc_floor
        call message(msg,level=0) 
!$OMP END CRITICAL(SFIRE_PHYS_CRIT)
    endif
    if(bad_wrf)then
        if(check_data.ge.3)then
            call crash('advance_moisture: invalid data passed from WRF')
        else
            call message('WARNING: advance_moisture: nonphysical input values replaced by floor',level=0)
        endif
    endif
endif

! one time step

rhmax=-huge(rhmax)
rhmin=huge(rhmin)
rainmax=-huge(rainmax)
rainmin= huge(rainmin)
do j=jts,jte
    do k=1,moisture_classes
        do i=its,ite
            ! old fuel moisture contents
            ! compute the rain intensity from the difference of accumulated rain
            rain_diff = ((rainc(i,j) + rainnc(i,j)) - rain_old(i,j))
            if(moisture_dt > 0.)then
                rain_int  = 3600. * rain_diff / moisture_dt 
            else
                rain_int  = 0.
            endif
            rainmax = max(rainmax,rain_int)
            rainmin = min(rainmin,rain_int)
            R = rain_int - rain_threshold(k)

            ! average the inputs for second order accuracy
            T = 0.5 * (t2_old(i,j) + t2(i,j))
            P = 0.5 * (psfc_old(i,j) + psfc(i,j))
            Q = 0.5 * (q2_old(i,j) + q2(i,j))

            ! replace nonphysical values by floor
            if(check_data .ge. 1)then
                T = max(T,t2_floor)
                P = max(P,psfc_floor)
                Q = max(Q,q2_floor)
        endif

            ! compute the relative humidity
            ! ES=610.78*exp(17.269*(T-273.161)/(T-35.861))
            ! QRS=0.622*ES/(P-0.378*ES)
            ! RH = Q/QRS
            ! function rh_from_q from Adam Kochanski following Murphy and Koop, Q.J.R. Meteorol. Soc (2005) 131 1539-1565 eq. (10)
            epsilon = 0.622 ! Molecular weight of water (18.02 g/mol) to molecular weight of dry air (28.97 g/mol)
            ! vapor pressure [Pa]
            Pw=q*P/(epsilon+(1-epsilon)*q); 
            ! saturation vapor pressure [Pa]
            Pws= exp( 54.842763 - 6763.22/T - 4.210 * log(T) + 0.000367*T + &
                tanh(0.0415*(T - 218.8)) * (53.878 - 1331.22/T - 9.44523 * log(T) + 0.014025*T))
            !realtive humidity [1]
            RH = Pw/Pws
            rh_fire(i,j)=RH
            rhmax=max(RH,rhmax)
            rhmin=min(RH,rhmin)

            deltaE = fmep(i,1,j)
            deltaS = fmep(i,2,j)

            if(.not.check_rh)then
                RH = min(RH,1.0)
            else
                if(RH < 0.0 .or. RH > 1.0 .or. RH .ne. RH )then
!$OMP CRITICAL(SFIRE_PHYS_CRIT)
                    write(msg,'(a,2i6,5(a,f10.2))')'At i,j ',i,j,' RH=',RH, &
                        ' from T=',T,' P=',P,' Q=',Q
                    call message(msg) 
                    call crash('Relative humidity must be between 0 and 1, saturated water contents must be >0')
!$OMP END CRITICAL(SFIRE_PHYS_CRIT)
                endif
            endif 
            !print *,'ADV_MOIST i=',i,' j=',j,' T=',T,' P=',P,' Q=',Q,' ES=',ES,' QRS=',QRS,' RH=',RH

            if (R > 0.) then
                select case(wetting_model(k))
                case(1) ! saturation_moisture=2.5 wetting_lag=14h saturation_rain=8 mm/h calibrated to VanWagner&Pickett 1985 per 24 hours
                    EMC_w=saturation_moisture(k) + deltaS
                    EMC_d=saturation_moisture(k) + deltaS
                    rlag=rec_wetting_lag_sec(k) * (1. - exp(-R/saturation_rain(k)))
                end select
            else ! not raining
                select case(drying_model(k))
                case(1) ! Van Wagner and Pickett (1972, 1985) per Viney (1991) eq (7) (8)
                    H = RH * 100.
                    d=0.942*H**0.679 + 0.4994e-4*exp(0.1*H) + 0.18*(21.1+273.15-T)*(1-exp(-0.115*H)) ! equilibrium moisture for drying
                    w=0.618*H**0.753 + 0.4540e-4*exp(0.1*H) + 0.18*(21.1+273.15-T)*(1-exp(-0.115*H)) ! equilibrium moisture for adsorbtion
                    if(d.ne.d.or.w.ne.w)call crash('equilibrium moisture calculation failed, result is NaN')
                    d = d*0.01
                    w = w*0.01
                    EMC_d = max(max(d,w)+deltaE,0.0)
                    EMC_w = max(min(d,w)+deltaE,0.0)
                    rlag=rec_drying_lag_sec(k)
                end select
            endif
            !*** MODELS THAT ARE NOT OF THE EXPONENTIAL TIME LAG KIND 
            ! ARE RESPONSIBLE FOR THEIR OWN LOGIC, THESE MODELS
            ! SHOULD COMPUTE fmc_gc(i,k,j) DIRECTLY AND SET TLAG < 0
            !
            if(rlag > 0.0)then

                if(.not.initialize .or. fmc_gc_initialization(k).eq.0)then ! take old from before, no initialization
                    fmc_old = fmc_gc(i,k,j)
                elseif(fmc_gc_initialization(k).eq.1)then ! from scalar fuelmc_g
                    fmc_old = fuelmc_g
                elseif(fmc_gc_initialization(k).eq.2)then ! from computed equilibrium
                    fmc_old=0.5*(EMC_d+EMC_w)
                else
                    call crash('bad value of fmc_gc_initialization(k), must be between 0 and 2')
                endif
                equi = max(min(fmc_old, EMC_d),EMC_w) ! take lower or upper equilibrium value 

                change = moisture_dt * rlag

                if(change  < tol)then
                     if(fire_print_msg.ge.3)call message('midpoint method')
                     fmc = fmc_old + (equi - fmc_old)*change*(1.0 - 0.5*change)  ! 2nd order Taylor
                else
                     if(fire_print_msg.ge.3)call message('exponential method')
                     fmc = fmc_old + (equi - fmc_old)*(1 - exp(-change))
                endif
                fmc_gc(i,k,j) = fmc

                ! diagnostics out
                fmc_equi(i,k,j)=equi
                fmc_lag(i,k,j)=1.0/(3600.0*rlag)

                ! diagnostic prints
                if(fire_print_msg.ge.3)then
!$OMP CRITICAL(SFIRE_PHYS_CRIT)
            write(msg,*)'i=',i,' j=',j,'EMC_w=',EMC_w,' EMC_d=',EMC_d
                    call message(msg)
                    write(msg,*)'fmc_old=',fmc,' equi=',equi,' change=',change,' fmc=',fmc
                    call message(msg)
!$OMP END CRITICAL(SFIRE_PHYS_CRIT)
                endif

            endif
      enddo
   enddo
enddo


! assimilated differences decay
do j=jts,jte
  do k=1,2
    do i=its,ite
      change = moisture_dt / (fmep_decay_tlag * 3600.)
      if(change < tol) then
        fmep(i,k,j) = fmep(i,k,j)*(1.0 - change * (1.0 - 0.5 * change))
      else
        fmep(i,k,j) = fmep(i,k,j)*exp(-change)
      endif
    enddo
  enddo
enddo


if(fire_print_msg.ge.2)then
!$OMP CRITICAL(SFIRE_PHYS_CRIT)
    write(msg,2)'Rain intensity    min',rainmin,  ' max',rainmax,' mm/h'
    call message(msg) 
    if(rainmin <0.)then
        call message('WARNING rain accumulation must increase')
    endif
    write(msg,2)'Relative humidity min',100*rhmin,' max',100*rhmax,'%'
    call message(msg) 
    if(.not.(rhmax<=1.0 .and. rhmin>=0))then
        call message('WARNING Relative humidity must be between 0 and 100%')
    endif
2   format(2(a,f10.2),a)
!$OMP END CRITICAL(SFIRE_PHYS_CRIT)
endif

call print_3d_stats_by_slice(its,ite,1,moisture_classes,jts,jte,ims,ime,1,nfmc,jms,jme,fmc_equi,'equilibrium fmc_equi')
call print_3d_stats_by_slice(its,ite,1,moisture_classes,jts,jte,ims,ime,1,nfmc,jms,jme,fmc_lag,'time lag')
call print_3d_stats_by_slice(its,ite,1,moisture_classes,jts,jte,ims,ime,1,nfmc,jms,jme,fmc_gc,'after advance fmc_gc')

call copy2old

return

contains

subroutine copy2old

do j=jts,jte
    do i=its,ite
        rain_old(i,j) = rainc(i,j) + rainnc(i,j)
        t2_old(i,j) = t2(i,j)
        q2_old(i,j) = q2(i,j)
        psfc_old(i,j) = psfc(i,j)
    enddo
enddo

end subroutine copy2old

subroutine get_equi_moist
end subroutine get_equi_moist

end subroutine advance_moisture


subroutine init_fuel_cats(init_fuel_moisture)
implicit none
!*** purpose: initialize fuel tables and variables by constants
!*** arguments: 
logical, intent(in)::init_fuel_moisture
logical, external:: wrf_dm_on_monitor
!$ integer, external:: OMP_GET_THREAD_NUM
!*** local
integer:: i,j,k,ii,iounit,ierr,kk
character(len=128):: msg
REAL   , DIMENSION( mfuelcats ) :: fwh, fz0
!*** executable

! read 
namelist /fuel_scalars/ cmbcnst,hfgl,fuelmc_g,fuelmc_c,nfuelcats,no_fuel_cat,no_fuel_cat2,fire_wind_height,ibeh, &
    ffwidth_thresh, balbi_msglevel
namelist /fuel_categories/ fuel_name,windrf,fgi,fueldepthm,savr, &
    fuelmce,fueldens,st,se,weight,fci_d,fct,ichap,itree,fwh,fz0,ffw,idynamic, &
    fmc_gl_ndwi_0,fmc_gl_ndwi_rate,fmc_gl_ndwi_stdev, fmc_gl_stdev, &
    adjr0,adjrw,adjrs,fmc_gw01,fmc_gw02,fmc_gw03,fmc_gw04,fmc_gw05,fmc_gw06, &
    savr_gc01,savr_gc02,savr_gc03,savr_gc04,savr_gc05,savr_gc06, &
    fuelload_gc01,fuelload_gc02,fuelload_gc03,fuelload_gc04,fuelload_gc05,fuelload_gc06, &
    is_dead, k_tc, tree_z_min, tree_z_max, tree_burn_time, layer_burn_time_b, layer_burn_time_c, layer_char_fraction, fuel_load_layers_hf

namelist /moisture/ moisture_classes,drying_lag,wetting_lag,saturation_moisture,saturation_rain,rain_threshold, &
    drying_model,wetting_model,moisture_class_name,fmc_gc_initialization,is_dead

!$  if (OMP_GET_THREAD_NUM() .ne. 0)then
!$     call crash('init_fuel_cats: must be called from master thread')
!$  endif
 
IF ( wrf_dm_on_monitor() ) THEN
    ! we are the master task
    
    ! copy in defaults
    fwh=fcwh
    fz0=fcz0
    
    ! read the file
    iounit=open_text_file('namelist.fire','read')
    read(iounit, nml=fuel_categories, iostat=ierr)
    if(ierr.ne.0)then
          call message('init_fuel_cats: error reading namelist fuel_categories in file namelist.fire')
          call error_namelist(iounit)
    endif
    rewind(iounit)
    read(iounit, nml=fuel_scalars, iostat=ierr)
    if(ierr.ne.0)then
          call message('init_fuel_cats: error reading namelist fuel_scalars in file namelist.fire',0)
          call error_namelist(iounit)
    endif
    rewind(iounit)
    if(init_fuel_moisture)then
        read(iounit, nml=moisture, iostat=ierr)
        if(ierr.ne.0)then
          call message('init_fuel_cats: error reading namelist moisture in file namelist.fire')
          call error_namelist(iounit)
        endif
        rewind(iounit)
    endif
    fmc_gw(1:mfuelcats,1)=fmc_gw01
    fmc_gw(1:mfuelcats,2)=fmc_gw02
    fmc_gw(1:mfuelcats,3)=fmc_gw03
    fmc_gw(1:mfuelcats,4)=fmc_gw04
    fmc_gw(1:mfuelcats,5)=fmc_gw05
    fmc_gw(1:mfuelcats,6)=fmc_gw06
    savr_gc(1:mfuelcats,1)=savr_gc01
    savr_gc(1:mfuelcats,2)=savr_gc02
    savr_gc(1:mfuelcats,3)=savr_gc03
    savr_gc(1:mfuelcats,4)=savr_gc04
    savr_gc(1:mfuelcats,5)=savr_gc05
    savr_gc(1:mfuelcats,6)=savr_gc06
    fuelload_gc(1:mfuelcats,1)=fuelload_gc01
    fuelload_gc(1:mfuelcats,2)=fuelload_gc02
    fuelload_gc(1:mfuelcats,3)=fuelload_gc03
    fuelload_gc(1:mfuelcats,4)=fuelload_gc04
    fuelload_gc(1:mfuelcats,5)=fuelload_gc05
    fuelload_gc(1:mfuelcats,6)=fuelload_gc06
    CLOSE(iounit)
    
    ! copy out to permanent names
    fcwh=fwh
    fcz0=fz0
    
    if (nfuelcats>mfuelcats) then
        write(msg,*)'nfuelcats=',nfuelcats,' too large, increase mfuelcats'
        call crash(msg)
    endif
    if (no_fuel_cat >= 1 .and. no_fuel_cat <= nfuelcats)then
        write(msg,*)'no_fuel_cat=',no_fuel_cat,' may not be between 1 and nfuelcats=',nfuelcats
        call crash(msg)
    endif
    if (no_fuel_cat > no_fuel_cat2)then
        write(msg,*)'no_fuel_cat=',no_fuel_cat,' must not be larger than no_fuel_cat2=',no_fuel_cat2
        call crash(msg)
    endif
ENDIF

! broadcast fuel tables
!$OMP CRITICAL(SFIRE_MODEL_CRIT)
write(*,*) "BCAST fuel data"
!$OMP END CRITICAL(SFIRE_MODEL_CRIT)

call wrf_dm_bcast_real(cmbcnst,1)
call wrf_dm_bcast_real(hfgl,1)
call wrf_dm_bcast_real(fuelmc_g,1)
call wrf_dm_bcast_real(fuelmc_c,1)
call wrf_dm_bcast_real(fire_wind_height,1)
call wrf_dm_bcast_real(ffwidth_thresh,1)
call wrf_dm_bcast_integer(nfuelcats,1)
call wrf_dm_bcast_integer(no_fuel_cat,1)
call wrf_dm_bcast_integer(no_fuel_cat2,1)
call wrf_dm_bcast_integer(ibeh,1)
call wrf_dm_bcast_real(windrf,    nfuelcats)
call wrf_dm_bcast_real(fgi,       nfuelcats)
call wrf_dm_bcast_real(fueldepthm,nfuelcats)
call wrf_dm_bcast_real(savr,      nfuelcats)
call wrf_dm_bcast_real(fuelmce,   nfuelcats)
call wrf_dm_bcast_real(fueldens,  nfuelcats)
call wrf_dm_bcast_real(st,        nfuelcats)
call wrf_dm_bcast_real(se,        nfuelcats)
call wrf_dm_bcast_real(weight,    nfuelcats)
call wrf_dm_bcast_real(fci_d,     nfuelcats)
call wrf_dm_bcast_real(fct,       nfuelcats)
call wrf_dm_bcast_integer(ichap,  nfuelcats)
call wrf_dm_bcast_integer(itree,  nfuelcats)
call wrf_dm_bcast_integer(idynamic,  nfuelcats)
call wrf_dm_bcast_real(savr_gc,     mfuelcats*max_moisture_classes)
call wrf_dm_bcast_real(fuelload_gc,     mfuelcats*max_moisture_classes)
call wrf_dm_bcast_real(fcwh,      nfuelcats)
call wrf_dm_bcast_real(fcz0,      nfuelcats)
call wrf_dm_bcast_real(ffw,       nfuelcats)
call wrf_dm_bcast_real(adjr0,     nfuelcats)
call wrf_dm_bcast_real(adjrw,     nfuelcats)
call wrf_dm_bcast_real(adjrs,     nfuelcats)
call wrf_dm_bcast_real(fmc_gl_ndwi_0,    nfuelcats)
call wrf_dm_bcast_real(fmc_gl_ndwi_rate, nfuelcats)
call wrf_dm_bcast_real(fmc_gl_ndwi_stdev,nfuelcats)
call wrf_dm_bcast_real(fmc_gl_stdev,     nfuelcats)
call wrf_dm_bcast_real(k_tc,   nfuelcats)
CALL wrf_dm_bcast_real(tree_z_max,   n_tree_type)
CALL wrf_dm_bcast_real(tree_z_min,   n_tree_type)
CALL wrf_dm_bcast_real(tree_burn_time,   n_tree_type)
! broadcast moisture tables
call wrf_dm_bcast_integer(moisture_classes,1)
call wrf_dm_bcast_real(drying_lag,     max_moisture_classes)
call wrf_dm_bcast_real(wetting_lag,     max_moisture_classes)
call wrf_dm_bcast_real(saturation_moisture,     max_moisture_classes)
call wrf_dm_bcast_real(saturation_rain,     max_moisture_classes)
call wrf_dm_bcast_real(rain_threshold,     max_moisture_classes)
call wrf_dm_bcast_integer(drying_model,     max_moisture_classes)
call wrf_dm_bcast_integer(wetting_model,     max_moisture_classes)
call wrf_dm_bcast_integer(fmc_gc_initialization,     max_moisture_classes)
call wrf_dm_bcast_integer(is_dead,     max_moisture_classes)
call wrf_dm_bcast_real(fmc_gw,     mfuelcats*max_moisture_classes)
call wrf_dm_bcast_real(fuel_load_layers_hf, n_fuel_cat_A13*n_layer_hf_surf)
call wrf_dm_bcast_real(layer_burn_time_b, n_layer_hf_surf)
call wrf_dm_bcast_real(layer_burn_time_c, n_layer_hf_surf)
call wrf_dm_bcast_real(layer_char_fraction, n_layer_hf_surf)

! moisture model derived scalars
do i=1,moisture_classes
  rec_drying_lag_sec(i)  = 1.0/(3600.0*drying_lag(i))
  rec_wetting_lag_sec(i) = 1.0/(3600.0*wetting_lag(i))
enddo

!-------------------------------- fuel model
! compute derived scalars

fuelheat = cmbcnst * 4.30e-04     ! convert J/kg to BTU/lb

! compute derived fuel category coefficients 

DO i = 1,nfuelcats
    fci(i) = (1.+fuelmc_c)*fci_d(i)
    if(fct(i) .ne.  0.)then
        fcbr(i) = fci_d(i)/fct(i) ! to avoid division by zero
    else
        fcbr(i) = 0
    endif
END DO

! prints

call message('**********************************************************')
call message('FUEL COEFFICIENTS')
write(msg,8)'cmbcnst    ',cmbcnst
call message(msg)
write(msg,8)'hfgl       ',hfgl
call message(msg)
write(msg,8)'fuelmc_g   ',fuelmc_g
call message(msg)
write(msg,8)'fuelmc_c   ',fuelmc_c
call message(msg)
write(msg,8)'fuelheat   ',fuelheat
call message(msg)
write(msg,7)'nfuelcats  ',nfuelcats
call message(msg)
write(msg,7)'no_fuel_cat',no_fuel_cat
call message(msg)
write(msg,7)'no_fuel_cat2',no_fuel_cat2
call message(msg)
if(init_fuel_moisture)then
    write(msg,7)'moisture_classes',moisture_classes
    call message(msg)
endif

j=1
7 format(a,5(i9,4x))
8 format(a,5(1x,g12.5e2))
9 format(a,5(1x,a))
10 format(a,i2.2,2x,5(1x,g12.5e2))
do i=1,nfuelcats,j
    k=min(i+j-1,nfuelcats)
    call message(' ')
    write(msg,7)'CATEGORY  ',(ii,ii=i,k)
    call message(msg)
    write(msg,9)'fuel name ',(trim(fuel_name(ii)),ii=i,k)
    call message(msg)
    write(msg,8)'fwh       ',(fcwh(ii),ii=i,k)
    call message(msg)
    write(msg,8)'fz0       ',(fcz0(ii),ii=i,k)
    call message(msg)
    write(msg,8)'windrf    ',(windrf(ii),ii=i,k)
    call message(msg)
    write(msg,8)'fgi       ',(fgi(ii),ii=i,k)
    call message(msg)
    write(msg,8)'fueldepthm',(fueldepthm(ii),ii=i,k)
    call message(msg)
    write(msg,8)'savr      ',(savr(ii),ii=i,k)
    call message(msg)
    write(msg,8)'fuelmce   ',(fuelmce(ii),ii=i,k)
    call message(msg)
    write(msg,8)'fueldens  ',(fueldens(ii),ii=i,k)
    call message(msg)
    write(msg,8)'st        ',(st(ii),ii=i,k)
    call message(msg)
    write(msg,8)'se        ',(se(ii),ii=i,k)
    call message(msg)
    write(msg,8)'weight    ',(weight(ii),ii=i,k)
    call message(msg)
    write(msg,8)'fci_d     ',(fci_d(ii),ii=i,k)
    call message(msg)
    write(msg,8)'fct       ',(fct(ii),ii=i,k)
    call message(msg)
    write(msg,7)'ichap     ',(ichap(ii),ii=i,k)
    call message(msg)
    write(msg,7)'itree     ',(itree(ii),ii=i,k)
    call message(msg)
    write(msg,7)'idynamic  ',(idynamic(ii),ii=i,k)
    call message(msg)
    write(msg,8)'fci       ',(fci(ii),ii=i,k)
    call message(msg)
    write(msg,8)'fcbr      ',(fcbr(ii),ii=i,k)
    call message(msg)
    write(msg,8)'ffw       ',(ffw(ii),ii=i,k)
    call message(msg)
    write(msg,8)'adjr0     ',(adjr0(ii),ii=i,k)
    call message(msg)
    write(msg,8)'adjrw     ',(adjrw(ii),ii=i,k)
    call message(msg)
    write(msg,8)'adjrs     ',(adjrs(ii),ii=i,k)
    call message(msg)
    do kk=1,moisture_classes 
        write(msg,10)'savr_gc',kk,(savr_gc(ii,kk),ii=i,k)
        call message(msg)
        write(msg,10)'fuelload_gc',kk,(fuelload_gc(ii,kk),ii=i,k)
        call message(msg)
    enddo
    if(init_fuel_moisture)then
        do kk=1,moisture_classes 
            write(msg,10)'fmc_gw',kk,(fmc_gw(ii,kk),ii=i,k)
            call message(msg)
        enddo
    endif
    if(kfmc_ndwi>0)then
        write(msg,8)'fmc_gl_stdev     ',(fmc_gl_stdev(ii),ii=i,k)
        call message(msg)
        write(msg,8)'fmc_gl_ndwi_0    ',(fmc_gl_ndwi_0(ii),ii=i,k)
        call message(msg)
        write(msg,8)'fmc_gl_ndwi_rate ',(fmc_gl_ndwi_rate(ii),ii=i,k)
        call message(msg)
        write(msg,8)'fmc_gl_ndwi_stdev',(fmc_gl_ndwi_stdev(ii),ii=i,k)
        call message(msg)
    endif
    write(msg,8)'k_tc   ',(k_tc(ii),ii=i,k)
    call message(msg)
enddo
call message(' ')
call message('**********************************************************')

if(init_fuel_moisture)then
    j=1
    do i=1,moisture_classes,j
        k=min(i+j-1,nfuelcats)
        call message(' ')
        write(msg,7)'FUEL MOISTURE CLASS',(ii,ii=i,k)
        call message(msg)
        write(msg,9)'moisture class name    ',(trim(moisture_class_name(ii)),ii=i,k)
        call message(msg)
        write(msg,7)'drying_model           ',(drying_model(ii),ii=i,k)
        call message(msg)
        write(msg,8)'drying_lag (h)         ',(drying_lag(ii),ii=i,k)
        call message(msg)
        write(msg,7)'wetting_model          ',(wetting_model(ii),ii=i,k)
        call message(msg)
        write(msg,7)'fmc_gc_initialization  ',(fmc_gc_initialization(ii),ii=i,k)
        call message(msg)
        write(msg,8)'wetting_lag (h)        ',(wetting_lag(ii),ii=i,k)
        call message(msg)    
        write(msg,8)'saturation_moisture (1)',(saturation_moisture(ii),ii=i,k)
        call message(msg)    
        write(msg,8)'saturation_rain (mm/h) ',(saturation_rain(ii),ii=i,k)
        call message(msg)    
        write(msg,8)'rain_threshold (mm/h)  ',(rain_threshold(ii),ii=i,k)
        call message(msg)    
        write(msg,7)'is_dead                ',(is_dead(ii),ii=i,k)
        call message(msg)
    enddo
    call message(' ')
    call message('**********************************************************')
    call message(' ')
endif
have_fuel_cats=.true.

! and print to file
IF ( wrf_dm_on_monitor() .and. call_write_fuels_m) THEN
  call write_fuels_m(61,30.,1.)
ENDIF


end subroutine init_fuel_cats


subroutine write_fuels_m(nsteps,maxwind,maxslope)
implicit none
integer, intent(in):: nsteps   ! number of steps for speed computation
real, intent(in):: maxwind,maxslope ! from zero to these

integer:: iounit,k,j,i,m,isave
type(fire_params)::fp
real, dimension(1:3,1:nsteps), target::vx,vy,zsf,dzdxf,dzdyf,bbb,phisc,phiwc,r_0,fgip,ischap,fmc_g,wind,nfuel_cat,rhof,tempf,ffwidth, cfb
real, dimension(1:3,max_moisture_classes,1:nsteps), target::fmc_gc_f
real, dimension(1:3,1:nsteps)::fuel_time,ros,fwh,fz0
real::ros_back,ros_wind,ros_slope,propx,propy,r,speed,tanphi
integer::ierrx
character(len=128)::msg

if(.not.have_fuel_cats)call crash('write_fuels_m: fuel categories not yet set')

fp%vx=>vx
fp%vy=>vy
fp%dzdxf=>dzdxf
fp%dzdyf=>dzdyf
fp%bbb=>bbb
fp%phisc=>phisc
fp%phiwc=>phiwc
fp%r_0=>r_0
fp%fgip=>fgip
fp%ischap=>ischap
fp%fmc_g=>fmc_g
fp%nfuel_cat=>nfuel_cat
fp%rhof=>rhof
fp%tempf=>tempf
fp%ffwidth=>ffwidth
fp%fmc_gc_f=>fmc_gc_f

rhof = 1.
tempf = 300.
ffwidth = 100.

iounit = open_text_file('fuels.m','write')


10 format('fuel(',i3,').',a,'=',"'",a,"'",';% ',a)
do k=1,nfuelcats
    write(iounit,10)k,'fuel_name',trim(fuel_name(k)),'FUEL MODEL NAME'
    call write_var(k,'ibeh',1.*ibeh,'Model: 0 = MacArthur, 1 = Rothermel, 2 = Balbi' )
    call write_var(k,'windrf',windrf(k),'WIND REDUCTION FACTOR FROM FCWH TO MIDFLAME HEIGHT' )
    call write_var(k,'fwh',fcwh(k),'WIND HEIGHT TO INTERPOLATE VERTICALLY TO (M)' )
    call write_var(k,'fz0',fcz0(k),'ROUGHNESS LENGTH FOR VERTICAL WIND LOG INTERPOLATION (M)' )
    call write_var(k,'fgi',fgi(k),'INITIAL TOTAL MASS OF SURFACE FUEL (KG/M**2)' )
    call write_var(k,'fueldepthm',fueldepthm(k),'FUEL DEPTH (M)')
    call write_var(k,'savr',savr(k),'FUEL PARTICLE SURFACE-AREA-TO-VOLUME RATIO, 1/FT')
    call write_var(k,'fuelmce',fuelmce(k),'MOISTURE CONTENT OF EXTINCTION')
    call write_var(k,'fueldens',fueldens(k),'OVENDRY PARTICLE DENSITY, LB/FT^3')
    call write_var(k,'st',st(k),'FUEL PARTICLE TOTAL MINERAL CONTENT')
    call write_var(k,'se',se(k),'FUEL PARTICLE EFFECTIVE MINERAL CONTENT')
    call write_var(k,'weight',weight(k),'WEIGHTING PARAMETER THAT DETERMINES THE SLOPE OF THE MASS LOSS CURVE')
    call write_var(k,'fci_d',fci_d(k),'INITIAL DRY MASS OF CANOPY FUEL')
    call write_var(k,'fct',fct(k),'BURN OUT TIME FOR CANOPY FUEL, AFTER DRY (S)')
    call write_var(k,'ichap',float(ichap(k)),'1 if chaparral, 0 if not')
    call write_var(k,'idynamic',float(idynamic(k)),'1 if dynamic, 0 if not')
    call write_var(k,'fci',fci(k),'INITIAL TOTAL MASS OF CANOPY FUEL')
    call write_var(k,'fcbr',fcbr(k),'FUEL CANOPY BURN RATE (KG/M**2/S)')
    call write_var(k,'adjr0',adjr0(k),'MULTIPLICATIVE ADJUSTMENT OF BACKING SPREAD RATE')
    call write_var(k,'adjrw',adjrw(k),'MULTIPLICATIVE ADJUSTMENT OF WIND CONTRIBUTION TO SPREAD RATE')
    call write_var(k,'adjrs',adjrs(k),'MULTIPLICATIVE ADJUSTMENT OF SLOPE CONTRIBUTION TO SPREAD RATE')
    call write_var(k,'ffw',ffw(k),'FUEL FRACTION CONSUMED IN THE FLAMING ZONE')
    call write_var(k,'hfgl',hfgl,'SURFACE FIRE HEAT FLUX THRESHOLD TO IGNITE CANOPY (W/m^2)')
    call write_var(k,'cmbcnst',cmbcnst,'JOULES PER KG OF DRY FUEL')
    call write_var(k,'fuelheat',fuelheat,'FUEL PARTICLE LOW HEAT CONTENT, BTU/LB')
    call write_var(k,'fuelmc_g',fuelmc_g,'FUEL PARTICLE (SURFACE) MOISTURE CONTENT')
    call write_var(k,'fuelmc_c',fuelmc_c,'FUEL PARTICLE (CANOPY) MOISTURE CONTENT')
    call write_var(k,'k_tc',k_tc(k),'FUEL THERMAL CONDUCTIVITY (W /m /K)')
    call write_var(k,'savr_1',savr_gc01(k),'FUEL PARTICLE SURFACE-AREA-TO-VOLUME RATIO FROM DEAD 1H FUEL (1/FT)')
    call write_var(k,'savr_10',savr_gc02(k),'FUEL PARTICLE SURFACE-AREA-TO-VOLUME RATIO FROM DEAD 10H FUEL (1/FT)')
    call write_var(k,'savr_100',savr_gc03(k),'FUEL PARTICLE SURFACE-AREA-TO-VOLUME RATIO FROM DEAD 100H FUEL (1/FT)')
    call write_var(k,'savr_1000',savr_gc04(k),'FUEL PARTICLE SURFACE-AREA-TO-VOLUME RATIO FROM DEAD 1000H FUEL (1/FT)')
    call write_var(k,'savr_herba',savr_gc05(k),'FUEL PARTICLE SURFACE-AREA-TO-VOLUME RATIO FROM LIVE HERBACEOUS FUEL (1/FT)')
    call write_var(k,'savr_woody',savr_gc06(k),'FUEL PARTICLE SURFACE-AREA-TO-VOLUME RATIO FROM LIVE WOODY FUEL (1/FT)')
    call write_var(k,'fuelload_1',fuelload_gc01(k),'OVEN-DRY FUEL LOAD FROM DEAD 1H FUEL (LB/FT^2)')
    call write_var(k,'fuelload_10',fuelload_gc02(k),'OVEN-DRY FUEL LOAD FROM DEAD 10H FUEL (LB/FT^2)')
    call write_var(k,'fuelload_100',fuelload_gc03(k),'OVEN-DRY FUEL LOAD FROM DEAD 100H FUEL (LB/FT^2)')
    call write_var(k,'fuelload_1000',fuelload_gc04(k),'OVEN-DRY FUEL LOAD FROM DEAD 1000H FUEL (LB/FT^2)')
    call write_var(k,'fuelload_herba',fuelload_gc05(k),'OVEN-DRY FUEL LOAD FROM LIVE HERBACEOUS FUEL (LB/FT^2)')
    call write_var(k,'fuelload_woody',fuelload_gc06(k),'OVEN-DRY FUEL LOAD FROM LIVE WOODY FUEL (LB/FT^2)')
    ! set up fuel arrays 
    !subroutine set_fire_params( &
    !                       ifds,ifde,jfds,jfde, &
    !                       ifms,ifme,jfms,jfme, &
    !                       ifts,ifte,jfts,jfte, &
    !                       fdx,fdy,nfuel_cat0,  &
    !                       nfuel_cat,fuel_time, &
    !                       fp ) 
    nfuel_cat = k
    do j=1,nsteps ! set moisture - must be before set_fire_params
       fmc_g(1,j)=fuelmc_g
       fmc_g(2,j)=fuelmc_g
       fmc_g(3,j)=(fuelmce(k)*(j-1))/(nsteps-2)
       do m=1,max_moisture_classes
           fmc_gc_f(1,m,j)=fuelmc_g
           fmc_gc_f(2,m,j)=fuelmc_g
           fmc_gc_f(3,m,j)=(fuelmce(k)*(j-1))/(nsteps-2)
       enddo
    enddo
    isave=fire_fmc_read
    fire_fmc_read=0
    call set_fire_params( &
                           1,3,1,nsteps, &
                           1,3,1,nsteps, &
                           1,3,1,nsteps, &
                           0.,0.,k,  &
                           nfuel_cat,fuel_time, &
                           fp ) 
    fire_fmc_read=isave
    ! set up windspeed slope moisture table
    propx=1.
    propy=0.
    do j=1,nsteps
       r=float(j-1)/float(nsteps-1)
       ! line 1 varies windspeed (in x direction), zero slope
       wind(1,j)=maxwind*r
       vx(1,j)=wind(1,j)*windrf(k)
       vy(1,j)=0.
       dzdxf(1,j)=0.
       dzdyf(1,j)=0.
       ! line 2 varies slope (in x direction), zero slope
       vx(2,j)=0.
       vy(2,j)=0.
       dzdxf(2,j)=maxslope*r
       dzdyf(2,j)=0.
       ! line 3 varies moisture, zero slope and wind
       vx(3,j)=0.
       vy(3,j)=0.
       dzdxf(3,j)=0.
       dzdyf(3,j)=0.
    enddo
    do j=1,nsteps
       do i=1,3
          call fire_ros(ros(i,j), &
                 propx,propy,i,j,fp,ierrx,msg)
       enddo
       write(iounit,13)k,'wind',j,wind(1,j),'wind speed at 6.1m'
       write(iounit,13)k,'ros_wind',j,ros(1,j),'rate of spread for the wind speed at 6.1m'
       write(iounit,13)k,'slope',j,dzdxf(2,j),'slope'
       write(iounit,13)k,'ros_slope',j,ros(2,j),'rate of spread for the slope'
       write(iounit,13)k,'fmc_g',j,fmc_g(3,j),'fuel moisture content'
       write(iounit,13)k,'ros_fmc_g',j,ros(3,j),'rate of spread for the fuel moisture content'
    enddo
enddo

do k = 1, n_tree_type
    call write_var(k,'tree_z_min',tree_z_min(itree(k)),'CANOPY BOTTOM HEIGHT (m)')
    call write_var(k,'tree_z_max',tree_z_max(itree(k)),'CANOPY TOP HEIGHT (m)')
    call write_var(k,'tree_burn_time',tree_burn_time(itree(k)),'CANOPY TOP HEIGHT (m)')
enddo

13 format('fuel(',i3,').',a,'(',i3,')=',g12.5e2,';% ',a)
 
close(iounit)
! stop

contains

subroutine write_var(k,name,value,descr)
! write entry for one variable
integer, intent(in)::k
character(len=*), intent(in)::name,descr
real, intent(in)::value
write(iounit,11)k,name,value
write(iounit,12)k,name,descr
11 format('fuel(',i3,').',a,'=',g12.5e2,  ';')
12 format('fuel(',i3,').',a,"_descr='",a,"';")
end subroutine write_var

end subroutine write_fuels_m

!
!*******************
!

subroutine set_fire_params( &
                           ifds,ifde,jfds,jfde, &
                           ifms,ifme,jfms,jfme, &
                           ifts,ifte,jfts,jfte, &
                           fdx,fdy,nfuel_cat0,  &
                           nfuel_cat,fuel_time, &
                           fp ) 

implicit none

!*** purpose: Set all fire model params arrays, constant values.

!*** arguments
integer, intent(in)::ifds,ifde,jfds,jfde                        ! fire domain bounds
integer, intent(in)::ifts,ifte,jfts,jfte                        ! fire tile bounds
integer, intent(in)::ifms,ifme,jfms,jfme                        ! memory array bounds
real, intent(in):: fdx,fdy                                      ! fire mesh spacing
integer,intent(in)::nfuel_cat0                                  ! default fuel category, if nfuel_cat=0
real, intent(inout),dimension(ifms:ifme, jfms:jfme)::nfuel_cat  ! fuel data
real, intent(out), dimension(ifms:ifme, jfms:jfme)::fuel_time   ! fire params arrays
type(fire_params),intent(inout)::fp

!*** local

real::  fuelload, fueldepth, rtemp1, rtemp2, &
        qig, epsilon, rhob, wn, betaop, e, c, &
        xifr, etas, etam, a, gammax, gamma, ratio, ir, &
        fuelloadm,fdxinv,fdyinv,betafl, bmst
integer:: i,j,k
integer::nerr
character(len=128)::msg

!*** executable

if(.not.have_fuel_cats)call crash('set_fire_params: fuel categories not yet set')

call print_2d_stats(ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme,fp%fmc_g,'set_fire_params: fmc_g')

nerr=0
do j=jfts,jfte
   do i=ifts,ifte
     ! fuel category 
     k=int( nfuel_cat(i,j) )
     if(k <= 0)then   ! no fuel 
        fp%fgip(i,j)=0.            ! no mass 
        fp%ischap(i,j)=0.
        fp%phisc(i,j)=0.          ! 
        fp%bbb(i,j)=1.             !
        fuel_time(i,j)=7./0.85  ! does not matter, just what was there before
        fp%phiwc(i,j)=0.
        fp%r_0(i,j)=0.             ! no fuel, no spread.
     else
        ! if(k.eq.0.and.nfuel_cat0.ge.1.and.nfuel_cat0.le.nfuelcats)then
        !     ! replace k=0 by default
        !     k=nfuel_cat0
        !     nerr=nerr+1
        ! endif
   
        ! crash if fuel_cat >= no_fuel_cat (ie max_fuel_cat)
        if(k > nfuelcats)then
!$OMP CRITICAL(SFIRE_PHYS_CRIT)
            write(msg,'(3(a,i5))')'nfuel_cat(', i ,',',j,')=',k
            call message(msg,level=0)
            if(k.eq.0)then
                call message('Possibly nfuel_cat is uninitialized on input',level=0)
            endif
            write(msg,'(a,i5)')'nfuel_cat must between 1 and nfuel_cat=',nfuelcats
            call message(msg,level=0)
            write(msg,'(a,i5,a,i10,a)')'or set it to 0 for no fuel'
            call message(msg,level=0)
!$OMP END CRITICAL(SFIRE_PHYS_CRIT)
            call crash('set_fire_params: fuel category out of bounds')
        endif
        

        fuel_time(i,j)=weight(k)/0.85 ! cell based
        
        ! do not understand calculations of stime in binit.m4
        ! set fuel time constant: weight=1000=>40% decrease over 10 min
        ! fuel decreases as exp(-t/fuel_time) 
        ! exp(-600*0.85/1000) = approx 0.6 

        fp%ischap(i,j)=ichap(k)
        fp%fgip(i,j)=fgi(k)
        if(fire_fmc_read.eq.1)then
           fp%fmc_g(i,j)=fuelmc_g
        endif
        ! print *,'fmc_g:',fire_fmc_read,i,j,fp%fmc_g(i,j)

        ! end jm addition

        !
        !*** rest copied from wf2_janice/fire_startup.m4 with minimal changes
        !

        !     ...Settings of fire spread parameters from Rothermel follows. These
        !        don't need to be recalculated later.
        
        bmst     = fp%fmc_g(i,j) / (1.+fp%fmc_g(i,j))
        fuelloadm= (1.-bmst) * fgi(k)  !  fuelload without moisture
        fuelload = fuelloadm * (.3048)**2 * 2.205    ! to lb/ft^2
        fueldepth = fueldepthm(k)/0.3048               ! to ft
        betafl = fuelload/(fueldepth * fueldens(k))! packing ratio
        betaop = 3.348 * savr(k)**(-0.8189)     ! optimum packing ratio
        qig = 250. + 1116.*fp%fmc_g(i,j)            ! heat of preignition, btu/lb
        epsilon = exp(-138./savr(k) )    ! effective heating number
        rhob = fuelload/fueldepth    ! ovendry bulk density, lb/ft^3

        c = 7.47 * exp( -0.133 * savr(k)**0.55)    ! const in wind coef
        fp%bbb(i,j) = 0.02526 * savr(k)**0.54      ! const in wind coef
        !if(fire_wind_log_interp .eq. 4 .or. fire_use_windrf .eq. 1) then
        !     c = c * windrf(k)**fp%bbb(i,j)             ! jm: multiply wind by reduction factor
        !endif
        e = 0.715 * exp( -3.59e-4 * savr(k))       ! const in wind coef
        fp%phiwc(i,j) = c * (betafl/betaop)**(-e)

        ! phis = 5.275 *(fp%betafl(i,j))**(-0.3) *tanphim**2   ! slope factor
        fp%phisc(i,j) = 5.275 *(betafl)**(-0.3)    ! const in slope coeff

        rtemp2 = savr(k)**1.5
        gammax = rtemp2/(495. + 0.0594*rtemp2)              ! maximum rxn vel, 1/min
        a = 1./(4.774 * savr(k)**0.1 - 7.27)   ! coef for optimum rxn vel
        ratio = betafl/betaop
        gamma = gammax *(ratio**a) *exp(a*(1.-ratio)) !optimum rxn vel, 1/min

        wn = fuelload/(1 + st(k))       ! net fuel loading, lb/ft^2
        rtemp1 = fp%fmc_g(i,j)/fuelmce(k)
        etam = 1.-2.59*rtemp1 +5.11*rtemp1**2 -3.52*rtemp1**3  !moist damp coef
        etam = max(etam,0.0)
        etas = 0.174* se(k)**(-0.19)                ! mineral damping coef
        ir = gamma * wn * fuelheat * etam * etas  !rxn intensity,btu/ft^2 min
        ! jm irm = ir * 1055./( 0.3048**2 * 60.) * 1.e-6     !for mw/m^2
        ! jm: irm set but never used??

        xifr = exp( (0.792 + 0.681*savr(k)**0.5) &
            * (betafl+0.1)) /(192. + 0.2595*savr(k)) ! propagating flux ratio

!        ... r_0 is the spread rate for a fire on flat ground with no wind.

        fp%r_0(i,j) = ir*xifr/(rhob * epsilon *qig)    ! default spread rate in ft/min
        fp%r_0(i,j) = fp%r_0(i,j) * .00508             ! convert to m/s
        fp%phiwc(i,j) = fp%phiwc(i,j) * fp%r_0(i,j)              ! premultiply wind coefficient so it can be used additively
        fp%phisc(i,j) = fp%phisc(i,j) * fp%r_0(i,j)              ! premultiply wind coefficient so it can be used additively

        ! apply adjustments
        fp%r_0(i,j)   = fp%r_0(i,j)   * adjr0(k)
        fp%phiwc(i,j) = fp%phiwc(i,j) * adjrw(k)
        fp%phisc(i,j) = fp%phisc(i,j) * adjrs(k)

        ! test fmc
        if(fp%r_0(i,j) > 1e-6 .and. fp%fmc_g(i,j) > fuelmce(k))then
!$OMP CRITICAL(SFIRE_PHYS_CRIT)
             write(msg,'(a,2i6,3(a,e14.6))') 'set_fire_params: at ',i,j,' base rate of spread',fp%r_0(i,j),' moisture ',fp%fmc_g(i,j),'> extinction =',fuelmce(k)
             call message(msg,level=0)
             write(msg,'(5(a,e14.6))')'rtemp1=',rtemp1,' etam=',etam
             call message(msg,level=0)
!$OMP END CRITICAL(SFIRE_PHYS_CRIT)
             call warning('propagation above extinction moisture',level=0)
        endif
     endif
  enddo
enddo

if(nerr.gt.1)then
!$OMP CRITICAL(SFIRE_PHYS_CRIT)
    write(msg,'(a,i6,a)')'set_fire_params: WARNING: fuel category 0 replaced in',nerr,' cells'
!$OMP END CRITICAL(SFIRE_PHYS_CRIT)
    call message(msg)
endif

end subroutine set_fire_params

!
!*******************
!

!------------------------------------------------------------------------------------------------------------------
! SUBROUTINE: heat_fluxes
!>   Compute heat fluxes on the fire grid
!!
!!  Compute ground and canopy heat fluxes  
!------------------------------------------------------------------------------------------------------------------
subroutine heat_fluxes(dt,fp,               &
  ifms, ifme, jfms, jfme,                   &
  ifts, ifte, jfts, jfte,                   &
  iffs, iffe, jffs, jffe,                   &
  hf_scheme_selection,                      &
  fgip, fuel_frac_burnt,                    & 
  grnhft, grnqft,                           &
  fire_area, canhft, canqft, can_fuel_frac  )

  use module_fr_sfire_fields, only: cfb, can_bulk_density, can_top, can_bottom, can_fuel_burnt

  implicit none

  ! 0. Declaration
  !    ===========
  
  type(fire_params),  intent(in   )  :: fp !< fire params
  real,               intent(in   )  :: dt !< the fire time step [s]

  integer,  intent(in   )  :: ifts  !< fire tile start index x axis
  integer,  intent(in   )  :: ifte  !< fire tile end   index x axis
  integer,  intent(in   )  :: jfts  !< fire tile start index y axis
  integer,  intent(in   )  :: jfte  !< fire tile end   index y axis

  integer,  intent(in   )  :: ifms  !< fire mesh start index x axis
  integer,  intent(in   )  :: ifme  !< fire mesh end   index x axis
  integer,  intent(in   )  :: jfms  !< fire mesh start index y axis
  integer,  intent(in   )  :: jfme  !< fire mesh end   index y axis

  integer,  intent(in   )  :: iffs  !< fire fuel frac burnt start index x axis
  integer,  intent(in   )  :: iffe  !< fire fuel frac burnt end   index x axis
  integer,  intent(in   )  :: jffs  !< fire fuel frac burnt start index y axis
  integer,  intent(in   )  :: jffe  !< fire fuel frac burnt end   index y axis

  integer,  intent(in   )  :: hf_scheme_selection   !< Choice of HF scheme (0: BURNUP, 1: Multilayers BioSmoke, 2: FOFEM)
  
  real, intent(in   ),  dimension(ifms:ifme,jfms:jfme)  :: fgip             !< init mass from fuel model [kg m-2]
  real, intent(in   ),  dimension(iffs:iffe,jffs:jffe)  :: fuel_frac_burnt  !< fraction burned this call [-]
  
  real, intent(  out),  dimension(ifms:ifme, jfms:jfme)             :: grnhft !< surface sensible heat flux [W m-2]   
  real, intent(  out),  dimension(ifms:ifme, jfms:jfme),  optional  :: grnqft !< latent sensible heat flux [W m-2]
  
  real, intent(in   ),  dimension(ifms:ifme, jfms:jfme),  optional    :: fire_area      !< proportion of the fire cell burning [-]
  real, intent(  out),  dimension(ifms:ifme, jfms:jfme),  optional    :: canhft         !< canopy sensible heat flux [W m-2]
  real, intent(  out),  dimension(ifms:ifme, jfms:jfme),  optional    :: canqft         !< canopy latent heat flux [W m-2]
  real, intent(inout),  dimension(ifms:ifme, jfms:jfme),  optional    :: can_fuel_frac  !< canopy fuel fraction [-]
  
  !*** local
  integer :: i, j
  real    :: dmass, bmst
  
  integer :: i_tree_cat !< tree class
  
  logical :: latent, canopy_sensible, canopy_latent
  
  real    :: can_fuel_cons      !< Fuel consumption for canopy [s-1]
  real    :: mlr_canopy         !< Mass Loss Rate for canopy [s-1]
  real    :: can_load_loss_rate !< Load loss rate for canopy [kg m-2 s-1]
  real    :: loc_fire_area      !< local fire area if present [-]
  
  ! Compute surface heat fluxes
  call heat_fluxes_surf(dt, fp, &
  ifms, ifme, jfms, jfme,     &
  ifts, ifte, jfts, jfte,     &
  iffs, iffe, jffs, jffe,     &
  hf_scheme_selection,        &
  fgip, fuel_frac_burnt,      & 
  grnhft, grnqft, fire_area   )
  
  ! Compute canopy heat fluxes
  ! A minima sensible heat fluxes
  if (present(fire_area) .and. present(canhft) .and. present(can_fuel_frac)) then
    call heat_fluxes_canopy(dt,fp,               &
    ifms, ifme, jfms, jfme,                   &
    ifts, ifte, jfts, jfte,                   &
    iffs, iffe, jffs, jffe,                   &
    fgip, fuel_frac_burnt,                    & 
    fire_area, canhft, can_fuel_frac, canqft  )
endif

end subroutine heat_fluxes


!------------------------------------------------------------------------------------------------------------------
! SUBROUTINE: heat_fluxes_surf
!>   Compute heat fluxes on the fire grid for surface only
!!
!!  Compute ground sensible and latent heat fluxes  
!------------------------------------------------------------------------------------------------------------------
subroutine heat_fluxes_surf(dt,fp,          &
    ifms, ifme, jfms, jfme,                   &
    ifts, ifte, jfts, jfte,                   &
    iffs, iffe, jffs, jffe,                   &
    hf_scheme_selection,                      &
    fgip, fuel_frac_burnt,                    & 
    grnhft, grnqft, fire_area                 )
    
    use module_fr_sfire_fields, only: surf_mass_frac_b, surf_mass_frac_c, surf_moisture_frac
    
    implicit none
    
    ! 0. Declaration
    !    ===========
    
    type(fire_params),  intent(in   )  :: fp !< fire params
    real,               intent(in   )  :: dt !< the fire time step [s]
    
    integer,  intent(in   )  :: ifts  !< fire tile start index x axis
    integer,  intent(in   )  :: ifte  !< fire tile end   index x axis
    integer,  intent(in   )  :: jfts  !< fire tile start index y axis
    integer,  intent(in   )  :: jfte  !< fire tile end   index y axis
    
    integer,  intent(in   )  :: ifms  !< fire mesh start index x axis
    integer,  intent(in   )  :: ifme  !< fire mesh end   index x axis
    integer,  intent(in   )  :: jfms  !< fire mesh start index y axis
    integer,  intent(in   )  :: jfme  !< fire mesh end   index y axis
    
    integer,  intent(in   )  :: iffs  !< fire fuel frac burnt start index x axis
    integer,  intent(in   )  :: iffe  !< fire fuel frac burnt end   index x axis
    integer,  intent(in   )  :: jffs  !< fire fuel frac burnt start index y axis
    integer,  intent(in   )  :: jffe  !< fire fuel frac burnt end   index y axis
    
    integer,  intent(in   )  :: hf_scheme_selection   !< Choice of HF scheme (0: BURNUP, 1: Multilayers BioSmoke, 2: FOFEM)

    real, intent(in   ),  dimension(ifms:ifme,jfms:jfme)  :: fgip             !< init mass from fuel model [kg m-2]
    real, intent(in   ),  dimension(iffs:iffe,jffs:jffe)  :: fuel_frac_burnt  !< fraction burned this call [-]
    
    real, intent(  out),  dimension(ifms:ifme, jfms:jfme)             :: grnhft     !< surface sensible heat flux [W m-2]   
    real, intent(  out),  dimension(ifms:ifme, jfms:jfme),  optional  :: grnqft     !< latent sensible heat flux [W m-2]
  real, intent(in   ),  dimension(ifms:ifme, jfms:jfme),  optional  :: fire_area  !< proportion of the fire cell burning [-]
  
  !*** local
  integer :: i, j, k, l
  real    :: dmass, bmst
  
  logical :: latent
  real, parameter :: water_from_combustion_prop_b = 0.56 !< water fraction produced by combustion from biomass [-]
  real, parameter :: water_from_combustion_prop_c = 0.56 !< water fraction produced by combustion from char[-]
  
  ! multi layer constants
  real, parameter :: delta_H_comb_b = 17.433e+06  !< biomass combustion enthalpy [J kg-1]
  real, parameter :: delta_H_comb_c = 17.433e+06  !< char combustion enthalpy [J kg-1]
  real, parameter :: moisture_evap_rate = 0.05    !< evaporation rate [s-1]
  
  ! multi layer local variables
  real :: mflr_b      !< biomass mass fraction loss rate [s-1]
  real :: mflr_c      !< char mass fraction loss rate [s-1]
  real :: moflr       !< moisture fraction loss rate [s-1]
  real :: tmp_flux_h  !< local sensible heat flux [W m-2]
  real :: tmp_flux_w  !< local latent heat flux [W m-2]

  !*** executable        
  latent = present(grnqft)
  select case(hf_scheme_selection)
  case(0) ! BURNUP
    ! default exeponential scheme
    do j = jfts, jfte
      do i = ifts, ifte

        ! Surface fluxes
        dmass = fgip(i, j) * fuel_frac_burnt(i, j) 
        bmst  = fp%fmc_g(i, j) / (1. + fp%fmc_g(i, j))
        
        ! Surface grounf heat flux
        grnhft(i, j) = (dmass / dt) * (1. - bmst) * cmbcnst
        
        if(latent) grnqft(i, j) = (bmst + (1. - bmst) * .56) * (dmass / dt) * xlv
        ! bmst = relative water contents; 0.56 = est. ratio of water from burning
        ! xlv = nominal specific latent heat of water J/kg (dependence on temperature ignored)
        ! xlv is  defined in module_model_constants 
      enddo
    enddo
  case(1) ! MLHF
    ! Crash if fire area not in inputs
    if (.not. present(fire_area)) then
      call crash("heat_fluxes_surf: fire_area not in inputs for MLHF scheme. Check call of heat_fluxes")
    end if

    do j = jfts, jfte
      do i = ifts, ifte

        ! get fuel class
        k = int(fp%nfuel_cat(i, j))

        ! init fluxes
        tmp_flux_h = 0.
        tmp_flux_w = 0.

        if (latent) then
          ! compute moisture evaporation rate (constant for now)
          moflr = min(fire_area(i, j) * moisture_evap_rate, surf_moisture_frac(i, j) / dt)
          
          ! update fuel fraction
          surf_moisture_frac(i, j) = surf_moisture_frac(i, j) - moflr * dt

          tmp_flux_w = fgip(i, j) * fp%fmc_g(i, j) * moflr * xlv
        endif

        ! loop over layers
        do l = 1, n_layer_hf_surf
          ! compute mass fraction loss rate
          mflr_b = min(fire_area(i, j) * surf_mass_frac_b(i, l, j) / layer_burn_time_b(l), surf_mass_frac_b(i, l, j) / dt)
          mflr_c = min(fire_area(i, j) * surf_mass_frac_c(i, l, j) / layer_burn_time_c(l), surf_mass_frac_c(i, l, j) / dt)

          ! update fuel fraction
          surf_mass_frac_b(i, l, j) = surf_mass_frac_b(i, l, j) - mflr_b * dt
          surf_mass_frac_c(i, l, j) = surf_mass_frac_c(i, l, j) - mflr_c * dt

          ! compute fluxes
          tmp_flux_h = tmp_flux_h + fgip(i, j) * (mflr_b * delta_H_comb_b + mflr_c * delta_H_comb_c)
          ! add fuel burning contribution to latent HF
          if(latent) tmp_flux_w = tmp_flux_w + fgip(i, j) * (mflr_b * water_from_combustion_prop_b + mflr_c * water_from_combustion_prop_c) * xlv
        enddo

        ! save fluxes
        grnhft(i, j) = tmp_flux_h
        if (latent) grnqft(i, j) = tmp_flux_w

      enddo
    enddo
  end select

end subroutine heat_fluxes_surf


!------------------------------------------------------------------------------------------------------------------
! SUBROUTINE: heat_fluxes_canopy
!>   Compute heat fluxes on the fire grid for canopy only
!!
!!  Compute canopy sensible and latent heat fluxes  
!------------------------------------------------------------------------------------------------------------------
subroutine heat_fluxes_canopy(dt,fp,        &
  ifms, ifme, jfms, jfme,                   &
  ifts, ifte, jfts, jfte,                   &
  iffs, iffe, jffs, jffe,                   &
  fgip, fuel_frac_burnt,                    & 
  fire_area, canhft, can_fuel_frac, canqft  )

  use module_fr_sfire_fields, only: cfb, can_bulk_density, can_top, can_bottom, can_fuel_burnt

  implicit none

  ! 0. Declaration
  !    ===========
  
  type(fire_params),  intent(in   )  :: fp !< fire params
  real,               intent(in   )  :: dt !< the fire time step [s]

  integer,  intent(in   )  :: ifts  !< fire tile start index x axis
  integer,  intent(in   )  :: ifte  !< fire tile end   index x axis
  integer,  intent(in   )  :: jfts  !< fire tile start index y axis
  integer,  intent(in   )  :: jfte  !< fire tile end   index y axis

  integer,  intent(in   )  :: ifms  !< fire mesh start index x axis
  integer,  intent(in   )  :: ifme  !< fire mesh end   index x axis
  integer,  intent(in   )  :: jfms  !< fire mesh start index y axis
  integer,  intent(in   )  :: jfme  !< fire mesh end   index y axis

  integer,  intent(in   )  :: iffs  !< fire fuel frac burnt start index x axis
  integer,  intent(in   )  :: iffe  !< fire fuel frac burnt end   index x axis
  integer,  intent(in   )  :: jffs  !< fire fuel frac burnt start index y axis
  integer,  intent(in   )  :: jffe  !< fire fuel frac burnt end   index y axis

  real, intent(in   ),  dimension(ifms:ifme,jfms:jfme)  :: fgip             !< init mass from fuel model [kg m-2]
  real, intent(in   ),  dimension(iffs:iffe,jffs:jffe)  :: fuel_frac_burnt  !< fraction burned this call [-]

  real, intent(in   ),  dimension(ifms:ifme, jfms:jfme)               :: fire_area      !< proportion of the fire cell burning [-]
  real, intent(  out),  dimension(ifms:ifme, jfms:jfme)               :: canhft         !< canopy sensible heat flux [W m-2]
  real, intent(inout),  dimension(ifms:ifme, jfms:jfme)               :: can_fuel_frac  !< canopy fuel fraction [-]
  real, intent(  out),  dimension(ifms:ifme, jfms:jfme),  optional    :: canqft         !< canopy latent heat flux [W m-2]

  !*** local
  integer :: i, j
  real    :: dmass, bmst

  integer :: i_tree_cat !< tree class

  logical :: canopy_latent

  real    :: can_fuel_cons      !< Fuel consumption for canopy [s-1]
  real    :: mlr_canopy         !< Mass Loss Rate for canopy [s-1]
  real    :: can_load_loss_rate !< Load loss rate for canopy [kg m-2 s-1]

  ! flag for canopy flux computation
  canopy_latent = present(canqft)

  do j = jfts, jfte
    do i = ifts, ifte
      ! get tree class
      i_tree_cat = itree(fp%nfuel_cat(i, j))
              
      ! compute mass loss rate for canopy
      if (i_tree_cat == 0) then
        can_fuel_cons = 0.
      else
        can_fuel_cons = can_fuel_frac(i, j) / tree_burn_time(i_tree_cat)
      endif

      mlr_canopy = min(fire_area(i, j) * cfb(i, j) * can_fuel_cons, can_fuel_frac(i, j) / dt)

      ! update canopy fuel fraction
      can_fuel_frac(i, j) = can_fuel_frac(i, j) - dt * mlr_canopy

      ! Canopy load loss rate (dot(sigma))
      can_load_loss_rate = can_bulk_density(i, j) * (can_top(i, j) - can_bottom(i, j)) * mlr_canopy

      ! save the load of fuel burnt 
      can_fuel_burnt(i, j) = dt * can_load_loss_rate

      ! compute canopy fluxes
      canhft(i, j) = cmbcnst * can_load_loss_rate
      if (canopy_latent)  canqft(i, j) = xlv * can_load_loss_rate * (.56 + fuelmc_c)
    enddo
  enddo

end subroutine heat_fluxes_canopy    


subroutine set_nfuel_cat(   &
    ifms,ifme,jfms,jfme,               &
    ifts,ifte,jfts,jfte,               &
    ifuelread,nfuel_cat0,zsf,nfuel_cat)

implicit none

! set fuel distributions for testing
integer, intent(in)::   ifts,ifte,jfts,jfte,               &
                        ifms,ifme,jfms,jfme               

integer, intent(in)::ifuelread,nfuel_cat0
real, intent(in), dimension(ifms:ifme, jfms:jfme)::zsf
real, intent(out), dimension(ifms:ifme, jfms:jfme)::nfuel_cat

!*** local

! parameters to control execution
integer:: i,j,iu1
real:: t1
character(len=128)msg

!$OMP CRITICAL(SFIRE_PHYS_CRIT)
    write(msg,'(a,i3)')'set_nfuel_cat: ifuelread=',ifuelread 
!$OMP END CRITICAL(SFIRE_PHYS_CRIT)
    call message(msg)

if (ifuelread .eq. -1 .or. ifuelread .eq. 2) then
!$OMP CRITICAL(SFIRE_PHYS_CRIT)
    call message('set_nfuel_cat: assuming nfuel_cat initialized already') 
    call message(msg)
!$OMP END CRITICAL(SFIRE_PHYS_CRIT)
else if (ifuelread .eq. 0) then
!
    do j=jfts,jfte
        do  i=ifts,ifte
            nfuel_cat(i,j)=real(nfuel_cat0)
        enddo
    enddo
!$OMP CRITICAL(SFIRE_PHYS_CRIT)
    write(msg,'(a,i3)')'set_nfuel_cat: fuel initialized with category',nfuel_cat0
!$OMP END CRITICAL(SFIRE_PHYS_CRIT)
    call message(msg)
         
else if (ifuelread .eq. 1) then
!
!         make dependent on altitude (co mountains/forest vs. plains)
!          2000 m : 6562 ft   ;    1600 m: 5249 ft

!        ... user defines fuel category spatial variability ! param!
    do j=jfts,jfte
        do  i=ifts,ifte
            ! nfuel_cat(i,j)= 2     ! grass with understory ! jm does nothing
            !jm t1=zsf(i,j)*slngth/100.
            t1 = zsf(i,j)  ! this is in m
            if(t1.le.1524.)then   !  up to 5000 ft
                nfuel_cat(i,j)= 3  ! tall grass
            else if(t1.ge.1524. .and. t1.le.2073.)then  ! 5.0-6.8 kft.
                nfuel_cat(i,j)= 2  ! grass with understory
            else if(t1.ge.2073..and.t1.le.2438.)then  ! 6.8-8.0 kft.
                nfuel_cat(i,j)= 8  ! timber litter - 10 (ponderosa)
            else if(t1.gt.2438. .and. t1.le. 3354.) then ! 8.0-11.0 kft.
!                 ... could also be mixed conifer.
                nfuel_cat(i,j)= 10 ! timber litter - 8 (lodgepole)
            else if(t1.gt.3354. .and. t1.le. 3658.) then ! 11.0-12.0 kft
                nfuel_cat(i,j)= 1  ! alpine meadow - 1
            else if(t1.gt.3658. ) then  ! > 12.0 kft
                nfuel_cat(i,j)= 14 ! no fuel.
            endif
        enddo
    enddo

    call message('set_nfuel_cat: fuel initialized by altitude')
else

    call crash('set_nfuel_cat: bad ifuelread')
endif
!     .............end  load fuel categories (or constant) here.

end subroutine set_nfuel_cat            


!------------------------------------------------------------------------------------------------------------------
! SUBROUTINE: set_canopy_data
!>   Fill canopy field variable for top height, bottom height and fuel load
!!
!!   The following fields are initialized by this routine (see description in module_fr_sfire_fields):
!!     - can_top
!!     - can_bottom
!!     - can_bulk_density
!!  
!!   Namelist input fire_can_fuel_read sets the source of information for setting the canopy data:
!!     - fire_can_fuel_read == -1 : default value meaning no initialization should be performed
!!     - fire_can_fuel_read == 0  : The canopy data are set according to surface fuel category.
!!     - fire_can_fuel_read == 1  : The canopy data are given by WRF input file and should not be modified.
!------------------------------------------------------------------------------------------------------------------
subroutine set_canopy_data(       &
  ifms, ifme, jfms, jfme,         &
  ifts, ifte, jfts, jfte,         &
  i_fire_can_fuel_read, nfuel_cat )

  use module_fr_sfire_fields, only: can_top, can_bottom, can_bulk_density

  implicit none

  ! 0. Declarations
  !    ------------

  ! 0.1 Arguments
  integer,  intent(in   ) :: ifms                 !< fire mesh index start on x axis
  integer,  intent(in   ) :: ifme                 !< fire mesh index end   on x axis
  integer,  intent(in   ) :: jfms                 !< fire mesh index start on y axis
  integer,  intent(in   ) :: jfme                 !< fire mesh index end   on y axis
  integer,  intent(in   ) :: ifts                 !< fire tile index start on x axis
  integer,  intent(in   ) :: ifte                 !< fire tile index end   on x axis
  integer,  intent(in   ) :: jfts                 !< fire tile index start on y axis
  integer,  intent(in   ) :: jfte                 !< fire tile index end   on y axis
  integer,  intent(in   ) :: i_fire_can_fuel_read !< source of data for canopy data

  real, dimension(ifms:ifme, jfms:jfme),  intent(in   )  :: nfuel_cat            !< surface fuel class

  ! 0.2 Locals
  integer :: surf_fuel_cat  !< local fuel category
  integer :: i, j           !< loop indexes

  

  ! 1. Fill the canopy data variable
  !    -----------------------------

  select case (i_fire_can_fuel_read)
  case (-1)
    ! The canopy data shouldbe as default value. Nothing to do
    can_top(:, :)           = 0.
    can_bottom(:, :)        = 0.
    can_bulk_density(:, :)  = 0.

  case (0)
    ! The canopy data are set according to surface fuel category.
    do j = jfts, jfte
      do  i = ifts, ifte
        ! get local fuel cat and convert to integer
        surf_fuel_cat = int(nfuel_cat(i, j))
        
        can_top(i, j)           = tree_z_max(itree(surf_fuel_cat))
        can_bottom(i, j)        = tree_z_min(itree(surf_fuel_cat))
        can_bulk_density(i, j)  = fci_d(surf_fuel_cat) / (can_top(i, j) - can_bottom(i, j))
      enddo
    enddo
    
  case (1)
    ! The canopy data should have been initialized in the WRF input file. Nothing to do
    return

  case default
    ! A non standard value of fire_can_fuel_read has been given.
    ! Don't know what the user wants to do = crash
    call crash("set_canopy_data: fire_can_fuel_read value is not valid. Please check namelist.input file.")

  end select

end subroutine set_canopy_data


!------------------------------------------------------------------------------------------------------------------
! SUBROUTINE: initialize_sfire_fields
!>   Initialize fields defined in sfire_fields module
!!
!!  Set a default value for each field in the module  
!------------------------------------------------------------------------------------------------------------------
subroutine initialize_sfire_fields( &
    ifms, ifme, jfms, jfme,         &
    ifts, ifte, jfts, jfte,         &
    nfuel_cat                       )

  use module_fr_sfire_fields, only:                                         &
    cfb, can_top, can_bottom, can_bulk_density, can_fx_z0, can_fuel_burnt,  &
    surf_mass_frac_b, surf_mass_frac_c, surf_moisture_frac
  
  implicit none

  ! 0. Declarations
  !    ============

  ! 0.1 Arguments
  integer,  intent(in   ) :: ifms                 !< fire mesh index start on x axis
  integer,  intent(in   ) :: ifme                 !< fire mesh index end   on x axis
  integer,  intent(in   ) :: jfms                 !< fire mesh index start on y axis
  integer,  intent(in   ) :: jfme                 !< fire mesh index end   on y axis
  integer,  intent(in   ) :: ifts                 !< fire tile index start on x axis
  integer,  intent(in   ) :: ifte                 !< fire tile index end   on x axis
  integer,  intent(in   ) :: jfts                 !< fire tile index start on y axis
  integer,  intent(in   ) :: jfte                 !< fire tile index end   on y axis

  real, dimension(ifms:ifme, jfms:jfme),  intent(in   )  :: nfuel_cat            !< surface fuel class

  ! 0.2 Locals
  integer :: surf_fuel_cat  !< local fuel category
  integer :: i, j, l        !< loop indexes

  ! 1. Initialize value on array
  !    =========================

  cfb(:, :)                 = 0.
  can_top(:, :)             = 0.
  can_bottom(:, :)          = 0.
  can_bulk_density(:, :)    = 0.
  can_fuel_burnt(:, :)      = 0.
  can_fx_z0(:, :)           = 0.
  surf_moisture_frac(:, :)  = 1.

  ! Initialize mass fraction for biomass and char
  do j = jfts, jfte
    do l = 1, n_layer_hf_surf
      do  i = ifts, ifte
        ! get local fuel cat and convert to integer
        surf_fuel_cat = int(nfuel_cat(i, j))
        
        ! ash fraction is considered to be effective mineral content se
        surf_mass_frac_b(i, l, j)   = 1. - layer_char_fraction(l) - se(surf_fuel_cat)
        surf_mass_frac_c(i, l, j)   = layer_char_fraction(l)
      enddo
    enddo
  enddo

end subroutine initialize_sfire_fields


!
!**********************
!            

real function fire_rate_of_spread(propx, propy, i,j,fp, save_cfb)
! compute rate of spread at grid node (i,j) in the direction (propx,propy)
! winds and everything else on the fire grid is passed in type fp
use module_fr_sfire_fields, only: cfb
implicit none
!***arguments
real, intent(in)::propx, propy! direction, need not be normalized
integer, intent(in)::i,j         ! node mesh coordinates
type(fire_params),intent(in)::fp
logical,  optional, intent(in)  :: save_cfb

!*** local
real:: ros_back,ros_wind,ros_slope,nvx,nvy,scale,rr
integer::ierrx     ! number of errors
character(len=128)::msg ! error message
real:: speed, tanphi, ros

!*** executable
!! Compute surface fire ros
call fire_ros(rr,propx,propy,i,j,fp,ierrx,msg)

!! Compute crown fire ros
call fire_crown_ros(fp, propx, propy, i, j, rr, ierrx, msg, save_cfb)

fire_rate_of_spread = rr

end function fire_rate_of_spread



subroutine fire_ros(ros, propx, propy, i, j, fp, ierrx, msg, surf_fuel_cat, wind_speed)

implicit none

! compute the wind speed and slope normal to the fireline and call fire_ros_cawfe or fire_ros_balbi

!*** arguments
real, intent(out)::ros ! rate of spread: backing, due to wind, due to slope
real, intent(in)::propx,propy  ! direction, not normalized
integer, intent(in)::i,j         ! node mesh coordinates
type(fire_params),intent(in)::fp ! everything else on the fire grid, indexed by i,j
integer, intent(out)::ierrx
character(len=*), intent(out)::msg
integer,    intent(in), optional  ::  surf_fuel_cat !< impose fuel category for surface fuel
real,       intent(in), optional  ::  wind_speed    !< impose wind speed

!*** local
real::ros_back,ros_wind,ros_slope, rr ! rate of spread: backing, due to wind, due to slope
real:: speed, tanphi ! windspeed and slope in the directino normal to the fireline
real::cor_wind,cor_slope,nvx,nvy,scale
integer::isave
integer::k

!*** executable
! Here we are checking which fuel category is being used
! If the fuel category is greater than or equal to the unburnable fuel, we set
! the ROS to 0
if (present(surf_fuel_cat)) then
    k = surf_fuel_cat
else
    k = int(fp%nfuel_cat(i, j))
endif

if (k == 0) then
  ! Unburnable fuel
  ros = 0. 

else if (k > 0 .and. k <= nfuelcats) then
  ! available fuel category     
        ! make normal direction length 1 unless 0
        ! scale=sqrt(propx*propx+propy*propy)
        ! if(scale > 0.)then
        !      nvx = propx/scale
        !      nvx = propx/scale
        ! else
        !     nvx = 0.
        !     nvy = 0.
        ! endif

        ! normalize direction 
        scale=sqrt(propx*propx+propy*propy + epsilon(0.))
        nvx=propx/scale
        nvy=propy/scale

        if (present(wind_speed)) then
            speed = wind_speed
        else
            speed = fp%vx(i, j) * nvx + fp%vy(i, j) * nvy
        endif

        ierrx=0
        if (ibeh .eq. 2) then
            call fire_ros_balbi(ros, nvx, nvy, i, j, fp, ierrx, msg, k, speed)
            ! if (ierrx.ne.0)then
            !     balbi_nc_count = balbi_nc_count + 1
            !     if(balbi_nc_count > balbi_nc_max) call crash('fire_ros_balbi failed limit exceeded')
            !     if ( balbi_msglevel < 4)then
            !         call message('running fire_ros_balbi again with more detailed messages')
            !         isave = balbi_msglevel
            !         balbi_msglevel=4
            !         call fire_ros_balbi(ros,nvx,nvy,i,j,fp,ierrx,msg)
            !         balbi_msglevel=isave
            !     endif
            !     call warning(msg)
            ! endif
        else if (ibeh .eq. 3) then
#ifdef USE_STARFIRE
            call fire_ros_wirc(ros, nvx, nvy, i, j, fp, ierrx, msg, k, speed)
#else
            call crash("fire_ros: fire_ros_wirc requires StarFIRE.")
#endif
        else if (ibeh .eq. 4) then
            call fire_ros_behave(ros,nvx,nvy,i,j,fp,ierrx,msg, k, speed)
        else
            if (fire_advection.ne.0) then ! from flags in module_fr_sfire_util
                ! wind speed is total speed
                if (present(wind_speed)) then
                    speed = wind_speed
                else
                    speed =  sqrt(fp%vx(i,j)*fp%vx(i,j)+ fp%vy(i,j)*fp%vy(i,j))+tiny(speed)
                endif
                ! slope is total slope
                tanphi = sqrt(fp%dzdxf(i,j)*fp%dzdxf(i,j) + fp%dzdyf(i,j)*fp%dzdyf(i,j))+tiny(tanphi)
                ! cos of wind and spread, if >0
                cor_wind =  max(0.,(fp%vx(i,j)*nvx + fp%vy(i,j)*nvy)/speed)
                ! cos of slope and spread, if >0
                cor_slope = max(0., (fp%dzdxf(i,j)*nvx + fp%dzdyf(i,j)*nvy)/tanphi)
            else
                ! wind speed in spread direction
                if (present(wind_speed)) then
                    speed = wind_speed
                else
                    speed =  fp%vx(i,j)*nvx + fp%vy(i,j)*nvy
                endif
                ! slope in spread direction
                tanphi = fp%dzdxf(i,j)*nvx + fp%dzdyf(i,j)*nvy
                cor_wind=1.
                cor_slope=1.
            endif
            call fire_ros_cawfe(ros_back,ros_wind,ros_slope, &
            speed,tanphi,cor_wind,cor_slope,i,j,fp,ierrx,msg)
            if (ierrx.ne.0)then
                call warning(msg)
            endif
            ros = ros_back+ros_wind+ros_slope
            if(compare_rothermel)then
                call fire_ros_rothermel(rr, k, speed, tanphi, fp%fmc_g(i,j))
                if (abs(ros-rr)>1e-4)then
                    write(msg,'(2(a,i6),2(a,f10.6),a,g15.4)')'i=',i,' j=',j,  &
                        ' fire_ros_cawfe:',ros,' rothermel:',rr,' diff:',ros-rr
                    call warning(msg)
                endif
            endif
        endif
else 
  ! fuel category not supported
  write(msg,*) "fire_ros: fuel categories ", k, " at ", i, j, " not available"
  call crash(msg)
endif
if(fire_grows_only.gt.0)then
ros=max(ros,0.)
endif
end subroutine fire_ros


!------------------------------------------------------------------------------------------------------------------
! SUBROUTINE: fire_crown_ros
!>   Compute the crown fire rate of spread
!!
!!   Note Method is from Scott, J. H. (2001). Assessing crown fire potential by linking models of surface and crown
!!   fire behavior. 
!!   US Department of Agriculture, Forest Service, Rocky Mountain Research Station.
!!
!!   todo:   - Use previous time step value as starting values for w_init and w_act computation
!!  
!------------------------------------------------------------------------------------------------------------------
subroutine fire_crown_ros(fp, propx, propy, i, j, ros, ierrx, msg, save_cfb)
    
  use module_fr_sfire_util, only: fire_use_windrf
  use module_fr_sfire_fields, only: cfb, can_top, can_bottom, can_bulk_density
  
  implicit none

  ! 0. Declaration
  !    -----------

  ! 0.1 Arguments
  !     ----------

  type(fire_params),  intent(in   ) :: fp       !< fire properties
  real,               intent(in   ) :: propx    !< x direction, not normalized
  real,               intent(in   ) :: propy    !< y direction, not normalized
  integer,            intent(in   ) :: i        !< position on the fire grid on x axis
  integer,            intent(in   ) :: j        !< position on the fire grid on y axis
  real,               intent(inout) :: ros      !< rate of spread with crowning
  integer,            intent(  out) :: ierrx    !< error flag
  character(len=*),   intent(  out) :: msg      !< error message
  logical,  optional, intent(in   ) :: save_cfb !< store cfb in global cfb variable

  ! 0.2 Locals
  !     ------
  
  integer, parameter  :: max_ite = 20         !< max iteration for convergence
  integer, parameter  :: crown_fuel_cat = 10  !< fuel category used for canopy
  real, parameter     :: tol = 1e-5           !< tolerance for convergence

  integer             :: n_ite                !< iteration counter for passive crowning bounds
  integer             :: surf_fuel_cat        !< fuel category used for surface
  
  real                :: f_n, f_n_1           !< score function for passive crowning bounds
  real                :: inc                  !< wind increment for passive crowning bounds
  real                :: nvx, nvy             !< normal direction
  real                :: ros_init_prime       !< ros that ingnite canopy
  real                :: ros_act_prime        !< active ros threshold
  real                :: ros_active_crown     !< active crown fire ros 
  real                :: scale                !< normal norm
  real                :: w_n                  !< current wind value for passive crowning bounds
  real                :: w_n_1                !< next iteration for wind value for passive crowning bounds
  real                :: wind_act             !< wind speed for active crowning initiation
  real                :: wind_init            !< wind speed for passive crowning initiation
  real                :: wind_speed           !< wind speed in the spread direction
  real                :: wind_speed_no_wrf    !< wind speed in the spread direction without wrf

  ! tmp
  real    :: cfb_local
  logical :: save_cfb_local


  ! 1. Initialisation
  !    --------------

  ! save crown fraction burnt flag
  save_cfb_local = .false.
  if (present(save_cfb)) save_cfb_local = save_cfb

  ! get local fuel class of surface fuel
  surf_fuel_cat = fp%nfuel_cat(i, j)

  ! get normal direction and wind
  scale = sqrt(propx ** 2 + propy ** 2 + epsilon(0.))
  nvx = propx / scale
  nvy = propy / scale
  wind_speed =  fp%vx(i, j) * nvx + fp%vy(i, j) * nvy

  ! remove the wind reduction factor from the projection
  if (fire_use_windrf == 1) then
    wind_speed_no_wrf = wind_speed / windrf(surf_fuel_cat)
  else
    wind_speed_no_wrf = wind_speed
  endif

  ! 2. Compute initiation ROS threshold
  !    --------------------------------

  ! R_init_prime = 60 * ((canopy_base_height * (460 + 25.9 * canopy_foliar_moisture_content))) ** 1.5 / (delta_H * sigma_d)
  ros_init_prime = 60. * ((can_bottom(i, j) * (460. + 25.9 * fuelmc_c))) ** 1.5 / (cmbcnst * fgi(surf_fuel_cat))

  ! surface ros is below initiation threshold, no passive of active crowning -> do not change ros and return
  if (ros <= ros_init_prime .or. can_bulk_density(i, j) == 0.) then
    if (save_cfb_local) cfb(i, j) = 0.
    return
  endif

  ! 3. Compute active ROS threshold
  !    ----------------------------

  ! active ros threshold
  ros_act_prime = 0.05 / can_bulk_density(i, j)

  ! compute active ros
  ! R_active = 3.34 * surface_ros(fuel=10, wind= 0.4 * 10m wind)
  call fire_ros(ros_active_crown, propx, propy, i, j, fp, ierrx, msg, surf_fuel_cat=crown_fuel_cat, wind_speed=0.4 * wind_speed_no_wrf)
  ros_active_crown = 3.34 * ros_active_crown

  ! active ros is above active threshold, return active crowning
  if (ros_active_crown >= ros_act_prime) then
    ros = ros_active_crown
    if (save_cfb_local) cfb(i, j) = 1.
    return
  endif

  ! 4. Compute passive ROS
  !    -------------------

  ! 4.1 get the wind that leads to surface_ros == R_init_prime

  ! Init values for wind and ros
  w_n = 0.
  w_n_1 = wind_speed
  ! f_n = ros(surface_fuel, w_n) - ros_init_prime
  call fire_ros(f_n, propx, propy, i, j, fp, ierrx, msg, surf_fuel_cat=surf_fuel_cat, wind_speed=w_n)
  f_n = f_n - ros_init_prime
  ! f_n_1 = ros - ros_init_prime (as surface ros with input wind has already been computed previously)
  f_n_1 = ros - ros_init_prime
  
  ! secante method to cenverge
  n_ite = 1
  do while (abs(f_n) > tol .AND. n_ite <= max_ite)
    ! update increment
    inc = (w_n - w_n_1) / (f_n - f_n_1)
    
    ! update n-1 values
    w_n_1 = w_n
    f_n_1 = f_n
    
    ! update n values
    w_n = min(30., max(0., w_n - inc * f_n))
    call fire_ros(f_n, propx, propy, i, j, fp, ierrx, msg, surf_fuel_cat=surf_fuel_cat, wind_speed=w_n)
    f_n = f_n - ros_init_prime
    
    ! update counter
    n_ite = n_ite + 1
  enddo
  
  ! save value of wind
  wind_init = w_n
  
  ! 4.2 get the wind that leads to active_ros == R_act_prime
  
  ! Init values for wind and ros
  w_n = 30.
  w_n_1 = wind_speed
  ! f_n = 3.34 * ros(crown_fuel, 0.4 * w_n) - ros_init_prime
  call fire_ros(f_n, propx, propy, i, j, fp, ierrx, msg, surf_fuel_cat=crown_fuel_cat, wind_speed=0.4 * w_n)
  f_n = 3.34 * f_n - ros_act_prime
  f_n_1 = ros_active_crown - ros_act_prime

  ! secante method to cenverge
  n_ite = 1
  do while (abs(f_n) > tol .AND. n_ite <= max_ite)
    ! update increment
    inc = (w_n - w_n_1) / (f_n - f_n_1)
    
    ! update n-1 values
    w_n_1 = w_n
    f_n_1 = f_n
    
    ! update n values
    w_n = min(30., max(0., w_n - inc * f_n))
    call fire_ros(f_n, propx, propy, i, j, fp, ierrx, msg, surf_fuel_cat=crown_fuel_cat, wind_speed=0.4 * w_n)
    f_n = 3.34 * f_n - ros_act_prime
    
    ! update counter
    n_ite = n_ite + 1

  enddo
  
  ! save value of wind
  wind_act = w_n
  
  ! 4.3 Compute crown fraction burnt

  cfb_local = min(max(0., (wind_speed - wind_init) / (wind_act - wind_init)), 1.)

  ! 4.4 Combine surface and active ros to get passive ros
  ros = min(6., ros + cfb_local * (ros_active_crown - ros))

  if (save_cfb_local) cfb(i, j) = cfb_local

end subroutine fire_crown_ros

!
!***
!
subroutine fire_ros_rothermel(ros,nfuel,speed,tanphi,fmc_g)
! in
!       nfuel   number of fuel description
!       speed   wind speed
!       tanphi  slope
!       fmc_g   fuel moisture content (0 to 1)
! out
!       ros     rate of spread

implicit none
!*** arguments

integer, intent(in) :: nfuel
real, intent(in):: speed,tanphi,fmc_g
real, intent(out):: ros

!*** local
real:: fuelheat, fuelloadm, fuelload, fueldepth, betafl, betaop, qig, epsilon, &
rhob, c, bbb, e, phiwc, rtemp2, gammax, a, ratio, gamma, wn, rtemp1, etam, etas, &
ir, irm, xifr, r_0,  spdms, umidm, umid, phiw, phis
integer::k

!*** executable

k = nfuel                        ! fuel index


fuelheat = cmbcnst * 4.30e-04             ! convert J/kg to BTU/lb
fuelloadm= fgi(k)                ! fuelload without moisture (fgi is oven-dry fuel load)
fuelload = fuelloadm * (.3048)**2 * 2.205  ! to lb/ft**2
fueldepth= fueldepthm(k)/0.3048              ! to ft
betafl   = fuelload/(fueldepth * fueldens(k))! packing ratio  jm: lb/ft**2/(ft * lb*ft**3) = 1
betaop   = 3.348 * savr(k)**(-0.8189)         ! optimum packing ratio jm: units??  
qig      = 250. + 1116.*fmc_g          ! heat of preignition, btu/lb
epsilon  = exp(-138./savr(k) )               ! effective heating number
rhob     = fuelload/fueldepth             ! ovendry bulk density, lb/ft**3
c        = 7.47 * exp(-0.133 * savr(k)**0.55) ! const in wind coef
bbb      = 0.02526 * savr(k)**0.54            ! const in wind coef
!c        = c * windrf**bbb                 ! jm: wind reduction from 20ft per Baughman&Albini(1980)
e        = 0.715 * exp( -3.59e-4 * savr(k))  ! const in wind coef
phiwc    = c * (betafl/betaop)**(-e) 
rtemp2   = savr(k)**1.5
gammax   = rtemp2/(495. + 0.0594*rtemp2)  ! maximum rxn vel, 1/min
a        = 1./(4.774 * savr(k)**0.1 - 7.27)   ! coef for optimum rxn vel
ratio    = betafl/betaop   
gamma    = gammax*(ratio**a)*exp(a*(1.-ratio)) ! optimum rxn vel, 1/min
wn       = fuelload/(1 + st(k))              ! net fuel loading, lb/ft**2
rtemp1   = fmc_g/fuelmce(k)
etam     = 1.-2.59*rtemp1 +5.11*rtemp1**2 -3.52*rtemp1**3  ! moist damp coef
etas     = 0.174* se(k)**(-0.19)              ! mineral damping coef
ir       = gamma * wn * fuelheat * etam * etas ! rxn intensity,btu/ft**2 min
irm      = ir * 1055./( 0.3048**2 * 60.) * 1.e-6! for mw/m**2 (set but not used)
xifr     = exp( (0.792 + 0.681*savr(k)**0.5) &
            * (betafl+0.1)) /(192. + 0.2595*savr(k)) ! propagating flux ratio
!        ... r_0 is the spread rate for a fire on flat ground with no wind.
r_0      = ir*xifr/(rhob * epsilon *qig)  ! default spread rate in ft/min

if (.not.ichap(k)>0.)then
    !       ... if wind is 0 or into fireline, phiw = 0, &this reduces to backing ros.
    spdms = max(speed,0.)
    umidm = min(spdms,30.)                    ! max input wind spd is 30 m/s   !param!
    umid = umidm * 196.850                    ! m/s to ft/min
    !  eqn.: phiw = c * umid**bbb * (betafl/betaop)**(-e) ! wind coef
    phiw = umid**bbb * phiwc                   ! wind coef
    phis = 5.275 * betafl**(-0.3) *max(0.,tanphi)**2   ! slope factor
    ros = r_0*(1. + phiw + phis)  * .00508 ! spread rate, m/s
else  ! chapparal
    !        .... spread rate has no dependency on fuel character, only windspeed.
    spdms = max(speed,0.)
    ros = max(.03333,1.2974 * spdms**1.41)       ! spread rate, m/s
endif

ros=min(ros,6.);

end subroutine fire_ros_rothermel

!
!***
!


subroutine fire_ros_balbi(ros, nvx, nvy, i, j, fp, ierrx, msg, surf_fuel_cat, wind_speed)
! in
!       nvx,nvy direction of propagation, normalized 
!       i,j     indices in the fire grid
!       fp      fire grid arrays 
! out
!       ros     rate of spread
!       ierrx   error code
!       msg     error message
! global
!       balbi_msglevel 0=nothing 1=errors 2=warnings 3=diagnostics 4=iterations
!       fuel categories from namelist.fire

!*** arguments

implicit none
real, intent(in):: nvx,nvy
integer, intent(in):: i,j
type(fire_params):: fp
real, intent(out):: ros
integer, intent(out)::ierrx 
character(len=*), intent(out)::msg
integer,    intent(in)  ::  surf_fuel_cat !< impose fuel category for surface fuel
real,       intent(in)  ::  wind_speed    !< impose wind speed

! Constants
real, parameter:: tol_R = 1e-5          ! tolerance to compute R
integer, parameter:: maxit_R = 20       ! max iterations to compute R
real, parameter:: pi=3.141592653589793

!*** local
real:: vx, vy, dzdxf, dzdyf, speed, tanphi, fmc_g, nfuel_cat, tempf, rhof, ffwidth, fgip
real:: dead_fract
real:: B_, chi_0, r_00, delta_h, tau_0, rho_a, s_t, T_i, g, C_pa, U, alpha, W_0
real:: m, sigma, sigma_t, h, s, rho_v, delta_H_, C_p
real:: beta, beta_t, S_, S_t_, nu, T_a, delta_T, q, T_f, u_0
real:: R_b, A, gamma, H_f, l_f
real:: W_00, a_up, a_lat, K_, b
real:: R_old, R_r, R_c, R, change
real:: C_pw, T_vap
integer:: k, mc, it
! nfuel_cat is real only because real array was easier than integer in WPS long ago

!*** executable

! initialize error flag
ierrx = 0

! retrieve balbi inputs
vx = fp%vx(i,j)
vy = fp%vy(i,j)
dzdxf = fp%dzdxf(i,j)
dzdyf = fp%dzdyf(i,j)
speed = wind_speed
tanphi = dzdxf*nvx + dzdyf*nvy      ! slope in spread the direction
fmc_g  = fp%fmc_g(i,j)              ! fuel moisture content (water / dry)
nfuel_cat = surf_fuel_cat           ! surface fuel categories
tempf = fp%tempf(i,j)               ! temperature interpolated to the fire grid
rhof  = fp%rhof(i,j)                ! air density interpolated to the fire grid
ffwidth = fp%ffwidth(i,j)           ! fire front width = 1/fireline curvature (?)
fgip = fp%fgip(i,j)                 ! fuel load

! retrieve static values by indexing from fuel categories, as fgip = fp%fgi(i,j) instead of above
! do not delete any fields from fp they are used by Rothermel model

if(balbi_msglevel.ge.3)then  ! print the inputs
!$OMP CRITICAL(SFIRE_PHYS_CRIT)
    write(msg,*)'fire_ros_balbi called with balbi_msglevel=',balbi_msglevel
    call message(msg)
    call message('0=nothing 1=errors 2=warnings 3=diagnostics 4=iterations')
    call message('inputs:')
    write(msg,*)'wind speed projected on propagation direction speed=',speed
    call message(msg)
    write(msg,*)'slope projected on propagation direction tanphi=',tanphi
    call message(msg)
    write(msg,*)'fuel moisture content (water / dry) fmc_g=',fmc_g
    call message(msg)
    write(msg,*)'nfuel_cat=',nfuel_cat
    call message(msg)
    write(msg,*)'temperature tempf=',tempf
    call message(msg)
    write(msg,*)'surface air density rhof=',rhof
    call message(msg)
    write(msg,*)'fire front width=',ffwidth
    call message(msg)
    write(msg,*)'total fuel load fgip=',fgip
    call message(msg)
!$OMP END CRITICAL(SFIRE_PHYS_CRIT)
endif
if (.not.ffwidth > 0.) call crash('fire_ros_balbi: ffwidth must be >0')

! not allowing negative for now

!*** continue here
speed=max(speed,0.)
tanphi=max(tanphi,0.)

k = nfuel_cat           ! converting fuel index to integer

! --- Constants and model parameters (Balbi 20220622 slide 9) ---
B_ = 5.670373e-8                ! Stefan-Boltzman constant              [W/m^2/K^4]
chi_0 = 0.3                     ! Radiative factor
r_00 = 2.5e-5                   ! Model parameter
delta_h = 2.3e6                 ! Heat of latent evaporation            [J/kg]
tau_0 = 75591.                  ! Anderson's residence time coefficient
rho_a = rhof                    ! Air density                           [kg/m3]
s_t = 17.                       ! Stoichiometric coefficient
T_i = 600.                      ! Ignition temperature                  [K]
g = 9.81                        ! Gravitational acceleration            [m/s^2]
C_pa = 1150                     ! Specific heat of air                  [J/kg/K]

U = speed                       ! Normal component of wind speed        [m/s]
alpha = atan(tanphi)            ! Local terrain slope angle             [rad]
W_0 = 50.                       ! Length of the ignition line           [m]


! --- Input fuel data (Balbi 20220622 slide 2) ---
m = fmc_g                       ! Dead fuel moisture content
! sigma = fgip * dead_fract     ! Dead fine fuel load                   [kg/m^2]
! fgip is the total oven dry fuel load. The total fuel load should take liquid water content into account.
sigma = fgip * (1 + m)          ! Dead fine fuel load                   [kg/m^2]
sigma_t = fgip * (1 + m)        ! Total fine fuel load                  [kg/m^2]
h = fueldepthm(k)               ! Fuel height                           [m]
s = savr(k)/0.3048              ! Surface area to volume ratio  (conversion from 1/ft to 1/m) [1/m]
rho_v = fueldens(k)*16.015721869155822   ! Fuel density (conversion from lb/ft^3 to kg/m^3) [kg/m^3]
delta_H_ = cmbcnst              ! Heat of combustion of pyrolysis gases [J/kg]
C_p = 1200.                     ! Specific heat of fuel                         [J/kg/K]
T_vap = 373.15                  ! Water evaporation temperature [K]
C_pw = 4180.                    ! Specific heat of liquid water                         [J/kg/K]


! --- Preliminary computations (Balbi 20220622 slide 3) ---
! Packing ratios - dead, total
beta = sigma / (h * rho_v)
beta_t = sigma_t / (h * rho_v)

! Leaf areas - dead, total
S_ = s * h * beta
S_t_ = s * h * beta_t

! Absorbtion coefficient
nu = S_t_

! Temperature difference
T_a = tempf
delta_T = T_i - T_a

! Activation energy
q = C_p * delta_T + m * (delta_h + C_pw * (T_vap - T_a))

! Flame temperature
T_f = T_a + (delta_H_ * (1 - chi_0)) / ((s_t + 1) * C_pa)

! Upward gas velocity
u_0 = 2. * (s_t + 1.0) * rho_v * T_f * min(S_, 2 * pi * S_ / S_t_) / (tau_0 * rho_a * T_a)

! --- Radiative contributions (slide 4) ---
! Flame base radiation contribution to ROS
R_b = min(S_t_ / (2. * pi), 1.) * (beta / beta_t)**2 * (B_ * T_f**4) / (beta * rho_v * q)

! Radiative coefficient
A = min(S_t_ / (2. * pi), S_ / S_t_) * chi_0 * delta_H_ / (4. * q)

! Flame tilt angle
gamma = max(0., atan(tanphi + U / u_0))

! Flame height
H_f = (u_0**2) / (g * (T_f / T_a - 1) * (cos(alpha))**2)

! Flame length
l_f = H_f / (cos(gamma - alpha))


! --- Convective contribution (slide 5) ---
! Heat loss coefficients
! For most fires, the fire front width (W_0) and corresponding threshold for
! convective cooling (W_00) can be ignored; typical values for W_00 are in the
! 33-100 m range, depending on the fuel.
! W_00 = ffwidth_thresh
! if (W_00 == 0.) then    ! force min(W_0 / W_00, 1.) to 1 => no convective cooling
!     W_0 = 1.
!     W_00 = 1.
! endif
a_up = 0.025 * s * beta_t * sqrt(h)
! a_lat = min(W_0 / W_00, 1.)
a_lat = 1.

! Drag force coefficient
! NOTE: denominator is identical to a_lat above; perhaps slide meant to
! specified a_lat instead?
! K_ = (beta_t / min(W_0 / W_00, 1.))
K_ = beta_t

! Convective coefficient
b = (a_up * a_lat * rho_a * delta_H_ * T_a) / (2 * q * (s_t + 1) * rho_v * beta_t * T_f)


! --- Algorithm (slide 7) ---
! Iteratively solve for the total ROS.  Seed the search using the flame
! base radiation contribution, which is invariant.
R_old = R_b
do it=1,maxit_R
    ! Update the iterative ROS contributions, then recompute the total ROS.
    R_r = A * R_old * ((1. + sin(gamma) - cos(gamma)) / (1. + ((R_old * cos(gamma)) / (s * r_00))))
    R_c = b * (min(2. * pi * S_ / S_t_, S_) * tanphi * (1 + s_t) * rho_v * T_f / (rho_a * T_a * tau_0)  + (U * exp(-K_ * R_old)))
    R = R_b + R_r + R_c
    change = R - R_old

    ! TODO: Review if this non-negative check is still needed
    if (R<0.)then
        if (balbi_msglevel .ge. 3)then
!$OMP CRITICAL(SFIRE_PHYS_CRIT)
            write(msg,1)it,i,j,R,change
            call message(msg)
            write(msg,'(3(a,f8.3,1x))') 'R_b=',R_b,'R_r=',R_r,'R_c=',R_c
            call message(msg)
            call message('fire_ros_balbi R<0 resetting to R=0')
!$OMP END CRITICAL(SFIRE_PHYS_CRIT)
        endif
        R=0.
    endif

    if (R .ne. R)call crash('fire_ros_balbi: R is NaN')

    if (balbi_msglevel .ge. 4)then 
!$OMP CRITICAL(SFIRE_PHYS_CRIT)
        write(msg,1)it,i,j,R,change
        call message(msg)
!$OMP END CRITICAL(SFIRE_PHYS_CRIT)
    endif

    if (abs(change) < tol_R)then
        if (balbi_msglevel .ge. 3)then
!$OMP CRITICAL(SFIRE_PHYS_CRIT)
            write(msg,1)it,i,j,R,change,' converged'
            call message(msg)
!$OMP END CRITICAL(SFIRE_PHYS_CRIT)
        endif
        ros = R
        return 
    endif
    R_old = R
enddo

1 format('fire_ros_balbi: iteration',i3,' at ij=',2i5,' R=',f6.3,' diff=',g12.3,a)
!$OMP CRITICAL(SFIRE_PHYS_CRIT)
write(msg,1)it,i,j,R,change,' did not converge'
!$OMP END CRITICAL(SFIRE_PHYS_CRIT)
ierrx=1
ros=R
end subroutine fire_ros_balbi




!
!***
!

subroutine fire_ros_behave(ros, nvx, nvy, i, j, fp, ierrx, msg, surf_fuel_cat, wind_speed)
! in
!       nvx,nvy direction of propagation, normalized 
!       i,j     indices in the fire grid
!       fp      fire grid arrays 
! out
!       ros     rate of spread
!       ierrx   error code
!       msg     error message
! global
!       fuel categories from namelist.fire
!
! Implementation of fire spread model in pages 16-38 of paper:
! Andrews, Patricia L.. “The Rothermel surface fire spread model and associated
! developments: A comprehensive explanation.” (2018).

!*** arguments

implicit none
real, intent(in):: nvx,nvy
integer, intent(in):: i,j
type(fire_params):: fp
integer,    intent(in), optional  ::  surf_fuel_cat !< impose fuel category for surface fuel
real,       intent(in), optional  ::  wind_speed    !< impose wind speed
real, intent(out):: ros
integer, intent(out)::ierrx 
character(len=*), intent(out)::msg

!*** local
real::fmc_1, fmc_10, fmc_100, fmc_1000, fmc_herba_dead, fmc_herba, fmc_woody, fmc_dead, fmc_live
real::fuelload_1, fuelload_10, fuelload_100, fuelload_1000, fuelload_herba_dead, fuelload_herba, fuelload_woody
real::savr_1, savr_10, savr_100, savr_1000, savr_herba_dead, savr_herba, savr_woody, savr_dead, savr_live
real::A_1, A_10, A_100, A_herba_dead, A_herba, A_woody, A_dead, A_live
real::f_1, f_10, f_100, f_herba_dead, f_herba, f_woody, f_dead, f_live
real::wn_1, wn_10, wn_100, wn_herba_dead, wn_herba, wn_woody, wn_dead, wn_live
real::exp_1, exp_10, exp_100, exp_herba_dead, exp_herba, exp_woody, exp_dead, exp_live
real::qig_1, qig_10, qig_100, qig_herba_dead, qig_herba, qig_woody
real::fuelheat_dead, fuelheat_live, se_dead, se_live, rtemp1_dead, rtemp1_live
real::etam_dead, etam_live, etas_dead, etas_live, heat_sink_dead, heat_sink_live
real::epsilon_1, epsilon_10, epsilon_100, epsilon_herba_dead, epsilon_herba, epsilon_woody
real::fueldepth, fuelload, savr_t, AT, W_prima, M_prima_f, heat_sink, fuelmce_live
real::a, rhob, betafl, betaop, rtemp2, gammax, ratio, gamma
real::xifr, ir, r_0, umid, c, bbb, e, phiwc, phiw, phis
real::dzdxf, dzdyf, vx, vy, speed, tanphi, nfuel_cat
integer::k

!*** executable
! retrieve inputs
vx = fp%vx(i,j)
vy = fp%vy(i,j)
dzdxf = fp%dzdxf(i,j)
dzdyf = fp%dzdyf(i,j)
speed  = wind_speed                 ! wind speed in the spread direction
tanphi = dzdxf*nvx + dzdyf*nvy      ! slope in spread the direction
nfuel_cat = fp%nfuel_cat(i,j)       ! fuel categories
k = surf_fuel_cat                  ! fuel index
! fuel moisture content (water / dry) classes
fmc_1          = fp%fmc_gc_f(i,1,j)
fmc_10         = fp%fmc_gc_f(i,2,j)
fmc_100        = fp%fmc_gc_f(i,3,j)
fmc_1000       = fp%fmc_gc_f(i,4,j)
fmc_herba      = fp%fmc_gc_f(i,5,j)
fmc_woody      = fp%fmc_gc_f(i,6,j)
fmc_herba_dead = fmc_1
! print the inputs
if(behave_msglevel.ge.3)then
!$OMP CRITICAL(SFIRE_PHYS_CRIT)
    write(msg,*)'fire_ros_behave called with behave_msglevel=',behave_msglevel
    call message(msg)
    call message('0=nothing 1=errors 2=warnings 3=diagnostics')
    call message('inputs:')
    write(msg,*)'wind speed projected on propagation direction speed=',speed
    call message(msg)
    write(msg,*)'slope projected on propagation direction tanphi=',tanphi
    call message(msg)
    write(msg,*)'fuel moisture content (water / dry) fmc_1=',fmc_1
    call message(msg)
    write(msg,*)'fuel moisture content (water / dry) fmc_10=',fmc_10
    call message(msg)
    write(msg,*)'fuel moisture content (water / dry) fmc_100=',fmc_100
    call message(msg)
    write(msg,*)'fuel moisture content (water / dry) fmc_1000=',fmc_1000
    call message(msg)
    write(msg,*)'fuel moisture content (water / dry) fmc_herba=',fmc_herba
    call message(msg)
    write(msg,*)'fuel moisture content (water / dry) fmc_woody=',fmc_woody
    call message(msg)
    write(msg,*)'nfuel_cat=',nfuel_cat
    call message(msg)
!$OMP END CRITICAL(SFIRE_PHYS_CRIT)
endif
! fuel load (lb/ft^2) classes
fuelload_1          = fuelload_gc(k,1) !fp%fuelload(i,k,1,j)
fuelload_10         = fuelload_gc(k,2) !fp%fuelload(i,k,2,j)
fuelload_100        = fuelload_gc(k,3) !fp%fuelload(i,k,3,j)
fuelload_1000       = fuelload_gc(k,4) !fp%fuelload(i,k,4,j)
fuelload_herba      = fuelload_gc(k,5) !fp%fuelload(i,k,5,j)
fuelload_woody      = fuelload_gc(k,6) !fp%fuelload(i,k,6,j)
fuelload_herba_dead = 0
! savr classes
savr_1          = savr_gc(k,1)
savr_10         = savr_gc(k,2)
savr_100        = savr_gc(k,3)
savr_1000       = savr_gc(k,4)
savr_herba      = savr_gc(k,5)
savr_woody      = savr_gc(k,6)
savr_herba_dead = savr_herba
! translating units
fueldepth= fueldepthm(k)/0.3048     ! from m to ft
! dynamic load transfer from live to dead
if (idynamic(k)>0.) then
    if (fmc_herba<0.3) then
        fuelload_herba_dead = fuelload_herba_dead + fuelload_herba
        fuelload_herba = 0
    else if (fmc_herba < 1.2) then
        fuelload_herba_dead = fuelload_herba_dead * (1.333 - 1.11 * fmc_herba)
        fuelload_herba = fuelload_herba - fuelload_herba_dead
    end if
end if
! total fuel load
fuelload     = sum(fuelload_gc(k,:)) ! sum(fp%fuelload(i,k,:,j))
! mean total surface area per unit fuel cell of each size class within each category
A_1          = savr_1 * fuelload_1 / fueldens(k)
A_10         = savr_10 * fuelload_10 / fueldens(k)
A_100        = savr_100 * fuelload_100 / fueldens(k)
A_herba      = savr_herba * fuelload_herba / fueldens(k)
A_woody      = savr_woody * fuelload_woody / fueldens(k)
A_herba_dead = savr_herba_dead * fuelload_herba_dead / fueldens(k)
! mean total surface area per unit fuel cell of the dead and live categories
A_dead       = A_1 + A_10 + A_100 + A_herba_dead
A_live       = A_herba + A_woody
! mean total surface area of the fuel
AT           = A_dead + A_live
! weighting factors of each size class within each category
if (A_dead.eq.0.) then
    f_1           = 0
    f_10          = 0
    f_100         = 0
    f_herba_dead  = 0
else
    f_1           = A_1 / A_dead
    f_10          = A_10 / A_dead
    f_100         = A_100 / A_dead
    f_herba_dead  = A_herba_dead / A_dead
end if
if (A_live.eq.0.) then
    f_herba       = 0
    f_woody       = 0
else
    f_herba       = A_herba / A_live
    f_woody       = A_woody / A_live
end if
! weighting factors of the dead and live categories
f_dead        = A_dead / AT
f_live        = A_live / AT
! fuel particle surface-to-volume ratio of the dead and live categories
savr_dead     = f_1 * savr_1 + f_10 * savr_10 + f_100 * savr_100 + f_herba_dead * savr_herba_dead
savr_live     = f_herba * savr_herba + f_woody * savr_woody
savr_t        = f_dead * savr_dead + f_live * savr_live
! net fuel loading of each size class within each category
wn_1          = fuelload_1 * (1 - st(k))
wn_10         = fuelload_10 * (1 - st(k))
wn_100        = fuelload_100 * (1 - st(k))
wn_herba      = fuelload_herba * (1 - st(k))
wn_woody      = fuelload_woody * (1 - st(k))
wn_herba_dead = fuelload_herba_dead * (1 - st(k))
! net fuel loading of the dead and live categories
wn_dead       = (f_1 + f_herba_dead) * (wn_1 + wn_herba_dead) + f_10 * wn_10 + f_100 * wn_100
wn_live       = f_herba * wn_herba + f_woody * wn_woody
! live fuel moisture of extinction
exp_1          = fuelload_1 * exp(-138.0 / savr_1)
exp_10         = fuelload_10 * exp(-138.0 / savr_10)
exp_100        = fuelload_100 * exp(-138.0 / savr_100)
exp_herba      = fuelload_herba * exp(-500 / savr_herba)
exp_woody      = fuelload_woody * exp(-500 / savr_woody)
exp_herba_dead = fuelload_herba_dead * exp(-138.0 / savr_herba_dead)
exp_dead       = exp_1 + exp_10 + exp_100 + exp_herba_dead
exp_live       = exp_herba + exp_woody
if (exp_live.eq.0.) then
    W_prima    = 0.
else
    W_prima    = exp_dead / exp_live
end if
M_prima_f      = (fmc_1 * exp_1 + fmc_10 * exp_10 + fmc_100 * exp_100 + fmc_herba_dead * exp_herba_dead) / exp_dead
fuelmce_live   = 2.9 * W_prima * (1 - M_prima_f / fuelmce(k)) - 0.226;
if (fuelmce_live < fuelmce(k)) then
    fuelmce_live = fuelmce(k)
end if
! heat content of the dead and live categories
fuelheat_dead = (f_1 + f_10 + f_100 + f_herba_dead) * fuelheat
fuelheat_live = (f_herba + f_woody) * fuelheat
! fuel particle effective mineral content of the dead and live categories
se_dead       = (f_1 + f_10 + f_100 + f_herba_dead) * se(k)
se_live       = (f_herba + f_woody) * se(k)
! fuel moisture content of the dead and live categories
fmc_dead      = f_1 * fmc_1 + f_10 * fmc_10 + f_100 * fmc_100 + f_herba_dead * fmc_herba_dead
fmc_live      = f_herba * fmc_herba + f_woody * fmc_woody
! moisture damping coefficient of the dead and live categories
rtemp1_dead   = min(1.,fmc_dead/fuelmce(k));
rtemp1_live   = min(1.,fmc_live/fuelmce_live);
etam_dead     = 1. - 2.59*rtemp1_dead + 5.11*rtemp1_dead**2 - 3.52*rtemp1_dead**3;
etam_live     = 1. - 2.59*rtemp1_live + 5.11*rtemp1_live**2 - 3.52*rtemp1_live**3;
! mineral damping coefficient of the dead and live categories
if (se_dead.eq.0.) then
    etas_dead = 0.
else
    etas_dead = 0.174 * se_dead**(-0.19)
end if
if (se_live.eq.0.) then
    etas_live = 0.
else
    etas_live = 0.174 * se_live**(-0.19)
end if
! packing ratio  jm: lb/ft**2/(ft * lb*ft**3) = 1
betafl   = fuelload/(fueldepth * fueldens(k))
! optimum packing ratio
betaop   = 3.348 * savr_t**(-0.8189)   
! heat of preignition, btu/lb 
qig_1         = 250. + 1116.*fmc_1
qig_10        = 250. + 1116.*fmc_10
qig_100       = 250. + 1116.*fmc_100
qig_herba_dead= 250. + 1116.*fmc_herba_dead
qig_herba     = 250. + 1116.*fmc_herba
qig_woody     = 250. + 1116.*fmc_woody   
! effective heating number
if (savr_1.eq.0.) then
    epsilon_1 = 1.
else
    epsilon_1 = exp(-138. / savr_1)
end if
if (savr_10.eq.0.) then
    epsilon_10 = 1.
else
    epsilon_10 = exp(-138. / savr_10)
end if
if (savr_100.eq.0.) then
    epsilon_100 = 1.
else
    epsilon_100 = exp(-138. / savr_100)
end if
if (savr_herba_dead.eq.0.) then
    epsilon_herba_dead = 1.
else
    epsilon_herba_dead = exp(-138. / savr_herba_dead)
end if
if (savr_herba.eq.0.) then
    epsilon_herba = 1.
else
    epsilon_herba = exp(-138. / savr_herba)
end if
if (savr_woody.eq.0.) then
    epsilon_woody = 1.
else
    epsilon_woody = exp(-138. / savr_woody)
end if
! ovendry bulk density, lb/ft**3
rhob     = fuelload/fueldepth
if (savr_t.eq.0.) then
    a      = 0.
    rtemp2 = 0.
    c      = 7.47
    bbb    = 0.
else
    a      = 133. * savr_t**(-0.7913);
    rtemp2 = savr_t**1.5
    c      = 7.47 * exp(-0.133 * savr_t**0.55)
    bbb    = 0.02526 * savr_t**0.54
end if
!c        = c * windrf**bbb                 ! jm: wind reduction from 20ft per Baughman&Albini(1980)
e        = 0.715 * exp(-3.59e-4 * savr_t)  ! const in wind coef
phiwc    = c * (betafl/betaop)**(-e) 
! maximum rxn vel, 1/min
gammax   = rtemp2/(495. + 0.0594*rtemp2)
ratio    = betafl/betaop   
! optimum rxn vel, 1/min
if (ratio.eq.0) then
    gamma = 0.
else
    gamma = gammax*(ratio**a)*exp(a*(1.-ratio)) 
end if
! rxn intensity,btu/ft**2 min
ir    = gamma * (wn_dead * fuelheat_dead * etam_dead * etas_dead + &
            wn_live * fuelheat_live * etam_live * etas_live);
! propagating flux ratio
if (savr_t.eq.0.) then
    xifr = exp(0.792 * (betafl + 0.1)) / 192.;
else
    xifr = exp( (0.792 + 0.681*savr_t**0.5) &
               * (betafl+0.1) ) /(192. + 0.2595*savr_t);
end if
! heat sink of the dead and live categories
heat_sink_dead = rhob * f_dead * ( &
                    f_1 * epsilon_1 * qig_1 + &
                    f_10 * epsilon_10 * qig_10 + &
                    f_100 * epsilon_100 * qig_100 + &
                    f_herba_dead * epsilon_herba_dead * qig_herba_dead)
heat_sink_live = rhob * f_live * ( &
                    f_herba * epsilon_herba * qig_herba + &
                    f_woody * epsilon_woody * qig_woody)
heat_sink = heat_sink_dead + heat_sink_live
! r_0 is the spread rate for a fire on flat ground with no wind (ft/min)
if (heat_sink.eq.0.) then
    r_0 = 0.
else
    r_0 = ir*xifr / heat_sink
end if
! if wind is 0 or into fireline, phiw = 0, this reduces to backing ros
umid = min(max(0., speed * 196.850), 0.9 * ir)
! wind coef
phiw = umid**bbb * phiwc                   
! slope factor
phis = 5.275 * betafl**(-0.3) *max(0.,tanphi)**2  
! spread rate, m/s
ros = r_0*(1. + phiw + phis)  * .00508
! print the outputs
if(behave_msglevel.ge.3)then
!$OMP CRITICAL(SFIRE_PHYS_CRIT)
    write(msg,*)'fuelheat_dead=',fuelheat_dead
    call message(msg)
    write(msg,*)'fuelheat_live=',fuelheat_live
    call message(msg)
    write(msg,*)'etas_dead=',etas_dead
    call message(msg)
    write(msg,*)'etas_live=',etas_live
    call message(msg)
    write(msg,*)'etam_dead=',etam_dead
    call message(msg)
    write(msg,*)'etam_live=',etam_live
    call message(msg)
    write(msg,*)'wn_dead=',wn_dead
    call message(msg)
    write(msg,*)'wn_live=',wn_live
    call message(msg)
    write(msg,*)'gamma=',gamma
    call message(msg)
    write(msg,*)'ir=',ir
    call message(msg)
    write(msg,*)'xifr=',xifr
    call message(msg)
    write(msg,*)'heat_sink=',heat_sink
    call message(msg)
    write(msg,*)'phiw=',phiw
    call message(msg)
    write(msg,*)'phis=',phis
    call message(msg)
    write(msg,*)'ros=',ros
    call message(msg)
!$OMP END CRITICAL(SFIRE_PHYS_CRIT)
endif

end subroutine fire_ros_behave

!
!***
!

subroutine fire_ros_cawfe(ros_back,ros_wind,ros_slope, &
speed,tanphi,cor_wind,cor_slope,i,j,fp,ierrx,msg)

implicit none

! find rate of spread from wind speed and slope
! copied from wf2_janice 
! with the following changes ONLY: 
!   0.5*(speed + abs(speed)) -> max(speed,0.)
!   index l -> j 
!   took out some prints
!   argument fuelloadm never used??
!   not using nfuel_cat here - cell info was put into arrays passed as arguments
!       in include file to avoid transcription errors when used elsewhere
!   betaop is absorbed in phiwc, see module_fr_sfire_model/fire_startup
!   return the backing, wind, and slope contributions to the rate of spread separately
!       because they may be needed to take advantage of known wind and slope vectors.
!       They should add up to get the total rate of spread.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     ... calculates fire spread rate with mcarthur formula or Rothermel
!           using fuel type of fuel cell
!
!      
!         m/s =(ft/min) *.3048/60. =(ft/min) * .00508   ! conversion rate
!         ft/min = m/s * 2.2369 * 88. = m/s *  196.850 ! conversion rate
!      
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!*** arguments
real, intent(out)::ros_back,ros_wind,ros_slope ! rate of spread: backing, due to wind, due to slope
real, intent(in)::speed,tanphi,cor_wind,cor_slope
integer, intent(in)::i,j         ! node mesh coordinates
type(fire_params),intent(in)::fp
integer, intent(out)::ierrx
character(len=*), intent(out)::msg

!*** local
real:: umid, phis, phiw, spdms, umidm, excess, tanphim,ros
real, parameter::ros_max=6.
integer::k

!*** executable

ierrx = 0

if(ibeh.eq.2)then
        call crash('fire_ros_cawfe: ibeh=2 Balbi not supported here')
endif

if (.not. fp%ischap(i,j) > 0.) then            ! if not chaparral, do not test for .eq. 0 for speed
    if (ibeh .eq. 1) then                ! use Rothermel formula
!       ... if wind is 0 or into fireline, phiw = 0, &this reduces to backing ros.
        spdms = max(speed,0.)          ! 
        umidm = min(spdms,30.)       ! max input wind spd is 30 m/s   !param!
        umid = umidm * 196.850                    ! m/s to ft/min
        !  eqn.: phiw = c * umid**bbb(i,j) * (fp%betafl(i,j)/betaop)**(-e) ! wind coef
        ros_wind = fp%phiwc(i,j) * (umid**fp%bbb(i,j)) ! wind coef
        tanphim=max(tanphi,0.0)
        tanphim=min(tanphim,5.0)  ! jm
        ! phis = 5.275 *(fp%betafl(i,j))**(-0.3) *tanphim**2   ! slope factor
        ros_slope = fp%phisc(i,j) *tanphim**2   ! slope factor
        ! rosm = fp%r_0(i,j)*(1. + phiw + phis)  * .00508 ! spread rate, m/s
        ros_back = fp%r_0(i,j)
    elseif(ibeh.eq.98)then ! for testing only, spread rate = wind but not < 0
        ros_back = 0.
        ros_wind = max(speed,0.)
        ros_slope= 0.
    elseif(ibeh.eq.99)then ! for testing only, spread rate = wind speed, no mask <0
        ros_back = 0.
        ros_wind = speed
        ros_slope= 0.
    elseif(ibeh.eq.0)then                                  ! MacArthur formula (Australian)
        ! rosm = 0.18*exp(0.8424*max(speed,0.))
        ros_back = 0.18
        ros_wind = 0.18*exp(0.8424*max(speed,0.)) - ros_back
        ros_slope =0.
        ! note: ros = ros_back + ros_wind + ros_slope
    else ! error, but prevent unintialized variables
        ! just so that there is something there
        ros_back=-999.
        ros_wind=-999.
        ros_slope=-999.
    endif
    k = int(fp%nfuel_cat(i,j))
    ros=ros_back+ros_wind+ros_slope
    if(ros > 1e-6 .and. fp%fmc_g(i,j) > fuelmce(k))then
!$OMP CRITICAL(SFIRE_PHYS_CRIT)
          write(msg,'(a,2i6,3(a,e13.5))') 'fire_ros_cawfe: at ',i,j,' rate of spread',ros,' moisture ',fp%fmc_g(i,j),'> extinction =',fuelmce(k)
!$OMP END CRITICAL(SFIRE_PHYS_CRIT)
          ! call warning(msg)
          ierrx = 1
    endif
!
else   ! chaparral model from Clark et al 2004
!        .... spread rate has no dependency on fuel character, only windspeed.
    spdms = max(speed,0.)      
    ! rosm = 1.2974 * spdms**1.41       ! spread rate, m/s
    ! note: backing ros is 0 for chaparral without setting nozero value below
    !sp_n=.03333  
    ! chaparral backing fire spread rate 0.033 m/s   ! param!
    !rosm= max(rosm, sp_n)   ! no less than backing r.o.s.

    ros_back=.03333    ! chaparral backing fire spread rate 0.033 m/s   ! param!
    ros_wind = 1.2974 * spdms**1.41       ! spread rate, m/s
    ros_wind = max(ros_wind, ros_back)-ros_back
    ros_slope =0.

endif
!  multiply by the correction factors (from angle calculations)
ros_wind=ros_wind*cor_wind
ros_slope=ros_slope*cor_slope
!
!     ----------note!  put an 6 m/s cap on max spread rate -----------
! rosm= min(rosm, 6.)         ! no faster than this cap   ! param !

excess = ros_back + ros_wind + ros_slope - ros_max

if (excess > 0.)then
    ! take it out of wind and slope in proportion
    ros_wind = ros_wind - excess*ros_wind/(ros_wind+ros_slope)
    ros_slope = ros_slope - excess*ros_slope/(ros_wind+ros_slope)
endif


!     ... to rescale to veloc. carried by model, mult x (svel*snorm(1,3))= .1
!jm: huh ???
!     fire_ros = 0.1*rosm
!
!write(msg,*)i,j,' speed=',speed,' tanphi',tanphi,' ros=',ros_back,ros_wind,ros_slope
!call message(msg)

return

contains
real function nrm2(u,v)
real, intent(in)::u,v
nrm2=sqrt(u*u+v*v)
end function nrm2

end subroutine fire_ros_cawfe

subroutine fire_risk(fp,               &
        ifms,ifme,jfms,jfme,                      &  ! memory dims
        ifts,ifte,jfts,jfte,                      &  ! tile dims
        nfuel_cat,                                &
        f_ros0,f_rosx,f_rosy,f_ros,               &  ! fire spread diagnostic variables
        f_int,f_lineint,f_lineint2)                  ! fireline intensities for danger rating

!*** arguments
type(fire_params), intent(in)::fp
integer, intent(in)::                             &
        ifms,ifme,jfms,jfme,                      &  ! memory dims
        ifts,ifte,jfts,jfte                          ! tile dims
real, intent(in), dimension(ifms:ifme,jfms:jfme) :: nfuel_cat 
real, intent(out), dimension(ifms:ifme,jfms:jfme) :: &
        f_ros0,f_rosx,f_rosy,f_ros,               &  ! fire spread diagnostic variables
        f_int,f_lineint,f_lineint2                   ! fire intensities for danger rating

!*** local
integer:: i,j,k, ierrx
real::r
character(len=128)::msg
real,parameter::s=1.0/1.41421356237,o=1.,z=0.
integer,parameter::nv=8
real, dimension(nv)::nvx=(/o,s,z,-s,-o,-s,z,s/),  &
                     nvy=(/z,s,o,s,z,-s,-o,-s/)

!*** executable

do j=jfts,jfte
  do i=ifts,ifte
!   compute the fire spread rate and vector

    f_ros0(i,j) = 0.      ! direction-less spread rate component

    r = 0.
    do k=1,nv 
       r = max(r,fire_rate_of_spread(nvx(k),nvy(k),i,j,fp))
    enddo
    
    ! max fire front speed in this location (m/s)
    f_ros(i,j) = r 

  enddo
enddo

call fire_intensity(fp,                       &  ! fuel properties
        ifms,ifme,jfms,jfme,                      &  ! memory dims
        ifts,ifte,jfts,jfte,                      &  ! tile dims
        ifms,ifme,jfms,jfme,                      &  ! f_ros dims
        f_ros,nfuel_cat,                          & !in
        f_lineint,f_lineint2,f_int)                 ! fireline intensities out

end subroutine fire_risk

!
!***
!

subroutine fire_intensity(fp,                     &  ! fuel params
        ifms,ifme,jfms,jfme,                      &  ! memory dims
        ifts,ifte,jfts,jfte,                      &  ! tile dims
        irms,irme,jrms,jrme,                      &  ! memory dims for ros
        ros,nfuel_cat,                            &  ! rate of spread in
        fibyram,filimit,f_int)                       ! intensities out

!*** arguments
type(fire_params), intent(in)::fp
integer, intent(in)::                             &
        ifms,ifme,jfms,jfme,                      &  ! memory dims
        ifts,ifte,jfts,jfte,                      &  ! tile dims
        irms,irme,jrms,jrme                         ! memory dims for ros
real, intent(in), dimension(irms:irme,jrms:jrme) :: ros ! in rate of spread
real, intent(in), dimension(ifms:ifme,jfms:jfme) :: nfuel_cat 
real, intent(out), dimension(ifms:ifme,jfms:jfme) :: &
        fibyram,filimit                             ! out fireline intensities
real, intent(out), dimension(ifms:ifme,jfms:jfme), optional :: f_int ! fire intensity (J/m^2/s)

!*** local
integer:: i,j,k
real, dimension(ifts:ifte,jfts:jfte):: rate_frac
real:: dt_fake=1.

!*** executable

call heat_fluxes(dt_fake,fp,                       &
        ifms,ifme,jfms,jfme,                      &  ! memory dims
        ifts,ifte,jfts,jfte,                      &  ! tile dims
        irms,irme,jrms,jrme,                      &  ! ros dims
        0,                                        & ! use burnup
        fp%fgip,ros,                              & !in
        fibyram)                                    !out

! multiply by fuel fraction consumed in flaming zone for the category
do j=jfts,jfte
    do i=ifts,ifte
       k=int( nfuel_cat(i,j) )
       fibyram(i,j)=fibyram(i,j)*ffw(k)
    enddo
enddo
    ! fuel fraction loss per fire front unit length traveled per unit time
    ! burn_rate(i,j) = 0.5 * front_speed / fp%fuel_time(i,j) 

    ! fireline element of length ds moves in time dt by front_speed * dt covering area ds * dt * front_speed (m^2)
    ! after time dt the fuel fraction decrease is 0 at the leading edge and dt/fuel_time at the trailing edge 
    ! so the average fuel consumption over this zone is 0.5 * dt/fuel_time  (1)
    ! fuel load is fgip (kg/m^2)
    ! and the amount of fuel burned fireline length ds travels over time dt is 0.5 fgip * ds * dt^2 * front_speed/fuel_time (kg)
    ! note: dt because 1. it travels more 2. it burns longer
    ! fgip*burn_rate_frac(i,j) = (kg/m^2) * (m/s^2) = kg/m/s^2

    ! http://www.forestencyclopedia.net/p/p487
    ! H = I*w*r  = (J/kg) * (kg/m^2) * (m/s) = J/m/s
    ! fuel fraction loss per fire front unit length traveled 

do j=jfts,jfte
    do i=ifts,ifte
       rate_frac(i,j)=0.5*ros(i,j)/fp%fuel_time(i,j)
    enddo
enddo

! multiply by heat contents * fuel load
call heat_fluxes(dt_fake,fp,                       &
        ifms,ifme,jfms,jfme,                      &  ! memory dims
        ifts,ifte,jfts,jfte,                      &  ! tile dims
        ifts,ifte,jfts,jfte,                      &  ! rate_frac dims
        0,                                        & ! use burnup
        fp%fgip,rate_frac,                        & !in
        filimit)                                    !out

if(present(f_int))then
  do j=jfts,jfte
    do i=ifts,ifte
       k=int( nfuel_cat(i,j) )
       ! in time tr the fraction of fuel consumed is ffw=1-exp(-tr/fuel_time)
       ! giving flame residence time tr = - log(1-ffw)*fuel_time
       ! fire intensity is heat contents * fuel load * fraction consumed / flame residence time
       ! J/kg * kg/m^2 * 1 / s = J/m^2/s
       ! this is averaged over the flame residence time, 
       rate_frac(i,j)=ffw(k)/(fp%fuel_time(i,j)*(-log(1.-ffw(k))))
    enddo
  enddo

  ! multiply by heat contents * fuel load J/m^2

  call heat_fluxes(dt_fake,fp,                      &
        ifms,ifme,jfms,jfme,                      &  ! memory dims
        ifts,ifte,jfts,jfte,                      &  ! tile dims
        ifts,ifte,jfts,jfte,                      &  ! rate_frac dims
        0,                                        & ! use burnup
        fp%fgip,rate_frac,                        & !in
        f_int)                                    !out

endif

end subroutine fire_intensity

!*** executable

end module module_fr_sfire_phys
