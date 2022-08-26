module module_domain
! define the main datastructure passed to the model
! follows WRFV3/frame/module_domain.F but only fields needed here are defined
   integer, parameter:: max_tiles=10

   type domain

      ! inputs to the fire model

      ! changes every time step
      real, pointer, dimension(:,:)::uf => null()              ! winds on the fire mesh 
      real, pointer, dimension(:,:)::vf => null()              ! winds on the fire mesh

      ! constant arrays - set once at initialization
      real, pointer, dimension(:,:)::zsf => null()             ! terrain height
      real, pointer, dimension(:,:)::dzdxf => null()           ! terrain gradient 
      real, pointer, dimension(:,:)::dzdyf => null()           ! terrain gradient 
      real, pointer, dimension(:,:)::fxlong=>null()            ! fire mesh nodes coordinates
      real, pointer, dimension(:,:)::fxlat =>null()            ! fire mesh nodes coordinates
      real, pointer, dimension(:,:)::nfuel_cat => null()       ! index to fuel data table
      ! nfuel_cat is expanded to the coefficient arrays below from the table in namelist.fire
      real, pointer, dimension(:,:)::bbb => null()             ! fire spread formula coefficients
      real, pointer, dimension(:,:)::betafl => null()
      real, pointer, dimension(:,:)::phiwc => null()
      real, pointer, dimension(:,:)::r_0 => null()
      real, pointer, dimension(:,:)::fgip  => null()
      real, pointer, dimension(:,:)::ischap => null()          ! fire spread formula coefficients
      real, pointer, dimension(:,:)::fuel_time => null()       ! fuel coefficient  
      real, pointer, dimension(:,:)::u_2 => null()
      real, pointer, dimension(:,:)::v_2 => null()
      real, pointer, dimension(:,:)::w_2 => null()             ! vertical wind, not used
      real, pointer, dimension(:,:)::ph_2 => null()
      real, pointer, dimension(:,:)::phb => null()
      real, pointer, dimension(:,:)::z0 => null()
      real, pointer, dimension(:,:)::ht => null()
      real, pointer, dimension(:,:)::xlong => null()
      real, pointer, dimension(:,:)::xlat => null()
      real, pointer, dimension(:,:)::avg_fuel_frac => null()
      real, pointer, dimension(:,:)::grnhfx => null()
      real, pointer, dimension(:,:)::grnqfx => null()
      real, pointer, dimension(:,:)::canhfx => null()
      real, pointer, dimension(:,:)::canqfx => null()
      real, pointer, dimension(:,:)::uah => null()
      real, pointer, dimension(:,:)::vah => null()
      real, pointer, dimension(:,:)::rainc => null()
      real, pointer, dimension(:,:)::rainnc => null()
      real, pointer, dimension(:,:)::t2 => null()
      real, pointer, dimension(:,:)::q2 => null()
      real, pointer, dimension(:,:)::psfc => null()
      real, pointer, dimension(:,:)::rain_old => null()
      real, pointer, dimension(:,:)::t2_old => null()
      real, pointer, dimension(:,:)::q2_old => null()
      real, pointer, dimension(:,:)::psfc_old => null()
      real, pointer, dimension(:,:)::rh_fire => null()
      real, pointer, dimension(:,:)::h_fire => null()
      real, pointer, dimension(:,:)::ndwi => null()
      real, pointer, dimension(:,:)::fndwi => null()
      real, pointer, dimension(:,:,:)::fmc_lag => null()
      real, pointer, dimension(:,:,:)::fmc_gc => null()
      real, pointer, dimension(:,:,:)::fmc_equi => null()
      real, pointer, dimension(:,:,:)::fmep => null()
      real, pointer, dimension(:,:,:,:)::chem => null()
      real, pointer, dimension(:,:,:,:)::tracer => null()

      ! scalars
      ! mesh and housekeeping

      real:: sr_x, sr_y                                        ! atmosphere/fire refinement ratio
      integer::itimestep                                       ! time step number
      real:: xtime                                             ! minutes since simulation start
      integer:: num_tiles                                      ! parallel execution tiles info
      integer, dimension(max_tiles):: i_start,i_end,j_start,j_end
      real:: dt,dx,dy                                          ! time and space stepsizes (atmosphere scale)
      real:: u_frame,v_frame                                   ! must be zero, for compatibility
      real:: unit_fxlong,unit_fxlat
      real:: fmoist_nexttime,fmoist_lasttime

      ! fire model state arrays 
       
      ! because of WRF architecture, the fire model may not maintain any state internally
      ! all state arrays must be passed in every call
      real, pointer, dimension(:,:)::lfn => null()             ! level set function
      real, pointer, dimension(:,:)::tign_g => null()          ! time of ignintion
      real, pointer, dimension(:,:)::tign_in => null()         ! given fire arrival time
      real, pointer, dimension(:,:)::fuel_frac => null()       ! remaining fuel fraction
      real, pointer, dimension(:,:)::fuel_frac_burnt => null() ! fuel fraction burned in the time step
      real, pointer, dimension(:,:)::lfn_out  => null()        ! a copy of the level set function
      real, pointer, dimension(:,:)::fire_area  => null()      ! fraction of cell area burning

      ! fire model outputs

      real, pointer, dimension(:,:)::fgrnhfx  => null()        ! ground sensible heat flux
      real, pointer, dimension(:,:)::fgrnqfx  => null()        ! ground latent heat flux
      real, pointer, dimension(:,:)::fcanhfx  => null()        ! canopy sensible heat flux
      real, pointer, dimension(:,:)::fcanqfx  => null()        ! canopy latent heat flux

      ! other fire model arrays
      real, pointer, dimension(:,:)::ros  => null()            ! rate of spread - for diagnostics only
      real, pointer, dimension(:,:)::phisc  => null()          ! spread formula coeff
      real, pointer, dimension(:,:)::fz0  => null()            ! 
      real, pointer, dimension(:,:)::fwh  => null()            ! 
      real, pointer, dimension(:,:)::fmc_g  => null()          ! fuel moisture, ground
      real, pointer, dimension(:,:)::flineint  => null()       !
      real, pointer, dimension(:,:)::flineint2  => null()      ! 
      real, pointer, dimension(:,:)::f_lineint  => null()      !
      real, pointer, dimension(:,:)::f_lineint2  => null()     ! 
      real, pointer, dimension(:,:)::f_int  => null()          ! 
      real, pointer, dimension(:,:)::f_ros  => null()          ! 
      real, pointer, dimension(:,:)::f_ros0  => null()         ! 
      real, pointer, dimension(:,:)::f_rosx  => null()         ! 
      real, pointer, dimension(:,:)::f_rosy  => null()         ! 
      real, pointer, dimension(:,:)::fire_hfx  => null()       ! 
      real, pointer, dimension(:,:)::f_ros11  => null()            ! rate of spread - for diagnostics only
      real, pointer, dimension(:,:)::f_ros12  => null()            ! rate of spread - for diagnostics only
      real, pointer, dimension(:,:)::f_ros13  => null()            ! rate of spread - for diagnostics only
      real, pointer, dimension(:,:)::f_ros21  => null()            ! rate of spread - for diagnostics only
      real, pointer, dimension(:,:)::f_ros23  => null()            ! rate of spread - for diagnostics only
      real, pointer, dimension(:,:)::f_ros31  => null()            ! rate of spread - for diagnostics only
      real, pointer, dimension(:,:)::f_ros32  => null()            ! rate of spread - for diagnostics only
      real, pointer, dimension(:,:)::f_ros33  => null()            ! rate of spread - for diagnostics only


   end type domain

end module module_domain
