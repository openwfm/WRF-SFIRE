! SFIRE - Spread fire model in WRF-Fire
!
!*** Jan Mandel August 2007 - March 2011 
!*** email: Jan.Mandel@gmail.com

! For support please subscribe to the wrf-fire mailing list at NCAR at
! http://mailman.ucar.edu/mailman/listinfo/wrf-fire
! or go to http://www.openwfm.org/wiki/WRF-Fire_user_support 

! Current drafts of the technical documentation and
! user's guide can be found at

! http://www.openwfm.org/wiki/WRF-Fire_documentation
! http://www.openwfm.org/wiki/WRF-Fire_publications

! This module is the only entry point from WRF-ARW to the wildland 
! fire model. The call to sfire_driver advances the fire model by 
! one timestep. The fire model inputs the wind and outputs 
! temperature and humidity tendencies. The fire model also inputs a 
! number of constant arrays (fuel data, topography). Additional 
! arguments are model state (for data assimilation) and constant arrays 
! the model gives to WRF for safekeeping because it is not allowed 
! to save anything.

! This code as of mid-2011 is described in [1]. If you use this code, 
! please acknowledge our work by citing [1].
! Thank you.

! Acknowledgements
!
! The fire physics code is adapted from the CAWFE code [2].
! The coupling with WRF is adapted from a code by Ned Patton, 
! coupling a Fortran 90 port of the CAWFE fire module to WRF [3].
! Support of refined fire grids in WRF was provided by John Michalakes.
! Jonathan D. Beezley has set up and maintained the WRF build and
! execution environment, provided software engineering infrastructure 
! including synchronization with the WRF repository, and was responsibe
! for all aspects of WRF modification. UCD students Minjeong Kim and
! Volodymyr Kondratenko have contributed to the implementation of the
! fire propagation by the level set method.

! Refefences
!
! [1] Jan Mandel, Jonathan D. Beezley, and Adam K. Kochanski, "Coupled
! atmosphere-wildland fire modeling with WRF 3.3 and SFIRE 2011, 
! Geoscientific Model Development (GMD) 4, 591-610, 2011. 
! doi:10.5194/gmd-4-591-2011
!
! [2] T. L. Clark, J. Coen, and D. Latham, Description of a coupled 
! atmosphere-fire model, Intl. J. Wildland Fire, vol. 13, pp. 49-64, 
! 2004
!
! [3] Edward G. Patton and Janice L. Coen, WRF-Fire: A Coupled 
! Atmosphere-Fire Module for WRF, Preprints of Joint MM5/Weather 
! Research and Forecasting Model Users' Workshop, Boulder, CO, 
! June 22-25, 2004, pp. 221-223, NCAR
!
! ---------------------------------------------
!
! CURRENT ACTIVITY
! 
! For current activity and development trends please check out
! http://ccm.ucdenver.edu/wiki/User:Jmandel/blog
! http://www.openwfm.org/wiki/WRF-Fire_development_notes
! 

module module_fr_sfire_driver
! use this module for standalone call, you only need to provide some mock-up wrf modules  

use module_fr_sfire_model, only: sfire_model
use module_fr_sfire_phys, only: fire_params, init_fuel_cats, fuel_moisture, &
   advance_moisture, moisture_classes, &
   fire_rate_of_spread
use module_fr_sfire_atm, only: apply_windrf,interpolate_wind2fire_height,interpolate_atm2fire, &
   interpolate_z2fire,setup_wind_log_interpolation,find_trees_fmesh,massman_fwh
use module_fr_sfire_util
!use module_fr_sfire_util, only: lines_type,fire_max_lines
! Driver layer modules
#ifdef DM_PARALLEL
    USE module_dm        , ONLY : ntasks_x,ntasks_y,local_communicator,mytask,ntasks,wrf_dm_sum_reals
    USE module_comm_dm , ONLY : halo_fire_fuel1_sub, halo_fire_tign_sub, halo_fire_wind_f_sub, &
HALO_SFIRE_WIND_A_sub, halo_fire_ph_sub, halo_fire_zsf_sub, halo_fire_longlat_sub, &
halo_fire_phb_sub, halo_fire_z0_sub, halo_fire_lfn1_sub, HALO_FIRE_LFN_OUT_sub, &
HALO_FIRE_MAG_sub, HALO_FIRE_MFG_sub, halo_fire_ndwi_sub
#endif
use module_fr_sfire_atm, only: read_emissions_table, add_fire_emissions


! WRF dependencies
USE module_domain, only: domain
USE module_configure, only: grid_config_rec_type
use module_model_constants, only: reradius,    & ! 1/earth radiusw
                                  pi2            ! 2*pi

implicit none


private
public sfire_driver_em,use_atm_vars,set_flags, &
       set_fp_from_grid, fire_ignition_convert
public ifun_beg, ifun_step, ifun_end

logical:: use_atm_vars=.true.   !  interpolate wind from atm mesh and average output fluxes back
logical:: interpolate_long_lat=.true. ! get fxlong fxlat by interpolation

logical:: fmoist_run, fmoist_interp, fire_run  ! which kind of model to run overall

integer, parameter:: ifun_beg=1, ifun_step=3, ifun_end=6  

TYPE(lines_type):: ignition, hfx

contains

! to write debugging information
#define DEBUG_OUT

subroutine sfire_driver_em ( grid , config_flags                    & 
            ,time_step_start,dt                                     &
            ,fire_ifun_start,fire_ifun_end,tsteps                   &
            ,ids,ide, kds,kde, jds,jde                              &
            ,ims,ime, kms,kme, jms,jme                              &
            ,ips,ipe, kps,kpe, jps,jpe                              &
            ,ifds,ifde, jfds,jfde                                   &
            ,ifms,ifme, jfms,jfme                                   &
            ,ifps,ifpe, jfps,jfpe                                   &
            ,rho,z_at_w,dz8w                                        &
)

!*** purpose: driver from grid structure



    implicit none
!*** arguments
    TYPE(domain) , TARGET :: grid                             ! state 
    TYPE (grid_config_rec_type) , INTENT(IN)  :: config_flags ! namelist
    real, intent(in):: time_step_start, dt
    integer, intent(in)::     fire_ifun_start,fire_ifun_end,tsteps ! driver cycle controls
    integer, intent(in):: &
             ids,ide, kds,kde, jds,jde,                              &
             ims,ime, kms,kme, jms,jme,                              &
             ips,ipe, kps,kpe, jps,jpe,                              &
             ifds,ifde, jfds,jfde,                                   &
             ifms,ifme, jfms,jfme,                                   &
             ifps,ifpe, jfps,jfpe 
    real,dimension(ims:ime, kms:kme, jms:jme),intent(in), optional::rho,z_at_w,dz8w

!*** local
    integer::fire_ifun,ir,jr,istep,itimestep,i,ipe1,kpe1,jpe1,j
    logical::restart,replay
    real:: corner_ll,corner_ul,corner_ur,corner_lr, max_u, max_v, max_w, max_rho, min_rho
    character(len=128)msg, msg2
    type(fire_params)::fp
    real:: moisture_time

    logical:: run_advance_moisture,run_fuel_moisture, moisture_initializing
    real::    dt_moisture



!*** executable

    call sfire_debug_hook(config_flags%fire_debug_hook_sec)
    call time_start
    
!   if(fire_ifun_start.le.1)call print_id  ! print id only once, during initialization


! **** THE FOLLOWING REALLY SHOULD BE DONE ONCE NOT EVERY TIMESTEP

! populate our structures from wrf

    call set_fp_from_grid(grid,fp)

! copy configuration flags to sfire internal structures
    call set_flags(config_flags)

    write(msg,'(a,i3)')' domain', grid%id
    call message(msg,level=1)

! prevent memory errors because of variables missing in registry packages
    call check_grid_alloc(grid,config_flags)
         
#ifndef SFIRE_STANDALONE

!   see what we got from wrf
!! need to replace ipe by min(ide-1,ipe) and similarly jpe

    if(fire_print_msg.ge.2 .and. fire_ifun_start .gt. 1)then

      ipe1=min(ide-1,ipe)
      jpe1=min(jde-1,jpe)
      kpe1=kpe-1

      max_u=fun_real(REAL_AMAX,  &
        ims,ime,kms,kme,jms,jme, &                ! memory dims
        ids,ide,kds,kde,jds,jde, &                ! domain dims
        ips,ipe1,kps,kpe1,jps,jpe1, &                ! patch or tile dims
        1,0,0,       &                            ! staggering
        grid%u_2,grid%u_2)

      max_v=fun_real(REAL_AMAX,  &
        ims,ime,kms,kme,jms,jme, &                ! memory dims
        ids,ide,kds,kde,jds,jde, &                ! domain dims
        ips,ipe1,kps,kpe1,jps,jpe1, &                ! patch or tile dims
        0,0,1,       &                            ! staggering
        grid%v_2,grid%v_2)

      max_w=fun_real(REAL_AMAX,  &
        ims,ime,kms,kme,jms,jme, &                ! memory dims
        ids,ide,kds,kde,jds,jde, &                ! domain dims
        ips,ipe1,kps,kpe1,jps,jpe1, &                ! patch or tile dims
        0,1,0,       &                            ! staggering
        grid%w_2,grid%w_2)

      !write(msg,93)time_step_start,'Maximal u v w  wind',max_u,max_v,max_w,'m/s'
      !call message(msg,0)
      !write(msg,92)time_step_start,'Min and max rho    ',min_rho,max_rho,'kg/m^3'
      !call message(msg,0)

      write(msg,91)time_step_start,'Maximal u wind      ',max_u,'m/s'
      call message(msg,0)
      write(msg,91)time_step_start,'Maximal v wind      ',max_v,'m/s'
      call message(msg,0)
      write(msg,91)time_step_start,'Maximal w wind      ',max_w,'m/s'
      call message(msg,0)

      if (present(rho)) then

        max_rho=fun_real(REAL_MAX,  &
          ims,ime,kms,kme,jms,jme, &                ! memory dims
          ids,ide,kds,kde,jds,jde, &                ! domain dims
          ips,ipe1,kps,kpe1,jps,jpe1, &                ! patch or tile dims
          0,0,0,       &                            ! staggering
          rho,rho)
  
        min_rho=fun_real(REAL_MIN,  &
          ims,ime,kms,kme,jms,jme, &                ! memory dims
          ids,ide,kds,kde,jds,jde, &                ! domain dims
          ips,ipe1,kps,kpe1,jps,jpe1, &                ! patch or tile dims
          0,0,0,       &                            ! staggering
          rho,rho)
  
  
        write(msg,91)time_step_start,'Minimal rho         ',min_rho,'kg/m^3'
        call message(msg,0)
        write(msg,91)time_step_start,'Maximal rho         ',max_rho,'kg/m^3'
        call message(msg,0)
 
      endif


93    format('Time ',f11.3,' s ',a,3e12.3,1x,a)
92    format('Time ',f11.3,' s ',a,2e12.3,1x,a)
91    format('Time ',f11.3,' s ',a,e12.3,1x,a)


    endif
#endif


    ! refinement r
    ir=grid%sr_x ! refinement ratio
    jr=grid%sr_y
    write(msg,'(a,2i4)')'fire mesh refinement ratios', ir,jr
    call message(msg)
    if(ir.le.0.or.jr.le.0)then
        call crash('fire mesh refinement ratio must be positive')
    endif

call print_2d_stats(ifps,min(ifpe,ifde-ir),jfps,min(jfpe,jfde-jr),ifms,ifme,jfms,jfme, grid%tign_g,'sfire_driver_em: grid%tign_g')
call print_2d_stats(ifps,min(ifpe,ifde-ir),jfps,min(jfpe,jfde-jr),ifms,ifme,jfms,jfme, grid%nfuel_cat,'sfire_driver_em: grid%nfuel_cat')

    itimestep=grid%itimestep
    restart=config_flags%restart .or. config_flags%cycling .or. config_flags%fire_restart ! skip state initialization
    replay= time_step_start+dt .le.  config_flags%fire_perimeter_time
94  format('Time step',i11,' from',f11.3,' to',f11.3,' perimeter_time',f11.3,' setting replay ',l1)
    write(msg,94)itimestep,time_step_start,time_step_start+dt,config_flags%fire_perimeter_time,replay
    call message(msg)
95  format('namelist.input restart ',l1,' cycling ',l1,' fire_restart ',l1,' setting restart ',l1)
    write(msg,95)config_flags%restart,config_flags%cycling,config_flags%fire_restart,restart
    call message(msg)
    
    
    ! **** moisture model

    ! decide what to run - moisture, interpolation, or fire model itself
    fmoist_run    = config_flags%fmoist_run
    fmoist_interp = config_flags%fmoist_interp 
    if(fire_fmc_read.ne.0.and.fmoist_run)call crash('fmoist_run=T requires fire_fmc_read=0')
    fire_run = .not. config_flags%fmoist_only

    !decide what to run
    moisture_time = time_step_start
    run_advance_moisture = .false. ! default
    run_fuel_moisture = .false. ! default
    moisture_initializing = fire_ifun_start < 3
    
    
    
    if(fmoist_run)then
        if(moisture_initializing)then
            if(fire_ifun_end>2)call crash('initialization must be run separately')
            grid%fmoist_lasttime=moisture_time ! initialize the last time the model has run to start of run
            grid%fmoist_nexttime=moisture_time 
            call message('moisture initialization')
            run_advance_moisture = .true.
        else ! regular timestep
            if(config_flags%fmoist_freq > 0)then  ! regular timestep. go by multiples?
                if(mod(grid%itimestep,config_flags%fmoist_freq) .eq. 0)then
                    write(msg,'(a,i10,a,i10)')'moisture model runs because timestep ',grid%itimestep,' is a multiple of ',config_flags%fmoist_freq
                    call message(msg)
                    run_advance_moisture = .true.
                endif
            else
                if(.not. moisture_time  < grid%fmoist_nexttime) then ! no, by time interval
                    write(msg,'(a,f12.2,a)')'moisture model runs because time ',grid%fmoist_nexttime,'s has arrived'
                    call message(msg)
                    run_advance_moisture = .true.
                endif
            endif
            if(run_advance_moisture)then ! decide on timing
                dt_moisture  = moisture_time - grid%fmoist_lasttime  ! Time since moisture model run the last time. Should be long.
                grid%fmoist_lasttime = moisture_time
                if(config_flags%fmoist_freq > 0)then
                    write(msg,'(a,f12.2,a,i10,a)')'moisture time step is ',dt_moisture,'s running every ',config_flags%fmoist_freq,' steps'
                    call message(msg)
                else
                    grid%fmoist_nexttime = moisture_time + config_flags%fmoist_dt
                    write(msg,'(a,f12.2,a,f12.2,a)')'moisture time step is ',dt_moisture,'s next run at ',grid%fmoist_nexttime,'s'
                    call message(msg)
                endif
                if(fmoist_interp)then
                    call message('moisture interpolation to fuels will run because moisture model does')
                    run_fuel_moisture=.true.
                endif
            endif
        endif
    elseif(itimestep.eq.1.and.fmoist_interp)then
            call message('initializing, moisture interpolation to fuels will run from input data')
            run_fuel_moisture=.true.
    endif

!$OMP CRITICAL(SFIRE_DRIVER_CRIT)
    write(msg,'(a,i1,a,i1,a,l1,a,l1)') &
       'sfire_driver_em: ifun from ',fire_ifun_start,' to ',fire_ifun_end, &
           ' restart=',restart,' replay=',replay
!$OMP END CRITICAL(SFIRE_DRIVER_CRIT)
    call message(msg)

    do istep=0,tsteps ! istep >0 is for testing only, exit after the first call
      itimestep = grid%itimestep + istep ! in the first call, do fire_test_steps steps of the fire model

      do fire_ifun=fire_ifun_start,fire_ifun_end

        write(msg,'(a,i1,a)')'*** stage ',fire_ifun,' ***'
        call message(msg)

        ! 1 = moisture_initialize run pass 1: interpolate height to zsf=terrain
        ! 2 = initialize run pass 2: set fuel data, terrain gradient
        ! 3 = initialize timestep: interpolate winds, check for ignition, time step on moisture model
        ! 4 = do one timestep
        ! 5 = copy timestep output to input
        ! 6 = compute output fluxes

#ifdef DM_PARALLEL

       if(fire_run)then

        if(fire_ifun.eq.1)then

!       halo exchange on topography
#include "HALO_FIRE_LONGLAT.inc"
!!            if(fire_topo_from_atm.eq.1)then
!!#include "HALO_FIRE_HT.inc"
!!            endif 
! base geopotential and roughness
#include "HALO_FIRE_PHB.inc"
#include "HALO_FIRE_Z0.inc"
        if(kfmc_ndwi > 0 .and. fndwi_from_ndwi .eq.1)then
#include "HALO_FIRE_NDWI.inc"
        endif

        elseif(fire_ifun.eq.2)then
!           halo exchange on zsf width 2
#include "HALO_FIRE_ZSF.inc"

            if(replay)then
                call message('replay, halo exchange on lfn and tign')
#include "HALO_FIRE_LFN1.inc"
#include "HALO_FIRE_TIGN.inc"
            endif

        elseif(fire_ifun.eq.3)then
!           halo exchange on atm winds and geopotential for interpolation
#include "HALO_SFIRE_WIND_A.inc"
#include "HALO_FIRE_PH.inc"

        elseif(fire_ifun.eq.4)then
!           halo exchange on fire winds width 2 for a 2-step RK method
#include "HALO_FIRE_WIND_F.inc"

            if(run_fuel_moisture)then
            ! have interpolated to the fire grid
#include "HALO_FIRE_MFG.inc"
            endif
 
        elseif(fire_ifun.eq.5)then
#include "HALO_FIRE_LFN_OUT.inc"

        elseif(fire_ifun.eq.6)then
!           computing fuel_left needs ignition time from neighbors
            call message('halo exchange on lfn width 2 and tign')
#include "HALO_FIRE_TIGN.inc"
#include "HALO_FIRE_LFN1.inc"

        endif ! fire_ifun
       endif ! fire_run
#endif

        if(fire_ifun.eq.2)then
            write(msg,'(a,i4,a,i4)') 'chem_opt=', config_flags%chem_opt, &
               ' tracer_opt=',config_flags%tracer_opt
            call message(msg)
            if(config_flags%chem_opt>0 .or. config_flags%tracer_opt > 0)then
                ! need reading fuel categories first
                call read_emissions_table(config_flags%chem_opt,config_flags%tracer_opt)
            endif ! config_flags
        endif  ! fire_ifun
        
        ! print *,'dt: ',dt,grid%dt,' diff ', dt-grid%dt
        ! need domain by 1 smaller, in last row.col winds are not set properly
        call sfire_driver_phys ( &
            fire_ifun,                  &
            config_flags,                                           &
            ids,ide-1, kds,kde, jds,jde-1,                          &
            ims,ime, kms,kme, jms,jme,                          &
            ips,min(ipe,ide-1), kps,kpe, jps,min(jpe,jde-1),                          & 
            ifds,ifde-ir, jfds,jfde-jr,                    &
            ifms,ifme, jfms,jfme,                    &
            ifps,min(ifpe,ifde-ir), jfps,min(jfpe,jfde-jr),      &
            ir,jr,                                      & ! atm/fire grid ratio
            grid%num_tiles,                             & ! atm grid tiling
            grid%i_start,min(grid%i_end,ide-1),                    &
            grid%j_start,min(grid%j_end,jde-1),                    &                 
            itimestep,restart,replay,config_flags%fire_fuel_read,config_flags%fire_fuel_cat, &  ! in scalars
            time_step_start,dt,grid%dx,grid%dy,                    &
            grid%u_frame,grid%v_frame,                  &
            config_flags%fire_ext_grnd,config_flags%fire_ext_crwn,config_flags%fire_crwn_hgt, &
            ignition,hfx,              &  ! lines
            grid%u_2,grid%v_2,           &          ! atm arrays in
            grid%ph_2,grid%phb,               & ! geopotential
            grid%z0,                        & ! roughness height
            grid%ht,                        &                         ! terrain height
            grid%xlong,grid%xlat,                         & ! coordinates of atm grid centers, for ignition location           
            grid%tign_in,                                 &
            grid%lfn,grid%tign_g,grid%fuel_frac,          & ! state arrays, fire grid
            grid%fire_area,                               & ! redundant, for display, fire grid
            grid%fuel_frac_burnt,                         &
            grid%lfn_out,                                 & ! work - one timestep    
            grid%avg_fuel_frac,                           & ! out redundant arrays, atm grid
            grid%grnhfx,grid%grnqfx,grid%canhfx,grid%canqfx, & ! out redundant arrays, atm grid
            grid%uah,grid%vah,                            &
            grid%fgrnhfx,grid%fgrnqfx,grid%fcanhfx,grid%fcanqfx, & ! out redundant arrays, atm grid
            grid%ros,grid%flineint,grid%flineint2,         & ! diagnostic variables
            grid%f_ros0,grid%f_rosx,grid%f_rosy,grid%f_ros,& ! fire risk spread 
            grid%f_int,grid%f_lineint,grid%f_lineint2,     & ! fire risk intensities 
            grid%f_ros11,grid%f_ros12,grid%f_ros13,grid%f_ros21,  & ! fire spread in nodal directions
            grid%f_ros23,grid%f_ros31,grid%f_ros32,grid%f_ros33,  & ! fire spread in nodal directions
            grid%fxlong,grid%fxlat,                           &       
            grid%fire_hfx,                                & !
            grid%nfuel_cat,                               & ! input, or internal for safekeeping
            grid%fuel_time,                      & 
            grid%wz0,                            &          !wrf roughness length on fire mesh
            grid%fz0, grid%fwh,grid%can_top,grid%cuf,grid%cvf, &
            fp,                                    & ! structure with pointers passed to spread rate calculation
            config_flags%nfmc,         & ! moisture model variables start
            run_advance_moisture,run_fuel_moisture,dt_moisture,     &    ! moisture model control
            config_flags%fmep_decay_tlag,                               & ! moisture extended model assim. diffs decay time lag
            grid%rainc, grid%rainnc,                       & ! accumulated rain from different sources
            grid%t2, grid%q2, grid%psfc,               & ! temperature (K), vapor contents (kg/kg), pressure (Pa) at the surface
            grid%rain_old,                   & ! previous value of accumulated rain
            grid%t2_old, grid%q2_old, grid%psfc_old,   & ! previous values of the atmospheric state at surface
            grid%rh_fire,                    & ! relative humidity, diagnostics
            grid%fmc_gc,                      & ! fuel moisture fields updated, by class, assumed set to something reasonable
            grid%fmep,                      & ! fuel moisture extended model parameters
            grid%fmc_equi,                      & ! fuel moisture fields updated, by class, equilibrium diagnostic 
            grid%fmc_lag,                      & ! fuel moisture fields updated, by class, tendency diagnostic
            fp%fmc_g, &                             !  write-only alias. need to exit before using fp again
            grid%ndwi, &
            grid%fndwi)
 

#ifdef DM_PARALLEL
          if(fire_run)then
              if(fire_ifun.eq.2)then
!                 halo exchange on all fuel data width 2
#include "HALO_FIRE_FUEL1.inc"
!                 fire state was initialized
                  call message('halo exchange on lfn width 2')
#include "HALO_FIRE_LFN1.inc"
              endif
              if(run_fuel_moisture)then
                  if(fire_ifun.eq.3)then
                     ! prepare for interpolation to the fire grid
#include "HALO_FIRE_MAG.inc"
                  endif
              endif
          endif
#endif

        if(fire_ifun.eq.2)then 
!           store mesh units computed in fire_ignition_convert in stage 2
            grid%unit_fxlong = ignition%unit_fxlong
            grid%unit_fxlat = ignition%unit_fxlat
        endif  ! fire_ifun

     if(fire_ifun.eq.6)then
         if(config_flags%chem_opt>0 .or. config_flags%tracer_opt>0)then
             if(.not.(present(rho).and.present(dz8w)))then
                 call crash('sfire_driver_em: must have rho and dz8w to call add_fire_emissions')
             endif
             call add_fire_emissions( &
                 config_flags%chem_opt,config_flags%tracer_opt,dt,grid%dx,grid%dy,      &
                 ifms,ifme,jfms,jfme,    &
                 ifps,ifpe,jfps,jfpe,    &               ! use patch instead of tile
                 ids,ide,kds,kde,jds,jde,          &
                 ims,ime,kms,kme,jms,jme,          &
                 ips,ipe,kps,kpe,jps,jpe,          &
                 rho,dz8w,                         & ! from atmosphere state
                 grid%fgip, grid%fuel_frac_burnt, grid%nfuel_cat, & ! from fire state
                 grid%chem,grid%tracer)              ! update/output
          endif
     endif
                

      enddo
    enddo

    if(tsteps>0)call crash('sfire_driver_em: test run of uncoupled fire model completed')
    call time_end('sfire')

end subroutine sfire_driver_em

!
!*******************
!

subroutine sfire_driver_phys (ifun,    &
    config_flags,                                 & ! namelist.input
    ids,ide, kds,kde, jds,jde,                    & ! atm grid dimensions
    ims,ime, kms,kme, jms,jme,                    &
    ips,ipe, kps,kpe, jps,jpe,                    &
    ifds, ifde, jfds, jfde,                       & ! fire grid dimensions
    ifms, ifme, jfms, jfme,                       &
    ifps, ifpe, jfps, jfpe,                       & ! fire patch in - will use smaller
    ir,jr,                                        & ! atm/fire grid ratio
    num_tiles,i_start,i_end,j_start,j_end,        & ! atm grid tiling
    itimestep,restart,replay,ifuelread,nfuel_cat0,       & ! in scalars
    time_step_start,dt,dx,dy,                     & ! in scalars
    u_frame,v_frame,                              &
    fire_ext_grnd,fire_ext_crwn,fire_crwn_hgt,    &
    ignition,hfx,                                 & ! lines
    u,v,                                          & ! in arrays, atm grid
    ph,phb,                                       &
    z0,zs,                                        & 
    xlong,xlat,                                   &
    tign_in,                                      &
    lfn,tign,fuel_frac,                           & ! state arrays, fire grid
    fire_area,                                    & ! redundant state, fire grid
    fuel_frac_burnt,                              & 
    lfn_out,                                      & ! out level set function    
    avg_fuel_frac,                                &
    grnhfx,grnqfx,canhfx,canqfx,                  & ! out redundant arrays, atm grid  
    uah,vah,                                      & ! out atm grid
    fgrnhfx,fgrnqfx,fcanhfx,fcanqfx,              & ! out redundant arrays, fire grid
    ros,flineint,flineint2,                       & ! diagnostic variables
    f_ros0,f_rosx,f_rosy,f_ros,                   & ! fire risk spread 
    f_int,f_lineint,f_lineint2,                   & ! fire risk intensities 
    f_ros11,f_ros12,f_ros13,f_ros21,  & ! fire spread in nodal directions
    f_ros23,f_ros31,f_ros32,f_ros33,  & ! fire spread in nodal directions
    fxlong,fxlat,                                 & !  
    fire_hfx,                                     & ! 
    nfuel_cat,                                    & ! in array, data, fire grid, or constant internal
    fuel_time,                                    & ! save constant internal data, fire grid
    wz0,                                          &
    fz0,fwh,can_top,cuf,cvf,                  &
    fp,                                           & ! fire params
    nfmc,                                     & ! number of fuel moisture classes
    run_advance_moisture,run_fuel_moisture,dt_moisture,& ! moisture model control
    fmep_decay_tlag,                              & ! moist. extended model assim. diffs time lag
    rainc,rainnc,               & ! accumulated rain from different sources
    t2, q2, psfc,               & ! temperature (K), vapor contents (kg/kg), pressure (Pa) at the surface
    rain_old,                   & ! previous value of accumulated rain
    t2_old, q2_old, psfc_old,   & ! previous values of the atmospheric state at surface
    rh_fire,                    & ! relative humidity, diagnostics
    fmc_gc,                     &  ! fuel moisture fields updated, by class, assumed set to something reasonable
    fmep,                       &  ! fuel moisture extended model parameters
    fmc_equi,                   &  ! fuel moisture fields updated, by class equilibrium diagnostic
    fmc_lag,                    &  ! fuel moisture fields updated, by class tendency diagnostic
    fmc_g,                      &  ! fuel moisture, alias of fp%fmc_g
    ndwi,                       &
    fndwi)                         ! ndwi on fire grid


implicit none

!*** arguments

TYPE (grid_config_rec_type) , INTENT(IN)  :: config_flags ! namelist, needed for ignion 
integer, intent(in)::ifun,                        &
    ids,ide, kds,kde, jds,jde,                    & ! atm domain bounds
    ims,ime, kms,kme, jms,jme,                    & ! atm memory bounds 
    ips,ipe, kps,kpe, jps,jpe,                    & ! atm patch bounds
    ifds, ifde, jfds, jfde,                       & ! fire domain bounds
    ifms, ifme, jfms, jfme,                       & ! fire memory bounds
    ifps, ifpe, jfps, jfpe,                       & ! fire patch bounds
    ir,jr,                                        & ! atm/fire grid refinement ratio
    nfmc,                                     & ! number of fuel moisture classes
    itimestep,                                    & ! number of this timestep
    ifuelread,                                    & ! how to initialize nfuel_cat:
                                                       ! -1=not at all, done outside 
                                                       ! 0=from nfuel_cat0
                                                       ! 1=from altitude
                                                       ! 2=from file
    nfuel_cat0,                                   & ! fuel category to initialize everything to
    num_tiles                                       ! number of tiles

logical, intent(in)::restart,replay
    
integer,dimension(num_tiles),intent(in) :: i_start,i_end,j_start,j_end  ! atm grid tiling

real, intent(in):: &
    time_step_start,                              & ! time step start
    dt,                                           & ! time step length
    dx,dy,                                        & ! atm grid step
    u_frame,v_frame,                              & ! velocity offset
    fire_crwn_hgt,                                & ! lowest height crown fire heat is released (m)
    fire_ext_grnd,                                & ! extinction depth of ground fire heat (m)
    fire_ext_crwn                                   ! and for the canopy (m) 


TYPE (lines_type), intent(inout):: ignition,hfx

real,intent(in),dimension(ims:ime,kms:kme,jms:jme)::u,v, & ! wind velocity (m/s) (staggered atm grid) 
                              ph, phb                      ! geopotential (w-points atm grid)
real,intent(in),dimension(ims:ime, jms:jme)::   z0, &    ! roughness height
                                                zs       ! terrain height  
real,intent(out),dimension(ims:ime,jms:jme)::&
    uah,                                           & ! atm wind at fire_wind_height, diagnostics
    vah                                              ! atm wind at fire_wind_height, diagnostics

real, dimension(ims:ime, jms:jme), intent(inout)::xlong, xlat, ndwi ! inout because of extension at bdry
    
real, intent(inout), dimension(ifms:ifme,jfms:jfme):: & ! fuel data; can be also set inside (cell based, fire grid)
    fz0,fwh,can_top,cuf,cvf,wz0,                            &  
    nfuel_cat,fndwi                                      

real, intent(inout), dimension(ifms:ifme, jfms:jfme)::     &
    tign_in,                                   &
    lfn,tign,fuel_frac,                        &     ! state: level function, ign time, fuel left
    lfn_out                                    ! fire wind velocities

real, intent(out), dimension(ifms:ifme, jfms:jfme)::  &
    fire_area, &                               ! fraction of each cell burning
    fuel_frac_burnt

real, intent(out), dimension(ims:ime, jms:jme):: &  ! redundant arrays, for display purposes only (atm grid)
    avg_fuel_frac,                               &  ! average fuel fraction
    grnhfx,                                      &  ! heat flux from ground fire (W/m^2) 
    grnqfx,                                      &  ! moisture flux from ground fire (W/m^2) 
    canhfx,                                      &  ! heat flux from crown fire (W/m^2) 
    canqfx                                         ! moisture flux from crown fire (W/m^2) 

real, intent(out), dimension(ifms:ifme, jfms:jfme):: &  ! redundant arrays, for display only, fire grid
    fgrnhfx,                                      &  ! heat flux from ground fire (W/m^2) 
    fgrnqfx,                                      &  ! moisture flux from ground fire (W/m^2) 
    fcanhfx,                                      &  ! heat flux from crown fire (W/m^2) 
    fcanqfx,                                      &  ! moisture flux from crown fire (W/m^2) 
    ros,flineint,flineint2,                       & ! diagnostic variables
    f_ros0,f_rosx,f_rosy,f_ros,                   & ! fire risk spread 
    f_int,f_lineint,f_lineint2,                   &  ! fire risk intensities 
    f_ros11,f_ros12,f_ros13,f_ros21,  & ! fire spread in nodal directions
    f_ros23,f_ros31,f_ros32,f_ros33     ! fire spread in nodal directions

    
! moisture model arguments
logical, intent(in)::run_advance_moisture,run_fuel_moisture
real, intent(in)::dt_moisture
real, intent(in)::fmep_decay_tlag
real, intent(in), dimension(ims:ime,jms:jme):: t2, q2, psfc, rainc, rainnc
real, intent(inout), dimension(ims:ime,jms:jme):: t2_old, q2_old, psfc_old, rain_old 
real, intent(out),dimension(ims:ime,jms:jme):: rh_fire
real, intent(inout), dimension(ims:ime,nfmc,jms:jme):: fmc_gc
real, intent(inout), dimension(ims:ime,2,jms:jme):: fmep
real, intent(out), dimension(ims:ime,nfmc,jms:jme):: fmc_equi,fmc_lag
real, intent(inout), dimension(ifms:ifme,jfms:jfme):: fmc_g



!  ***** data (constant in time) *****

real, dimension(ifms:ifme, jfms:jfme), intent(inout)::fxlong,fxlat, & ! fire mesh coordinates
    fire_hfx
real, intent(out), dimension(ifms:ifme, jfms:jfme)::fuel_time   ! fire params arrays

type(fire_params),intent(inout)::fp
    
!*** local
real :: dxf,dyf,time_start,latm, s
integer :: its,ite,jts,jte,kts,kte, &            ! tile
    ij,i,j,k,id,pid,ipe1,jpe1,ite1,jte1, &
    ii,jj,                                 &
    ifts,ifte,jfts,jfte                          ! fire tile
character(len=128)::msg
character(len=3)::kk
real, parameter::zero=0.

!*** executable

! time - assume dt does not change
! time_start = (itimestep-1) * dt     ! timestep 1 starts at 0
! print *,'time_start: ',time_start,time_step_start,' diff ', time_start-time_step_start
time_start = time_step_start ! use the time passed from wrf

! fire mesh step
dxf=dx/ir
dyf=dy/jr


!$OMP CRITICAL(SFIRE_DRIVER_CRIT)
write(msg,'(a,i5)')'sfire_driver_phys stage ',ifun
call message(msg)
write(msg,'(a,2f15.6)')'atmosphere mesh step:',dx,dy
call message(msg)
write(msg,'(a,2f15.6)')'fire mesh step:      ',dxf,dyf
call message(msg)
write(msg,7001)'atm domain      ','ids',ids,ide,jds,jde
call message(msg)                    
write(msg,7001)'atm memory      ','ims',ims,ime,jms,jme
call message(msg)                    
write(msg,7001)'atm patch       ','ips',ips,ipe,jps,jpe
call message(msg)                    
write(msg,7001)'fire domain     ','ifds',ifds,ifde,jfds,jfde
call message(msg)                    
write(msg,7001)'fire memory     ','ifms',ifms,ifme,jfms,jfme
call message(msg)                    
write(msg,7001)'fire patch      ','ifps',ifps,ifpe,jfps,jfpe
call message(msg)                    
!$OMP END CRITICAL(SFIRE_DRIVER_CRIT)

! check mesh dimensions
call check_fmesh(ids,ide,ifds,ifde,ir,'id')           ! check if atm and fire grids line up
call check_fmesh(jds,jde,jfds,jfde,jr,'jd')
call check_fmesh(ips,ipe,ifps,ifpe,ir,'ip')
call check_fmesh(jps,jpe,jfps,jfpe,jr,'jp')
call check_mesh_2dim(ips,ipe,jps,jpe,ims,ime,jms,jme)        ! check if atm patch fits in atm array
call check_mesh_2dim(ifps,ifpe,jfps,jfpe,ifms,ifme,jfms,jfme) ! check if fire patch fits in fire array
call check_mesh_2dim(ips,ipe,jps,jpe,ids,ide,jds,jde)        ! check if atm patch fits in atm domain
call check_mesh_2dim(ifps,ifpe,jfps,jfpe,ifds,ifde,jfds,jfde) ! check if fire patch fits in fire domain

pid=0
if(fire_print_file.gt.0)then
    if(itimestep.le.fire_print_file.or.mod(itimestep,fire_print_file).eq.0)pid=itimestep ! print 1-fire_print_file then every fire_print_file-th
endif


if(ifun.eq.1)then
         call init_fuel_cats(fmoist_run .or. fmoist_interp) ! properties of fuel categories and moisture classes from namelist.fire
endif

if(ifun.eq.2)then
 call print_chsum(itimestep,ims,ime,kms,kme,jms,jme,ids,ide,kds,kde,jds,jde,ips,ipe,kps,kpe,jps,jpe,1,0,0,u,'u')
 call print_chsum(itimestep,ims,ime,kms,kme,jms,jme,ids,ide,kds,kde,jds,jde,ips,ipe,kps,kpe,jps,jpe,0,0,1,v,'v')
 call print_chsum(itimestep,ims,ime,kms,kme,jms,jme,ids,ide,kds,kde,jds,jde,ips,ipe,kps,kpe,jps,jpe,0,1,0,ph,'ph')
endif ! ifun.eq.3

call print_chsum(itimestep,ifms,ifme,1,1,jfms,jfme,ifds,ifde,1,1,jfds,jfde,ifps,ifpe,1,1,jfps,jfpe,0,0,0,lfn,'lfn')
call print_chsum(itimestep,ifms,ifme,1,1,jfms,jfme,ifds,ifde,1,1,jfds,jfde,ifps,ifpe,1,1,jfps,jfpe,0,0,0,tign,'tign')


! fake atm tile bounds
kts=kps
kte=kpe

! staggered atm patch bounds
ipe1=ifval(ipe.eq.ide,ipe+1,ipe)
jpe1=ifval(jpe.eq.jde,jpe+1,jpe)

! set up fire tiles & interpolate to fire grid
!$OMP PARALLEL DO PRIVATE(ij,its,ite,jts,jte,ite1,jte1,ifts,ifte,jfts,jfte,msg,id) &
!$OMP SCHEDULE(STATIC)
do ij=1,num_tiles

    id = ifval(pid.ne.0,pid+(ij-1)*10000,0) ! for print

    ! set up tile bounds    
    its = i_start(ij)  ! start atmospheric tile in i
    ite = i_end(ij)    ! end atmospheric tile in i
    jts = j_start(ij)  ! start atmospheric tile in j
    jte = j_end(ij)    ! end atmospheric tile in j
    ifts= (its-ids)*ir+ifds       ! start fire tile in i
    ifte= (ite-ids+1)*ir+ifds-1   ! end fire tile in i
    jfts= (jts-jds)*jr+jfds       ! start fire tile in j
    jfte= (jte-jds+1)*jr+jfds-1   ! end fire tile in j
        
! staggered atm tile bounds
    ite1=ifval(ite.eq.ide,ite+1,ite)
    jte1=ifval(jte.eq.jde,jte+1,jte)

!$OMP CRITICAL(SFIRE_DRIVER_CRIT)
    write(msg,'(a,i3,1x,a,i7,1x,a,i3)')'tile=',ij,' id=',id,' ifun=',ifun
    call message(msg)
    write(msg,7001)'atm tile   ','its',its,ite,jts,jte
    call message(msg)                   
    write(msg,7001)'fire tile  ','ifts',ifts,ifte,jfts,jfte
    call message(msg)                    
!$OMP END CRITICAL(SFIRE_DRIVER_CRIT)

    ! check the tiles
    call check_mesh_2dim(its,ite,jts,jte,ips,ipe,jps,jpe)                 ! check if atm tile fits in atm patch
    call check_mesh_2dim(ifts,ifte,jfts,jfte,ifps,ifpe,jfps,jfpe)         ! check if fire tile fits in fire patch
    call check_mesh_2dim(ifts-2,ifte+2,jfts-2,jfte+2,ifms,ifme,jfms,jfme)! check if fire node tile fits in memory


!$OMP CRITICAL(SFIRE_DRIVER_CRIT)
    write(msg,'(a,i6,a,2(f15.6,a))')'time step',itimestep,' at',time_start,' duration',dt,'s'
    call message(msg)
    7001 format(a,' dimensions ',a4,':',i6,' to ',i6,' by ',i6,' to ',i6)
    write(msg,'(a,2i9)')'refinement ratio:',ir,jr
!$OMP END CRITICAL(SFIRE_DRIVER_CRIT)

    if(run_advance_moisture)then
      if(ifun.eq.3)then 
      
      ! one timestep of the moisture model
          call message('advance_moisture start')
          call advance_moisture(    &
              itimestep.eq.1,             & ! initialize?
              ims,ime,  jms,jme,          & ! memory dimensions
              its,ite,  jts,jte,          & ! tile dimensions
              nfmc,                       & ! number of moisture fields
              dt_moisture,                & ! moisture model time step
              fmep_decay_tlag,            & ! moisture extended model assim. diffs decay tlag
              rainc, rainnc,              & ! accumulated rain 
              t2, q2, psfc,               & ! temperature (K), vapor contents (kg/kg), pressure (Pa) at the surface
              rain_old,                   & ! previous value of accumulated rain
              t2_old, q2_old, psfc_old,   & ! previous values of the atmospheric state at surface
              rh_fire,                    & ! relative humidity, diagnostics
              fmc_gc,                     & ! fuel moisture fields updated, by class, assumed set to something reasonable
              fmep,                       & ! fuel moisture extended model parameters
              fmc_equi,                   & ! fuel moisture fields updated, by class equilibrium diagnostic
              fmc_lag                     & ! fuel moisture fields updated, by class tendency diagnostic
          )
          call message('advance_moisture end')
      endif
    endif
        
    if(fire_run)then

     if(ifun.eq.2)then   

      if(restart)then
          
          call message('restart - interpolation skipped')

      else ! restart
          call message ('Interpolating static data from atmosphere to fire mesh')
        if(kfmc_ndwi > 0 .and. fndwi_from_ndwi .eq.1)then
            call print_2d_stats(ips,ipe,jps,jpe,ims,ime,jms,jme,fndwi,'driver:ndwi')
            call interpolate_z2fire(id,0,                 & ! for debug output, <= 0 no output, extend strip
                ids,ide,  jds,jde,                    & ! atm grid dimensions
                ims,ime,  jms,jme,                    &
                ips,ipe,jps,jpe,                              &
                its,ite,jts,jte,                              &
                ifds, ifde, jfds, jfde,                       & ! fire grid dimensions
                ifms, ifme, jfms, jfme,                       &
                ifts,ifte,jfts,jfte,                          &
                ir,jr,                                        & ! atm/fire grid ratio
                ndwi,                                       & ! atm grid arrays in
                fndwi)                                      ! fire grid arrays out
            call print_2d_stats(ifps,ifpe,jfps,jfpe,ifms,ifme,jfms,jfme,fndwi,'driver:fndwi')
         endif
!
!        call print_2d_stats(ips,ipe,jps,jpe,ims,ime,jms,jme,zs,'driver:zs')
!    
!        ! interpolate terrain height
!        if(fire_topo_from_atm.eq.1)then
!            call interpolate_z2fire(id,1,                 & ! for debug output, <= 0 no output
!                ids,ide,  jds,jde,                    & ! atm grid dimensions
!                ims,ime,  jms,jme,                    &
!                ips,ipe,jps,jpe,                              &
!                its,ite,jts,jte,                              &
!                ifds, ifde, jfds, jfde,                       & ! fire grid dimensions
!                ifms, ifme, jfms, jfme,                       &
!                ifts,ifte,jfts,jfte,                          &
!                ir,jr,                                        & ! atm/fire grid ratio
!                zs,                                       & ! atm grid arrays in
!                fp%zsf)                                      ! fire grid arrays out
!        else
!!$OMP CRITICAL(SFIRE_DRIVER_CRIT)
!           write(msg,'(a,i3,a)')'fire_topo_from_atm=',fire_topo_from_atm,' assuming ZSF set, interpolation skipped'
!!$OMP END CRITICAL(SFIRE_DRIVER_CRIT)
!        endif

#ifdef DEBUG_OUT
         call write_array_m(its,ite,jts,jte,ims,ime,jms,jme,xlat,'xlat',id)
         call write_array_m(its,ite,jts,jte,ims,ime,jms,jme,xlong,'xlong',id)
#endif
        if (interpolate_long_lat)then
         call message('Interpolating node longitude and latitude to fire mesh')
         call interpolate_z2fire(id,1,                 & ! for debug output, <= 0 no output
            ids,ide,  jds,jde,                    & ! atm grid dimensions
            ims,ime,  jms,jme,                    &
            ips,ipe,jps,jpe,                              &
            its,ite,jts,jte,                              &
            ifds, ifde, jfds, jfde,                       & ! fire grid dimensions
            ifms, ifme, jfms, jfme,                       &
            ifts,ifte,jfts,jfte,                          &
            ir,jr,                                        & ! atm/fire grid ratio
            xlat,                                       & ! atm grid arrays in
            fxlat)                                      ! fire grid arrays out

         call interpolate_z2fire(id,1,                 & ! for debug output, <= 0 no output
            ids,ide,  jds,jde,                    & ! atm grid dimensions
            ims,ime,  jms,jme,                    &
            ips,ipe,jps,jpe,                              &
            its,ite,jts,jte,                              &
            ifds, ifde, jfds, jfde,                       & ! fire grid dimensions
            ifms, ifme, jfms, jfme,                       &
            ifts,ifte,jfts,jfte,                          &
            ir,jr,                                        & ! atm/fire grid ratio
            xlong,                                       & ! atm grid arrays in
            fxlong)                                      ! fire grid arrays out
        endif !interpolate_long_lat 

        ! cannot initialize moisture model because T2 Q2 PSFC are not set yet

        ! not using zsf any more but anyway
        call print_2d_stats(ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme,fp%zsf,'driver_phys:zsf')        
 
     endif !restart

    elseif(ifun.eq.3)then  ! interpolate winds to the fire grid
 
      if(use_atm_vars)then                                  
     
        call write_array_m(its,ite,jts,jte,ims,ime,jms,jme,z0,'z0',id)
        call write_array_m3(its,ite1,kts,kde-1,jts,jte,ims,ime,kms,kme,jms,jme,u,'u_2',id)
        call write_array_m3(its,ite,kts,kde-1,jts,jte1,ims,ime,kms,kme,jms,jme,v,'v_2',id)
        call write_array_m3(its,ite,kts,kde,jts,jte,ims,ime,kms,kme,jms,jme,ph,'ph_2',id)
        call write_array_m3(its,ite,kts,kde,jts,jte,ims,ime,kms,kme,jms,jme,phb,'phb',id)
        
        if(fire_wind_log_interp.eq.4)then
          call print_2d_stats(its,ite,jts,jte,ims,ime,jms,jme,z0,'driver_phys:z0')
          call interpolate_atm2fire(id,                     & ! flag for debug output
            ids,ide, kds,kde, jds,jde,                    & ! atm grid dimensions
            ims,ime, kms,kme, jms,jme,                    &
            ips,ipe, jps,jpe,                             &
            its,ite,jts,jte,                              &                    
            ifds, ifde, jfds, jfde,                       & ! fire grid dimensions
            ifms, ifme, jfms, jfme,                       &
            ifps, ifpe, jfps, jfpe,                       & ! fire patch bounds
            ifts, ifte, jfts, jfte,                       &
            ir,jr,                                        & ! atm/fire grid ratio
            u_frame, v_frame,                             & ! velocity frame correction
            u,v,                                          & ! 3D atm grid arrays in
            ph,phb,                                       &
            z0,zs,                                        & ! 2D atm grid arrays in
            uah,vah,                                      & ! 2D atm grid out
            fp%vx,fp%vy)                                    ! fire grid arrays out

          call apply_windrf(                        &
            ifms,ifme,jfms,jfme,                    &
            ifts,ifte,jfts,jfte,                    &
            nfuel_cat,fp%vx,fp%vy)

            wz0 = 0             !Set WRF roughness length to 0 since we dont use this for option 4.

        else
          call print_2d_stats(ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme,fz0,'driver_phys:fz0')        
          call interpolate_wind2fire_height(id,       & ! to identify debugging prints and files if needed
            ids,ide, kds,kde, jds,jde,                    & ! atm grid dimensions
            ims,ime, kms,kme, jms,jme,                    &
            ips,ipe,jps,jpe,                              &
            its,ite,jts,jte,                              &
            ifds, ifde, jfds, jfde,                       & ! fire grid dimensions
            ifms, ifme, jfms, jfme,                       &
            ifps, ifpe, jfps, jfpe,                       & ! fire patch bounds
            ifts,ifte,jfts,jfte,                          &
            ir,jr,                                        & ! atm/fire grid ratio
            u_frame, v_frame,                             & ! velocity frame correction
            u,v,ph,phb,                                   & ! input atmospheric arrays
            fz0,fwh,                                      & ! input fire arrays
            fp%vx,fp%vy)                                          ! output fire arrays

          if(fire_use_windrf.eq.1)then
            call apply_windrf(                      &
            ifms,ifme,jfms,jfme,                    &
            ifts,ifte,jfts,jfte,                    &
            nfuel_cat,fp%vx,fp%vy)
          endif


          !If we are deploying the Massman parameterization, we need to interpolate our winds
          !down to the canopy height instead of the FWH. Here we will need to call the interpolate_wind2can
          !subroutine, for our specified canopy height assigned in the namelist (not yet implemented, 
          !this is just going to be hard coded for now). Once this is done, we can use the massman 
          !parameterization to scale the winds down to the FWH.

          if(fire_use_windrf.eq.4)then
          
            call message('fire_use_windrf option 4 has been selected....!')
            write(msg, '(A20,I1)') 'fire_can_top_read = ', fire_can_top_read
            call message(msg)

            !Lets set the canopy height, using either a look up table, or just explictly setting it 
            !using WRF input data from LANDFIRE.

            if(fire_can_top_read == 0)then

                can_top = 0                             ! Initialize canopy top array on fire mesh 

                call find_trees_fmesh(id,             & ! to identify debugging prints and files if needed
                ids,ide, kds,kde, jds,jde,                    & ! atm grid dimensions
                ims,ime, kms,kme, jms,jme,                    &
                ips,ipe,jps,jpe,                              &
                its,ite,jts,jte,                              &
                ifds, ifde, jfds, jfde,                       & ! fire grid dimensions
                ifms, ifme, jfms, jfme,                       &
                ifps, ifpe, jfps, jfpe,                       & ! fire patch bounds
                ifts,ifte,jfts,jfte,                          &
                ir,jr,                                        & ! atm/fire grid ratio
                can_top, nfuel_cat)                             ! tree canopy height info 
            
            end if 


            !Compute canopy top winds. Save this variable to cuf and cvf, which is like
            !vf and uf in terms of dimensionality, except that its canopy top winds.

            call interpolate_wind2fire_height(id,         & ! to identify debugging prints and files if needed
            ids,ide, kds,kde, jds,jde,                    & ! atm grid dimensions
            ims,ime, kms,kme, jms,jme,                    &
            ips,ipe,jps,jpe,                              &
            its,ite,jts,jte,                              &
            ifds, ifde, jfds, jfde,                       & ! fire grid dimensions
            ifms, ifme, jfms, jfme,                       &
            ifps, ifpe, jfps, jfpe,                       & ! fire patch bounds
            ifts,ifte,jfts,jfte,                          &
            ir,jr,                                        & ! atm/fire grid ratio
            u_frame, v_frame,                             & ! velocity frame correction
            u,v,ph,phb,                                   & ! input atmospheric arrays
            wz0,can_top,                                  & ! input fire arrays
            cuf,cvf)                                        ! output fire arrays


            !Lets call our massman parameterization and scale winds winds within 'forests'
            !within 'forests' using this scheme. uf and vf values that fall within 
            !forests are modifed accordingly, which will change the fire ROS.

            call massman_fwh(id,                          & ! to identify debugging prints and files if needed   
            ids,ide, kds,kde, jds,jde,                    & ! atm grid dimensions
            ims,ime, kms,kme, jms,jme,                    &
            ips,ipe,jps,jpe,                              &
            its,ite,jts,jte,                              &
            ifds, ifde, jfds, jfde,                       & ! fire grid dimensions
            ifms, ifme, jfms, jfme,                       &
            ifps, ifpe, jfps, jfpe,                       & ! fire patch bounds
            ifts,ifte,jfts,jfte,                          &
            ir,jr,                                        & ! atm/fire grid ratio
            cuf,cvf,                                      & ! canopy top winds (input)
            nfuel_cat,                                    & ! fire fuel category (input)
            can_top,                                      & ! canopy top (input)
            fp%vx,fp%vy)

          endif
        endif
      endif

    elseif(ifun.eq.4)then
    
      ! interpolate and compute weighted average to get the fuel moisture
      !! print *,'ifun=4, run_fuel_moisture=',run_fuel_moisture
      if(run_fuel_moisture)then
        call message('fuel_moisture start')
        call fuel_moisture(                &
        id,                                     & ! for prints and maybe file names
        nfmc,                                &
        ids,ide, jds,jde,               & ! atm grid dimensions
        ims,ime, jms,jme,           &
        ips,ipe,jps,jpe,                &
        its,ite,jts,jte,                     &
        ifds, ifde, jfds, jfde,         & ! fire grid dimensions
        ifms, ifme, jfms, jfme,     &
        ifts,ifte,jfts,jfte,                 &
        ir,jr,                                   & ! atm/fire grid ratio
        nfuel_cat,                         & ! fuel data
        fndwi,                             & ! satellite sensing interpolated on fire grid
        fmc_gc,                             & ! moisture contents by class on atmospheric grid
        fmc_g                                & ! weighted fuel moisture contents on fire grid
        )
        call message('fuel_moisture end')
      endif

      call print_2d_stats(ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme,fmc_g,'driver_phys:fmc_g')        
    endif

      

    call sfire_model (id,ifun,restart,replay,  &
        run_fuel_moisture,                      & ! if fuel moisture needs to be updated
        ifuelread,nfuel_cat0,                   & ! initialize fuel categories
        ifds,ifde,jfds,jfde,                    & ! fire domain dims
        ifms,ifme,jfms,jfme,                    & ! fire memory dims
        ifps,ifpe,jfps,jfpe,                    &
        ifts,ifte,jfts,jfte,                    & ! fire patch dims
        time_start,dt,                          & ! time and increment
        dxf,dyf,                                & ! fire mesh spacing
        ignition,hfx,                           & ! description of ignition lines
        fxlong,fxlat,                           & ! fire mesh coordinates
        fire_hfx,                               & ! given heat flux
        tign_in,                                & ! given igntion time
        lfn,lfn_out,tign,fuel_frac,             & ! state: level function, ign time, fuel left
        fire_area,                              & ! output: fraction of cell burning
        fuel_frac_burnt,                        & ! output: fuel fraction burned in this step
        fgrnhfx,fgrnqfx,                        & ! output: heat fluxes
        ros,flineint,flineint2,                 & ! diagnostic variables
        f_ros0,f_rosx,f_rosy,f_ros,             & ! fire risk spread 
        f_int,f_lineint,f_lineint2,             & ! fire risk intensities 
        nfuel_cat,                              & ! fuel data per point 
        fuel_time,fwh,fz0,                      & ! save derived internal data
        fp                                      & ! fire coefficients
    )

     if(ifun.eq.2)then
        call setup_wind_log_interpolation(           &
            ids,ide,  jds,jde,                    & ! atm grid dimensions
            ims,ime,  jms,jme,                    &
            ips,ipe,jps,jpe,                              &
            its,ite,jts,jte,                              &
            ifds, ifde, jfds, jfde,                       & ! fire grid dimensions
            ifms, ifme, jfms, jfme,                       &
            ifts,ifte,jfts,jfte,                          &
            ir,jr,                                        & ! atm/fire grid ratio
            z0,                                           & ! atm grid arrays in
            nfuel_cat,                                    & ! fuel data per point 
            fz0,fwh,wz0)                                    ! fire arrays out

    
    elseif(ifun.eq.6)then 


    ! populate the rate of spread in the 8 directions
    do j=jfts,jfte
       do i=ifts,ifte
          f_ros11(i,j)=fire_rate_of_spread( dxf*(1-2), dyf*(1-2), i,j,fp)
          f_ros12(i,j)=fire_rate_of_spread( dxf*(1-2), dyf*(2-2), i,j,fp)
          f_ros13(i,j)=fire_rate_of_spread( dxf*(1-2), dyf*(3-2), i,j,fp)
          f_ros21(i,j)=fire_rate_of_spread( dxf*(2-2), dyf*(1-2), i,j,fp)
          f_ros23(i,j)=fire_rate_of_spread( dxf*(2-2), dyf*(3-2), i,j,fp)
          f_ros31(i,j)=fire_rate_of_spread( dxf*(3-2), dyf*(1-2), i,j,fp)
          f_ros32(i,j)=fire_rate_of_spread( dxf*(3-2), dyf*(2-2), i,j,fp)
          f_ros33(i,j)=fire_rate_of_spread( dxf*(3-2), dyf*(3-2), i,j,fp)
       enddo
    enddo
    call print_2d_stats(ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme,f_ros11,'driver_phys:f_ros11')        
    call print_2d_stats(ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme,f_ros12,'driver_phys:f_ros12')        
    call print_2d_stats(ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme,f_ros13,'driver_phys:f_ros13')        
    call print_2d_stats(ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme,f_ros21,'driver_phys:f_ros21')        
    call print_2d_stats(ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme,f_ros23,'driver_phys:f_ros23')        
    call print_2d_stats(ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme,f_ros31,'driver_phys:f_ros31')        
    call print_2d_stats(ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme,f_ros32,'driver_phys:f_ros32')        
    call print_2d_stats(ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme,f_ros33,'driver_phys:f_ros33')        
    
    ! heat fluxes into the atmosphere    

        call print_2d_stats(ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme,ros,'sfire_driver:ros')        
        call print_2d_stats(ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme,fgrnhfx,'sfire_driver:fgrnhfx')
        call print_2d_stats(ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme,fgrnqfx,'sfire_driver:fgrnqfx')
    
        ! sum the fluxes over atm cells
        if(use_atm_vars)then                                  
          call sum_2d_cells(        &
            ifms,ifme,jfms,jfme,  &
            ifts,ifte,jfts,jfte,  &
            fuel_frac,              &
            ims, ime, jms, jme,   &
            its,ite,jts,jte,      &
            avg_fuel_frac)
          call sum_2d_cells(        &
            ifms,ifme,jfms,jfme,  &
            ifts,ifte,jfts,jfte,  &
            fgrnhfx,              &
            ims, ime, jms, jme,   &
            its,ite,jts,jte,      &
            grnhfx)
!comment out the next call to get results as before commit 55fd92051196b796891b60cb7ec1c4bdb8800078
          call sum_2d_cells(        &
            ifms,ifme,jfms,jfme,  &
            ifts,ifte,jfts,jfte,  &
            fgrnqfx,              &
            ims, ime, jms, jme,   &
            its,ite,jts,jte,      &
            grnqfx)

!$OMP CRITICAL(SFIRE_DRIVER_CRIT)
          write(msg,'(a,f6.3)')'fire-atmosphere feedback scaling ',fire_atm_feedback
!$OMP end CRITICAL(SFIRE_DRIVER_CRIT)
	  call message(msg)
          s = 1./(ir*jr)
          do j=jts,jte
            do i=its,ite
                ! scale ground fluxes to get the averages
                avg_fuel_frac(i,j)=avg_fuel_frac(i,j)*s
                grnhfx(i,j)=fire_atm_feedback*grnhfx(i,j)*s
                grnqfx(i,j)=fire_atm_feedback*grnqfx(i,j)*s
                ! we do not have canopy fluxes yet...
                canhfx(i,j)=0
                canqfx(i,j)=0
            enddo
          enddo

          call print_2d_stats(its,ite,jts,jte,ims,ime,jms,jme,grnhfx,'fire_driver:grnhfx')
          call print_2d_stats(its,ite,jts,jte,ims,ime,jms,jme,grnqfx,'fire_driver:grnqfx')
       endif

    endif ! ifun=6
  endif

enddo ! tiles
!$OMP END PARALLEL DO

#ifdef DEBUG_OUT
if(ifun.eq.1)then
    if(pid.ne.0)then
        call write_array_m(ips,ipe,jps,jpe,ims,ime,jms,jme,zs,'zs',pid)
        call write_array_m(ifps,ifpe,jfps,jfpe,ifms,ifme,jfms,jfme,fp%zsf,'zsf',pid)
    endif
endif
#endif

if(ifun.eq.2)then
! ****      THIS REALLY SHOULD BE DONE ONCE NOT EVERY TIMESTEP
!           we cannot store the ignition structure in grid
!           so we set global ignition instead
!           and we need to have fxlong fxlat for the corners (for checking)
!           which requires MPI calls 
!           must be run after the parallel do because of all that
            call message('Processing ignition lines')
            call fire_ignition_convert (config_flags,ignition,          &
                fxlong, fxlat,                                          &
                ifds,ifde, jfds,jfde,                                   &
                ifms,ifme, jfms,jfme,                                   &
                ifps,ifpe, jfps,jfpe )
            call fire_hfx_convert (config_flags,hfx) ! do we even use heat flux lines any more?
endif

if (ifun.eq.3)then
 call print_3d_stats_by_slice(ips,ipe,1,moisture_classes,jps,jpe,ims,ime,1,nfmc,jms,jme,fmc_gc,'fmc_gc')
 call print_chsum(itimestep,ims,ime,1,nfmc,jms,jme,ids,ide,1,moisture_classes,jds,jde,ips,ipe,1,moisture_classes,jps,jpe,0,0,0,fmc_gc,'fmc_gc')
endif

if (ifun.eq.4)then

 call print_chsum(itimestep,ifms,ifme,1,1,jfms,jfme,ifds,ifde,1,1,jfds,jfde,ifps,ifpe,1,1,jfps,jfpe,0,0,0,fmc_g,'fmc_g')
 !call print_chsum(itimestep,ims,ime,1,1,jms,jme,ids,ide,1,1,jds,jde,ips,ipe,1,1,jps,jpe,1,0,0,uah,'uah')
 !call print_chsum(itimestep,ims,ime,1,1,jms,jme,ids,ide,1,1,jds,jde,ips,ipe,1,1,jps,jpe,0,0,1,vah,'vah')
 call print_chsum(itimestep,ifms,ifme,1,1,jfms,jfme,ifds,ifde,1,1,jfds,jfde,ifps,ifpe,1,1,jfps,jfpe,0,0,0,fp%vx,'uf')
 call print_chsum(itimestep,ifms,ifme,1,1,jfms,jfme,ifds,ifde,1,1,jfds,jfde,ifps,ifpe,1,1,jfps,jfpe,0,0,0,fp%vy,'vf')
 call print_chsum(itimestep,ifms,ifme,1,1,jfms,jfme,ifds,ifde,1,1,jfds,jfde,ifps,ifpe,1,1,jfps,jfpe,0,0,0,lfn,'lfn')
 call print_chsum(itimestep,ifms,ifme,1,1,jfms,jfme,ifds,ifde,1,1,jfds,jfde,ifps,ifpe,1,1,jfps,jfpe,0,0,0,tign,'tign')

#ifdef DEBUG_OUT
    if(pid.gt.0)then
 !       call write_array_m(ips,ipe1,jps,jpe,ims,ime,jms,jme,uah,'uah',pid)
 !       call write_array_m(ips,ipe,jps,jpe1,ims,ime,jms,jme,vah,'vah',pid)
        call write_array_m(ips,ipe,jps,jpe,ims,ime,jms,jme,grnhfx,'grnhfx',pid)
        call write_array_m(ips,ipe,jps,jpe,ims,ime,jms,jme,grnqfx,'grnqfx',pid)
        call write_array_m3(ips,ipe1,kds,kde+1,jps,jpe,ims,ime,kms,kme,jms,jme,u,'u',pid)
        call write_array_m3(ips,ipe,kds,kde+1,jps,jpe1,ims,ime,kms,kme,jms,jme,v,'v',pid)
        call write_array_m(ifps,ifpe,jfps,jfpe,ifms,ifme,jfms,jfme,fp%vx,'uf',pid)
        call write_array_m(ifps,ifpe,jfps,jfpe,ifms,ifme,jfms,jfme,fp%vy,'vf',pid)
    endif
#endif
endif

if(ifun.eq.5)then
    call print_chsum(itimestep,ifms,ifme,1,1,jfms,jfme,ifds,ifde,1,1,jfds,jfde,ifps,ifpe,1,1,jfps,jfpe,0,0,0,lfn,'lfn')
    call print_chsum(itimestep,ifms,ifme,1,1,jfms,jfme,ifds,ifde,1,1,jfds,jfde,ifps,ifpe,1,1,jfps,jfpe,0,0,0,tign,'tign')
#ifdef DEBUG_OUT
    if(pid.gt.0)then
        call write_array_m(ifps,ifpe,jfps,jfpe,ifms,ifme,jfms,jfme,lfn,'lfn',pid)
        call write_array_m(ifps,ifpe,jfps,jfpe,ifms,ifme,jfms,jfme,tign,'tign',pid)
    endif
#endif
endif

if(ifun.eq.6)then
    call print_chsum(itimestep,ifms,ifme,1,1,jfms,jfme,ifds,ifde,1,1,jfds,jfde,ifps,ifpe,1,1,jfps,jfpe,0,0,0,fgrnhfx,'fgrnhfx')
    call print_chsum(itimestep,ifms,ifme,1,1,jfms,jfme,ifds,ifde,1,1,jfds,jfde,ifps,ifpe,1,1,jfps,jfpe,0,0,0,fgrnqfx,'fgrnqfx')
    call print_chsum(itimestep,ims,ime,1,1,jms,jme,ids,ide,1,1,jds,jde,ips,ipe,1,1,jps,jpe,0,0,0,grnhfx,'grnhfx')
    call print_chsum(itimestep,ims,ime,1,1,jms,jme,ids,ide,1,1,jds,jde,ips,ipe,1,1,jps,jpe,0,0,0,grnqfx,'grnqfx')
#ifdef DEBUG_OUT
    if(pid.gt.0)then
        call write_array_m(ips,ipe,jps,jpe,ims,ime,jms,jme,grnhfx,'grnhfx',pid)
        call write_array_m(ips,ipe,jps,jpe,ims,ime,jms,jme,grnqfx,'grnqfx',pid)
        call write_array_m(ifps,ifpe,jfps,jfpe,ifms,ifme,jfms,jfme,fuel_frac,'fuel_frac',pid)
        call write_array_m(ifps,ifpe,jfps,jfpe,ifms,ifme,jfms,jfme,fgrnhfx,'fgrnhfx',pid)
        call write_array_m(ifps,ifpe,jfps,jfpe,ifms,ifme,jfms,jfme,fgrnqfx,'fgrnqfx',pid)
    endif
#endif
endif

end subroutine sfire_driver_phys

!
!***
!

subroutine check_fmesh(ids,ide,ifds,ifde,ir,s)
!*** purpose: check if fire and atm meshes line up
implicit none
!*** arguments
integer, intent(in)::ids,ide,ifds,ifde,ir
character(len=*),intent(in)::s
!*** local
character(len=128)msg
!*** executable
if ((ide-ids+1)*ir.ne.(ifde-ifds+1))then
!$OMP CRITICAL(SFIRE_DRIVER_CRIT)
    write(msg,1)s,ids,ide,ifds,ifde,ir
1   format('module_fr_sfire_driver: incompatible bounds ',a,' atm ',i5,':',i5,' fire ',i5,':',i5,' ratio ',i3)    
!$OMP END CRITICAL(SFIRE_DRIVER_CRIT)
    call crash(msg)
endif
end subroutine check_fmesh

subroutine fire_ignition_convert (config_flags,ignition,             &
             fxlong, fxlat,                                          &
             ifds,ifde, jfds,jfde,                                   &
             ifms,ifme, jfms,jfme,                                   &
             ifps,ifpe, jfps,jfpe )
    implicit none
! create ignition arrays from scalar flags
!*** arguments
    TYPE (grid_config_rec_type) , INTENT(IN)          :: config_flags
    TYPE (lines_type), INTENT(OUT):: ignition ! any values from input discarded 
    integer::ifds,ifde, jfds,jfde,                                   &
             ifms,ifme, jfms,jfme,                                   &
             ifps,ifpe, jfps,jfpe 
    real, dimension(ifms:ifme,jfms:jfme):: fxlong,fxlat
!*** local
    integer::i,j,ii,jj
    logical:: real,ideal
    character(len=128)msg
    real:: corner_longlat(2,2,2), corner_longlat_1(8), corner_longlat_2(8),lon(2),lat(2)
    real, dimension(2,2):: corner_long,corner_lat ! coordinates of fire mesh corner cells

!*** executable


    ignition%max_lines=5 ! number of lines that have entries in the namelist
    ignition%num_lines=config_flags%fire_num_ignitions

    ! this is only until I figure out how to input arrays through the namelist...
    if(fire_max_lines.lt.ignition%max_lines)call crash('fire_max_lines too small')

    ! figure out which kind of coordinates from the first given
    ideal=config_flags%fire_ignition_start_x1 .ne.0. .or. config_flags%fire_ignition_start_y1 .ne. 0.
    real=config_flags%fire_ignition_start_lon1 .ne. 0. .or. config_flags%fire_ignition_start_lat1 .ne. 0.
    if(ideal)call message('Using ideal ignition coordinates, m from the lower left domain corner')
    if(real)call message('Using real ignition coordinates, longitude and latitude')
    if(ideal.and.real)call crash('Only one of the ideal or real coordinates may be given')

    ignition%longlat=0  ! default, if no ignition
    if(ideal)then
       ! use values from _x and _y variables
       ignition%longlat=0
       ignition%line(1)%start_x=config_flags%fire_ignition_start_x1
       ignition%line(1)%start_y=config_flags%fire_ignition_start_y1
       ignition%line(1)%end_x=config_flags%fire_ignition_end_x1
       ignition%line(1)%end_y=config_flags%fire_ignition_end_y1
       ignition%line(2)%start_x=config_flags%fire_ignition_start_x2
       ignition%line(2)%start_y=config_flags%fire_ignition_start_y2
       ignition%line(2)%end_x=config_flags%fire_ignition_end_x2
       ignition%line(2)%end_y=config_flags%fire_ignition_end_y2
       ignition%line(3)%start_x=config_flags%fire_ignition_start_x3
       ignition%line(3)%start_y=config_flags%fire_ignition_start_y3
       ignition%line(3)%end_x=config_flags%fire_ignition_end_x3
       ignition%line(3)%end_y=config_flags%fire_ignition_end_y3
       ignition%line(4)%start_x=config_flags%fire_ignition_start_x4
       ignition%line(4)%start_y=config_flags%fire_ignition_start_y4
       ignition%line(4)%end_x=config_flags%fire_ignition_end_x4
       ignition%line(4)%end_y=config_flags%fire_ignition_end_y4
       ignition%line(5)%start_x=config_flags%fire_ignition_start_x5
       ignition%line(5)%start_y=config_flags%fire_ignition_start_y5
       ignition%line(5)%end_x=config_flags%fire_ignition_end_x5
       ignition%line(5)%end_y=config_flags%fire_ignition_end_y5
    endif
    if(real)then
        ! use values from _long and _lat
       ignition%longlat=1
       ignition%line(1)%start_x=config_flags%fire_ignition_start_lon1
       ignition%line(1)%start_y=config_flags%fire_ignition_start_lat1
       ignition%line(1)%end_x=config_flags%fire_ignition_end_lon1
       ignition%line(1)%end_y=config_flags%fire_ignition_end_lat1
       ignition%line(2)%start_x=config_flags%fire_ignition_start_lon2
       ignition%line(2)%start_y=config_flags%fire_ignition_start_lat2
       ignition%line(2)%end_x=config_flags%fire_ignition_end_lon2
       ignition%line(2)%end_y=config_flags%fire_ignition_end_lat2
       ignition%line(3)%start_x=config_flags%fire_ignition_start_lon3
       ignition%line(3)%start_y=config_flags%fire_ignition_start_lat3
       ignition%line(3)%end_x=config_flags%fire_ignition_end_lon3
       ignition%line(3)%end_y=config_flags%fire_ignition_end_lat3
       ignition%line(4)%start_x=config_flags%fire_ignition_start_lon4
       ignition%line(4)%start_y=config_flags%fire_ignition_start_lat4
       ignition%line(4)%end_x=config_flags%fire_ignition_end_lon4
       ignition%line(4)%end_y=config_flags%fire_ignition_end_lat4
       ignition%line(5)%start_x=config_flags%fire_ignition_start_lon5
       ignition%line(5)%start_y=config_flags%fire_ignition_start_lat5
       ignition%line(5)%end_x=config_flags%fire_ignition_end_lon5
       ignition%line(5)%end_y=config_flags%fire_ignition_end_lat5
    endif
    ! common to both cases
       ignition%line(1)%ros=config_flags%fire_ignition_ros1 
       ignition%line(1)%radius=config_flags%fire_ignition_radius1 
       ignition%line(1)%start_time=config_flags%fire_ignition_start_time1 
       ignition%line(1)%end_time=config_flags%fire_ignition_end_time1 
       ignition%line(2)%ros=config_flags%fire_ignition_ros2 
       ignition%line(2)%radius=config_flags%fire_ignition_radius2 
       ignition%line(2)%start_time=config_flags%fire_ignition_start_time2 
       ignition%line(2)%end_time=config_flags%fire_ignition_end_time2 
       ignition%line(3)%ros=config_flags%fire_ignition_ros3 
       ignition%line(3)%radius=config_flags%fire_ignition_radius3 
       ignition%line(3)%start_time=config_flags%fire_ignition_start_time3 
       ignition%line(3)%end_time=config_flags%fire_ignition_end_time3 
       ignition%line(4)%ros=config_flags%fire_ignition_ros4 
       ignition%line(4)%radius=config_flags%fire_ignition_radius4 
       ignition%line(4)%start_time=config_flags%fire_ignition_start_time4 
       ignition%line(4)%end_time=config_flags%fire_ignition_end_time4 
       ignition%line(5)%ros=config_flags%fire_ignition_ros5 
       ignition%line(5)%radius=config_flags%fire_ignition_radius5 
       ignition%line(5)%start_time=config_flags%fire_ignition_start_time5
       ignition%line(5)%end_time=config_flags%fire_ignition_end_time5

       call postprocess_lines(ignition,'ros',config_flags)

!       get the coordinates of the corner cells
        corner_longlat=0.
        if(ifds.eq.ifps.and.jfds.eq.jfps)then     
            corner_longlat(1,1,1)=fxlong(ifps,jfps)
            corner_longlat(1,1,2)=fxlat(ifps,jfps)
        endif
        if(ifds.eq.ifps.and.jfde.eq.jfpe)then 
            corner_longlat(1,2,1)=fxlong(ifps,jfpe)
            corner_longlat(1,2,2)=fxlat(ifps,jfpe)
        endif
        if(ifde.eq.ifpe.and.jfds.eq.jfps)then 
            corner_longlat(2,1,1)=fxlong(ifpe,jfps)
            corner_longlat(2,1,2)=fxlat(ifpe,jfps)
        endif
        if(ifde.eq.ifpe.and.jfde.eq.jfpe)then 
            corner_longlat(2,2,1)=fxlong(ifpe,jfpe)
            corner_longlat(2,2,2)=fxlat(ifpe,jfpe)
        endif
        corner_longlat_1=reshape(corner_longlat,(/8/))
#ifdef DM_PARALLEL
        call wrf_dm_sum_reals(corner_longlat_1,corner_longlat_2)
#else
        corner_longlat_2=corner_longlat_1
#endif
        corner_longlat=reshape(corner_longlat_2,(/2,2,2/))
        corner_long=corner_longlat(1:2,1:2,1)
        corner_lat=corner_longlat(1:2,1:2,2)
        if(fire_print_msg.ge.2)then
            do i=1,2
                do j=1,2
                     write(msg,'(a,2i2,a,2f14.8)')'corner',i,j,' coordinates ',corner_long(i,j),corner_lat(i,j)
                     call message(msg)
                enddo
            enddo
        endif
	lon(1)=(corner_long(1,1)+corner_long(1,2))/2.
	lon(2)=(corner_long(2,1)+corner_long(2,2))/2.
	lat(1)=(corner_lat(1,1)+corner_lat(2,1))/2.
	lat(2)=(corner_lat(1,2)+corner_lat(2,2))/2.
        if(fire_print_msg.ge.2)then
            write(msg,'(4(a,f14.8))')'coordinates ',lon(1),':',lon(2),',',lat(1),':',lat(2)
            call message(msg)
        endif

       do i=1,ignition%num_lines
           call check_ignition_coordinate(ignition%line(i)%start_x,lon(1),lon(2))
           call check_ignition_coordinate(ignition%line(i)%start_y,lat(1),lat(2))
           call check_ignition_coordinate(ignition%line(i)%end_x,lon(1),lon(2))
           call check_ignition_coordinate(ignition%line(i)%end_y,lat(1),lat(2))
       enddo

       if (fire_ignition_clamp>0) then
       do i=1,ignition%num_lines
           call clamp_to_grid(ignition%line(i)%start_x,lon(1),lon(2),ifds,ifde,ignition%line(i)%start_x,ii)
           call clamp_to_grid(ignition%line(i)%start_y,lat(1),lat(2),jfds,jfde,ignition%line(i)%start_y,jj)
           call display_clamp
           call clamp_to_grid(ignition%line(i)%end_x,lon(1),lon(2),ifds,ifde,ignition%line(i)%end_x,ii)
           call clamp_to_grid(ignition%line(i)%end_y,lat(1),lat(2),jfds,jfde,ignition%line(i)%end_y,jj)
           call display_clamp
           ! for now, ii jj ignored. In future replace by fxlong(ii,jj), fxlat(ii,jj) to guard against rounding
       enddo
       endif
       contains
       subroutine display_clamp
        character(len=128)::msg
        real::d1,d2
           if(ii>=ifps.and.ii<=ifpe.and.jj>=jfps.and.jj<=jfpe)then
               write(msg,'(a,2f14.8,a,2i6)')'grid node ',fxlong(ii,jj),fxlat(ii,jj),' index',ii,jj
               call message(msg)
           endif
       end subroutine display_clamp
    end subroutine fire_ignition_convert

     subroutine check_ignition_coordinate(x,x1,x2)
!***	arguments
        real, intent(in)::x,x1,x2
        character(len=128)::msg
        if (.not.(x>x1 .and. x<x2))then
            write(msg,'(a,f14.8,a,2f14.8)')'ignition point coordinate',x,' must be inside the bounds',x1,x2
            call crash(msg)
        endif
     end subroutine check_ignition_coordinate

     subroutine clamp_to_grid(x,x1,x2,i1,i2,xout,iout)
!***	round point on uniform mesh to a mesh node
!***	arguments
        real, intent(in)::x,x1,x2
        integer, intent(in)::i1,i2
        real, intent(out)::xout
        integer, intent(out)::iout
!***	local
        character(len=128)::msg
        integer:: i
        real::r,dx,xr
!***    executable
        dx=(x2-x1)/(i2-i1)
        r=i1+(x-x1)/dx
        iout=nint(r)
        xr=x1+(iout-i1)*dx
        if(fire_print_msg.ge.2)then
            write(msg,'(a,f14.8,a,f14.8,a,i6)')'coordinate ',x,' clamped to ',xr,' index',iout
            call message(msg)
        endif
        xout=xr
    end subroutine clamp_to_grid
!
!***
!
     subroutine postprocess_lines(lines,value_name,config_flags)
!      count lines, fill endpoints, print stats
!***   arguments
        type(lines_type), intent(inout)::lines
        character(len=3), intent(in)::value_name  
        TYPE (grid_config_rec_type) , INTENT(IN)  :: config_flags ! namelist
!***   local
     integer::i,n
     real::value
     real::lat_ctr,lon_ctr
     character(len=128)msg,f2,f3
!***   executable

        n=lines%num_lines
        do i=1,n
            ! find the last line that has positive radius and reset num_lines if needed
            if(lines%line(i)%radius.gt.0.)lines%num_lines=i
            ! expand ignition data given as zero
            if(lines%line(i)%end_x.eq.0.)lines%line(i)%end_x=lines%line(i)%start_x
            if(lines%line(i)%end_y.eq.0.)lines%line(i)%end_y=lines%line(i)%start_y
            if(lines%line(i)%end_time.eq.0.)lines%line(i)%end_time=lines%line(i)%start_time
        enddo

    if(lines%longlat .eq. 0)then
       ! ideal
       !  ignition is in m
       lines%unit_fxlong=1.  
       lines%unit_fxlat=1.
       ! will set fire mesh coordinates to uniform mesh below
    else
       ! real
       lat_ctr=config_flags%cen_lat
       lon_ctr=config_flags%cen_lon
       ! 1 degree in m (approximate OK)
       lines%unit_fxlat=pi2/(360.*reradius)  ! earth circumference in m / 360 degrees
       lines%unit_fxlong=cos(lat_ctr*pi2/360.)*lines%unit_fxlat  ! latitude
    endif

    if(fire_print_msg.ge.2)then
!$OMP CRITICAL(SFIRE_DRIVER_CRIT)
       if(lines%longlat .eq. 0)then
         write(msg,1)'start x','start y','end x','end y','start t','end t',value_name,'radius'
       else
         write(msg,1)'start lon','start lat','end lon','end lat','start time','end time',value_name,'radius'
       endif
1      format(4a10,4a9)
       call message(msg)
       do i=1,lines%num_lines
         select case (value_name)
              case ('ros')
                  value = lines%line(i)%ros
                  f2='(4f10.3,4f9.2)'
                  f3='(4f10.5,4f9.2)'
              case ('hfx')
                  value = lines%line(i)%hfx_value
                  f2='(4f10.3,2f9.2,e9.2,f9.2)'
                  f3='(4f10.5,2f9.2,e9.2,f9.2)'
              case default
                  call crash('postprocess_lines: bad value_name '//value_name)
         end select
         if(lines%longlat .eq. 0)then
           write(msg,f2)lines%line(i)%start_x, lines%line(i)%start_y, &
                      lines%line(i)%end_x,  lines%line(i)%end_y,        &
                      lines%line(i)%start_time, lines%line(i)%end_time, &
                      value, lines%line(i)%radius
         else
           write(msg,f3)lines%line(i)%start_x, lines%line(i)%start_y, &
                      lines%line(i)%end_x,  lines%line(i)%end_y,        &
                      lines%line(i)%start_time, lines%line(i)%end_time, &
                      value, lines%line(i)%radius
         endif
         call message(msg)
         if(lines%line(i)%start_time > lines%line(i)%end_time)then
              call crash('start time may not be after end time')
         endif
       enddo
!$OMP END CRITICAL(SFIRE_DRIVER_CRIT)
     endif
end subroutine postprocess_lines


subroutine fire_hfx_convert (config_flags,hfx)
    implicit none
! create heat flux line(s) from scalar flags
!*** arguments
    TYPE (grid_config_rec_type) , INTENT(IN)          :: config_flags
    TYPE (lines_type), INTENT(OUT):: hfx ! any values from input discarded 
!*** local
    integer::i
    logical:: real,ideal
    real::lat_ctr,lon_ctr
    character(len=128)msg
!*** executable
    ! this is only until I figure out how to input arrays through the namelist...
    hfx%num_lines=config_flags%fire_hfx_num_lines
    if(fire_max_lines.lt.hfx%num_lines)call crash('fire_max_lines too small')

    ! figure out which kind of coordinates from the first given
    ideal=config_flags%fire_hfx_start_x1 .ne.0. .or. config_flags%fire_hfx_start_y1 .ne. 0.
    real=config_flags%fire_hfx_start_lon1 .ne. 0. .or. config_flags%fire_hfx_start_lat1 .ne. 0.
    if(ideal)call message('Using ideal heat flux line coordinates, m from the lower left domain corner')
    if(real)call message('Using real heat flux line coordinates, longitude and latitude')
    if(ideal.and.real)call crash('Only one of the ideal or real coordinates may be given')

    hfx%longlat=0  ! default, if no ignition
    if(ideal)then
       ! use values from _x and _y variables
       hfx%longlat=0
       hfx%line(1)%start_x=config_flags%fire_hfx_start_x1
       hfx%line(1)%start_y=config_flags%fire_hfx_start_y1
       hfx%line(1)%end_x=config_flags%fire_hfx_end_x1
       hfx%line(1)%end_y=config_flags%fire_hfx_end_y1
    endif
    if(real)then
        ! use values from _long and _lat
       hfx%longlat=1
       hfx%line(1)%start_x=config_flags%fire_hfx_start_lon1
       hfx%line(1)%start_y=config_flags%fire_hfx_start_lat1
       hfx%line(1)%end_x=config_flags%fire_hfx_end_lon1
       hfx%line(1)%end_y=config_flags%fire_hfx_end_lat1
    endif
    ! common to both cases
       hfx%line(1)%radius=config_flags%fire_hfx_radius1 
       hfx%line(1)%start_time=config_flags%fire_hfx_start_time1 
       hfx%line(1)%end_time=config_flags%fire_hfx_end_time1 
       hfx%line(1)%trans_time=config_flags%fire_hfx_trans_time1 
       hfx%line(1)%hfx_value=config_flags%fire_hfx_value1 

       call postprocess_lines(hfx,'hfx',config_flags)

end subroutine fire_hfx_convert

subroutine set_flags(config_flags)
USE module_configure
implicit none
TYPE (grid_config_rec_type) , INTENT(IN)          :: config_flags
! copy flags from wrf to module_fr_sfire_util
! for instructions how to add a flag see the top of module_fr_sfire_util.F

fire_perimeter_time     = config_flags%fire_perimeter_time
fire_tign_in_time       = config_flags%fire_tign_in_time
fire_print_msg          = config_flags%fire_print_msg
fire_print_file         = config_flags%fire_print_file
fuel_left_method        = config_flags%fire_fuel_left_method
fuel_left_irl           = config_flags%fire_fuel_left_irl
fuel_left_jrl           = config_flags%fire_fuel_left_jrl
fire_atm_feedback       = config_flags%fire_atm_feedback
fire_hfx_given          = config_flags%fire_hfx_given
fire_hfx_num_lines      = config_flags%fire_hfx_num_lines
fire_hfx_latent_part    = config_flags%fire_hfx_latent_part
fire_update_fuel_frac   = config_flags%fire_update_fuel_frac
boundary_guard          = config_flags%fire_boundary_guard
fire_back_weight        = config_flags%fire_back_weight
fire_grows_only         = config_flags%fire_grows_only
sfire_upwinding         = config_flags%sfire_upwinding
fire_viscosity          = config_flags%fire_viscosity 
fire_lfn_ext_up         = config_flags%fire_lfn_ext_up 
fire_test_steps         = config_flags%fire_test_steps 
!fire_topo_from_atm     = config_flags%fire_topo_from_atm
fire_advection          = config_flags%fire_advection
fire_wind_log_interp    = config_flags%fire_wind_log_interp
fire_use_windrf         = config_flags%fire_use_windrf
fire_fmc_read           = config_flags%fire_fmc_read
fire_ignition_clamp     = config_flags%fire_ignition_clamp
kfmc_ndwi               = config_flags%kfmc_ndwi
fndwi_from_ndwi         = config_flags%fndwi_from_ndwi 
fire_can_top_read       = config_flags%fire_can_top_read
end subroutine set_flags

!
!*****************************
!


subroutine set_fp_from_grid(grid,fp)
    implicit none
    type(domain),intent(in)::grid
    type(fire_params),intent(out)::fp

    ! pointers to be passed to fire spread formulas
    fp%vx => grid%uf         ! fire winds
    fp%vy => grid%vf         ! fire winds
    fp%zsf => grid%zsf       ! terrain height
    fp%dzdxf => grid%dzdxf   ! terrain grad
    fp%dzdyf => grid%dzdyf   ! terrain grad
    fp%bbb => grid%bbb       ! spread formula coeff
    fp%phisc => grid%phisc   ! spread formula coeff
    fp%phiwc => grid%phiwc   ! spread formula coeff
    fp%r_0 => grid%r_0       ! spread formula coeff
    fp%fgip => grid%fgip     ! spread formula coeff
    fp%ischap => grid%ischap ! spread formula coeff
    fp%fuel_time => grid%fuel_time ! time for fuel to burn to 1/e
    fp%fmc_g => grid%fmc_g   ! fuel moisture, ground
    fp%nfuel_cat => grid%nfuel_cat ! fuel category


end subroutine set_fp_from_grid

subroutine check_grid_alloc(grid,config_flags)

! check if needed arrays in grid are fully allocated
! they may not be if they are missing in the registry packages
! leading to memory corruption which is hard to debug

    TYPE(domain) , TARGET, INTENT(IN) :: grid                             ! state 
    TYPE (grid_config_rec_type) , INTENT(IN)  :: config_flags ! namelist

    call message('check_grid_alloc: checking if grid arrays used sfire_driver_phys are allocated',level=2)
    call check_size(size(grid%u_2),'u_2')
    call check_size(size(grid%v_2),'v_2')
    call check_size(size(grid%ph_2),'ph_2')
    call check_size(size(grid%phb),'phb')
    call check_size(size(grid%z0),'z0')
    call check_size(size(grid%xlong),'xlong')
    call check_size(size(grid%xlat),'xlat')
    call check_size(size(grid%tign_in),'tign_in')
    call check_size(size(grid%lfn),'lfn')
    call check_size(size(grid%tign_g),'tign_g')
    call check_size(size(grid%fire_area),'fire_area')
    call check_size(size(grid%fuel_frac),'fuel_frac')
    call check_size(size(grid%fuel_frac_burnt),'fuel_frac_burnt')
    call check_size(size(grid%lfn_out),'lfn_out')
    call check_size(size(grid%avg_fuel_frac),'avg_fuel_frac')
    call check_size(size(grid%grnhfx),'grnhfx')
    call check_size(size(grid%grnqfx),'grnqfx')
    call check_size(size(grid%canhfx),'canhfx')
    call check_size(size(grid%canqfx),'canqfx')
    call check_size(size(grid%fgrnhfx),'fgrnhfx')
    call check_size(size(grid%fgrnqfx),'fgrnqfx')
    call check_size(size(grid%fcanhfx),'fcanhfx')
    call check_size(size(grid%ros),'ros')
    call check_size(size(grid%flineint),'flineint')
    call check_size(size(grid%flineint2),'flineint2')
    call check_size(size(grid%f_ros0),'f_ros0')
    call check_size(size(grid%f_rosx),'f_rosx')
    call check_size(size(grid%f_rosy),'f_rosy')
    call check_size(size(grid%f_ros),'f_ros')
    call check_size(size(grid%f_int),'f_int')
    call check_size(size(grid%f_lineint),'f_lineint')
    call check_size(size(grid%f_lineint2),'f_lineint2')
    call check_size(size(grid%f_ros11),'f_ros11')
    call check_size(size(grid%f_ros12),'f_ros12')
    call check_size(size(grid%f_ros13),'f_ros13')
    call check_size(size(grid%f_ros21),'f_ros21')
    call check_size(size(grid%f_ros23),'f_ros23')
    call check_size(size(grid%f_ros31),'f_ros31')
    call check_size(size(grid%f_ros32),'f_ros32')
    call check_size(size(grid%f_ros33),'f_ros33')
    call check_size(size(grid%fxlong),'fxlong')
    call check_size(size(grid%fxlat),'fxlat')
    call check_size(size(grid%fire_hfx),'fire_hfx')
    call check_size(size(grid%nfuel_cat),'nfuel_cat')
    call check_size(size(grid%fuel_time),'fuel_time')
    call check_size(size(grid%fz0),'fz0')
    call check_size(size(grid%fwh),'fwh')
    call check_size(size(grid%uah),'uah')
    call check_size(size(grid%vah),'vah')
    call check_size(size(grid%rainc),'rainc')
    call check_size(size(grid%rainnc),'rainnc')
    call check_size(size(grid%t2),'t2')
    call check_size(size(grid%q2),'q2')
    call check_size(size(grid%psfc),'psfc')
    call check_size(size(grid%ndwi),'ndwi')
    call check_size(size(grid%fndwi),'fndwi')
    call message('check_grid_alloc: checking if grid arrays used in fire params are allocated',level=2)
    call check_size(size(grid%uf),'uf')
    call check_size(size(grid%vf),'vf')
    call check_size(size(grid%zsf),'zsf')
    call check_size(size(grid%dzdxf),'dzdxf')
    call check_size(size(grid%dzdyf),'dzdyf')
    call check_size(size(grid%bbb),'bbb')
    call check_size(size(grid%phiwc),'phiwc')
    call check_size(size(grid%r_0),'r_0')
    call check_size(size(grid%fgip),'fgip')
    call check_size(size(grid%ischap),'ischap')
    call check_size(size(grid%fmc_g),'fmc_g')
    call message('check_grid_alloc: checking if grid arrays used in fuel moisture are allocated',level=2)
    if (config_flags%fmoisti_run==1)then
        call check_size(size(grid%fmc_gc),'fmc_gc')
        call check_size(size(grid%fmep),'fmep')
        call check_size(size(grid%fmc_equi),'fmc_equi')
        call check_size(size(grid%fmc_lag),'fmc_lag')
        call check_size(size(grid%rain_old),'rain_old')
        call check_size(size(grid%t2_old),'t2_old')
        call check_size(size(grid%q2_old),'q2_old')
        call check_size(size(grid%rain_old),'rain_old')
        call check_size(size(grid%psfc_old),'psfc_old')
        call check_size(size(grid%rh_fire),'rh_fire')
        call check_size(size(grid%fmc_g),'fmc_g')
    endif
    if (config_flags%fmoisti_interp .eq. 1)then
        call check_size(size(grid%fmc_gc),'fmc_gc')
    endif
    call message('check_grid_alloc: checking if grid arrays used in emissions are allocated',level=2)
    call check_size(size(grid%chem),'chem',config_flags%chem_opt)
    call check_size(size(grid%tracer),'tracer',config_flags%tracer_opt)

end subroutine check_grid_alloc

subroutine check_size(sz,name,req)
integer, intent(in)::sz
integer, optional::req
character(len=*), intent(in)::name
character(len=128)::msg
logical::strict = .false.

write(msg,'(a,a10,a,i12)')'array ',name,' has size ',sz
call message(msg,level=3)
if(present(req)) strict=req > 0
if (sz < 2 .and. strict) then
    write(msg,'(a,a,a,i2,a)')'array ',name,' is size ',sz,': not allocated, add to registry package'
    call crash(trim(msg))
endif
end subroutine check_size

            
!subroutine print_id
!character(len=128)::id,msg
!include "sfire_id.inc"
!msg=id
!call message(msg,level=1)
!end subroutine print_id

end module module_fr_sfire_driver
