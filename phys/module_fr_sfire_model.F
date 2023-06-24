!
#define DEBUG_OUT

module module_fr_sfire_model

use module_fr_sfire_core
use module_fr_sfire_util
use module_fr_sfire_phys

implicit none

contains

subroutine sfire_model (                    &
    id,                                     & ! unique number for prints and debug
    ifun,                                   & ! what to do see below
    restart,replay,                         & ! use existing state, prescribe spread
    run_fuel_moisture,                      & ! if need update fuel moisture in pass 4
    ifuelread,nfuel_cat0,                   & ! initialize fuel categories
    ifds,ifde,jfds,jfde,                    & ! fire domain dims - the whole domain
    ifms,ifme,jfms,jfme,                    & ! fire memory dims - how declared
    ifps,ifpe,jfps,jfpe,                    & ! patch - nodes owned by this process
    ifts,ifte,jfts,jfte,                    & ! fire tile dims  - this thread
    time_start,dt,                          & ! time and increment
    fdx,fdy,                                & ! fire mesh spacing,
    ignition,hfx,                           & ! small array of ignition line descriptions
    coord_xf,coord_yf,                      & ! fire mesh coordinates
    fire_hfx,                               & ! input: given heat flux, or set inside
    tign_in,                                & ! ignition time, if given
    lfn,lfn_out,tign,fuel_frac,fire_area,   & ! state: level function, ign time, fuel left, area burning
    fuel_frac_burnt,                        & ! output: fuel fraction burnt in this timestep
    fgrnhfx,fgrnqfx,                          & ! output: heat fluxes
    ros,flineint,flineint2,                 & ! diagnostic variables
    f_ros0,f_rosx,f_rosy,f_ros,             & ! fire risk spread 
    f_int,f_lineint,f_lineint2,             & ! fire risk intensities 
    nfuel_cat,                              & ! fuel data per point 
    fuel_time,fwh,fz0,                      & ! save derived internal data
    fp &
) 

! This subroutine implements the fire spread model.
! All quantities are on the fire grid. It inputs
! winds given on the nodes of the fire grid
! and outputs the heat fluxes on the cells of the fire grid.
! This subroutine has no knowledge of any atmospheric model.
! This code was written to conform with the WRF parallelism model, however it
! does not depend on it. It can be called with domain equal to tile.
! Wind and height must be given on 1 more node beyond the domain bounds. 
! The subroutine changes only array entries of the arguments in the tile.
! Upon exit with ifun=2 (time step), lfn_out is to be copied into lfn by the caller.
! When this subroutine is used on separate tiles that make a domain the value, the
! it uses lfn on a strip of width 2 from neighboring tiles.
!
! All computation is done on one tile. 
!
! This subroutine is intended to be called in a loop like
!
! 
! do ifun=1,6 (if initizalize run, otherwise 3,6)
!   start parallel loop over tiles
!       if ifun=1, set z and fuel data
!       if ifun=3, set the wind arrays
!       call sfire_model(....)
!   end parallel loop over tiles
!
!   
!   if ifun=0
!       halo exchange on z width 2
!       halo exchange on fuel data width 1
!   endif
!   
!   if ifun=3, halo exchange on winds width 2
!    
! enddo

implicit none

!*** arguments

! control switches
integer, intent(in) :: id
integer, intent(in) :: ifun                 ! 1 = initialize run pass 1
                                            ! 2 = initialize run pass 2
                                            ! 3 = initialize timestep
                                            ! 4 = do one timestep 
                                            ! 5 = copy timestep output to input
                                            ! 6 = compute output fluxes
logical, intent(in):: restart               ! if true, use existing state
logical, intent(in):: replay                ! if true, use tign_g for level set
logical, intent(in)::run_fuel_moisture      ! 
! scalar data
integer, intent(in) :: ifuelread,nfuel_cat0 ! for set_fire_params
integer, intent(in) :: ifds,ifde,jfds,jfde,&  ! fire domain bounds
        ifps,ifpe,jfps,jfpe                ! patch - nodes owned by this process
integer, intent(in) :: ifts,ifte,jfts,jfte  ! fire tile bounds
integer, intent(in) :: ifms,ifme,jfms,jfme  ! fire memory array bounds
REAL,INTENT(in) :: time_start,dt            ! starting time, time step
REAL,INTENT(in) :: fdx,fdy                  ! spacing of the fire mesh
! array data
type(lines_type), intent(in):: ignition,hfx  ! descriptions of ignition lines and hfx lines
real, dimension(ifms:ifme, jfms:jfme), intent(in):: & 
    coord_xf,coord_yf                      ! node coordinates  
real, dimension(ifms:ifme, jfms:jfme), intent(inout):: & 
    fire_hfx                                ! given heat flux
    
real, dimension(ifms:ifme, jfms:jfme), intent(inout):: tign_in ! given ignition times
! state
REAL, INTENT(inout), dimension(ifms:ifme,jfms:jfme):: &
    lfn   , &                               ! level function: fire is where lfn<0 (node)
    tign  , &                               ! absolute time of ignition (node)
    fuel_frac                               ! fuel fraction (node), currently redundant

REAL, INTENT(out), dimension(ifms:ifme,jfms:jfme):: &
    fire_area, &                            ! fraction of each cell burning
    fuel_frac_burnt                         ! fuel fraction burned in this timestep
    
! output
REAL, INTENT(out), dimension(ifms:ifme,jfms:jfme):: &
    lfn_out, &                              !                              
    fgrnhfx,fgrnqfx, &                        ! heat fluxes J/m^2/s  (cell)             
    ros,flineint,flineint2,                 & ! diagnostic variables
    f_ros0,f_rosx,f_rosy,f_ros,             & ! fire risk spread 
    f_int,f_lineint,f_lineint2                ! fire risk intensities 

 
! constant arrays - set at initialization
real, intent(inout), dimension(ifms:ifme, jfms:jfme)::nfuel_cat ! cell based, data, constant
real,intent(inout),dimension(ifms:ifme,jfms:jfme):: fuel_time,fwh,fz0
type(fire_params),intent(inout)::fp

!*** local

integer :: xifms,xifme,xjfms,xjfme  ! memory bounds for pass-through arguments to normal spread
real, dimension(ifts:ifte,jfts:jfte)::fuel_frac_end
integer::ignited,ig,i,j,itso,iteo,jtso,jteo
real::tbound,err,erri,errj,maxgrad,grad,tfa,thf,mhf,tqf,mqf,aw,mw,t
character(len=128)::msg
logical:: freeze_fire
real::fireline_mask=0.

!*** executable

call check_mesh_2dim(ifts-1,ifte+1,jfts-1,jfte+1,ifms,ifme,jfms,jfme)

xifms=ifms  ! dimensions for the include file
xifme=ifme
xjfms=jfms
xjfme=jfme


! init flags
freeze_fire = fire_hfx_given .ne.  0


if(ifun.eq.1)then       ! do nothing, init pass 1 is outside only
! !$OMP SINGLE
!! done in driver now
!   call init_fuel_cats  ! initialize fuel subsystem
! !$OMP END SINGLE
        if(replay) then    
        ! when starting replay, recompute lfn and fuel_frac
          call message('replay, setting the level set function',level=1)
          do j=jfts,jfte
            do i=ifts,ifte
              lfn(i,j) = tign(i,j) - time_start ! <0 if burning at the end of the time step
            enddo
          enddo
          call check_lfn_tign('replay, lfn from tign',time_start,ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme,lfn,tign)
        endif
elseif(ifun.eq.2)then   
        ! initialize all arrays that the model will not change later

        ! assuming halo on zsf done
        ! extrapolate on 1 row of cells beyond the domain boundary
        ! including on the halo regions 

        call continue_at_boundary(1,1,0., & ! do x direction or y direction
            ifms,ifme,jfms,jfme,           &                ! memory dims
            ifds,ifde,jfds,jfde, &                     ! domain dims 
            ifps,ifpe,jfps,jfpe, &            ! patch dims - winds defined up to +1
            ifts,ifte,jfts,jfte, &                ! tile dims
            itso,iteo,jtso,jteo, &              ! where set now
            fp%zsf)                               ! array

!       compute the gradients once for all
        err=0.
        maxgrad=0.
        do j=jfts,jfte
            do i=ifts,ifte
                erri = fp%dzdxf(i,j) - (fp%zsf(i+1,j)-fp%zsf(i-1,j))/(2.*fdx)
                errj = fp%dzdyf(i,j) - (fp%zsf(i,j+1)-fp%zsf(i,j-1))/(2.*fdy)
                err=max(err,abs(erri),abs(errj))
                grad=sqrt(fp%dzdxf(i,j)**2+fp%dzdyf(i,j)**2)
                maxgrad=max(maxgrad,grad)
            enddo
        enddo
!$OMP CRITICAL(SFIRE_MODEL_CRIT)
        write(msg,*)'max gradient ',maxgrad,' max error against zsf',err
!$OMP END CRITICAL(SFIRE_MODEL_CRIT)
        call message(msg)

        call set_nfuel_cat(      &       ! also on restart
            ifms,ifme,jfms,jfme, &
            ifts,ifte,jfts,jfte, &
            ifuelread,nfuel_cat0,&
            fp%zsf,nfuel_cat)            ! better not use the extrapolated zsf!!

        ! uses nfuel_cat to set the other fuel data arrays
        ! needs zsf on halo width 1 to compute the terrain gradient
        call set_fire_params(    &       ! also on restart
            ifds,ifde,jfds,jfde, &
            ifms,ifme,jfms,jfme, &
            ifts,ifte,jfts,jfte, &
            fdx,fdy,nfuel_cat0,  &
            nfuel_cat,fuel_time, &
            fp) 

        ! initialize model state to no fire / set tign_in
        if(.not.restart)then
            call init_no_fire  ( &
            ifds,ifde,jfds,jfde, &
            ifms,ifme,jfms,jfme, &
            ifts,ifte,jfts,jfte, &
            fdx,fdy,time_start,dt,  &
            fuel_frac,fire_area,lfn,tign_in,tign)
            
        endif

        call check_lfn_tign('after initialization:',time_start,ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme,lfn,tign)

        if(replay) then
          call message('replay, recomputing fuel fraction')
          call fuel_left( &
            ifds,ifde,jfds,jfde, &
            ifms,ifme,jfms,jfme, &
            ifts,ifte,jfts,jfte, &
            ifms,ifme,jfms,jfme, &
            lfn,tign,fuel_time,time_start,fuel_frac,fire_area) !fuel_frac is global
          call check_lfn_tign('recomputed fuel fraction',time_start,ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme,lfn,tign)
        endif

elseif(ifun.eq.3)then   ! ignition if so specified

    
elseif (ifun.eq.4) then  ! do the timestep
     
    ! #35 call check_lfn_tign('time step start',time_start,ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme,lfn, tign)
    if(run_fuel_moisture)then
        ! uses nfuel_cat to set the other fuel data arrays
        ! needs zsf on halo width 1 to compute the terrain gradient
        call set_fire_params(    &       ! also on restart
            ifds,ifde,jfds,jfde, &
            ifms,ifme,jfms,jfme, &
            ifts,ifte,jfts,jfte, &
            fdx,fdy,nfuel_cat0,  &
            nfuel_cat,fuel_time, &
            fp) 
    endif

    if(fire_print_msg.ge.stat_lev)then
      aw=fun_real(RNRM_SUM,  &
        ifms,ifme,1,1,jfms,jfme, &                ! memory dims
        ifds,ifde,1,1,jfds,jfde, &                ! domain dims
        ifts,ifte,1,1,jfts,jfte, &                ! patch or tile dims
        0,0,0,       &                            ! staggering
        fp%vx,fp%vy)/((ifde-ifds+1)*(jfde-jfds+1))
      mw=fun_real(RNRM_MAX,  &
        ifms,ifme,1,1,jfms,jfme, &                ! memory dims
        ifds,ifde,1,1,jfds,jfde, &                ! domain dims
        ifts,ifte,1,1,jfts,jfte, &                ! patch or tile dims
        0,0,0,       &                            ! staggering
        fp%vx,fp%vy)
!$OMP MASTER 
      write(msg,91)time_start,'Average surface wind',aw,'m/s'
      call message(msg,stat_lev)
      write(msg,91)time_start,'Maximum surface wind',mw,'m/s'
      call message(msg,stat_lev)
!$OMP END MASTER 
    endif

!   compute fuel fraction at start
!    call fuel_left( &
!        ifms,ifme,jfms,jfme, &
!        ifts,ifte,jfts,jfte, &
!        ifms,ifme,jfms,jfme, &
!        lfn,tign,fuel_time,time_start,fuel_frac,fire_area) ! fuel frac is shared

    call print_2d_stats(ifts,ifte,jfts,jfte, &
                   ifms,ifme,jfms,jfme, &
                   fuel_frac,'model: fuel_frac start')

    ! advance the model from time_start to time_start+dt
    ! return the fuel fraction burnt this call in each fire cell
    ! will call module_fr_sfire_speed::normal_spread for propagation speed
    ! We cannot simply compute the spread rate here because that will change with the
    ! angle of the wind and the direction of propagation, thus it is done in subroutine
    ! normal_spread at each fire time step. Instead, we pass arguments that 
    ! the speed function may use as fp. 

!   propagate level set function in time
!   set lfn_out tign
!   lfn does not change, tign has no halos

    if(.not. freeze_fire)then
      if(.not.replay) then    
        call prop_ls(id,1,     &
        ifds,ifde,jfds,jfde,                      & ! fire domain dims - the whole domain
        ifms,ifme,jfms,jfme,                      &
        ifps,ifpe,jfps,jfpe, &                ! patch - nodes owned by this process
        ifts,ifte,jfts,jfte,                      &
        time_start,dt,fdx,fdy,tbound,  &
        lfn,lfn_out,tign,ros, fp &
        ) 
      else
        do j=jfts,jfte
          do i=ifts,ifte
            lfn_out(i,j) = tign(i,j) - (time_start + dt) ! <0 if burning at the end of the time step
          enddo
        enddo
      endif
    else
        call message('sfire_model: EXPERIMENTAL: skipping fireline propagation')

    endif

    call print_2d_stats(ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme, ros,'model: ros')

elseif (ifun.eq.5) then ! copy the result of timestep back to input
    ! this cannot be done in the time step itself because of race condition
    ! some thread may still be using lfn as input in their tile halo

    ! fix lfn_out if needed and compute tign
    if(.not. freeze_fire)then
      if(.not.replay) then    
        call prop_ls(id,2,     &
        ifds,ifde,jfds,jfde,                      & ! fire domain dims - the whole domain
        ifms,ifme,jfms,jfme,                      &
        ifps,ifpe,jfps,jfpe, &                ! patch - nodes owned by this process
        ifts,ifte,jfts,jfte,                      &
        time_start,dt,fdx,fdy,tbound,  &
        lfn,lfn_out,tign,ros, fp &
        ) 
      endif
    endif

    call check_lfn_tign('time step end',time_start+dt,ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme,lfn_out, tign)

    if(.not. freeze_fire)then
    do j=jfts,jfte
        do i=ifts,ifte
            lfn(i,j)=lfn_out(i,j)
            ! if want to try timestep again treat tign the same way here
            ! even if tign does not need a halo
        enddo
    enddo

    call check_lfn_tign('before ignition',time_start+dt,ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme,lfn, tign)

    ! check for ignitions
    if( fire_tign_in_time > 0.) then
         if (ignition%num_lines >  0) then
              call crash('ignition from lines and from tign_in are not compatible')
         endif
         call ignite_from_tign_in(                                   &
                        ifds,ifde,jfds,jfde,                      & ! fire domain dims - the whole domain
                        ifms,ifme,jfms,jfme,                      &
                        ifts,ifte,jfts,jfte,                      &
                        time_start,time_start+dt,                 &
                        tign_in,                                  &
                        lfn,tign,ignited)
         call check_lfn_tign('after ignite_from_tign_in',time_start+dt,ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme,lfn, tign)
    endif
    

    do ig = 1,ignition%num_lines
    
!  for now, check for ignition every time step...
!        if(ignition%line(ig)%end_time>=time_start.and.ignition%line(ig)%start_time<time_start+dt)then 
            call ignite_fire(                             &
                ifds,ifde,jfds,jfde,                      & ! fire domain dims - the whole domain
                ifms,ifme,jfms,jfme,                      &
                ifts,ifte,jfts,jfte,                      &
                ignition%line(ig),                        &
                time_start,time_start+dt,                 &
                coord_xf,coord_yf,ignition%unit_fxlong,ignition%unit_fxlat,        & 
                lfn,tign,ignited)
                
#ifdef DEBUG_OUT    
            call write_array_m(ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme,lfn,'lfn_ig',id)
            call write_array_m(ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme,coord_xf,'coord_xf_ig',id)
            call write_array_m(ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme,coord_yf,'coord_yf_ig',id)
#endif
!        endif
        
    enddo
    else
        call message('sfire_model: EXPERIMENTAL: skipping ignition')
    endif
            
    call print_2d_stats(ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme, &
                   lfn,'sfire_model: lfn out')

    call check_lfn_tign('after ignition',time_start+dt,ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme,lfn, tign)

elseif (ifun.eq.6) then ! timestep postprocessing
    ! have halo on lfn now

    ! diagnostics
    call fire_intensity(fp,                   &  ! fuel properties
        ifms,ifme,jfms,jfme,                      &  ! memory dims
        ifts,ifte,jfts,jfte,                      &  ! tile dims
        ifms,ifme,jfms,jfme,                      &  ! ros dims
        ros,nfuel_cat,                            & !in
        flineint,flineint2)                       ! fireline intensities out

    if(fireline_mask < 0.) then
      do j=jfts,jfte
        do i=ifts,ifte
          ! if the sign of lfn is the same as all of its neighbors or we are at domain boundary, we are not near the fireline
          if( (lfn(i-1,j-1)>0. .and. lfn(i-1,j)>0. .and. lfn(i,j-1)>0. .and. lfn(i,j)>0.  .and. &
            lfn(i+1,j+1)>0. .and. lfn(i+1,j)>0. .and. lfn(i,j+1)>0. ) .or.                   &
           (lfn(i-1,j-1)<0. .and. lfn(i-1,j)<0. .and. lfn(i,j-1)<0. .and. lfn(i,j)<0.  .and. &
            lfn(i+1,j+1)<0. .and. lfn(i+1,j)<0. .and. lfn(i,j+1)<0. ) .or.                   &
            i.eq.ifds .or. i .eq. ifde .or. j.eq.jfds .or. j.eq.jfde) then
               ros(i,j)=fireline_mask
               flineint(i,j)=fireline_mask
               flineint2(i,j)=fireline_mask
           endif
        enddo
      enddo
    endif
    
  call fire_risk(fp,                              &
        ifms,ifme,jfms,jfme,                      &  ! memory dims
        ifts,ifte,jfts,jfte,                      &  ! tile dims
        nfuel_cat,                                &  !
        f_ros0,f_rosx,f_rosy,f_ros,               &  ! fire spread 
        f_int,f_lineint,f_lineint2)                  ! fire intensities for danger rating


  select case(fire_hfx_given)
     
    case(0)   ! normal

    ! compute the heat fluxes from the fuel burned
    ! needs lfn and tign from neighbors so halo must be updated before
!$OMP CRITICAL(SFIRE_MODEL_CRIT)
    write(msg,*)'time_start=',time_start,' dt=',dt,' before fuel_left'
!$OMP END CRITICAL(SFIRE_MODEL_CRIT)
    call message(msg)
    call print_2d_stats(ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme,lfn,'model: lfn')
    call print_2d_stats(ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme,tign,'model: tign')
    call check_lfn_tign('before fuel_left',time_start+dt,ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme,lfn, tign)

    call print_2d_stats(ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme,fuel_time,'model: fuel_time')

    if(fire_update_fuel_frac.eq.1)then
        ! run model the normal way
        call fuel_left(&
            ifds,ifde,jfds,jfde, &
            ifms,ifme,jfms,jfme, &
            ifts,ifte,jfts,jfte, &
            ifts,ifte,jfts,jfte, &
            lfn,tign,fuel_time,time_start+dt,fuel_frac_end,fire_area) !fuel_frac_end is private and tile based
    
        call print_2d_stats(ifts,ifte,jfts,jfte,ifts,ifte,jfts,jfte,fuel_frac_end,'model: fuel_frac end')
        
        do j=jfts,jfte
            do i=ifts,ifte
                t = min(fuel_frac(i,j),fuel_frac_end(i,j)) ! do not allow fuel fraction to increase, in case of approximation error 
                fuel_frac_burnt(i,j)=fuel_frac(i,j)-t ! fuel lost this timestep
                fuel_frac(i,j)=t ! copy new value to state array
            enddo
        enddo
    elseif(fire_update_fuel_frac.eq.2)then
        do j=jfts,jfte
            do i=ifts,ifte
                if(lfn(i,j)<0.)then
                    fuel_frac_burnt(i,j) = dt / fuel_time(i,j)
                else
                    fuel_frac_burnt(i,j) = 0.
                endif
            enddo
        enddo
    else
        call crash('fire_update_fuel_frac value not supported')
    endif

    call print_2d_stats(ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme,fuel_frac_burnt,'model: fuel_frac burned')
        
    call heat_fluxes(dt,fp,                       &
        ifms,ifme,jfms,jfme,                      &
        ifts,ifte,jfts,jfte,                      &
        ifms,ifme,jfms,jfme,                      &  ! fuel_frac_burned has standard memory dimensions
        fp%fgip,                                     &
        fuel_frac_burnt,                          & !
        fgrnhfx,fgrnqfx)                              !out

    case (1, 2) 
!$OMP CRITICAL(SFIRE_MODEL_CRIT)
       write(msg,*)"model: expecting fire_hfx to be set in WRF, from wrfinput or wrfrst files"
       call message(msg)
!$OMP END CRITICAL(SFIRE_MODEL_CRIT)

       do j=jfts,jfte
           do i=ifts,ifte
              fgrnhfx(i,j) = (1. - fire_hfx_latent_part)*fire_hfx(i,j)
              fgrnqfx(i,j) =       fire_hfx_latent_part *fire_hfx(i,j)
           enddo
       enddo

    case (3)
    
       call message('artificial heat flux from parameters given in namelist.input')

       call param_hfx( time_start, &
                            ifms,ifme,jfms,jfme, &
                            ifts,ifte,jfts,jfte, &
                            coord_xf,coord_yf,   &
                            hfx,                 &
                            fire_area,fgrnhfx,fgrnqfx)

   case default
        call crash('bad fire_hfx_given')
   end select

   ! this should run in any case

    if(fire_print_msg.ge.stat_lev)then
      tfa=fun_real(REAL_SUM,  &
        ifms,ifme,1,1,jfms,jfme, &                ! memory dims
        ifds,ifde,1,1,jfds,jfde, &                ! domain dims
        ifts,ifte,1,1,jfts,jfte, &                ! patch or tile dims
        0,0,0,       &                            ! staggering
        fire_area,fire_area) * fdx * fdy
      thf=fun_real(REAL_SUM,  &
        ifms,ifme,1,1,jfms,jfme, &                ! memory dims
        ifds,ifde,1,1,jfds,jfde, &                ! domain dims
        ifts,ifte,1,1,jfts,jfte, &                ! patch or tile dims
        0,0,0,       &                            ! staggering
        fgrnhfx,fgrnhfx) * fdx * fdy
      mhf=fun_real(REAL_MAX,  &
        ifms,ifme,1,1,jfms,jfme, &                ! memory dims
        ifds,ifde,1,1,jfds,jfde, &                ! domain dims
        ifts,ifte,1,1,jfts,jfte, &                ! patch or tile dims
        0,0,0,       &                            ! staggering
        fgrnhfx,fgrnhfx) 
      tqf=fun_real(REAL_SUM,  &
        ifms,ifme,1,1,jfms,jfme, &                ! memory dims
        ifds,ifde,1,1,jfds,jfde, &                ! domain dims
        ifts,ifte,1,1,jfts,jfte, &                ! patch or tile dims
        0,0,0,       &                            ! staggering
        fgrnqfx,fgrnqfx) * fdx * fdy
      mqf=fun_real(REAL_MAX,  &
        ifms,ifme,1,1,jfms,jfme, &                ! memory dims
        ifds,ifde,1,1,jfds,jfde, &                ! domain dims
        ifts,ifte,1,1,jfts,jfte, &                ! patch or tile dims
        0,0,0,       &                            ! staggering
        fgrnqfx,fgrnqfx) 
!$OMP MASTER 
      write(msg,91)time_start,'Fire area           ',tfa,'m^2'
      call message(msg,stat_lev)
      write(msg,91)time_start,'Heat output         ',thf,'W'
      call message(msg,stat_lev)
      write(msg,91)time_start,'Max heat flux       ',mhf,'W/m^2'
      call message(msg,stat_lev)
      write(msg,91)time_start,'Latent heat output  ',tqf,'W'
      call message(msg,stat_lev)
      write(msg,91)time_start,'Max latent heat flux',mqf,'W/m^2'
      call message(msg,stat_lev)
!$OMP END MASTER
91  format('Time ',f11.3,' s ',a,e12.3,1x,a)
    endif
        

   call print_2d_stats(ifts,ifte,jfts,jfte, &
                   ifms,ifme,jfms,jfme, &
                   fgrnhfx,'model: heat flux(J/m^2/s)')

else
!$OMP CRITICAL(SFIRE_MODEL_CRIT)
    write(msg,*)'sfire_model: bad ifun=',ifun
!$OMP END CRITICAL(SFIRE_MODEL_CRIT)
    call crash(msg)
endif

end subroutine sfire_model

       subroutine param_hfx( time_now,&
                            ifms,ifme,jfms,jfme, &
                            ifts,ifte,jfts,jfte, &
                            coord_xf,coord_yf,   &
                            hfx,                 &
                            fire_area,fgrnhfx,fgrnqfx)
!***   generate artifical heat flux from a parametric description
!***   arguments
       real, intent(in)::time_now
       integer, intent(in):: & 
                            ifms,ifme,jfms,jfme, &
                            ifts,ifte,jfts,jfte  
       type(lines_type), intent(in)::hfx
       real, dimension(ifms:ifme,jfms:jfme), intent(in)::coord_xf,coord_yf ! nodal coordinates
       real, dimension(ifms:ifme,jfms:jfme), intent(out)::fire_area,fgrnhfx,fgrnqfx ! the imposed heat flux
       character(len=128):: msg
!***   local
       integer::i,j,k,nfa,ncells
       real:: d2,ax,ay,hfrac,fa,thfx,t,r,radius
       real, parameter:: sigma_mult=3.   ! 3 gives drop to 1% in trans. time from gaussian kernel 
       real:: maxspeed=100  ! max speed how the circle can move

       do j=jfts,jfte  ! zero out the outputs
           do i=ifts,ifte
              fire_area(i,j)=0
              fgrnhfx(i,j) = 0.
              fgrnqfx(i,j) = 0.
           enddo
       enddo

       do k=1,hfx%num_lines
          if(hfx%line(k)%radius > 0.)then
             ! processing heatflux line i
             ! find the time multiplier
             t = max(hfx%line(k)%start_time - time_now, time_now - hfx%line(k)%end_time)
             if(t > 0.)then  ! postitive distance from the time interval
                 r = t / hfx%line(k)%trans_time  ! position in the transition - 1 is at transition distance 
                 hfrac = exp(-(sigma_mult * r)**2/2.)  ! gaussian kernel
             else
                 hfrac = 1.0
             endif
             

             ! find the coordinates of the center of the heat flux circle now
             ax = safe_prop(time_now, &
                              hfx%line(k)%start_time,&
                              hfx%line(k)%end_time,&
                              hfx%line(k)%start_x,&
                              hfx%line(k)%end_x, &
                              hfx%unit_fxlong*maxspeed)
             ay = safe_prop(time_now,&
                              hfx%line(k)%start_time,&
                              hfx%line(k)%end_time,&
                              hfx%line(k)%start_y,&
                              hfx%line(k)%end_y, &
                              hfx%unit_fxlat*maxspeed)

             radius=hfx%line(k)%radius

!$OMP CRITICAL(SFIRE_MODEL_CRIT)
             write(msg,*)'hfx line ',i,' at ',time_now,'s ',hfrac,' of max ', hfx%line(k)%hfx_value
             call message(msg)
             write(msg,*)'center ',ax,ay,' radius ',radius
             call message(msg)
!$OMP END CRITICAL(SFIRE_MODEL_CRIT)

             nfa=0
             ncells=0
             do j=jfts,jfte
                 do i=ifts,ifte
                     ! distance squared from the center
                     d2 = (hfx%unit_fxlong*(ax - coord_xf(i,j)))**2 +  (hfx%unit_fxlat*(ay - coord_yf(i,j)))**2
                     if(d2 < radius**2)then
                         fa=1.
                     else
                         fa=0.
                     endif
                     ! set heat fluxes
                     thfx= hfx%line(k)%hfx_value * hfrac * fa  ! total heat flux at this point
                     fgrnhfx(i,j)= fgrnhfx(i,j) + (1.-fire_hfx_latent_part) * thfx
                     fgrnqfx(i,j)= fgrnqfx(i,j) + fire_hfx_latent_part * thfx
                     ! set fire area
                     fire_area(i,j) = max(fire_area(i,j),fa)
                     ! set stats
                     nfa=nfa+fa;
                     ncells=ncells+1
                 enddo
             enddo

!$OMP CRITICAL(SFIRE_MODEL_CRIT)
             write(msg,*)'Number of cells in fire area in this tile ',nfa,' ',(100.*nfa)/ncells,' %'
             call message(msg)
!$OMP END CRITICAL(SFIRE_MODEL_CRIT)

         endif
      enddo
 end subroutine param_hfx
                          

                 
real function safe_prop(t,t1,t2,x1,x2,ms)
!  return x between x1 and x2 in the same proportion as is t between t1 and t2, safe in the case when t1=t2 and x1=x2
!  safe_prop = x1 + (t-t1)*(x2-x1)/(t2-t1) 
!  future: abs((x2-x1)/(t2-t1)) capped at ms but still return x1 when t=t1 and x2 when t=t2
real, intent(in)::t,t1,t2,x1,x2,ms
real:: p,x
       if(t2 < t1)call crash('safe_prop: must have t2>t1')
       if(.not.ms>0.)call crash('safe_prop: must have ms>0')
       if(t1 .eq. t2)then
           if(x1.eq.x2)then
               x=x1
           else
               call crash('safe_prop: infinite speed')
           endif
       else
           p = (t - t1)/(t2 - t1)   ! 0 at t1, 1 at t2
           x = x1*(1.-p) + x2*p
       endif
       safe_prop=x
end function safe_prop 
!
!*****************
!
            
end module module_fr_sfire_model
