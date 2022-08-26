module module_fire_standalone

use module_fr_sfire_driver, only: set_flags, fire_ignition_convert, &
                                  set_fp_from_grid
use module_fr_sfire_util, only: message,crash, &
          lines_type, print_2d_stats
use module_fr_sfire_phys, only: fire_params, init_fuel_cats
use module_fr_sfire_model, only: sfire_model
use module_domain, only: domain
use module_configure, only: grid_config_rec_type,read_namelist
use wrf_netcdf, only : grid_info, read_info, &
                       create_output_file,write_vars, &
                       read_vars, debug_print
implicit none

contains 

subroutine sub_main

!*** purpose: standalone driver with compatible files to WRF-Fire

implicit none

!*** local

! arguments to SFIRE

type(domain)::grid          ! all: state+inputs+outputs, compatible with wrf
TYPE (grid_config_rec_type):: config_flags ! the namelist
integer::  &                ! fire mesh dimensions
    ifds,ifde,jfds,jfde, &  ! the physical domain
    ifps,ifpe,jfps,jfpe, &  ! patch - assigned to one process. Here the same as domain.
    ifts,ifte,jfts,jfte, &  ! memory allocated, needs a strip around the patch
    ifms,ifme,jfms,jfme     ! memory allocated, needs a strip around the patch

! I/O interface
character(len=*),parameter::inputfile='fire_input.nc'
character(len=*),parameter::outputfile='fire_output.nc'
real, pointer, dimension(:,:) ::  uf1, vf1, uf2, vf2, fmc_g1, fmc_g2   ! stored input fields

! other derived types
type(grid_info)::info                    ! dimensions, grid controls

! scalars
integer:: nsteps,itimestep,ifun_start,ifun_end,id,ifun,iframe,istep
integer::nhalo=5
double precision:: dt_double,dt_config,duration_s,frame_s  ! may need more accurate time computation to get the number of timesteps right
real:: time_start,dt,t
logical::do_ouput
TYPE(lines_type) :: ignition, hfx
type(fire_params)::fp
logical::restart=.false.,replay=.false.,uniform=.false.
integer::iframe_start,iframe_end
logical::run_fuel_moisture=.false.

!*** executable

call read_namelist(config_flags)           ! read flags from namelist.input
call set_flags(config_flags)               ! copy configuration flags to sfire internal structures

debug_print = config_flags%fire_print_msg.ge.2 ! if we write a lot
debug_print = config_flags%fire_print_msg.ge.1 ! if we write a lot

call read_info(inputfile,info)             ! get dimensions

! start empty NetCDF file with the dimensions
! not here, may want to overwrite something
! call create_output_file(outputfile,info)



! set dimensions
ifds=1
ifde=info%nfirex
jfds=1
jfde=info%nfirey
ifms=ifds-nhalo
ifme=ifde+nhalo
jfms=ifds-nhalo
jfme=ifde+nhalo
ifps=1
ifpe=ifde
jfps=1
jfpe=jfde
ifts=1
ifte=ifde
jfts=1
jfte=jfde

write(6,2)'fire domain dimensions       ',ifds,ifde,jfds,jfde
write(6,2)'fire memory dimensions       ',ifms,ifme,jfms,jfme
write(6,2)'fire patch  dimensions       ',ifps,ifpe,jfps,jfpe
write(6,2)'fire tile   dimensions       ',ifts,ifte,jfts,jfte
2 format(a,6i6)

! allocate

! inputs
call allocate2d(grid%uf,ifms,ifme,jfms,jfme,'uf')              ! fire winds
call allocate2d(grid%vf,ifms,ifme,jfms,jfme,'vf')              ! fire winds
call allocate2d(grid%zsf,ifms,ifme,jfms,jfme,'zsf')             ! terrain height
call allocate2d(grid%dzdxf,ifms,ifme,jfms,jfme,'dzdxf')           ! terrain grad
call allocate2d(grid%dzdyf,ifms,ifme,jfms,jfme,'dzdyf')           ! terrain grad
call allocate2d(grid%fxlong,ifms,ifme,jfms,jfme,'fxlong')          ! 
call allocate2d(grid%fxlat,ifms,ifme,jfms,jfme,'fxlat')           !
call allocate2d(grid%nfuel_cat,ifms,ifme,jfms,jfme,'nfuel_cat')          ! 
call allocate2d(grid%fmc_g,ifms,ifme,jfms,jfme,'fmc_g')          ! 

! state
call allocate2d(grid%bbb,ifms,ifme,jfms,jfme,'bbb')             ! spread formula coeff
call allocate2d(grid%betafl,ifms,ifme,jfms,jfme,'betafl')          ! spread formula coeff
call allocate2d(grid%phiwc,ifms,ifme,jfms,jfme,'phiwc')           ! spread formula coeff
call allocate2d(grid%phisc,ifms,ifme,jfms,jfme,'phisc')           ! spread formula coeff
call allocate2d(grid%r_0,ifms,ifme,jfms,jfme,'r_0')             ! spread formula coeff
call allocate2d(grid%fgip,ifms,ifme,jfms,jfme,'fgip')            ! spread formula coeff
call allocate2d(grid%ischap,ifms,ifme,jfms,jfme,'ischap')          ! spread formula coeff
call allocate2d(grid%fuel_time,ifms,ifme,jfms,jfme,'fuel_time')        ! 
call allocate2d(grid%lfn,ifms,ifme,jfms,jfme,'lfn') 
call allocate2d(grid%tign_g,ifms,ifme,jfms,jfme,'tign_g') 
call allocate2d(grid%fuel_frac,ifms,ifme,jfms,jfme,'fuel_frac') 
call allocate2d(grid%fuel_frac_burnt,ifms,ifme,jfms,jfme,'fuel_frac_burnt') 
call allocate2d(grid%lfn_out,ifms,ifme,jfms,jfme,'lfn_out') 

! outputs
call allocate2d(grid%fire_area,ifms,ifme,jfms,jfme,'fire_area') 
call allocate2d(grid%ros,ifms,ifme,jfms,jfme,'ros') 
call allocate2d(grid%flineint,ifms,ifme,jfms,jfme,'flineint') 
call allocate2d(grid%flineint2,ifms,ifme,jfms,jfme,'flineint2') 
call allocate2d(grid%fgrnhfx,ifms,ifme,jfms,jfme,'fgrnhfx')          ! 
call allocate2d(grid%fgrnqfx,ifms,ifme,jfms,jfme,'fgrnqfx')          ! 
call allocate2d(grid%fcanhfx,ifms,ifme,jfms,jfme,'fcanhfx')          ! 
call allocate2d(grid%fcanqfx,ifms,ifme,jfms,jfme,'fcanqfx')          ! 
call allocate2d(grid%f_ros,ifms,ifme,jfms,jfme,'f_ros')              ! 
call allocate2d(grid%f_ros0,ifms,ifme,jfms,jfme,'f_ros0')            ! 
call allocate2d(grid%f_rosx,ifms,ifme,jfms,jfme,'f_rosx')            ! 
call allocate2d(grid%f_rosy,ifms,ifme,jfms,jfme,'f_rosy')            ! 
call allocate2d(grid%f_lineint,ifms,ifme,jfms,jfme,'f_lineint')      ! 
call allocate2d(grid%f_lineint2,ifms,ifme,jfms,jfme,'f_lineint2')    ! 
call allocate2d(grid%f_int,ifms,ifme,jfms,jfme,'f_int')              ! 

! local
call allocate2d(uf1,ifms,ifme,jfms,jfme,'uf1')              ! fire winds
call allocate2d(vf1,ifms,ifme,jfms,jfme,'vf1')              ! fire winds
call allocate2d(uf2,ifms,ifme,jfms,jfme,'uf2')              ! fire winds
call allocate2d(vf2,ifms,ifme,jfms,jfme,'vf2')              ! fire winds
call allocate2d(fmc_g1,ifms,ifme,jfms,jfme,'fmc_g1')              ! moisture 
call allocate2d(fmc_g2,ifms,ifme,jfms,jfme,'fmc_g2')              ! moisture 

! copy pointers to grid fields, to pass to the spread rate calculation
call set_fp_from_grid(grid,fp)
call init_fuel_cats(.true.)

! time control
! NOTE: dt in the netcdf input file as returned in info%dt is WRONG !!
dt_config=config_flags%time_step
if(config_flags%time_step_fract_den.ne.0)then
  dt_config=dt_config+dble(config_flags%time_step_fract_num)/dble(config_flags%time_step_fract_den)
endif
duration_s = config_flags%run_seconds           &
           + 60d0*(config_flags%run_minutes     &
           + 60d0*(config_flags%run_hours       &
           + 24d0*(config_flags%run_days)))       

if(config_flags%history_interval.ne.0)config_flags%history_interval_m=config_flags%history_interval
frame_s = config_flags%history_interval_s           &
           + 60d0*(config_flags%history_interval_m     &
           + 60d0*(config_flags%history_interval_h       &
           + 24d0*(config_flags%history_interval_d)))       

write(*,'(a,f10.2,a,f10.6,a,f10.6)')'from namelist.input: history_interval=',frame_s, &
     ' time_step=',dt_config,'s time steps in history frame ', frame_s / dt_config 

nsteps = nint( frame_s / dt_config ) ! number of time steps for the duration
dt_double = frame_s / nsteps
dt = dt_double

write(*,'(a,f10.2,a,f8.4,a,i6)')'adjusted: ',frame_s,' time_step=',dt,'s time steps in history frame ',nsteps

! divide up for shared memory parallel execution
!!call set_tiles(1,1,ips,ipe,jps,jpe,grid%num_tiles,grid%i_start,grid%i_end,grid%j_start,grid%j_end)

! set the scalars in grid type
grid%dt = dt
grid%itimestep=0
grid%xtime=0.
grid%u_frame=0.
grid%v_frame=0.
info%dt = dt     ! dt may be different than it was in the input file

! start output file
call create_output_file(outputfile,info)


if(info%ntimes.lt.3)then
  !call crash('need at least 3 steps')
  uniform=.true.
  call read_vars(inputfile,info,1,grid)
  iframe_start=1
  iframe_end=int(duration_s/frame_s)
else
  uniform=.false.
  call read_vars(inputfile,info,2,grid)
  iframe_start=3
  write(*,'(a)')'starting from frame 3 because need to interpolate between frames and frame 1 may be invalid'
  iframe_end=info%ntimes
  uf1=grid%uf
  vf1=grid%vf
  fmc_g1=grid%fmc_g
endif

print *,'Fire mesh:'
print *,'fxlat  lower bounds:',lbound(grid%fxlat)
print *,'fxlat  upper bounds:',ubound(grid%fxlat)
print *,'fxlat(1,1)=',grid%fxlat(1,1),' fxlat(',ifpe,',',jfpe,')=',grid%fxlat(ifpe,jfpe)
print *,'fxlong lower bounds:',lbound(grid%fxlong)
print *,'fxlong upper bounds:',ubound(grid%fxlong)
print *,'fxlong(1,1)=',grid%fxlong(1,1),' fxlong(',ifpe,',',jfpe,')=',grid%fxlong(ifpe,jfpe)
call print_2d_stats(ifps,ifpe,jfps,jfpe,ifms,ifme,jfms,jfme,grid%fxlong,'fire:fxlong')
call print_2d_stats(ifps,ifpe,jfps,jfpe,ifms,ifme,jfms,jfme,grid%fxlat,'fire:fxlat')

! get ignition data - should have fxlong fxlat now
call fire_ignition_convert (config_flags,ignition,                   &
             grid%fxlong, grid%fxlat,                                &
             ifds,ifde, jfds,jfde,                                   &
             ifms,ifme, jfms,jfme,                                   &
             ifps,ifpe, jfps,jfpe )

write(*,'(a,i5)')'number of frames in the file: ntimes=',info%ntimes
write(*,'(2(a,i5))')'using frames from',iframe_start,' to',iframe_end

itimestep = 0
ifun_start=1
do iframe=iframe_start,iframe_end ! interval ending with iframe
  if(.not.uniform)then
    call read_vars(inputfile,info,iframe,grid)
    uf2=grid%uf
    vf2=grid%vf
    fmc_g2=grid%fmc_g
  endif
  do istep=1,nsteps
    itimestep=info%ntimes * (iframe - 1) + istep
    grid%itimestep = itimestep
    grid%xtime = itimestep * grid%dt / 60.
    id=itimestep
    ifun_end=6
    ! interpolate time
    time_start = dt_double * (nsteps * (iframe - 1) + istep - 1)
    ! interpolate wind
    if(.not.uniform)then
      t = (istep - 1.)/real(nsteps)
      write(*,'(a,i4,a,i3,a,i8,a,f10.3,a,f10.3)')'frame',iframe,' istep',istep, &
         ' itimestep',itimestep, &
         ' start at ',time_start,'s weight ',t
      grid%uf = (1. - t)*uf1 + t*uf2
      grid%vf = (1. - t)*vf1 + t*vf2
      grid%fmc_g = (1. - t)*fmc_g1 + t*fmc_g2
    endif
 
    do ifun=ifun_start,ifun_end
  
      if(ifun.eq.4)then
        call print_2d_stats(ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme,fp%fmc_g,'fire:fmc_g')
      endif

      call sfire_model (                    &
        id,                                     & ! unique number for prints and debug
        ifun,                                   & ! what to do see below
        restart,replay,                         & ! use existing state; prescribe fire arrival time
        run_fuel_moisture,                      & ! run the moisture model
        config_flags%fire_fuel_read,config_flags%fire_fuel_cat,   & ! legacy initial constant fuel category
        ifds,ifde,jfds,jfde,                    & ! fire domain dims - the whole domain
        ifms,ifme,jfms,jfme,                    & ! fire memory dims - how declared
        ifps,ifpe,jfps,jfpe,                    & ! patch - nodes owned by this process
        ifts,ifte,jfts,jfte,                    & ! fire tile dims  - this thread
        time_start,dt,                          & ! time and increment
        info%fdx,info%fdy,                                & ! fire mesh spacing,
        ignition,hfx,                               & ! small array of ignition line descriptions
        grid%fxlong,grid%fxlat,                      & ! fire mesh coordinates
        grid%fire_hfx,                          & ! given heat flux (experimental)
        grid%tign_in,                           & ! ignition time, if given
        grid%lfn,grid%lfn_out,grid%tign_g,grid%fuel_frac,grid%fire_area,   & ! state: level function, ign time, fuel left, area burning
        grid%fuel_frac_burnt,                   &
        grid%fgrnhfx,grid%fgrnqfx,                          & ! output: heat fluxes
        grid%ros,grid%flineint,grid%flineint2,                 & ! diagnostic variables
        grid%f_ros0,grid%f_rosx,grid%f_rosy,grid%f_ros,             & ! fire risk spread
        grid%f_int,grid%f_lineint,grid%f_lineint2,             & ! fire risk intensities
        grid%nfuel_cat,                              & ! fuel data per point
        grid%fuel_time,grid%fwh,grid%fz0,                      & ! save derived internal data
        fp &
      )
    
    enddo
    ifun_start=3
  enddo
  call write_vars(outputfile,grid,info,iframe)
  if(.not.uniform)then
    uf1=uf2
    vf1=vf2
    fmc_g1=fmc_g2
  endif
enddo

end subroutine sub_main


!subroutine model_driver(grid,config_flags)

!
!******************************
!

subroutine set_tiles(itiles,jtiles,ips,ipe,jps,jpe,num_tiles,i_start,i_end,j_start,j_end)
!*** set tiles for standalone/testing
implicit none
!*** arguments
integer,intent(in)::itiles,jtiles,ips,ipe,jps,jpe
integer,intent(out)::num_tiles
integer,intent(out),dimension(itiles*jtiles)::i_start,i_end,j_start,j_end
!*** local
integer::i,j,istep,jstep,ij
character(len=128)::msg
write(msg,1)'patch',ips,':',ipe,jps,':',jpe
1 format(a,5x,i6,a,2i6,a,i6)
call message(msg,level=-1)
!if(ips.ge.ipe.or.jps.ge.jpe)call crash('bad domain bounds')
!num_tiles=itiles*jtiles
!istep=(ipe-ips+itiles)/itiles
!jstep=(jpe-jps+jtiles)/jtiles
!do i=1,itiles
!    do j=1,jtiles
!        ij=j+(i-1)*jtiles
!        i_start(ij)=min(ipe,ips+(i-1)*istep)
!        i_end(ij)  =min(ipe,ips+(i  )*istep-1)
!        j_start(ij)=min(jpe,jps+(j-1)*jstep)
!        j_end(ij)  =min(jpe,jps+(j  )*jstep-1)
!    enddo
!enddo
!call check_tiles(ips,ipe,jps,jpe,num_tiles,i_start,i_end,j_start,j_end)
end subroutine set_tiles


subroutine check_tiles(ips,ipe,jps,jpe,num_tiles,i_start,i_end,j_start,j_end)
implicit none
!*** purpose: check if tiles fit
!*** arguments
integer,intent(in)::ips,ipe,jps,jpe,num_tiles
integer,intent(in),dimension(num_tiles)::i_start,i_end,j_start,j_end
!*** local
character(len=128)::msg
integer:: ij,ie
!*** executable
if(num_tiles.lt.1)call crash('check_tiles: need at least one tile')
ie=0
do ij=1,num_tiles
    if(i_start(ij).lt.ips.or.i_end(ij).gt.ipe &
    .or.j_start(ij).lt.jps.or.j_end(ij).gt.jpe)then
        write(msg,1)'patch',ips,':',ipe,jps,':',jpe
1       format(a,5x,i6,a,2i6,a,i6)
        call message(msg,level=-1)
        write(msg,2)'tile',ij,i_start(ij),':',i_end(ij),j_start(ij),':',j_end(ij)
2       format(a,2i6,a,2i6,a,i6)
        call message(msg,level=-1)
        call crash('bad tile bounds')
    endif
enddo
end subroutine check_tiles


subroutine allocate2d(p,ims,ime,jms,jme,s) 
!*** allocate a pointer with error checking and initialization
implicit none
!*** arguments
real, pointer, dimension(:,:)::p
integer, intent(in):: ims,ime,jms,jme
character(len=*),intent(in)::s
!*** local
integer::err
!*** executable
if(debug_print)write(6,1) ims,ime,jms,jme,trim(s)
if(associated(p))call crash(trim(s) // ' already allocated')
1 format('allocate2d',2(1x,i6,' :',i6),1x,a)
allocate(p(ims:ime,jms:jme),stat=err)
if(err.ne.0)then
   write(6,1)ims,ime,jms,jme,trim(s)
   call crash('memory allocation failed')
endif
p=0.
end subroutine allocate2d

subroutine allocate3d(p,ims,ime,jms,jme,kms,kme,s) 
!*** allocate a pointer with error checking and initialization
implicit none
!*** arguments
real, pointer, dimension(:,:,:)::p
integer, intent(in):: ims,ime,jms,jme,kms,kme
character(len=*),intent(in)::s
!*** local
integer::err
!*** executable
if(debug_print)write(6,1) ims,ime,jms,jme,kms,kme,trim(s)
1 format('allocate3d',3(1x,i6,' :',i6),1x,a)
if(associated(p))call crash('already allocated')
allocate(p(ims:ime,jms:jme,kms:kme),stat=err)
if(err.ne.0)then
   write(6,1)ims,ime,jms,jme,kms,kme,trim(s)
   call crash('memory allocation failed')
endif
p=0.
end subroutine allocate3d

end module module_fire_standalone

!
!******************************
!


program fire
use module_fire_standalone, only: sub_main
call  sub_main
end program fire
