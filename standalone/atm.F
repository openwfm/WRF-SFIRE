module module_model_test
use module_fr_sfire_util, only: message,crash
implicit none

contains 

subroutine sub_main

!*** purpose: standalone driver with compatible files to WRF-Fire

use module_fr_sfire_driver, only: sfire_driver_em, use_atm_vars
use module_domain, only: domain
use module_configure, only: grid_config_rec_type,read_namelist
use wrf_netcdf, only : grid_info, set_info_from_file, &
                       create_output_file,write_vars,output_vars, &
                       input_vars_fire,read_vars_fire
implicit none

!*** local

! arguments to SFIRE

type(domain)::grid          ! all: state+inputs+outputs, compatible with wrf
TYPE (grid_config_rec_type):: config_flags ! the namelist
integer::  &                ! atmosphere mesh dimensions, for compatibility
    ids,ide, kds,kde, jds,jde,& 
    ims,ime, kms,kme, jms,jme,&
    ips,ipe, kps,kpe, jps,jpe
integer::  &                ! fire mesh dimensions
    ifds,ifde,jfds,jfde, &  ! the physical domain
    ifps,ifpe,jfps,jfpe, &  ! patch - assigned to one process. Here the same as domain.
    ifms,ifme,jfms,jfme     ! memory allocated, needs a strip around the patch

! I/O interface
character(len=*),parameter::inputfile='fire_input.nc'
character(len=*),parameter::outputfile='fire_output.nc'
type(output_vars)::output                ! output arrays
type(input_vars_fire)::input                  ! input arrays_fire

! other derived types
type(grid_info)::info                    ! dimensions, grid controls

! scalars
integer:: nx,ny,nz,nfx,nfy,nfz,nsteps,itimestep,sr_x,sr_y
integer::nhalo=5
double precision:: dt, duration_s  ! may need more accurate time computation to get the number of timesteps right
real:: time,time_step_start, dts
logical::do_ouput

!*** executable

call read_namelist(config_flags)           ! read flags from namelist.input
call set_info_from_file(inputfile,info)     ! get dimensions

! set dimensions
!nx=info%natmx
!ny=info%natmy
!nz=info%natmz
nfx=info%nfirex
nfy=info%nfirey
sr_x=info%sr_x
sr_y=info%sr_y

!write(6,2)' atmospheric domain size               ',config_flags%e_we,config_flags%e_sn,config_flags%e_vert
!if(config_flags%e_we.ne.nx+1 .or. config_flags%e_sn.ne.ny+1 .or. config_flags%e_vert.ne.nz+1)then
!  write(6,*)'dimensions in input file are ',nx,ny,nz,' must be 1 less than in namelist.input'
!  call crash('inconsistent input files')
!endif


! set atm domain size
!ids=1
!ide=nx+1
!jds=1
!jde=ny+1
!kds=1
!kde=nz+1
!ims=ids-nhalo
!ime=ide+nhalo
!jms=jds-nhalo
!jme=jde+nhalo
!kms=kds
!kme=kde
!ips=ids
!ipe=ide
!jps=jds
!jpe=jde
!kps=kds
!kpe=kde
!
!write(6,2)'fire_mesh refinement ratio  ',info%sr_x,info%sr_y

write(6,2)'fire domain size            ',nfx,nfy

!if(nfx.ne.sr_x*nx .or. nfy.ne.sr_y*ny)then
!   call crash('fire domain size is not multiple of atmospheric domain size')
!endif


! set fire domain size
ifds=1
ifde=nfx+sr_x
jfds=1
jfde=nfy+sr_y
ifms=ifds-sr_x*nhalo
ifme=ifde+sr_x*nhalo
jfms=ifds-sr_y*nhalo
jfme=ifde+sr_y*nhalo
ifps=1
ifpe=nfx
jfps=1
jfpe=nfy

!write(6,2)'atmospheric domain dimensions',ids,ide,jds,jde,kds,kde
!write(6,2)'atmospheric memory dimensions',ims,ime,jms,jme,kms,kme
write(6,2)'fire domain dimensions       ',ifds,ifde,jfds,jfde
write(6,2)'fire memory dimensions       ',ifms,ifme,jfms,jfme
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

! state
call allocate2d(grid%bbb,ifms,ifme,jfms,jfme,'bbb')             ! spread formula coeff
call allocate2d(grid%betafl,ifms,ifme,jfms,jfme,'betafl')          ! spread formula coeff
call allocate2d(grid%phiwc,ifms,ifme,jfms,jfme,'phiwc')           ! spread formula coeff
call allocate2d(grid%r_0,ifms,ifme,jfms,jfme,'r_0')             ! spread formula coeff
call allocate2d(grid%fgip,ifms,ifme,jfms,jfme,'fgip')            ! spread formula coeff
call allocate2d(grid%ischap,ifms,ifme,jfms,jfme,'ischap')          ! spread formula coeff
call allocate2d(grid%fuel_time,ifms,ifme,jfms,jfme,'fuel_time')        ! 
call allocate2d(grid%lfn,ifms,ifme,jfms,jfme,'lfn') 
call allocate2d(grid%tign_g,ifms,ifme,jfms,jfme,'tign') 
call allocate2d(grid%fuel_frac,ifms,ifme,jfms,jfme,'fuel_frac') 
call allocate2d(grid%fuel_frac_burnt,ifms,ifme,jfms,jfme,'fuel_frac_burnt') 
call allocate2d(grid%fire_area,ifms,ifme,jfms,jfme,'fire_area') 
call allocate2d(grid%lfn_out,ifms,ifme,jfms,jfme,'lfn_out') 
call allocate2d(grid%ros,ifms,ifme,jfms,jfme,'ros') 

! output
call allocate2d(grid%fgrnhfx,ifms,ifme,jfms,jfme,'fgrnhfx')          ! 
call allocate2d(grid%fgrnqfx,ifms,ifme,jfms,jfme,'fgrnqfx')          ! 
call allocate2d(grid%fcanhfx,ifms,ifme,jfms,jfme,'fcanhfx')          ! 
call allocate2d(grid%fcanqfx,ifms,ifme,jfms,jfme,'fcanqfx')          ! 

! atmosphere input/output compatibility arrays
!call allocate2d(grid%avg_fuel_frac,ims,ime,jms,jme,'avg_fuel_frac') 
!call allocate2d(grid%xlong,ims,ime,jms,jme,'xlong') 
!call allocate2d(grid%xlat,ims,ime,jms,jme,'xlat') 
!call allocate2d(grid%ht,ims,ime,jms,jme,'ht') 
!call allocate2d(grid%z0,ims,ime,jms,jme,'z0') 
!call allocate2d(grid%grnhfx,ims,ime,jms,jme,'grnhfx') 
!call allocate2d(grid%grnqfx,ims,ime,jms,jme,'qrnqfx') 
!call allocate2d(grid%canhfx,ims,ime,jms,jme,'canhfx')
!call allocate2d(grid%canqfx,ims,ime,jms,jme,'canqfx') 
!call allocate3d(grid%ph_2,ims,ime,kms,kme,jms,jme,'ph_2') 
!call allocate3d(grid%phb,ims,ime,kms,kme,jms,jme,'phb') 
!call allocate3d(grid%u_2,ims,ime,kms,kme,jms,jme,'u_2') 
!call allocate3d(grid%v_2,ims,ime,kms,kme,jms,jme,'v_2') 
!call allocate2d(grid%uah,ims,ime,jms,jme,'uah')                  ! 
!call allocate2d(grid%vah,ims,ime,jms,jme,'vah')                  ! 

! prepare reading input file
input%nfuel_cat=>grid%nfuel_cat(1:nfx,1:nfy)
input%zsf      =>grid%zsf(1:nfx,1:nfy)
input%dzdxf    =>grid%dzdxf(1:nfx,1:nfy)
input%dzdyf    =>grid%dzdyf(1:nfx,1:nfy)
input%vf       =>grid%uf(1:nfx,1:nfy)
input%uf       =>grid%vf(1:nfx,1:nfy)
input%lfn      =>grid%lfn(1:nfx,1:nfy)
input%tign     =>grid%tign_g(1:nfx,1:nfy)
input%fuel_frac   =>grid%fuel_frac(1:nfx,1:nfy)
!call allocate3d(input%ph  ,1,nx,1,ny,1,nz+1,'ph') 
!call allocate3d(input%phb, 1,nx,1,ny,1,nz+1,'phb') 
!call allocate3d(input%u,   1,nx+1,1,ny,1,nz,'u') 
!call allocate3d(input%v,   1,nx,1,ny+1,1,nz,'v') 

call read_vars_fire(inputfile,info,input)

! sfire uses wrf ordering of dimensions
!grid%ph_2(1:nx,1:nz+1,1:ny)=reshape(input%ph, (/nx,nz+1,ny/),order=(/1,3,2/))
!grid%phb (1:nx,1:nz+1,1:ny)=reshape(input%phb,(/nx,nz+1,ny/),order=(/1,3,2/))
!deallocate(input%ph)
!deallocate(input%phb)
!grid%u_2(1:nx+1,1:nz,1:ny)=reshape(input%u,(/nx+1,nz,ny/),order=(/1,3,2/))
!grid%v_2(1:nx,1:nz,1:ny+1)=reshape(input%v,(/nx,nz,ny+1/),order=(/1,3,2/))
!deallocate(input%u)
!deallocate(input%v)
!grid%z0(1:nx,1:ny)=input%z0
!deallocate(input%z0)

! time control
! NOTE: dt in the netcdf input file as returned in info%dt is WRONG !!
dt=config_flags%time_step
if(config_flags%time_step_fract_den.ne.0)then
  dt=dt+dble(config_flags%time_step_fract_num)/dble(config_flags%time_step_fract_den)
endif
duration_s = config_flags%run_seconds           &
           + 60d0*(config_flags%run_minutes     &
           + 60d0*(config_flags%run_hours       &
           + 24d0*(config_flags%run_days)))       
nsteps = nint( duration_s / dt ) ! number of time steps

! divide up for shared memory parallel execution
call set_tiles(1,1,ips,ipe,jps,jpe,grid%num_tiles,grid%i_start,grid%i_end,grid%j_start,grid%j_end)

! set the scalars in grid type
grid%sr_x = sr_x
grid%sr_y = sr_y
grid%dt = dt
grid%dx = info%dx
grid%dy = info%dy
grid%itimestep=0
grid%u_frame=0.
grid%v_frame=0.

! start output file
call create_output_file(outputfile,info)

! initialize model
   time_step_start  =0.
   dts=dt
call sfire_driver_em ( grid , config_flags                          &
            ,time_step_start,dts                                    &
            ,1,2,0                                                 &
            ,ids,ide, kds,kde, jds,jde                              &
            ,ims,ime, kms,kme, jms,jme                              &
            ,ips,ipe, kps,kpe, jps,jpe                              &
            ,ifds,ifde, jfds,jfde                                   &
            ,ifms,ifme, jfms,jfme                                   &
            ,ifps,ifpe, jfps,jfpe )

do itimestep=1,nsteps
! run  model
   grid%itimestep = itimestep
   time_step_start = itimestep*dt
   dts=dt
   call sfire_driver_em ( grid , config_flags                          &
            ,time_step_start,dts                                    &
            ,3,6,0                                                  &
            ,ids,ide, kds,kde, jds,jde                              &
            ,ims,ime, kms,kme, jms,jme                              &
            ,ips,ipe, kps,kpe, jps,jpe                              &
            ,ifds,ifde, jfds,jfde                                   &
            ,ifms,ifme, jfms,jfme                                   &
            ,ifps,ifpe, jfps,jfpe )
    
    if(itimestep.le.10.or.mod(itimestep,10).eq.0)then
        time = dt*itimestep
        output%lfn=>grid%lfn(ifps:ifpe,jfps:jfpe)
        output%tign=>grid%tign_g(ifps:ifpe,jfps:jfpe)
        output%fgrnhfx=>grid%fgrnhfx(ifps:ifpe,jfps:jfpe)
        call write_vars(outputfile,output,time)
    endif
enddo

end subroutine sub_main

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
if(ips.ge.ipe.or.jps.ge.jpe)call crash('bad domain bounds')
num_tiles=itiles*jtiles
istep=(ipe-ips+itiles)/itiles
jstep=(jpe-jps+jtiles)/jtiles
do i=1,itiles
    do j=1,jtiles
        ij=j+(i-1)*jtiles
        i_start(ij)=min(ipe,ips+(i-1)*istep)
        i_end(ij)  =min(ipe,ips+(i  )*istep-1)
        j_start(ij)=min(jpe,jps+(j-1)*jstep)
        j_end(ij)  =min(jpe,jps+(j  )*jstep-1)
    enddo
enddo
call check_tiles(ips,ipe,jps,jpe,num_tiles,i_start,i_end,j_start,j_end)
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
real, pointer, intent(out), dimension(:,:)::p
integer, intent(in):: ims,ime,jms,jme
character(len=*),intent(in)::s
!*** local
integer::err
!*** executable
write(6,1) ims,ime,jms,jme,trim(s)
if(associated(p))call crash('already allocated')
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
real, pointer, intent(out), dimension(:,:,:)::p
integer, intent(in):: ims,ime,jms,jme,kms,kme
character(len=*),intent(in)::s
!*** local
integer::err
!*** executable
write(6,1) ims,ime,jms,jme,kms,kme,trim(s)
1 format('allocate3d',3(1x,i6,' :',i6),1x,a)
if(associated(p))call crash('already allocated')
allocate(p(ims:ime,jms:jme,kms:kme),stat=err)
if(err.ne.0)then
   write(6,1)ims,ime,jms,jme,kms,kme,trim(s)
   call crash('memory allocation failed')
endif
p=0.
end subroutine allocate3d


end module module_model_test

!
!******************************
!


program model_test_main
use module_model_test, only: sub_main
call  sub_main
end program model_test_main
