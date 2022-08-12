program fire_ros_main

use module_fr_sfire_phys
use module_fr_sfire_util

type(fire_params)::fp

real,allocatable,dimension(:,:):: propx,propy          ! direction of propagation 
real,allocatable,dimension(:,:):: vx,vy                ! wind velocity (m/s)
real,allocatable,dimension(:,:):: dzdxf,dzdyf          ! terrain grad (1)
real,allocatable,dimension(:,:):: fmc_g                ! fuel moisture contents, ground (1)
real,allocatable,dimension(:,:):: nfuel_cat            ! fuel category (integer values)
real,allocatable,dimension(:,:):: t2f                  ! temperature at 2m (K) (Balbi model only)
real,allocatable,dimension(:,:):: psfcf                ! surface pressure (Pa) (Balbi model only)
real,allocatable,dimension(:,:):: ffw                  ! fire front width (m) (Balbi model only)

integer:: n=10, iounit=9, ierr
character(len=50)::file_in

! read namelist.fire into fuels
! also create file fuels.m read and visualize the results as described in 
! https://wiki.openwfm.org/wiki/How_to_diagnose_fuel_properties_in_WRF-SFIRE

call init_fuel_cats(init_fuel_moisture=.false.) 

namelist /size/n
namelist /inputs/propx,propy,vx,vy,dzdxf,dzdyf,fmc_g,nfuel_cat,t2f,psfcf,ffw 

OPEN(iounit, FILE=trim(file_in),FORM='FORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=ierr)
if(ierr.ne.0)call crash('cannot open file '//trim(file_in)//' for reading')
read(iounit,size,iostat=ierr)
if(ierr.ne.0)call crash('cannot read namelist size in file '//trim(filein))
print *,'n=',n
allocate(propx(1,n),propy(1,n),vx(1,n),vy(1,n),dzdxf(1,n),dzdyf(1,n),fmc_g(1,n),nfuel_cat(1,n),t2f(1,n),psfcf(1,n),ffw(1,n))
read(iounit,inputs,iostat=ierr)
if(ierr.ne.0)call crash('cannot read namelist inputs in file '//trim(filein))
close(iounit,iostat=ierr)
if(ierr.ne.0)icall crash('cannot close file '//trim(file_in))
print *,'propx=',propx
print *,'propy=',propy
print *,'vx=',vx
print *,'vy=',vy
print *,'dzdxf=',dzdxf
print *,'dzdyf=',dzdyf
print *,'fmc_g=',fmc_g
print *,'nfuel_cat=',nfuel_cat
print *,'t2t=',t2t
print *,'psfc=',psfc
print *,'ffw=',ffw
fp.propx => propx
fp.propy => propy
fp.vx => vx
fp.vy => vy
fp.dzdxf => dzdxf
fp.dzdyf => dzdyf
fp.fmc_g => fmc_g
fp.t2f => t2f
fp.psfcf => psfcf
fp.ffw => ffw

    ! set up fuel arrays in fp
    !subroutine set_fire_params( &
    !                       ifds,ifde,jfds,jfde, &
    !                       ifms,ifme,jfms,jfme, &
    !                       ifts,ifte,jfts,jfte, &
    !                       fdx,fdy,nfuel_cat0,  &
    !                       nfuel_cat,fuel_time, &
    !                       fp )
call set_fire_params( &
                           1,1,1,n, &
                           1,1,1,n, &
                           1,1,1,n, &
                           0.,0.,0,  &
                           nfuel_cat,fuel_time, &
                           fp )

end program fire_ros_main
