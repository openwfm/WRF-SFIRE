program fire_ros_main

use module_fr_sfire_phys
use module_fr_sfire_util

implicit none
type(fire_params)::fp

real,pointer,dimension(:,:):: propx,propy          ! direction of propagation 
real,pointer,dimension(:,:):: vx,vy                ! wind velocity (m/s)
real,pointer,dimension(:,:):: dzdxf,dzdyf          ! terrain grad (1)
real,pointer,dimension(:,:):: fmc_g                ! fuel moisture contents, ground (1)
real,pointer,dimension(:,:):: nfuel_cat            ! fuel category (integer values)
real,pointer,dimension(:,:):: t2f                  ! temperature at 2m (K) (Balbi model only)
real,pointer,dimension(:,:):: psfcf                ! surface pressure (Pa) (Balbi model only)
real,pointer,dimension(:,:):: ffw                  ! fire front width (m) (Balbi model only)
real,pointer,dimension(:,:):: ros                  ! fire rate of spread (m/s) 


integer:: n=10, iounit=9, ierr, i, j
character(len=50)::file_in="namelist_standalone.input",file_out="namelist_standalone.output"
character(len=128)::msg

! read namelist.fire into fuels
! also create file fuels.m read and visualize the results as described in 
! https://wiki.openwfm.org/wiki/How_to_diagnose_fuel_properties_in_WRF-SFIRE

call init_fuel_cats(init_fuel_moisture=.false.) 

namelist /size/n
namelist /data/propx,propy,vx,vy,dzdxf,dzdyf,fmc_g,nfuel_cat,t2f,psfcf,ffw,ros

OPEN(iounit, FILE=trim(file_out),FORM='FORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=ierr)
if(ierr.ne.0)call crash('cannot open file '//trim(file_in)//' for reading')
read(iounit,size,iostat=ierr)
if(ierr.ne.0)call crash('cannot read namelist &size in file '//trim(file_in))
allocate(propx(1,n),propy(1,n),vx(1,n),vy(1,n),dzdxf(1,n),dzdyf(1,n),fmc_g(1,n),nfuel_cat(1,n),t2f(1,n),psfcf(1,n),ffw(1,n),ros(1,n))
read(iounit,data,iostat=ierr)
if(ierr.ne.0)call crash('cannot read namelist &data in file '//trim(file_in))
close(iounit,iostat=ierr)
if(ierr.ne.0)call crash('cannot close file '//trim(file_in))
print *,'propx=',propx
print *,'propy=',propy
print *,'vx=',vx
print *,'vy=',vy
print *,'dzdxf=',dzdxf
print *,'dzdyf=',dzdyf
print *,'fmc_g=',fmc_g
print *,'nfuel_cat=',nfuel_cat
print *,'t2f=',t2f
print *,'psfcf=',psfcf
print *,'ffw=',ffw
fp%vx => vx
fp%vy => vy
fp%dzdxf => dzdxf
fp%dzdyf => dzdyf
fp%fmc_g => fmc_g
fp%t2f => t2f
fp%psfcf => psfcf
fp%ffw => ffw
ros=0.

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
                           nfuel_cat,fp%fuel_time, &
                           fp )
i=1
do j=1,n
         call fire_ros(ros(i,j), &
                 propx(i,j),propy(i,j),i,j,fp,ierr,msg)
         if(ierr.ne.0)call crash(msg)
enddo

OPEN(iounit, FILE=trim(file_in),FORM='FORMATTED',STATUS='UNKNOWN',ACTION='WRITE',IOSTAT=ierr)
if(ierr.ne.0)call crash('cannot open file '//trim(file_out)//' for writing')
write(iounit,size,iostat=ierr)
if(ierr.ne.0)call crash('cannot write namelist &size in file '//trim(file_out))
write(iounit,data,iostat=ierr)
if(ierr.ne.0)call crash('cannot write namelist &data in file '//trim(file_out))
close(iounit,iostat=ierr)
if(ierr.ne.0)call crash('cannot close file '//trim(file_out))


end program fire_ros_main
