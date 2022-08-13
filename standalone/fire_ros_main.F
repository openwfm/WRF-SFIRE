program fire_ros_main

use module_fr_sfire_phys
use module_fr_sfire_util

!*** variables

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
logical:: have_data=.true.
character(len=50)::file_in="namelist_standalone.input",file_out="namelist_standalone.output"
character(len=128)::msg

!*** executable

fire_print_msg = -1  ! suppress messages and warnings from WRF-SFIRE

! read namelist.fire into fuels
! also create file fuels.m read and visualize the results as described in 
! https://wiki.openwfm.org/wiki/How_to_diagnose_fuel_properties_in_WRF-SFIRE

print *,'Initializing fuel categories from namelist.fire, writing fuels.m'
call init_fuel_cats(init_fuel_moisture=.false.) 

namelist /size/n
namelist /data/propx,propy,vx,vy,dzdxf,dzdyf,fmc_g,nfuel_cat,t2f,psfcf,ffw,ros

print *,'Reading file '//trim(file_in)
OPEN(iounit, FILE=trim(file_in),FORM='FORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=ierr)
if(ierr.ne.0)then 
     print *,'cannot open file '//trim(file_in)//' for reading'
     n=1
     have_data=.false.
     goto 100
endif
read(iounit,size,iostat=ierr)
if(ierr.ne.0)then
     n=1
     print *,'cannot read namelist &size in file '//trim(file_in)
     have_data=.false.
     goto 100
endif

100 continue
print *,'Using size n=',n

! allocate our arrays used in namelist
allocate(propx(1,n),propy(1,n),vx(1,n),vy(1,n),dzdxf(1,n),dzdyf(1,n),fmc_g(1,n),nfuel_cat(1,n),t2f(1,n),psfcf(1,n),ffw(1,n),ros(1,n))

! allocate the rest of the arrays in fp 
allocate(fp%fgip(1,n),fp%ischap(1,n),fp%phisc(1,n),fp%bbb(1,n),fp%fuel_time(1,n),fp%phiwc(1,n),fp%r_0(1,n))

! initialize our arrays to reasonable values that won't crash the code
propx = 0.
propy = 0.
vx = 0.
vy = 0.
dzdxf = 0.
dzdyf = 0.
fmc_g = 0.
nfuel_cat=14  ! no fuel
t2f = 273.    ! temp K 
psfcf = 1e5   ! pressure Pa
ffw = 10.     ! fire front width m
ros = 0.

if(have_data)then
     read(iounit,data,iostat=ierr)
     if(ierr.ne.0)then
          print *,'cannot read namelist &data in file '//trim(file_in)
          have_data=.false.
     endif
endif

close(iounit,iostat=ierr)

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
fp%nfuel_cat => nfuel_cat


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

print *,'writing file '//trim(file_out)
OPEN(iounit, FILE=trim(file_out),FORM='FORMATTED',STATUS='UNKNOWN',ACTION='WRITE',IOSTAT=ierr)
if(ierr.ne.0)call crash('cannot open file '//trim(file_out)//' for writing')
write(iounit,size,iostat=ierr)
if(have_data)then
     write(*,*)       'ros computed from inputs  file '//trim(file_in)
     write(iounit,*)'! ros computed from inputs  file '//trim(file_in)
else
     write(*,*)       'ros computed from defaults, use as a template ONLY:'
     write(iounit,*)'! ros computed from defaults, use as a template ONLY:'
     write(*,*)       'cp '//trim(file_out)//' '//trim(file_in)//' and edit'
     write(iounit,*)'! cp '//trim(file_out)//' '//trim(file_in)//' and edit'
endif
write(iounit,data,iostat=ierr)
close(iounit,iostat=ierr)

end program fire_ros_main
