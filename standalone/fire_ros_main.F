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
real,pointer,dimension(:,:):: tempf                ! temperature at 2m (K) (Balbi model only)
real,pointer,dimension(:,:):: rhof                 ! surface air density (kg/m^3) (Balbi model only)
real,pointer,dimension(:,:):: ffwidth              ! fire front width (m) (Balbi model only)
real,pointer,dimension(:,:):: ros                  ! fire rate of spread (m/s) 

integer:: n=10, iounit=9, ierr, i, j
logical:: have_data=.true.
character(len=50)::file_in="namelist_standalone.input",file_out="namelist_standalone.output"
character(len=128)::msg
character(len=9), dimension(2)::model=(/'Rothermel','Balbi    '/)

namelist /control/n,fire_print_msg,balbi_msglevel,call_write_fuels_m,compare_rothermel
namelist /data/propx,propy,vx,vy,dzdxf,dzdyf,fmc_g,nfuel_cat,tempf,rhof,ffwidth,ros

!*** executable

! defaults, can set in &control
fire_print_msg = -1  ! suppress messages and warnings from WRF-SFIRE
balbi_msglevel = 1           ! prints in fire_ros_balbi
call_write_fuels_m = .false. ! create fuels.m for graphics 

! read namelist.fire into fuels
! also create file fuels.m read and visualize the results as described in 
! https://wiki.openwfm.org/wiki/How_to_diagnose_fuel_properties_in_WRF-SFIRE

print *,'Reading file '//trim(file_in)
OPEN(iounit, FILE=trim(file_in),FORM='FORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=ierr)
if(ierr.ne.0)then 
     print *,'Cannot open file '//trim(file_in)//' for reading,using defaults'
     n=1
     have_data=.false.
     goto 100
endif
read(iounit,control,iostat=ierr)
if(ierr.ne.0)then
     n=1
     print *,'Cannot read namelist &control in file '//trim(file_in)//', using defaults'
     have_data=.false.
     goto 100
endif

100 continue
print *,'fire_ros using control parameters:'
print *,'number of cases in namelist &data                n=',n
print *,'WRF-SFIRE print level               fire_print_msg=',fire_print_msg
print *,'print level in fire_ros_balbi       balbi_msglevel=',balbi_msglevel
print *,'write file fuels.m for graphics call_write_fuels_m=',call_write_fuels_m
print *,'compare codes (ibeh=1 only)      compare_rothermel=',compare_rothermel
write(*,control)

print *,'Initializing fuel categories from file namelist.fire'
if(call_write_fuels_m)then
     print *,'and writing file fuels.m. For a guide how to use fuels.m, see' 
     print *,'https://wiki.openwfm.org/wiki/How_to_diagnose_fuel_properties_in_WRF-SFIRE'
endif
call init_fuel_cats(init_fuel_moisture=.false.) 
print *,'Fire ROS model from file namelist.fire is ibeh=',ibeh,' ',model(ibeh)
 
! allocate our arrays used in namelist
allocate(propx(1,n),propy(1,n),vx(1,n),vy(1,n),dzdxf(1,n),dzdyf(1,n),fmc_g(1,n),nfuel_cat(1,n),tempf(1,n),rhof(1,n),ffwidth(1,n),ros(1,n))

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
nfuel_cat=1   ! some fuel
tempf = 273.    ! temp K 
rhof = 1.2    ! air density (kg/m^3)
ffwidth = 10. ! fire front width m
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
print *,'tempf=',tempf
print *,'rhof=',rhof
print *,'ffwidth=',ffwidth
fp%vx => vx
fp%vy => vy
fp%dzdxf => dzdxf
fp%dzdyf => dzdyf
fp%fmc_g => fmc_g
fp%tempf => tempf
fp%rhof => rhof
fp%ffwidth => ffwidth
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

print *,'computed ros=',ros

print *,'writing file '//trim(file_out)
OPEN(iounit, FILE=trim(file_out),FORM='FORMATTED',STATUS='UNKNOWN',ACTION='WRITE',IOSTAT=ierr)
if(ierr.ne.0)call crash('cannot open file '//trim(file_out)//' for writing')
write(iounit,control,iostat=ierr)
if(have_data)then
     write(*,*)       'ros computed from inputs in file '//trim(file_in)
     write(iounit,*)'! ros computed from inputs in file '//trim(file_in)
else
     write(*,*)       'ros computed from defaults, use as a template ONLY:'
     write(iounit,*)'! ros computed from defaults, use as a template ONLY:'
     write(*,*)       'cp '//trim(file_out)//' '//trim(file_in)//' and edit'
     write(iounit,*)'! cp '//trim(file_out)//' '//trim(file_in)//' and edit'
endif
write(iounit,data,iostat=ierr)
close(iounit,iostat=ierr)

end program fire_ros_main
