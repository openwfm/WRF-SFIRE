! A simple program to initialize a test case for wrf-fire standalone driver.  
! It will create a netcdf file called fire_input.nc containing initial
! conditions for the fire simulation, somewhat compatible with wrfinput_d01
! for WRF, created by ideal.exe or real.exe

program model_test_init

use module_configure
use wrf_netcdf
use module_fr_sfire_util, only : crash
implicit none

type(domain)::grid
type(grid_info)::info
integer::ierr,i,j,iu
!real,parameter::wind=1.,alpha=0.0
real::windx,windy,slopex,slopey
type(grid_config_rec_type)::config_flags ! for the namelist

namelist /uniform_init/ windx,windy,slopex,slopey

iu=97
windx=10
windy=0
slopex=0
slopey=0

print*,'Reading ',trim(nmlfile)
open(iu,file=nmlfile,err=99)
read(iu,uniform_init,err=99,end=99)
close(iu,err=99)
goto 100

99 continue
print*,'Error in reading namelist uniform_init, using defaults'
100 continue

! read namelist parameters
call read_namelist(config_flags)

! get grid information
info%nstagx=config_flags%e_we-config_flags%s_we+1
info%nstagy=config_flags%e_sn-config_flags%s_sn+1
!info%natmz=config_flags%e_vert-config_flags%s_vert+1
info%sr_x=config_flags%sr_x
info%sr_y=config_flags%sr_y
info%nfirex=(info%nstagx)*config_flags%sr_x
info%nfirey=(info%nstagy)*config_flags%sr_y
info%io_nfirex=info%nfirex
info%io_nfirey=info%nfirey
info%len_time_string=19
info%fdx=config_flags%dx / info%sr_x
info%fdy=config_flags%dy / info%sr_y
info%dt=real(config_flags%time_step)+real(config_flags%time_step_fract_num)/real(config_flags%time_step_fract_den)
info%times='1980-01-01_00:00:00'

! allocate space for input variables
allocate(grid%nfuel_cat(info%nfirex,info%nfirey))
allocate(grid%dzdxf(info%nfirex,info%nfirey)) 
allocate(grid%dzdyf(info%nfirex,info%nfirey))
allocate(grid%zsf(info%nfirex,info%nfirey))
allocate(grid%uf(info%nfirex,info%nfirey))
allocate(grid%vf(info%nfirex,info%nfirey))
allocate(grid%fxlat(info%nfirex,info%nfirey))
allocate(grid%fxlong(info%nfirex,info%nfirey))
allocate(grid%fmc_g(info%nfirex,info%nfirey))

! initialize input variables
do j=1,info%nfirey
  do i=1,info%nfirex
    grid%fxlat(i,j)=(j-1)*info%fdy
    grid%fxlong(i,j)=(i-1)*info%fdx
    grid%zsf(i,j)=grid%fxlong(i,j)*slopex + grid%fxlat(i,j)*slopey
    grid%dzdxf(i,j)=slopex
    grid%dzdyf(i,j)=slopey
    grid%nfuel_cat(i,j)=3
    grid%fmc_g(i,j)=.08
  enddo
enddo

grid%uf(:,:)=windx
grid%vf(:,:)=windy

! write variables to the file
call write_inputs(inputfile,grid,info)

! clean up and exit
deallocate(grid%nfuel_cat,grid%dzdxf,grid%dzdyf,grid%uf,grid%vf,grid%zsf,&
           grid%fxlat,grid%fxlong,grid%fmc_g)



end program model_test_init
