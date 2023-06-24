!*** A simple output module for wrf-fire standalone that produces output 
!    files and reads input files that are somewhat compatible with wrf files
!

module util_netcdf
use netcdf
implicit none

end type

contains

subroutine create_dim_attr(filename,grid)

! Create an empty NetCDF file with proper dimension names 
! (a.k.a. attributes) defined.

implicit none

!*** arguments
character(len=*),intent(in)::filename  ! filename to create
type(grid_info),intent(in)::grid       ! grid information structure

!*** local
integer::ncid,dimid,timeid,strlenid,varid

!*** executable

! create an empty file
call check(nf90_create(filename,nf90_clobber,ncid))

! define all dimensions

!!The function NF90_DEF_DIM adds a new dimension to an open netCDF dataset in 
!!define mode. It returns (as an argument) a dimension ID, given the netCDF ID,
!!the dimension name, and the dimension length. (From NetCDF docs)

call check(nf90_def_dim(ncid,dim_time,nf90_unlimited,dimid))
if(compat_fire_grid)then
  call check(nf90_def_dim(ncid,dim_fire_x,(grid%natmx+1)*grid%sr_x,dimid))
  call check(nf90_def_dim(ncid,dim_fire_y,(grid%natmy+1)*grid%sr_y,dimid))
else
  call check(nf90_def_dim(ncid,dim_fire_x,grid%nfirex,dimid))
  call check(nf90_def_dim(ncid,dim_fire_y,grid%nfirey,dimid))
endif
call check(nf90_def_dim(ncid,dim_atm_x,grid%natmx,dimid))
call check(nf90_def_dim(ncid,dim_atm_y,grid%natmy,dimid))
call check(nf90_def_dim(ncid,dim_atm_z,grid%natmz,dimid))
call check(nf90_def_dim(ncid,dim_atm_x_s,grid%natmx+1,dimid))
call check(nf90_def_dim(ncid,dim_atm_y_s,grid%natmy+1,dimid))
call check(nf90_def_dim(ncid,dim_atm_z_s,grid%natmz+1,dimid))
call check(nf90_def_dim(ncid,'DateStrLen',12,dimid))

! define global attributes
call check(nf90_put_att(ncid,nf90_global,'DX',grid%dx))
call check(nf90_put_att(ncid,nf90_global,'DY',grid%dy))
call check(nf90_put_att(ncid,nf90_global,'DT',grid%dt))
call check(nf90_put_att(ncid,nf90_global,'STANDALONE_DRIVER',1))

call check(nf90_enddef(ncid))
call check(nf90_close(ncid))

end subroutine create_dim_attr

!
!***
!

subroutine create_output_file(filename,grid)

!*** Create a file containing meta-data suitable for output of this simulation.
!    Only creates dimensions, variables, and attributes... does not populate
!    data.  Squash file if it already exists.

implicit none

!*** arguments
character(len=*),intent(in)::filename  ! filename to create
type(grid_info),intent(in)::grid       ! grid information structure

!*** local
integer::ncid,dimid,timeid,strlenid,varid

!*** executable

call create_dim_attr(filename,grid)

call check(nf90_open(filename,nf90_write,ncid))
call check(nf90_inq_dimid(ncid,'DateStrLen',strlenid))
call check(nf90_inq_dimid(ncid,dim_time,timeid))

! define a timekeeping variable
call check(nf90_redef(ncid))
call check(nf90_def_var(ncid,'Times',nf90_char,(/strlenid,timeid/),varid))
call check(nf90_enddef(ncid))
call check(nf90_close(ncid))

! create all of the output variables
call define_var(filename,var_lfn,dim_fire,unit_lfn,desc_lfn)
call define_var(filename,var_tign,dim_fire,unit_tign,desc_tign)
call define_var(filename,var_vx,dim_fire,unit_vx,desc_vx)
call define_var(filename,var_vy,dim_fire,unit_vy,desc_vy)
call define_var(filename,var_grnhfx,dim_fire,unit_grnhfx,desc_grnhfx)

end subroutine create_output_file

!
!***
!

subroutine define_var(filename,varname,dims,units,description)

!*** define a variable in a netcdf data set, the file is assumed to exist and
!    have valid meta-data (as created by create_output_file)

implicit none

!*** arguments
character(len=*),intent(in)::filename,varname  ! create variable varname in filename
character(len=*),dimension(:),intent(in)::dims ! the dimension names of the variable
character(len=*),intent(in)::units,description ! attributes created by wrf (not used at the moment)

!*** local
integer::ncid,i,varid
integer,dimension(4)::dimids
character(len=*),parameter::memorder='XYZ'
character(len=3)::stag

!*** executable
! open the file
call check(nf90_open(filename,nf90_write,ncid))

! get dimension id's
do i=1,size(dims)
  call check(nf90_inq_dimid(ncid,dims(i),dimids(i)))
enddo

! enter define mode and define the variable
call check(nf90_redef(ncid))
call check(nf90_def_var(ncid,varname,vartype,dimids(1:size(dims)),varid))

! add attributes
call check(nf90_put_att(ncid,varid,'FieldType',field_type))
call check(nf90_put_att(ncid,varid,'MemoryOrder',memorder(1:size(dims))))
call check(nf90_put_att(ncid,varid,'description',description))
call check(nf90_put_att(ncid,varid,'units',units))
if(size(dims).eq.3)then
  stag='Z'
elseif(trim(dims(1)).eq.dim_atm_x_s)then
  stag='X'
elseif(trim(dims(2)).eq.dim_atm_y_s)then
  stag='Y'
elseif(trim(dims(3)).eq.dim_atm_z_s)then
  stag='Z'
else
  stag=''
endif
call check(nf90_put_att(ncid,varid,'stagger',stag))
call check(nf90_put_att(ncid,varid,'coordinates','XLONG XLAT'))

! close file
call check(nf90_enddef(ncid))
call check(nf90_close(ncid))
end subroutine define_var

!
!***
!

subroutine write_inputs(filename,input,grid)
implicit none

! write simulation input file in a way compatible 
! with the wrf input file

!*** arguments
character(len=*),intent(in)::filename
type(input_vars),intent(in)::input
type(grid_info),intent(in)::grid

!*** local
integer::ncid,varid
integer,dimension(4)::s,c

!*** executable
call create_dim_attr(filename,grid)
call define_var(filename,var_nfuel_cat,dim_fire,'','')
call define_var(filename,var_ux,dim_atm_u,'','')
call define_var(filename,var_uy,dim_atm_v,'','')
call define_var(filename,var_dzdxf,dim_fire,'','')
call define_var(filename,var_dzdyf,dim_fire,'','')
call define_var(filename,var_zsf,dim_fire,'','')

call check(nf90_open(filename,nf90_write,ncid))

s(:)=1
c(:)=1
c(1)=grid%nfirex
c(2)=grid%nfirey

call check(nf90_inq_varid(ncid,var_nfuel_cat,varid))
call check(nf90_put_var(ncid,varid,input%nfuel_cat,start=s(1:3),count=c(1:3)))
call check(nf90_inq_varid(ncid,var_dzdxf,varid))
call check(nf90_put_var(ncid,varid,input%dzdxf,start=s(1:3),count=c(1:3)))
call check(nf90_inq_varid(ncid,var_dzdyf,varid))
call check(nf90_put_var(ncid,varid,input%dzdyf,start=s(1:3),count=c(1:3)))
call check(nf90_inq_varid(ncid,var_zsf,varid))
call check(nf90_put_var(ncid,varid,input%zsf,start=s(1:3),count=c(1:3)))

c(1)=grid%natmx+1
c(2)=grid%natmy

call check(nf90_inq_varid(ncid,var_ux,varid))
call check(nf90_put_var(ncid,varid,input%vx,start=s,count=c))

c(1)=grid%natmx
c(2)=grid%natmy+1

call check(nf90_inq_varid(ncid,var_uy,varid))
call check(nf90_put_var(ncid,varid,input%vy,start=s,count=c))

call check(nf90_close(ncid))

end subroutine write_inputs

!
!***
!

subroutine write_vars(filename,output,time)
!*** append variables to an output file (extending by the unlimited time
!    dimension)
implicit none

!*** arguments
character(len=*),intent(in)::filename
type(output_vars),intent(in)::output
real,intent(in)::time

!*** local
integer::ncid,tstep,dimid,varid
integer,dimension(3)::s,c
character(len=12)::timestr

!*** executable

! open the file
call check(nf90_open(filename,nf90_write,ncid))

! get the current number of timeslices already written
call check(nf90_inq_dimid(ncid,dim_time,dimid))
call check(nf90_inquire_dimension(ncid,dimid,len=tstep))

! write out the current simulation time
call check(nf90_inq_varid(ncid,'Times',varid))
write(timestr,'(E12.3)')time
call check(nf90_put_var(ncid,varid,timestr,start=(/1,tstep+1/),count=(/12,1/)))

! set up start and count variables for nf90_put_var so that it writes to the
! correct time slice and write data for each variable
s(:)=1
c(:)=0
s(3)=tstep+1
c(3)=1
c(1:2)=(/size(output%lfn,1),size(output%lfn,2)/)
call check(nf90_inq_varid(ncid,var_lfn,varid))
call check(nf90_put_var(ncid,varid,output%lfn,start=s,count=c))
c(1:2)=(/size(output%tign,1),size(output%tign,2)/)
call check(nf90_inq_varid(ncid,var_tign,varid))
call check(nf90_put_var(ncid,varid,output%tign,start=s,count=c))
c(1:2)=(/size(output%vx,1),size(output%vx,2)/)
call check(nf90_inq_varid(ncid,var_vx,varid))
call check(nf90_put_var(ncid,varid,output%vx,start=s,count=c))
c(1:2)=(/size(output%vy,1),size(output%vy,2)/)
call check(nf90_inq_varid(ncid,var_vy,varid))
call check(nf90_put_var(ncid,varid,output%vy,start=s,count=c))
c(1:2)=(/size(output%grnhfx,1),size(output%grnhfx,2)/)
call check(nf90_inq_varid(ncid,var_grnhfx,varid))
call check(nf90_put_var(ncid,varid,output%grnhfx,start=s,count=c))

! close file
call check(nf90_close(ncid))
end subroutine write_vars

!
!***
!

subroutine set_grid_from_file(filename,grid)
! get grid sizes from input file
implicit none

!*** arguments
character(len=*),intent(in)::filename
type(grid_info),intent(out)::grid

!*** local
integer::ncid,it
integer,dimension(4)::fdimid,adimid

!*** executable
call check(nf90_open(filename,nf90_nowrite,ncid))
call check(nf90_inq_dimid(ncid,dim_time,fdimid(3)))
adimid(4)=fdimid(3)
call check(nf90_inq_dimid(ncid,dim_fire_x,fdimid(1)))
call check(nf90_inq_dimid(ncid,dim_fire_y,fdimid(2)))
call check(nf90_inq_dimid(ncid,dim_atm_x, adimid(1)))
call check(nf90_inq_dimid(ncid,dim_atm_y, adimid(2)))
call check(nf90_inq_dimid(ncid,dim_atm_z, adimid(3)))

call check(nf90_inquire_dimension(ncid,fdimid(3),len=it))
if(it.ne.1)then
  call crash('invalid number of time steps in input file, must be 1')
endif
call check(nf90_inquire_dimension(ncid,fdimid(2),len=grid%nfirey))
call check(nf90_inquire_dimension(ncid,fdimid(1),len=grid%nfirex))
call check(nf90_inquire_dimension(ncid,adimid(2),len=grid%natmy))
call check(nf90_inquire_dimension(ncid,adimid(1),len=grid%natmx))
call check(nf90_inquire_dimension(ncid,adimid(3),len=grid%natmz))

call check(nf90_get_att(ncid,nf90_global,'DX',grid%dx))
call check(nf90_get_att(ncid,nf90_global,'DY',grid%dy))
call check(nf90_get_att(ncid,nf90_global,'DT',grid%dt))

if(compat_fire_grid)then
  grid%sr_x=grid%nfirex/(grid%natmx+1)
  grid%sr_y=grid%nfirey/(grid%natmy+1)
  if( (grid%natmx+1)*grid%sr_x .ne. grid%nfirex .or. &
      (grid%natmy+1)*grid%sr_y .ne. grid%nfirey)then
      call crash('invalid dimensions in input file or compat_fire_grid set incorrectly')
  endif
  grid%nfirex=grid%natmx*grid%sr_x
  grid%nfirey=grid%natmy*grid%sr_y
else
  grid%sr_x=grid%nfirex/grid%natmx
  grid%sr_y=grid%nfirey/grid%natmy
  if( grid%natmx*grid%sr_x .ne. grid%nfirex .or. &
      grid%natmy*grid%sr_y .ne. grid%nfirey)then
      call crash('invalid dimensions in input file or compat_fire_grid set incorrectly')
  endif
endif
call check(nf90_close(ncid))
end subroutine set_grid_from_file

!
!***
!

subroutine read_vars(filename,input,grid)
! read all variables from input file
implicit none

!*** arguments
character(len=*),intent(in)::filename
type(input_vars),intent(inout)::input
type(grid_info),intent(in)::grid

!*** local
integer::ncid,varid,ierr
integer,dimension(4)::s,c
real,dimension(:,:),allocatable::u,v

!*** executable
call check(nf90_open(filename,nf90_nowrite,ncid))
s(:)=1
c(:)=1
c(1)=grid%nfirex
c(2)=grid%nfirey
call check(nf90_inq_varid(ncid,var_nfuel_cat,varid))
call check(nf90_get_var(ncid,varid,input%nfuel_cat,start=s(1:3),count=c(1:3)))
call check(nf90_inq_varid(ncid,var_dzdxf,varid))
call check(nf90_get_var(ncid,varid,input%dzdxf,start=s(1:3),count=c(1:3)))
call check(nf90_inq_varid(ncid,var_dzdyf,varid))
call check(nf90_get_var(ncid,varid,input%dzdyf,start=s(1:3),count=c(1:3)))

allocate(u(0:grid%natmx+2,0:grid%natmy+2),v(0:grid%natmx+2,0:grid%natmy+2),stat=ierr)
if(ierr.ne.0)call crash('memory allocation error')

c(1)=grid%natmx+1
c(2)=grid%natmy
c(3)=1
c(4)=1

call check(nf90_inq_varid(ncid,var_ux,varid))
call check(nf90_get_var(ncid,varid,u(1:grid%natmx+1,1:grid%natmy),start=s,count=c))

c(1)=grid%natmx
c(2)=grid%natmy+1
c(3)=1
c(4)=1

call check(nf90_inq_varid(ncid,var_uy,varid))
call check(nf90_get_var(ncid,varid,v(1:grid%natmx,1:grid%natmy+1),start=s,count=c))

  call continue_at_boundary(1,0,0.,       &
          0,grid%natmx+2,0,grid%natmy+2,  &
          1,grid%natmx,1,grid%natmy+1,    &
          1,grid%natmx,1,grid%natmy+1,    &
          1,grid%natmx,1,grid%natmy+1,    &
          v)

  call continue_at_boundary(0,1,0.,       &
          0,grid%natmx+2,0,grid%natmy+2,  &
          1,grid%natmx+1,1,grid%natmy,    &
          1,grid%natmx+1,1,grid%natmy,    &
          1,grid%natmx+1,1,grid%natmy,    &
          u)


  call interpolate_2d(0,grid%natmx+2,0,grid%natmy+2,         &
                      1,grid%natmx+1,0,grid%natmy+1,         &
                      1,grid%nfirex,1,grid%nfirey,           &
                      1,grid%nfirex,1,grid%nfirey,           &
                      grid%sr_x,grid%sr_y,                   &
                      1.,1.,.5,1-.5+(grid%sr_y)*.5,          &
                      u,input%vx)

  call interpolate_2d(0,grid%natmx+2,0,grid%natmy+2,         &
                      0,grid%natmx+1,1,grid%natmy+1,         &
                      1,grid%nfirex,1,grid%nfirey,           &
                      1,grid%nfirex,1,grid%nfirey,           &
                      grid%sr_x,grid%sr_y,                   &
                      1.,1.,1-.5+(grid%sr_x)*.5,.5,          &
                      v,input%vy)

deallocate(u,v,stat=ierr)
if(ierr.ne.0)call crash('deallocation error')
call check(nf90_close(ncid))
end subroutine read_vars

subroutine check(ncerr)
implicit none
integer,intent(in)::ncerr
if(ncerr.ne.nf90_noerr)then
  print*,"Error calling NetCDF subroutine"
  print*,trim(nf90_strerror(ncerr))
  call crash("NETCDF ERROR")
endif
end subroutine 

end module wrf_netcdf
