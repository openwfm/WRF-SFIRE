!    files and reads input files that are somewhat compatible with wrf files
!

module wrf_netcdf
use netcdf
use module_fr_sfire_util , only : crash, interpolate_2d, continue_at_boundary
use module_domain, only: domain
implicit none

! default file names
character(len=*),parameter::inputfile='fire_input.nc'
character(len=*),parameter::outputfile='fire_output.nc'

! control whether the dimensions of the fire grid variables should have 
! the same (incorrect) sizes that wrf outputs and same dimension names 
! 
logical::compat_fire_grid=.true.
logical::debug_print=.false.
logical::debug_print_f=.true.
logical::read_check=.false.

! output variable type
integer,parameter::vartype=nf90_float,field_type=104

! max number of dimensions
integer, parameter:: mdims=4

! variable names
character(len=nf90_max_name)::                        &
                              var_xtime='XTIME',        &
                              unit_xtime='min',                &
                              desc_xtime='minutes since simulation start',                &
                              var_itimestep='ITIMESTEP',        &
                              unit_itimestep='',                &
                              desc_itimestep='',                &
                              var_nfuel_cat='NFUEL_CAT',        &
                              unit_nfuel_cat='',                &
                              desc_nfuel_cat='',                &
                              var_dzdxf='DZDXF',                &
                              unit_dzdxf='',                    &
                              desc_dzdxf='',                    &
                              var_dzdyf='DZDYF',                &
                              unit_dzdyf='',                    &
                              desc_dzdyf='',                    &
                              var_zsf='ZSF',                    &
                              unit_zsf='',                      &
                              desc_zsf='',                      &
                              var_lfn='LFN',                    &
                              unit_lfn='',                      &
                              desc_lfn='',                      &
                              var_tign_g='TIGN_G',              &
                              unit_tign_g='',                   &
                              desc_tign_g='',                   &
                              var_tign_in='TIGN_IN',            &
                              unit_tign_in='',                  &
                              desc_tign_in='',                  &
                              var_fmc_g='FMC_G',                &
                              unit_fmc_g='',                    &
                              desc_fmc_g='',                    &
                              var_fmc_gc='FMC_GC',              &
                              unit_fmc_gc='',                   &
                              desc_fmc_gc='',                   &
                              var_fxlong='FXLONG',              &
                              unit_fxlong='',                   &
                              desc_fxlong='',                   &
                              var_fxlat='FXLAT',                &
                              unit_fxlat='',                    &
                              desc_fxlat='',                    &
                              var_unit_fxlong='UNIT_FXLONG',    &
                              unit_unit_fxlong='',              &
                              desc_unit_fxlong='',              &
                              var_unit_fxlat='UNIT_FXLAT',      &
                              unit_unit_fxlat='',               &
                              desc_unit_fxlat='',               &
                              var_uf='UF',                      &
                              unit_uf='',                       &
                              desc_uf='',                       &
                              var_vf='VF',                      &
                              unit_vf='',                       &
                              desc_vf='',                       &
                              var_fuel_frac='FUEL_FRAC',        &
                              unit_fuel_frac='',                &
                              desc_fuel_frac='',                &
                              var_fire_area='FIRE_AREA',        &
                              unit_fire_area='',                &
                              desc_fire_area='',                &
                              var_fgrnhfx='FGRNHFX',            &
                              unit_fgrnhfx='',                  &
                              desc_fgrnhfx='',                  &
                              var_fgrnqfx='FGRNQFX',            &
                              unit_fgrnqfx='J/m^2',             &
                              desc_fgrnqfx='heat flux',         &
                              var_ros='ROS',                    &
                              unit_ros='m/s',                   &
                              desc_ros='ROS',                   &
                              var_r_0='R_0',                    &
                              unit_r_0='m/s',                   &
                              desc_r_0='R_0',                   &
                              var_flineint='FLINEINT',          &
                              unit_flineint='J/m/s',            &
                              desc_flineint='Byram fireline intensity',         &
                              var_flineint2='FLINEINT2',        & 
                              unit_flineint2='J/m/s^2',         & 
                              desc_flineint2='New fireline intensity' , &
                              var_f_ros0='F_ROS0',          &
                              unit_f_ros0='m/s',            &
                              desc_f_ros0='base rate of spread in all directions',  &
                              var_f_rosx='F_ROSX',          &
                              unit_f_rosx='m/s',            &
                              desc_f_rosx='potential rate of spread in direction X',  &
                              var_f_rosy='F_ROSY',          &
                              unit_f_rosy='m/s',            &
                              desc_f_rosy='potential rate of spread in direction Y',  &
                              var_f_ros='F_ROS',          &
                              unit_f_ros='m/s',            &
                              desc_f_ros='potential fire max spread rate in any direction',  &
                              var_f_int='F_INT',          &
                              unit_f_int='J/m^2/s',            &
                              desc_f_int='potential fire reaction intensity for risk rating',  &
                              var_f_lineint='F_LINEINT',          &
                              unit_f_lineint='J/m/s',            &
                              desc_f_lineint='potential Byram fireline intensity for risk rating',  &
                              var_f_lineint2='F_LINEINT2',          &
                              unit_f_lineint2='J/m/s^2',            &
                              desc_f_lineint2='potential alternative fireline intensity for risk rating'

character(len=nf90_max_name),parameter::var_times='Times' 

! grid information structure
integer, parameter::max_times_length=19
type grid_info
  integer::nfirex,nfirey,      & ! fire grid size
           ntimes                ! number of time frames
  real::fdx,fdy,dt               ! mesh and time resolution
  character(len=max_times_length)::times ! the time string from the file
  integer::len_time_string     ! length of the times string
  character(len=nf90_max_name)::      &
    dim_time_string='DateStrLen',     &
    dim_fire_x='west_east_subgrid',   &
    dim_fire_y='south_north_subgrid', &
    dim_time='Time'
  integer::io_nfirex,io_nfirey        ! fire grid size in files
  ! WRF compatibility
  integer :: sr_x,sr_y           ! refinement ratios
  integer :: nstagx,nstagy
  character(len=nf90_max_name)::      &
    dim_atm_x='west_east',            &
    dim_atm_y='south_north',          &
    dim_atm_z='bottom_top',           &
    dim_atm_x_s='west_east_stag',     &
    dim_atm_y_s='south_north_stag',   &
    dim_atm_z_s='bottom_top_stag'
end type

contains

subroutine create_dim_attr(filename,info)

! Create an empty NetCDF file with proper dimension names 
! (a.k.a. attributes) defined.

implicit none

!*** arguments
character(len=*),intent(in)::filename  ! filename to create
type(grid_info),intent(in)::info       ! grid information structure

!*** local
integer::ncid,dimid,timeid,strlenid,varid

!*** executable

! create an empty file
call check(nf90_create(filename,nf90_clobber,ncid))

! define all dimensions

call create_dim(ncid,info%dim_time,nf90_unlimited)
call create_dim(ncid,info%dim_fire_x,info%io_nfirex)
call create_dim(ncid,info%dim_fire_y,info%io_nfirey)
call create_dim(ncid,info%dim_time_string,info%len_time_string)
if (compat_fire_grid)then
  call create_dim(ncid,info%dim_atm_x_s,info%nstagx)
  call create_dim(ncid,info%dim_atm_y_s,info%nstagy)
endif

! define global attributes
call check(nf90_put_att(ncid,nf90_global,'DX',info%fdx*info%sr_x))
call check(nf90_put_att(ncid,nf90_global,'DY',info%fdy*info%sr_y))
call check(nf90_put_att(ncid,nf90_global,'FDX',info%fdx))
call check(nf90_put_att(ncid,nf90_global,'FDY',info%fdy))
call check(nf90_put_att(ncid,nf90_global,'DT',info%dt))
call check(nf90_put_att(ncid,nf90_global,'STANDALONE_DRIVER',1))

call check(nf90_enddef(ncid))
call check(nf90_close(ncid))

end subroutine create_dim_attr

subroutine create_dim(ncid,dimname,dimsize)
integer, intent(in):: ncid
character(len=*),intent(in)::dimname  ! dimension name to create
integer, intent(in):: dimsize
integer::dimid

!!The function NF90_DEF_DIM adds a new dimension to an open netCDF dataset in 
!!define mode. It returns (as an argument) a dimension ID, given the netCDF ID,
!!the dimension name, and the dimension length. (From NetCDF docs)

if(debug_print)write(*,'(3a,i10)')'creating dimension ',trim(dimname),' size',dimsize
call check(nf90_def_dim(ncid,dimname,dimsize,dimid))
end subroutine create_dim

!
!***
!

subroutine create_output_file(filename,info)

!*** Create a file containing meta-data suitable for output of this simulation.
!    Only creates dimensions, variables, and attributes... does not populate
!    data.  Squash file if it already exists.

implicit none

!*** arguments
character(len=*),intent(in)::filename  ! filename to create
type(grid_info),intent(in)::info       ! grid information structure

!*** local
integer::ncid,dimid,timeid,strlenid,varid
character(len=NF90_MAX_NAME), dimension(3)::dim_fire

!*** executable

call create_dim_attr(filename,info)

call ncopen(filename,nf90_write,ncid)
call check(nf90_inq_dimid(ncid,info%dim_time_string,strlenid))
call check(nf90_inq_dimid(ncid,info%dim_time,timeid))

call check(nf90_redef(ncid))
! define a timekeeping variable
call check(nf90_def_var(ncid,var_times,nf90_char,(/strlenid,timeid/),varid))
! define scalar variables
call check(nf90_def_var(ncid,var_unit_fxlong,nf90_float,(/timeid/),varid))
call check(nf90_def_var(ncid,var_unit_fxlat,nf90_float,(/timeid/),varid))
call check(nf90_def_var(ncid,var_xtime,nf90_float,(/timeid/),varid))
call check(nf90_def_var(ncid,var_itimestep,nf90_int,(/timeid/),varid))
call check(nf90_enddef(ncid))
call check(nf90_close(ncid))

! create all of the output variables
dim_fire=(/info%dim_fire_x,info%dim_fire_y,info%dim_time/)
call define_var(filename,info,var_zsf,3,dim_fire,unit_zsf,desc_zsf)
call define_var(filename,info,var_dzdxf,3,dim_fire,unit_dzdxf,desc_dzdxf)
call define_var(filename,info,var_dzdyf,3,dim_fire,unit_dzdyf,desc_dzdyf)
call define_var(filename,info,var_nfuel_cat,3,dim_fire,unit_nfuel_cat,desc_nfuel_cat)
call define_var(filename,info,var_uf,3,dim_fire,unit_uf,desc_uf)
call define_var(filename,info,var_vf,3,dim_fire,unit_vf,desc_vf)
call define_var(filename,info,var_fmc_g,3,dim_fire,unit_fmc_g,desc_fmc_g)
call define_var(filename,info,var_lfn,3,dim_fire,unit_lfn,desc_lfn)
call define_var(filename,info,var_tign_g,3,dim_fire,unit_tign_g,desc_tign_g)
call define_var(filename,info,var_tign_in,3,dim_fire,unit_tign_in,desc_tign_in)
call define_var(filename,info,var_fxlong,3,dim_fire,unit_fxlong,desc_fxlong)
call define_var(filename,info,var_fxlat,3,dim_fire,unit_fxlat,desc_fxlat)
call define_var(filename,info,var_fgrnhfx,3,dim_fire,unit_fgrnhfx,desc_fgrnhfx)
call define_var(filename,info,var_fgrnqfx,3,dim_fire,unit_fgrnqfx,desc_fgrnqfx)
call define_var(filename,info,var_fuel_frac,3,dim_fire,unit_fuel_frac,desc_fuel_frac)
call define_var(filename,info,var_fire_area,3,dim_fire,unit_fire_area,desc_fire_area)
call define_var(filename,info,var_flineint,3,dim_fire,unit_flineint,desc_flineint)
call define_var(filename,info,var_flineint2,3,dim_fire,unit_flineint2,desc_flineint2)
call define_var(filename,info,var_ros,3,dim_fire,unit_ros,desc_ros)
call define_var(filename,info,var_r_0,3,dim_fire,unit_r_0,desc_r_0)
call define_var(filename,info,var_f_ros0,3,dim_fire,unit_f_ros0,desc_f_ros0)
call define_var(filename,info,var_f_rosx,3,dim_fire,unit_f_rosx,desc_f_rosx)
call define_var(filename,info,var_f_rosy,3,dim_fire,unit_f_rosy,desc_f_rosy)
call define_var(filename,info,var_f_ros,3,dim_fire,unit_f_ros,desc_f_ros)
call define_var(filename,info,var_f_int,3,dim_fire,unit_f_int,desc_f_int)
call define_var(filename,info,var_f_lineint,3,dim_fire,unit_f_lineint,desc_f_lineint)
call define_var(filename,info,var_f_lineint2,3,dim_fire,unit_f_lineint2,desc_f_lineint2)

end subroutine create_output_file

!
!***
!

subroutine define_var(filename,info,varname,ndims,dims,units,description)

!*** define a variable in a netcdf data set, the file is assumed to exist and
!    have valid meta-data (as created by create_output_file)

implicit none

!*** arguments
character(len=*),intent(in)::filename,varname  ! create variable varname in filename
type(grid_info), intent(in)::info
integer, intent(in)::ndims
character(len=NF90_MAX_NAME),intent(in)::dims(ndims) ! the dimension names of the variable
character(len=*),intent(in) ::units,description ! attributes created by wrf 

!*** local
integer::ncid,i,varid
integer,dimension(ndims)::dimids
character(len=*),parameter::memorder='XYZ'
character(len=3)::stag

!*** executable
! open the file
call ncopen(filename,nf90_write,ncid)

! get dimension id's
do i=1,ndims
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
if(ndims.eq.3)then
  stag='Z'
elseif(trim(dims(1)).eq.info%dim_atm_x_s)then
  stag='X'
elseif(trim(dims(2)).eq.info%dim_atm_y_s)then
  stag='Y'
elseif(trim(dims(3)).eq.info%dim_atm_z_s)then
  stag='Z'
else
  stag=''
endif
call check(nf90_put_att(ncid,varid,'stagger',stag))
call check(nf90_put_att(ncid,varid,'coordinates','XLONG XLAT'))

! close file
call check(nf90_enddef(ncid))
call check(nf90_close(ncid))
call print_var_info(filename,varname)  ! check
end subroutine define_var

!
!***
!

subroutine write_vars(filename,grid,info,iframe)
!*** append variables to an output file (extending by the unlimited time
!    dimension)
implicit none

!*** arguments
character(len=*),intent(in)::filename
type(domain),intent(inout)::grid      ! the mother of all arrays
type(grid_info),intent(in)::info    ! dimensions
integer, intent(in)::iframe         ! frame to write into

!*** local
integer::ncid,dimid,varid

!*** executable

! open the file
call print_var_info(filename,var_lfn)
call ncopen(filename,nf90_write,ncid)

if(debug_print_f)write(*,'(3a,i4,2a)')'write_vars: writing file ',trim(filename), &
   ' frame ',iframe,' time ',info%times 

! write out the current simulation time
call check(nf90_inq_varid(ncid,var_times,varid),'cannot find '//trim(var_times))
call check(nf90_put_var(ncid,varid,info%times,start=(/1,iframe/), &
     count=(/info%len_time_string,1/)), 'error writing '//trim(var_times))

call write_integer(ncid,iframe,var_itimestep,grid%itimestep)
call write_real(ncid,iframe,var_xtime,grid%xtime)
call write_fire_var(ncid,info,iframe,var_lfn,grid%lfn)
call write_fire_var(ncid,info,iframe,var_tign_g,grid%tign_g)
call write_fire_var(ncid,info,iframe,var_fxlong,grid%fxlong)
call write_fire_var(ncid,info,iframe,var_fxlat,grid%fxlat)
call write_fire_var(ncid,info,iframe,var_zsf,grid%zsf)
call write_fire_var(ncid,info,iframe,var_fuel_frac,grid%fuel_frac)
call write_fire_var(ncid,info,iframe,var_fire_area,grid%fire_area)
call write_fire_var(ncid,info,iframe,var_fgrnhfx,grid%fgrnhfx)
call write_fire_var(ncid,info,iframe,var_fgrnqfx,grid%fgrnqfx)
call write_fire_var(ncid,info,iframe,var_ros,grid%ros)
call write_fire_var(ncid,info,iframe,var_r_0,grid%r_0)
call write_fire_var(ncid,info,iframe,var_flineint,grid%flineint)
call write_fire_var(ncid,info,iframe,var_flineint2,grid%flineint2)
call write_real(ncid,iframe,var_unit_fxlong,grid%unit_fxlong)
call write_real(ncid,iframe,var_unit_fxlat,grid%unit_fxlat)
call write_fire_var(ncid,info,iframe,var_uf,grid%uf)
call write_fire_var(ncid,info,iframe,var_vf,grid%vf)
call write_fire_var(ncid,info,iframe,var_f_ros0,grid%f_ros0)
call write_fire_var(ncid,info,iframe,var_f_rosx,grid%f_rosx)
call write_fire_var(ncid,info,iframe,var_f_rosy,grid%f_rosy)
call write_fire_var(ncid,info,iframe,var_f_ros,grid%f_ros)
call write_fire_var(ncid,info,iframe,var_f_int,grid%f_int)
call write_fire_var(ncid,info,iframe,var_f_lineint,grid%f_lineint)
call write_fire_var(ncid,info,iframe,var_f_lineint2,grid%f_lineint2)

! close file
call latch
call check(nf90_close(ncid))
call print_var_info(filename,var_lfn)
end subroutine write_vars

subroutine write_inputs(filename,grid,info)
implicit none
character(len=*),intent(in)::filename
type(domain),intent(in)::grid
type(grid_info),intent(in)::info

integer::ncid,dimid,varid

call create_output_file(filename,info)
call ncopen(filename,nf90_write,ncid)

! write out the current simulation time
call check(nf90_inq_varid(ncid,var_times,varid),'cannot find '//trim(var_times))
call check(nf90_put_var(ncid,varid,info%times,start=(/1,1/), &
     count=(/info%len_time_string,1/)), 'error writing '//trim(var_times))

call write_fire_var(ncid,info,1,var_nfuel_cat,grid%nfuel_cat)
call write_fire_var(ncid,info,1,var_dzdxf,grid%dzdxf)
call write_fire_var(ncid,info,1,var_dzdyf,grid%dzdyf)
call write_fire_var(ncid,info,1,var_zsf,grid%zsf)
call write_fire_var(ncid,info,1,var_uf,grid%uf)
call write_fire_var(ncid,info,1,var_vf,grid%vf)
call write_fire_var(ncid,info,1,var_fmc_g,grid%fmc_g)
call write_fire_var(ncid,info,1,var_fxlat,grid%fxlat)
call write_fire_var(ncid,info,1,var_fxlong,grid%fxlong)

call latch
call check(nf90_close(ncid))
end subroutine write_inputs

!
!*** READING 
!

subroutine read_info(filename,info)
! get fire grid sizes from input file
implicit none

!*** arguments
character(len=*),intent(in)::filename
type(grid_info),intent(inout)::info

!*** local
integer::ncid,it,ndims,idim,ratio,stagid,i
integer,dimension(mdims)::dims,sr,dimids,stag
real::dx,dy
character(len=128)::msg
character(len=NF90_MAX_NAME)::dimnames(mdims),dimname,stagname,stagnames(mdims)
character(len=*), parameter:: subgrid='_subgrid'

!*** executable

! get time string information
call read_var_info(filename,var_times,ndims,dims,dimnames)
if(ndims.ne.2)call crash(trim(var_times)//' must have 2 dimensions')
info%dim_time_string = dimnames(1)
info%len_time_string = dims(1)
if(dims(1).gt.max_times_length)call crash('time string too long')
  
! get dimensions
call read_var_info(filename,var_nfuel_cat,ndims,dims,dimnames)

! store in info what you can now
info%io_nfirex=dims(1)
info%io_nfirey=dims(2)
info%ntimes=dims(3)
info%dim_fire_x=dimnames(1)
info%dim_fire_y=dimnames(2)
info%dim_time=dimnames(3)

! attributes
call ncopen(filename,nf90_nowrite,ncid)
call att_read(ncid,'DX',info%fdx)
call att_read(ncid,'DY',info%fdy)
call att_read(ncid,'DT',info%dt)

sr=0
stag=0
stagnames=' '
if(compat_fire_grid)then
  do idim=1,2
    dimname=dimnames(idim)
    i=index(dimname,subgrid,.true.)
    if ( i + len(subgrid) .eq. len_trim(dimname)+1) then ! subgrid, need to fix
        stagname=dimname(1:i-1) // '_stag'
        call check(nf90_inq_dimid(ncid, stagname, stagid))
        call check(nf90_inquire_dimension(ncid, stagid, len=stag(idim)))
        sr(idim) = dims(idim)/stag(idim)
        if(debug_print)write(*,'(3a,i5,a,i5)')'dimension ',trim(stagname), &
           ' length ',stag(idim),' ratio ',sr(idim)
    endif
    stagnames(idim)=stagname
    dims(idim) = dims(idim) - sr(idim)
  enddo
  if((sr(1).ne.0.and.sr(2).eq.0).or.sr(1).lt.0.or.sr(2).lt.0)then
    write(msg,'(a,2i5)')'bad subgrid refinement ratios',sr(1),sr(2)
    call crash(msg)
  endif
endif

call check(nf90_close(ncid))

! store the rest in info
info%nfirex=dims(1)
info%nfirey=dims(2)
info%sr_x=sr(1)
info%sr_y=sr(2)
info%nstagx=stag(1)
info%nstagy=stag(2)
info%dim_atm_x_s=stagnames(1)
info%dim_atm_y_s=stagnames(2)
if(sr(1).ne.0)then
    info%fdx=info%fdx/sr(1)
    info%fdy=info%fdy/sr(2)
endif
	
if(debug_print)then
  write(*,'(4(a,1x))')'dimension names:',trim(info%dim_fire_x),trim(info%dim_fire_y), &
     trim(info%dim_time)
  write(*,'(a,2i6)')'fire grid dimensions:',info%nfirex,info%nfirey
  write(*,'(a,i6)')'number of time frames',info%ntimes
  write(*,'(3(a,f8.4,1x))')'stepsizes fdx=',info%fdx,'fdy=',info%fdy,'dt=',info%dt
endif

end subroutine read_info

subroutine dim_read(ncid,dim_name,dim_len)
character(len=*),intent(in)::dim_name
integer, intent(in)::ncid
integer,intent(out):: dim_len
integer:: dim_id
call check(nf90_inq_dimid(ncid,dim_name,dim_id))
call check(nf90_inquire_dimension(ncid,dim_id,len=dim_len))
if(debug_print)write(*,'(a,1x,a,i6)')'dimension',trim(dim_name),dim_len
end subroutine dim_read

subroutine att_read(ncid,att_name,att_val)
integer, intent(in)::ncid
character(len=*),intent(in)::att_name
real, intent(out)::att_val
call check(nf90_get_att(ncid,nf90_global,att_name,att_val))
if(debug_print)write(*,'(a,1x,a,g20.5)')'attribute',trim(att_name),att_val
end subroutine att_read

!
!***
!

subroutine print_var_info(filename,varname)

!*** arguments
character(len=*), intent(in)::filename,varname

!*** local
integer::ncid,ndims,dimlengths(mdims),dimids(mdims),i,type
character(len=NF90_MAX_NAME)::dimnames(mdims)

!*** executable
if(debug_print)write(*,'(4a)')'reading file ',filename,' dimensions of variable ',trim(varname)
call read_var_info(filename,varname,ndims,dimlengths,dimnames)
if(debug_print)write(*,'(3a,4(2a,i5,1x))')'variable ',trim(varname),' dimensions ', &
  (trim(dimnames(i)),'=',dimlengths(i),i=1,ndims)

end subroutine print_var_info

!
!***
!

subroutine read_var_info(filename,varname,ndims,dimlengths,dimnames,type)

! get variable dimensions from a file
implicit none

! arguments
character(len=*), intent(in):: filename ! variable name
character(len=*), intent(in):: varname ! variable name
integer, intent(out)::ndims            ! number of dimensions of this variable
integer, intent(out)::dimlengths(mdims)      ! the dimensions
character(len=NF90_MAX_NAME),intent(out)::dimnames(mdims) ! dimension names
integer, intent(out), optional :: type

! local
integer :: ncid             ! open netcdf file
integer::varid,dimid,idimids(mdims),idim,i,dimlen,xtype
character(len=NF90_MAX_NAME)::dimname

! executable
if(debug_print)write(*,'(4a)')'reading file ',trim(filename),' variable ',trim(varname)
call ncopen(filename,nf90_nowrite,ncid)
call check(nf90_inq_varid(ncid, varname, varid))
call check(nf90_inquire_variable(ncid,varid,ndims=ndims,xtype=xtype))
if(ndims.gt.mdims)then
    write(*,1)'variable ',trim(varname),' has ',ndims,' dimensions >',mdims
    call crash('variable has too many dimensions')
1 format(3a,i5,a,(4i6))
endif
call check(nf90_inquire_variable(ncid, varid, dimids=idimids))
do idim=1,ndims
    if(debug_print)write(*,'(a,i3,a,i6)')'inquiring dimension ',idim,' id',idimids(idim)
    call check(nf90_inquire_dimension(ncid, idimids(idim), dimname, dimlen))
    if(debug_print)write(*,'(3a,i5)')'got dimension name ',trim(dimname),' length',dimlen
    dimnames(idim)=dimname
    dimlengths(idim) = dimlen  
enddo
call check(nf90_close(ncid))
if(present(type)) type=xtype
2 format(3a,i4,a,(4i6))
if(debug_print)write(*,2)'variable ',trim(varname), &
   ' type',xtype,' dimensions ',(dimlengths(idim),idim=1,ndims)
end subroutine read_var_info

!
!***
!

subroutine write_fire_var(ncid,info,iframe,varname,v)
implicit none
!*** arguments
integer, intent(in)::ncid             ! id of netcdf file open in data mode
type(grid_info),intent(in)::info      ! dimensions
integer, intent(in):: iframe          ! number of frame in the file
character(len=*),intent(in)::varname  ! the variable name
real, pointer, intent(in):: v(:,:)            ! values
!*** local
integer::varid,nx,ny,i,j
real,pointer::v2(:,:)
real::vmin,vmax,err
!*** executable
nx=info%nfirex
ny=info%nfirey
allocate(v2(nx,ny))

vmin=huge(vmin)
vmax=-huge(vmax)
do j=1,ny
  do i=1,nx
    vmin=min(vmin,v(i,j))
    vmax=max(vmax,v(i,j))
  enddo
enddo

if(debug_print)write(*,'(3a,2i5,2(a,g19.6))')'writing variable ',trim(varname), &
    ' size ',nx,ny,' min',vmin,' max',vmax
call check(nf90_inq_varid(ncid,varname,varid),'cannot find '//trim(varname))
call check(nf90_put_var(ncid,varid,v(1:nx,1:ny),start=(/1,1,iframe/),count=(/nx,ny,1/)), &
     'error writing '//trim(varname))

if(read_check)then
  call read_fire_var(ncid,info,iframe,varname,v2) ! read back to check
  err=0
  do j=1,ny
    do i=1,nx
      err=max(err,v(i,j)-v2(i,j))
    enddo
  enddo
  write(*,'(a,g19.6)')'max write-read error',err
endif

end subroutine write_fire_var

!
!***
!

subroutine write_integer(ncid,iframe,varname,iv)
implicit none
!*** arguments
integer, intent(in)::ncid             ! id of netcdf file open in data mode
integer, intent(in):: iframe          ! number of frame in the file
character(len=*),intent(in)::varname  ! the variable name
integer, intent(in):: iv                  ! values
!*** local
integer::iv2,ierr
integer::varid
!*** executable

if(debug_print)write(*,'(3a,i5,a,i10)')'writing int ',trim(varname),'(',iframe,')=',iv
call check(nf90_inq_varid(ncid,varname,varid),'cannot find '//trim(varname))
call check(nf90_put_var(ncid,varid,iv,start=(/iframe/)), &
     'error writing '//trim(varname))
if(read_check)then
  call read_integer(ncid,iframe,varname,iv2) ! read back to check
  ierr=iv2-iv
  write(*,'(a,i6)')'write-read error',ierr
endif

end subroutine write_integer
!
!***
!

subroutine write_real(ncid,iframe,varname,v)
implicit none
!*** arguments
integer, intent(in)::ncid             ! id of netcdf file open in data mode
integer, intent(in):: iframe          ! number of frame in the file
character(len=*),intent(in)::varname  ! the variable name
real, intent(in):: v                  ! values
!*** local
real::v2,val(1),err
integer::varid
!*** executable

if(debug_print)write(*,'(3a,i5,a,g19.6)')'writing real ',trim(varname),'(',iframe,')=',v
call check(nf90_inq_varid(ncid,varname,varid),'cannot find '//trim(varname))
val(1)=v
call check(nf90_put_var(ncid,varid,val,start=(/iframe/),count=(/1/)), &
     'error writing '//trim(varname))
if(read_check)then
  call read_real(ncid,iframe,varname,v2) ! read back to check
  err=v2-v
  write(*,'(a,g19.6)')'write-read error',err
endif

end subroutine write_real

!
!***
!

subroutine read_fire_var(ncid,info,iframe,varname,v)
implicit none
!*** arguments
integer, intent(in)::ncid             ! id of netcdf file open in data mode
type(grid_info),intent(in)::info      ! dimensions
integer, intent(in):: iframe          ! number of frame in the file
character(len=*),intent(in)::varname  ! the variable name
real, pointer :: v(:,:)            ! values
!*** local
integer::varid,nx,ny,start(3),count(3)
!*** executable
nx=info%nfirex
ny=info%nfirey

if(debug_print)write(*,'(2a)')'reading variable ',trim(varname)
if(debug_print)write(*,'(a,4i10)')'lower bounds',lbound(v)
if(debug_print)write(*,'(a,4i10)')'upper bounds',ubound(v)
call check(nf90_inq_varid(ncid,varname,varid),'cannot find '//trim(varname))
call check(nf90_get_var(ncid,varid,v(1:nx,1:ny),start=(/1,1,iframe/),count=(/nx,ny,1/)), &
     'error reading '//trim(varname))
if(debug_print)write(*,'(a,2i5,a,1x,e15.5,3(1x,a,e15.5))')'dimensions',nx,ny, &
     'min',minval(v(1:nx,1:ny)),'max',maxval(v(1:nx,1:ny)),'(1,1)=',v(1,1),'end=',v(nx,ny)

end subroutine read_fire_var

!
!***
!

subroutine read_integer(ncid,iframe,varname,iv)
implicit none
!*** arguments
integer, intent(in)::ncid             ! id of netcdf file open in data mode
integer, intent(in):: iframe          ! number of frame in the file
character(len=*),intent(in)::varname  ! the variable name
integer :: iv                ! value
!*** local
integer::varid

!*** executable

if(debug_print)write(*,'(3a,i5)')'reading integer',trim(varname),' timestep ',iframe
call check(nf90_inq_varid(ncid,varname,varid),'cannot find '//trim(varname))
call check(nf90_get_var(ncid,varid,iv,start=(/iframe/)), &
    'error reading '//trim(varname))
if(debug_print)write(*,'(2a,i5,a,i8)')trim(varname),'(',iframe,')=',iv

end subroutine read_integer

!
!***
!

subroutine read_real(ncid,iframe,varname,v)
implicit none
!*** arguments
integer, intent(in)::ncid             ! id of netcdf file open in data mode
integer, intent(in):: iframe          ! number of frame in the file
character(len=*),intent(in)::varname  ! the variable name
real :: v                ! value
!*** local
integer::varid
real:: val(1)

!*** executable

if(debug_print)write(*,'(3a,i5)')'reading real ',trim(varname),' timestep ',iframe
call check(nf90_inq_varid(ncid,varname,varid),'cannot find '//trim(varname))
call check(nf90_get_var(ncid,varid,val,start=(/iframe/),count=(/1/)), &
    'error reading '//trim(varname))
v=val(1)
if(debug_print)write(*,'(2a,i5,a,g18.8)')trim(varname),'(',iframe,')=',v

end subroutine read_real

!
!***
!

subroutine read_vars(filename,info,iframe,grid)
! read all variables from input file
implicit none

!*** arguments
character(len=*),intent(in)::filename ! the input file
type(grid_info),intent(inout)::info   ! dimensions
integer, intent(in):: iframe          ! number of frame in the file
type(domain),intent(inout)::grid      ! the mother of all arrays

!*** local
integer::ncid,varid,ierr
integer,dimension(4)::s,c

!*** executable
if(debug_print_f)write(*,'(3a,i4)')'read_vars: reading file ',trim(filename),' frame ',iframe 
call ncopen(filename,nf90_nowrite,ncid)

call check(nf90_inq_varid(ncid,var_times,varid),'cannot find '//trim(var_times))
call check(nf90_get_var(ncid,varid,info%times,start=(/1,iframe/),count=(/info%len_time_string,1/)), &
     'error reading '//trim(var_times))
if(debug_print_f)write(*,'(2a)')'Time ',info%times

call read_real(ncid,iframe,var_unit_fxlong,grid%unit_fxlong)
call read_real(ncid,iframe,var_unit_fxlat,grid%unit_fxlat)
call read_fire_var(ncid,info,iframe,var_nfuel_cat,grid%nfuel_cat) 
call read_fire_var(ncid,info,iframe,var_dzdxf,grid%dzdxf) 
call read_fire_var(ncid,info,iframe,var_dzdyf,grid%dzdyf) 
call read_fire_var(ncid,info,iframe,var_zsf,grid%zsf) 
call read_fire_var(ncid,info,iframe,var_fxlong,grid%fxlong) 
call read_fire_var(ncid,info,iframe,var_fxlat,grid%fxlat) 
call read_fire_var(ncid,info,iframe,var_fmc_g,grid%fmc_g) 
!call read_fire_var(ncid,info,iframe,var_unit_fxlong,grid%unit_fxlong) 
!call read_fire_var(ncid,info,iframe,var_unit_fxlat,grid%unit_fxlat) 
call read_fire_var(ncid,info,iframe,var_uf,grid%uf) 
call read_fire_var(ncid,info,iframe,var_vf,grid%vf) 

call check(nf90_close(ncid))
end subroutine read_vars

!
!***
!

subroutine check(ncerr,msg,cont)
implicit none
integer,intent(in)::ncerr
character(len=*), optional, intent(in)::msg
logical, intent(in), optional::cont
character(len=128)::message
if(ncerr.ne.nf90_noerr)then
  write(6,'(2a)')"NetCDF error: ",trim(nf90_strerror(ncerr))
  if(present(msg))then
    message=msg
  else
    message="NETCDF ERROR"
  endif
  if(present(cont))then
    if(cont)return
  endif
  call crash(message)
endif
end subroutine check

!
!***
!

subroutine ncopen(filename,mode,ncid)
!*** purpose: open netcdf file with an informative error message 
implicit none
!*** arguments
character(len=*), intent(in):: filename
integer, intent(in)::mode
integer, intent(out):: ncid
!*** executable
call check(nf90_open(filename,mode,ncid),"Cannot open file "//trim(filename))
end subroutine ncopen

end module wrf_netcdf
