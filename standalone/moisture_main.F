!
!*** moisture_main.F
!
! This is the main driver for the standalone moisture model.  To compile,
! create a make.inc file in this directory from one of the example 
! make.inc.* files.  Make sure to set NETCDF to the prefix path of your
! netcdf installation.  Then, type `make moisture` to create moisture.exe.
!
! This program requires one or more wrfout files to run.  It relies on 
! surface variables T2, Q2, and PSFC to be computed from the WRF run,
! which generally requires a nonzero sfclay option.  The command line
! syntax for this program is:
!
!   ./moisture.exe [-r] wrfout1 [wrfout2 [wrfout3 [...]]]
!
! -r, when present indicates a restart run.  The program will open up a 
! previous output file named moisture.nc and use the last time step
! present in that file to initialize a new run.
!
! The wrfout files specified on the command line must be in sequential
! order, but they can contain overlapping time steps.  This is useful
! for outputs from operational reanalysis runs.  When two files contain
! data for the same time step, this program will always use the second
! file (with time step closer to the reanalysis).  
!
! For every time step computed by this program, it will output surface
! fields to moisture.nc, a netcdf file with similar conventions as a
! wrfout file.

module moisture_util
use module_fr_sfire_util, only : crash
use esmf_mod
use netcdf
implicit none

! size parameters that depend on WRF
integer,parameter::TIMESTRLEN=19,nfmc=5

! if true, the program will only output moisture fields
! but the outputs cannot be used for restarting
logical,parameter::smalloutput=.false.

! model array allocation size for compatibility with WRF
integer,save:: ims,ime,jms,jme

! model array domain size
integer,save::its,ite,jts,jte

! initially set to false to indicate that the sizes above
! have not been initialized
logical,private,save::initialized=.false.

! a simple data structure to group together all variables used
! by the moisture model at a particular time step
type ncvars
    ! model variables
    real,pointer,dimension(:,:,:)::fmc_gc,fmc_equi,fmc_tend
    real,pointer,dimension(:,:)::t2,q2,psfc,rainc,rainnc, &
               t2_old,q2_old,psfc_old,rain_old,rh_fire
    ! time stamp
    type(ESMF_Time)::time

    ! rain accumulation variable
    real,pointer,dimension(:,:)::rain_accum,rain_zero
end type

! a data structure that contains metadata information about a 
! wrfout file
type ncfile
    ! number of time steps in the file
    integer::ntime

    ! array of time stamps contained in the file
    type(ESMF_Time),dimension(:),pointer::times

    ! the start time of the simulation
    type(ESMF_Time)::starttime

    ! file handle
    integer::ncid
end type

contains

subroutine initialize_and_check(filename,wrffile)
! open a file and initialize an ncfile type
! if module is not initialized, then initialize
! otherwise check that sizes are consistent

! the name of the file to open
character(len=*),intent(in)::filename

! ncfile type initialized on output
type(ncfile),intent(out)::wrffile

! local variables dealing with netcdf i/o
integer::dimidx,dimidy,dimidtime,dimidfuel
integer::varidtime
integer::nx,ny,nfuel,i,rc
character(len=TIMESTRLEN),dimension(:),allocatable::timestr
character(len=TIMESTRLEN)::stimestr

! open the file and collect metadata
call check(nf90_open(filename,nf90_nowrite,wrffile%ncid))
call check(nf90_inq_dimid(wrffile%ncid,'west_east',  dimidx))
call check(nf90_inq_dimid(wrffile%ncid,'south_north',dimidy))
call check(nf90_inq_dimid(wrffile%ncid,'Time',dimidtime))
call check(nf90_inquire_dimension(wrffile%ncid,dimidx,len=nx))
call check(nf90_inquire_dimension(wrffile%ncid,dimidy,len=ny))
call check(nf90_inquire_dimension(wrffile%ncid,dimidtime,len=wrffile%ntime))

if(.not.initialized)then
    ! the ESMF time modules need to be initialized on first call
    call ESMF_Initialize(defaultCalKind=ESMF_CAL_GREGORIAN)
    
    ! on first call set module level variables containing array sizes
    ims=1
    ime=nx
    jms=1
    jme=ny
    its=ims
    ite=ime
    jts=jms
    jte=jme
    initialized=.true.
else
    ! on later calls check that the file sizes are consistant with the sizes
    ! as stored
    if( (nx.ne.ime-ims+1) .or. &
        (ny.ne.jme-jms+1) ) then
        call crash('Incompatible file: '//trim(filename))
    endif
endif

! allocate memory for ncfile data type
allocate(wrffile%times(wrffile%ntime))
allocate(timestr(wrffile%ntime))

call check(nf90_inq_varid(wrffile%ncid,'Times',varidtime))
call check(nf90_get_var(wrffile%ncid,varidtime,timestr))

! go through all time stamps in the file and store them as ESMF types
do i=1,wrffile%ntime
    call parse_wrf_time(timestr(i),wrffile%times(i))
enddo

! get start time
call check(nf90_get_att(wrffile%ncid,NF90_GLOBAL,'START_DATE',stimestr))
call parse_wrf_time(stimestr,wrffile%starttime)

! free local memory
deallocate(timestr)

end subroutine initialize_and_check

subroutine alldone()
implicit none
! called at then end of execution to free allocated memory
call ESMF_Finalize()
end subroutine alldone

subroutine destroy_file(wrffile)
implicit none
! free memory allocated in a ncfile data structure
type(ncfile),intent(inout)::wrffile
call check(nf90_close(wrffile%ncid))
deallocate(wrffile%times)
wrffile%ntime=0
end subroutine destroy_file

subroutine initialize_vars(vars)
implicit none
! Allocate memory for the model arrays
type(ncvars),intent(inout)::vars
if(.not.initialized)then
  call crash('Module not initialized')
endif
allocate(vars%fmc_gc(ims:ime,1:nfmc,jms:jme))
allocate(vars%fmc_equi(ims:ime,1:nfmc,jms:jme))
allocate(vars%fmc_tend(ims:ime,1:nfmc,jms:jme))
allocate(vars%t2(ims:ime,jms:jme))
allocate(vars%q2(ims:ime,jms:jme))
allocate(vars%psfc(ims:ime,jms:jme))
allocate(vars%rainc(ims:ime,jms:jme))
allocate(vars%rainnc(ims:ime,jms:jme))
allocate(vars%t2_old(ims:ime,jms:jme))
allocate(vars%q2_old(ims:ime,jms:jme))
allocate(vars%psfc_old(ims:ime,jms:jme))
allocate(vars%rain_old(ims:ime,jms:jme))
allocate(vars%rh_fire(ims:ime,jms:jme))
allocate(vars%rain_accum(ims:ime,jms:jme))
allocate(vars%rain_zero(ims:ime,jms:jme))
vars%rain_accum(:,:)=0.
vars%rain_zero(:,:)=0.
end subroutine initialize_vars

subroutine destroy_vars(vars)
implicit none
! free memory from model arrays
type(ncvars),intent(inout)::vars
deallocate(vars%fmc_gc)
deallocate(vars%fmc_equi)
deallocate(vars%fmc_tend)
deallocate(vars%t2)
deallocate(vars%q2)
deallocate(vars%psfc)
deallocate(vars%rainc)
deallocate(vars%rainnc)
deallocate(vars%t2_old)
deallocate(vars%q2_old)
deallocate(vars%psfc_old)
deallocate(vars%rain_old)
deallocate(vars%rh_fire)
deallocate(vars%rain_accum)
deallocate(vars%rain_zero)
end subroutine destroy_vars

subroutine create_output(filename)
implicit none
! create an output netcdf file with metadata initialized
! for the model being run.... does not write any data
character(len=*),intent(in)::filename

integer::ncid,nxid,nyid,timeid,fuelid,dateid
integer::tmp

! create dataset
call check(nf90_create(filename,nf90_clobber,ncid))

! create dimensions
call check(nf90_def_dim(ncid,'west_east',ime-ims+1,nxid))
call check(nf90_def_dim(ncid,'south_north',jme-jms+1,nyid))
call check(nf90_def_dim(ncid,'Time',nf90_unlimited,timeid))
call check(nf90_def_dim(ncid,'fuel_moisture_classes_stag',nfmc,fuelid))
call check(nf90_def_dim(ncid,'DateStrLen',TIMESTRLEN,dateid))

! create variables
call check(nf90_def_var(ncid,'Times',NF90_CHAR,(/dateid,timeid/),tmp))
call check(nf90_def_var(ncid,'FMC_GC',NF90_REAL,(/nxid,nyid,fuelid,timeid/),tmp))
if(.not.smalloutput)then
    call check(nf90_def_var(ncid,'FMC_EQUI',NF90_REAL,(/nxid,nyid,fuelid,timeid/),tmp))
    call check(nf90_def_var(ncid,'FMC_TEND',NF90_REAL,(/nxid,nyid,fuelid,timeid/),tmp))
    call check(nf90_def_var(ncid,'T2',NF90_REAL,(/nxid,nyid,timeid/),tmp))
    call check(nf90_def_var(ncid,'Q2',NF90_REAL,(/nxid,nyid,timeid/),tmp))
    call check(nf90_def_var(ncid,'PSFC',NF90_REAL,(/nxid,nyid,timeid/),tmp))
    call check(nf90_def_var(ncid,'RAIN',NF90_REAL,(/nxid,nyid,timeid/),tmp))
    call check(nf90_def_var(ncid,'RH_FIRE',NF90_REAL,(/nxid,nyid,timeid/),tmp))
endif
call check(nf90_close(ncid))
end subroutine create_output

subroutine write_output(filename,vars)
implicit none
! write arrays from the model into a preinitialized output file
character(len=*),intent(in)::filename
type(ncvars),intent(in)::vars

integer::ncid,varid,dimid,itime
character(len=TIMESTRLEN)::timestr
real,dimension(:,:,:),allocatable::tmp

! allocate temporary memory for transposing 3D arrays
allocate(tmp(its:ite,jts:jte,1:nfmc))

print*,'writing output to '//trim(filename)
call check(nf90_open(filename,NF90_WRITE,ncid))

! get the number of time frames already present and set the
! output time frame to the end of the file
call check(nf90_inq_dimid(ncid,'Time',dimid))
call check(nf90_inquire_dimension(ncid,dimid,len=itime))
itime=itime+1

! write time string
call check(nf90_inq_varid(ncid,'Times',varid))
call ESMF_TimeGet(vars%time,timeString=timestr)
call check(nf90_put_var(ncid,varid,timestr, &
             start=(/1,itime/),count=(/TIMESTRLEN,1/)))

! write FMC_GC
call check(nf90_inq_varid(ncid,'FMC_GC',varid))
call transpose_var(vars%fmc_gc,tmp)
call check(nf90_put_var(ncid,varid,tmp, &
             start=(/its,jts,1,itime/),count=(/ite,jte,nfmc,1/)))

! the following variables are diagnostic only, but useful for restarting
if(.not.smalloutput)then

    call check(nf90_inq_varid(ncid,'FMC_EQUI',varid))
    call transpose_var(vars%fmc_equi,tmp)
    call check(nf90_put_var(ncid,varid,tmp, &
                 start=(/its,jts,1,itime/),count=(/ite,jte,nfmc,1/)))

    call check(nf90_inq_varid(ncid,'FMC_TEND',varid))
    call transpose_var(vars%fmc_tend,tmp)
    call check(nf90_put_var(ncid,varid,tmp, &
                 start=(/its,jts,1,itime/),count=(/ite,jte,nfmc,1/)))

    call check(nf90_inq_varid(ncid,'T2',varid))
    call check(nf90_put_var(ncid,varid,vars%t2, &
                 start=(/its,jts,itime/),count=(/ite,jte,1/)))

    call check(nf90_inq_varid(ncid,'Q2',varid))
    call check(nf90_put_var(ncid,varid,vars%q2, &
                 start=(/its,jts,itime/),count=(/ite,jte,1/)))

    call check(nf90_inq_varid(ncid,'PSFC',varid))
    call check(nf90_put_var(ncid,varid,vars%psfc, &
                 start=(/its,jts,itime/),count=(/ite,jte,1/)))

    call check(nf90_inq_varid(ncid,'RAIN',varid))
    call check(nf90_put_var(ncid,varid,vars%rain_accum, &
                 start=(/its,jts,itime/),count=(/ite,jte,1/)))

    call check(nf90_inq_varid(ncid,'RH_FIRE',varid))
    call check(nf90_put_var(ncid,varid,vars%rh_fire, &
                 start=(/its,jts,itime/),count=(/ite,jte,1/)))

endif

! free allocated memory
call check(nf90_close(ncid))
deallocate(tmp)
end subroutine write_output

subroutine read_file(wrffile,atime,vars,ierr,restart)
implicit none
! read array data for the model from a file 
type(ncfile),intent(in)::wrffile ! initialized ncfile object
integer,intent(in)::atime        ! the time step in the file to read
type(ncvars),intent(inout)::vars ! on return contains data from the file
integer,intent(out)::ierr        ! error flag, 0 on success, otherwise read error
logical,optional,intent(in)::restart ! if this a restart, we will read fmc_gc

integer::varid,i,j,k,itime
logical::r
real,dimension(:,:,:),allocatable::tmp
real,dimension(ims:ime,jms:jme)::rainold

r=.false.
if(present(restart))then
    r=restart
endif

! if the time step is unspecified (<= 0) then
! just read the last time step
if(atime.gt.0)then
    itime=atime
else
    itime=size(wrffile%times)
endif
vars%time=wrffile%times(itime)

! read FMC_GC if it is present (generally, this test will
! only pass when we are reading for a restart)
if(r.and.nf90_inq_varid(wrffile%ncid,'FMC_GC',varid).eq.0)then
    allocate(tmp(its:ite,jts:jte,nfmc))
    call check(nf90_get_var(wrffile%ncid,varid,tmp, &
            start=(/its,jts,1,itime/), count=(/ite,jte,nfmc,1/)))
    call transpose_var(tmp,vars%fmc_gc)
    deallocate(tmp)
endif

! read variables that should always be present
call check(nf90_inq_varid(wrffile%ncid,'T2',varid))
call check(nf90_get_var(wrffile%ncid,varid,vars%t2, &
        start=(/its,jts,itime/), count=(/ite,jte,1/)))

call check(nf90_inq_varid(wrffile%ncid,'Q2',varid))
call check(nf90_get_var(wrffile%ncid,varid,vars%q2, &
        start=(/its,jts,itime/), count=(/ite,jte,1/)))

call check(nf90_inq_varid(wrffile%ncid,'PSFC',varid))
call check(nf90_get_var(wrffile%ncid,varid,vars%psfc, &
        start=(/its,jts,itime/), count=(/ite,jte,1/)))

if(nf90_inq_varid(wrffile%ncid,'RAINC',varid).eq.0)then
    call check(nf90_get_var(wrffile%ncid,varid,vars%rainc, &
            start=(/its,jts,itime/), count=(/ite,jte,1/)))
    
    if(itime.gt.1)then
        call check(nf90_get_var(wrffile%ncid,varid,rainold, &
               start=(/its,jts,itime-1/),count=(/ite,jte,1/)))
        vars%rainc(:,:)=vars%rainc(:,:)-rainold(:,:)
    endif

    call check(nf90_inq_varid(wrffile%ncid,'RAINNC',varid))
    call check(nf90_get_var(wrffile%ncid,varid,vars%rainnc, &
            start=(/its,jts,itime/), count=(/ite,jte,1/)))
    
    if(itime.gt.1)then
        call check(nf90_get_var(wrffile%ncid,varid,rainold, &
               start=(/its,jts,itime-1/),count=(/ite,jte,1/)))
        vars%rainnc(:,:)=vars%rainnc(:,:)-rainold(:,:)
    endif
else
    call check(nf90_inq_varid(wrffile%ncid,'RAIN',varid))
    call check(nf90_get_var(wrffile%ncid,varid,vars%rainc, &
            start=(/its,jts,itime/), count=(/ite,jte,1/)))
    vars%rainnc(:,:)=0.
    
    if(itime.gt.1)then
        call check(nf90_get_var(wrffile%ncid,varid,rainold, &
               start=(/its,jts,itime-1/),count=(/ite,jte,1/)))
        vars%rainc(:,:)=vars%rainc(:,:)-rainold(:,:)
    endif
endif

! calculate total rain accumulation for this moisture run
vars%rain_accum(:,:)=vars%rain_accum(:,:) + vars%rainc(:,:) + vars%rainnc(:,:)

! here we preempt any model errors by detecting if the data is correct
! if not, the main code will skip this time step and move on to the
! next
if(any(vars%t2.le.0).or.any(vars%q2.lt.0).or.any(vars%psfc.le.0))then
    ierr=1
else
    ierr=0
endif
end subroutine read_file

subroutine transpose_var(A,B)
implicit none
! do a transpose on a 3D array to account for different
! storage conventions between wrf output and runtime
! B(i,k,j)=A(i,j,k)
real,dimension(:,:,:),intent(in)::A
real,dimension(:,:,:),intent(out)::B

integer::i,j,k
if( size(A,1).ne.size(B,1) .or. &
    size(A,2).ne.size(B,3) .or. &
    size(A,3).ne.size(B,2) )then
  call crash('Invalid array sizes in transpose')
endif
do j=1,size(A,2)
    do k=1,size(A,3)
        do i=1,size(A,1)
            B(i,k,j)=A(i,j,k)
        enddo
    enddo
enddo
end subroutine transpose_var

subroutine parse_wrf_time(wrfstr,time)
! parse a wrf time string into an esmf time class
! example wrf time string: "2011-01-01_00:00:00"
implicit none
character(len=TIMESTRLEN),intent(in)::wrfstr
type(ESMF_Time),intent(out)::time

integer::year,month,day,hour,minute,second
character(len=1)::a1,a2,a3,a4,a5

read(wrfstr,'(I04,A1,I02,A1,I02,A1,I02,A1,I02,A1,I02)') &
       year,a1,month,a2,day,a3,hour,a4,minute,a5,second
call ESMF_TimeSet(time,YY=year,MM=month,DD=day,H=hour,M=minute,S=second)
end subroutine parse_wrf_time

subroutine get_next_timestep(files,ifile,istep)
! choose a file and time step within the file to use
! for the next call to the moisture model
! 
! This assumes files are in sequential order, and 
! allows for reanalysis runs where there may be
! time step overlaps from one file to the next.
! Here, we prefer to move on to the next file
! rather than the next time step in the old file
! to use model arrays closer to reanalysis.
! If there is no overlap, this routine will
! increment from one file to the next as expected.
!
! This subroutine will return istep=ifile=0, when
! there are no more steps to process.
implicit none
type(ncfile),dimension(:),intent(in)::files
integer,intent(inout)::ifile,istep

integer::nfile,nstep
type(ESMF_Time)::time1,time2,time3
nfile=size(files)
if(ifile.lt.1)then
    ! for initialization only
    ifile=1
    istep=1
    return
endif
nstep=size(files(ifile)%times)
time1=files(ifile)%times(istep)

if(istep+1.le.size(files(ifile)%times))then
    time2=files(ifile)%times(istep+1)
else
    ifile=ifile+1
    istep=1
    goto 990
endif

if(ifile+1.le.nfile)then
    ! if there is another file after the current one
    time3=files(ifile+1)%times(1)
else
    ! if there is not another file pick next time step in current file
    istep=istep+1
    goto 990
endif

if(time2.ge.time3)then
    ! if next file has a time step closer to the current time
    ! or if there are no more time steps in the current file,
    ! use the next file
    ifile=ifile+1
    istep=1
else
    ! otherwise use the next step in the current file
    istep=istep+1
endif

 990 continue
if(ifile.gt.nfile.or.istep.gt.size(files(ifile)%times))then
    ! make sure the current time step is valid
    ! end if not
    istep=0
    ifile=0
endif
end subroutine get_next_timestep

subroutine check(status)
integer,intent(in)::status
! check a netcdf status code for error
if(status.ne.nf90_noerr) then
    call crash(trim(nf90_strerror(status)))
endif
end subroutine check

end module moisture_util

program moisture_main
use moisture_util
use module_fr_sfire_phys, only: advance_moisture, init_fuel_cats, moisture_classes
use module_fr_sfire_util, only: fire_print_msg 
implicit none

logical::initialize
integer::fmoist_init=2
real::moisture_dt,dt
integer::istep,i,iglobal
logical::restart

integer :: numarg,ifile,ii,nfile,s,sn,sd,ierr
character(len=120) :: wrfoutfname
character(len=*),parameter :: outputfname='moisture.nc'
character(len=64)::timediffstr

type(ncfile),dimension(:),allocatable::wrffile
type(ncfile)::wrfoutfile
type(ncvars)::vars
type(ESMF_Time)::oldtime,timenow,starttime
type(ESMF_TimeInterval)::timedt

integer,parameter::iprint=100,jprint=100
integer::k
real::hours,t,q,p,rain
real,dimension(5)::fmc_equi,fmc_tend,fmc_gc
logical,parameter::debugprint=.false.
integer,parameter::debugfile=8

#ifndef GFORTRAN_IARGC
integer,external::iargc
#endif

1 format(a4,a8,  a6,  a9,  a6  ,a8, a6, 5(a5,i1,a8,i1,a5,i1))
2 format(a4,a8,  a6,  a9,  a6  ,a8, a6, 5(a6,a9,a6))
3 format(i4,f8.2,f6.1,f9.1,f6.3,f8.1,f6.3,5(f6.3,f9.1,f6.3))

if(debugprint)then
    open(debugfile,file='moisture_output.txt',form='formatted',status='unknown')
    write(debugfile,1)'Step','Time', 'T','P','Q','RAIN','RH',('EQUI',k,'TLAG',k,'FMC',k,k=1,moisture_classes)
    write(debugfile,2)' ','hours','K','Pa','kg/kg','mm','1',('kg/kg','hours','kg/kg',k=1,moisture_classes)
endif
hours=0

! parse commandline argument for wrfout file name
numarg=iargc()

if(numarg .lt. 1)then
    print*, 'usage: ./moisture.exe [-r] wrfout1 [wrfout2 [wrfout3 [...]]]'
    call abort()
endif

! check for restart flag
call getarg(1,wrfoutfname)

! allocate local storage for wrfout metadata
if(trim(wrfoutfname).eq.'-r')then
    allocate(wrffile(numarg-1))
else
    allocate(wrffile(numarg))
endif

! open files and populate ncfile objects for each
! check that the sizes are consistent
nfile=0
restart=.false.
do i=1,numarg
    call getarg(i,wrfoutfname)
    if(trim(wrfoutfname).eq.'-r')then
        restart=.true.
        call initialize_and_check(outputfname,wrfoutfile)
    else
        print*,'using: '//trim(wrfoutfname)
        nfile=nfile+1
        call initialize_and_check(wrfoutfname,wrffile(nfile))
    endif
enddo

! allocate memore for model arrays
call initialize_vars(vars)

oldtime=wrffile(1)%times(1)
timenow=oldtime

! if this isn't a restart run, we need to create
! a new output file, clobber anything that might
! be there, otherwise we append data to the restart
if(.not.restart) &
call create_output(outputfname)

! model initialization calls
fire_print_msg = 2
call init_fuel_cats(.true.)
print *,moisture_classes,' moisture classes'

! initialize various looping indices
ifile=0
ii=0
istep=0

! get the first file/time step we will read
call get_next_timestep(wrffile,ifile,ii)
if(restart)then
    starttime=wrfoutfile%times(size(wrfoutfile%times))
else
    starttime=wrffile(1)%times(1)
endif

! for restarts, we have to find the first time step in the 
! wrfout's after the the time in the restart
do while(restart.and.wrffile(ifile)%times(ii).le.starttime)
    call get_next_timestep(wrffile,ifile,ii)
    if(ifile.eq.0)goto 1000
enddo
 1000 continue

! after all the initialization code, we crash here if there are no
! new time steps to read
if(ifile.eq.0)then
    call crash('no new time steps available')
endif

! let's tell the user what we are doing
print*,'starting at file ',ifile,' step ',ii
call ESMF_TimeGet(wrffile(ifile)%times(ii),TimeString=timediffstr)
print*,'wrffile time=',trim(timediffstr)
call ESMF_TimeGet(starttime,TimeString=timediffstr)
print*,'(re)start time=',trim(timediffstr)

! get the next target time step
call get_next_timestep(wrffile,ifile,ii)

iglobal=0

! main program loop, go until get_next_timestep returns ifile==0
do while(ifile.gt.0)

    ! keep track of the number of loops we have made
    istep=istep+1
    iglobal=iglobal+1
    
    ! read initialization data
    if(restart)then
        call read_file(wrfoutfile,-1,vars,ierr,restart=.true.)
        if(ierr.ne.0)then
            call crash('invalid restart file')
        endif
        fmoist_init=3
    else
        call read_file(wrffile(ifile),ii,vars,ierr)
        if(ierr.ne.0)then
            print*,'invalid data: skipping file ',ifile,' step ',ii
            goto 10001
        endif
    endif
    
    ! increment time keeping
    oldtime=timenow
    timenow=wrffile(ifile)%times(ii)

    ! let's tell the user what we are doing
    print*,'reading file ',ifile,' of ',nfile
    print*,'time step ',ii,' of ',wrffile(ifile)%ntime
    call ESMF_TimeGet(timenow,timeString=timediffstr)
    print*,'time: ',trim(timediffstr)

    if(iglobal.eq.1)then
        dt=0.  ! required on initialization
        initialize=.true.
    else
        if(timenow.lt.oldtime)then
            call crash('Negative time step detected.')
        endif
        ! get the number of seconds that the model needs to advance
        timedt=timenow-oldtime
        call ESMF_TimeIntervalGet(timedt,S=s,sn=sn,sd=sd,TimeString=timediffstr)
        if(sd.gt.0)then
            dt=dble(s)+dble(sn)/dble(sd)
        else
            dt=dble(s)
        endif
        initialize=.false.
    endif


    call advance_moisture(    &
        initialize,                 & ! initialize timestepping. true on the first call at time 0, then false
        ims,ime,  jms,jme,          & ! memory dimensions
        its,ite,  jts,jte,          & ! tile dimensions
        nfmc,                       & ! dimension of moisture fields
        dt,                & ! timestep = time step time elapsed from the last call
        vars%rain_accum, vars%rain_zero,              & ! accumulated rain
        vars%t2, vars%q2, vars%psfc,               & ! temperature (K), vapor contents (kg/kg), pressure (Pa) at the surface
        vars%rain_old,                   & ! previous value of accumulated rain
        vars%t2_old, vars%q2_old, vars%psfc_old,   & ! previous values of the atmospheric state at surface
        vars%rh_fire,                    & ! relative humidity diagnostics
        vars%fmc_gc,                     & ! fuel moisture by class, updated
        vars%fmc_equi,                   & ! fuel moisture equilibrium by class, for diagnostics only
        vars%fmc_tend                    & ! fuel moisture tendency by classe, for diagnostics only
        )

    if(debugprint)then
        hours=hours+dt/3600.
        t=vars%t2(iprint,jprint)
        q=vars%q2(iprint,jprint)
        p=vars%psfc(iprint,jprint)
        rain=vars%rain_accum(iprint,jprint) + vars%rain_zero(iprint,jprint)
        fmc_equi(:)=vars%fmc_equi(iprint,:,jprint)
        fmc_tend(:)=vars%fmc_tend(iprint,:,jprint)
        fmc_gc(:)=vars%fmc_gc(iprint,:,jprint)
        write(debugfile,3)iglobal,hours,t,p,q,rain, &
                  vars%rh_fire(iprint,jprint),   &
               (fmc_equi(k),fmc_tend(k),fmc_gc(k), &
                k=1,moisture_classes)
    endif
    ! write the results to the output file
    call write_output(outputfname,vars)

    10001 continue ! we skip to here on invalid data

    ! after the first pass through the loop, nothing special needs to be done
    ! for restart runs, so turn of the restart flag here
    if(restart)then
        restart=.false.
    else
        call get_next_timestep(wrffile,ifile,ii)
    endif
enddo ! end main loop

! free up allocated memory
call destroy_vars(vars)
do i=1,nfile
    call destroy_file(wrffile(i))
enddo
call alldone()
end program moisture_main
