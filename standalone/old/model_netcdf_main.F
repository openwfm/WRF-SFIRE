! the main program is at the end because of fortran limitations
#ifndef NF_TYPE
#define NF_TYPE NF90_FLOAT
#endif

module jb_utils
contains

subroutine check(status)
use netcdf
  integer, intent ( in) :: status
  
  if(status /= nf90_noerr) then 
    print *, nf90_strerror(status)
    stop "Stopped"
  end if
end subroutine check 

end module jb_utils


module module_model_test
use module_fr_sfire_util
use module_fr_sfire_model
use module_fr_sfire_phys
use jb_utils
use netcdf

contains

subroutine main_sub
implicit none

!*** declarations

integer:: nx,ny, msteps  ! problem dimension, in cells, number of steps
real:: fdx,fdy, dt,wind,alpha       ! fire mesh spacing (m), time step (s), wind azimuth
integer::rfx,rfy ! refinement, convenience only

nx=400
ny=400
msteps=200
msteps=6
! msteps=100
fdx=6
fdy=6

rfx=1
rfy=1

nx=nx*rfx
ny=ny*rfy
fdx=fdx/rfx
fdy=fdy/rfy

dt=0.5
wind=10
alpha=0.0

print *,'nx=',nx,' ny=',ny,' msteps=',msteps,' fdx=',fdx,' fdy=',fdy,' dt=',dt,' wind=',wind,' alpha=',alpha
call model_test(1,nx,1,ny, &
!   -1,nx+3,-2,ny+2, &
   0,nx+1,0,ny+1, &
   1,nx,1,ny, &
   fdx,fdy,wind,alpha,dt,msteps)

end subroutine main_sub

!
!****************
!


subroutine model_test(  &
    ifds,ifde,jfds,jfde, &
    ifms,ifme,jfms,jfme, &
    ifps,ifpe,jfps,jfpe, &
    fdx,fdy,wind,alpha,dt,msteps)
implicit none

!*** arguments
integer, intent(in):: &
    ifds,ifde,jfds,jfde, &
    ifps,ifpe,jfps,jfpe, &
    ifms,ifme,jfms,jfme, msteps
real, intent(in)::fdx,fdy,dt,wind,alpha

!*** local
real, dimension(ifms:ifme,jfms:jfme):: zsf,     &
                 lfn,tign,fuel_frac,                    &
                 grnhfx,grnqfx,lfn_out
integer::  num_ignitions,i,j,ifuelread,istep,nfuel_cat0,ifun,ifun_start
logical::need_lfn_update
real:: t0,time_start,sm,sn
integer, dimension(ifms:ifme,jfms:jfme)::nfuel_cat,ischap
real, dimension(ifms:ifme,jfms:jfme)::fuel_time,vx,vy,dzfsdx,dzfsdy,bbb,betafl,phiwc,r_0,fgip
integer, parameter::max_tiles=10
integer::num_tiles,ij,ifts,ifte,jfts,jfte
integer, dimension(max_tiles)::if_start,if_end,jf_start,jf_end
integer, parameter :: max_ignitions=10
real, dimension(max_ignitions) :: ignition_start_x,ignition_start_y, &
    ignition_end_x,ignition_end_y,ignition_radius,ignition_time
character (len=NF90_MAX_NAME), parameter :: iofile = "fire_io.nc"
integer :: ncid,xid,yid,tid,lfnid,vfxid,vfyid,tignid,grnhftid,itmp,tdim,tvalsid
integer :: xsid,ysid  ! staggered dims
integer :: ostat,stepnum,grnqftid,fuel_fracid,zsfid,fueltimeid
real :: tvals
integer,parameter :: ndims = 3
integer,dimension(ndims):: dimids,count,start,atmp,counts,starts
integer::sstep

!*** executable

call set_tiles(2,2,ifps,ifpe,jfps,jfpe,num_tiles,if_start,if_end,jf_start,jf_end)

print *,'mesh size in cells: ',ifps,ifpe,jfps,jfpe
print *,'array allocation:   ',ifms,ifme,jfms,jfme

t0=0                                        ! starting time, arbitrary
ostat=nf90_open(iofile,NF90_WRITE,ncid)
print*,ostat,nf90_noerr
if(ostat.ne.NF90_NOERR)then


! populate the arrays somehow
do j=jfps,jfpe
    do i=ifps,ifpe
        zsf(i,j)=1000   ! flat ground
        vx(i,j)=wind*cos(alpha)    ! constant wind
        vy(i,j)=wind*sin(alpha)
    enddo
enddo
call check( nf90_create(iofile,NF90_NOCLOBBER,ncid) )
call check( nf90_def_dim(ncid,"x",ifpe-ifps+1,xid) )
call check( nf90_def_dim(ncid,"y",jfpe-jfps+1,yid) )
call check( nf90_def_dim(ncid,"x_stag",ifpe-ifps+2,xsid) )
call check( nf90_def_dim(ncid,"y_stag",jfpe-jfps+2,ysid) )
call check( nf90_def_dim(ncid,"Time",NF90_UNLIMITED,tid) )
call check( nf90_put_att(ncid,NF90_GLOBAL,"DX",fdx) )
call check( nf90_put_att(ncid,NF90_GLOBAL,"DY",fdy) )
stepnum=0
call check( nf90_put_att(ncid,NF90_GLOBAL,"stepnum",stepnum) )

dimids=(/xsid,ysid,tid/)
call check( nf90_def_var(ncid,"lfn",NF_TYPE,dimids,lfnid) )
call check( nf90_def_var(ncid,"vx",NF_TYPE,dimids,vfxid) )
call check( nf90_def_var(ncid,"vy",NF_TYPE,dimids,vfyid) )
call check( nf90_def_var(ncid,"Times",NF_TYPE,(/tid/),tvalsid) )

dimids=(/xid,yid,tid/)
call check( nf90_def_var(ncid,"tign",NF_TYPE,dimids,tignid) )
call check( nf90_def_var(ncid,"grnhfx",NF_TYPE,dimids,grnhftid) )
call check( nf90_def_var(ncid,"grnqfx",NF_TYPE,dimids,grnqftid) )
call check( nf90_def_var(ncid,"fuelfrac",NF_TYPE,dimids,fuel_fracid) )
call check( nf90_def_var(ncid,"fueltime",NF_TYPE,dimids,fueltimeid) )
call check( nf90_enddef(ncid) )
print*,'askjdfksdjf'
tvals=0.0
tdim=0

else

call check( nf90_inq_dimid(ncid,"x_stag",xsid) )
call check( nf90_inquire_dimension(ncid,xsid,len=itmp) )
if(itmp.ne.ifpe-ifps+2)call check(NF90_EVARSIZE)
call check( nf90_inq_dimid(ncid,"y_stag",ysid) )
call check( nf90_inquire_dimension(ncid,ysid,len=itmp) )
if(itmp.ne.jfpe-jfps+2)call check(NF90_EVARSIZE)

call check( nf90_inq_dimid(ncid,"x",xid) )
call check( nf90_inquire_dimension(ncid,xid,len=itmp) )
if(itmp.ne.ifpe-ifps+1)call check(NF90_EVARSIZE)
call check( nf90_inq_dimid(ncid,"y",yid) )
call check( nf90_inquire_dimension(ncid,yid,len=itmp) )
if(itmp.ne.jfpe-jfps+1)call check(NF90_EVARSIZE)
call check( nf90_inq_dimid(ncid,"Time",tid) )
call check( nf90_inquire(ncid,unlimitedDimID=itmp) )
if(itmp.ne.tid)call check(NF90_EVARSIZE)

call check( nf90_inq_varid(ncid,"lfn",lfnid) )
call check( nf90_inq_varid(ncid,"vx",vfxid) )
call check( nf90_inq_varid(ncid,"vy",vfyid) )
call check( nf90_inq_varid(ncid,"tign",tignid) )
call check( nf90_inq_varid(ncid,"grnhfx",grnhftid) )
call check( nf90_inq_varid(ncid,"grnqfx",grnqftid) )
call check( nf90_inq_varid(ncid,"fuelfrac",fuel_fracid) )
call check( nf90_inq_varid(ncid,"fueltime",fueltimeid) )
call check( nf90_inquire_dimension(ncid,tid,len=tdim) )
call check( nf90_inq_varid(ncid,"Times",tvalsid) )
call check( nf90_get_var(ncid,tvalsid,atmp,start=(/tdim/),count=(/1/)) )
tvals=atmp(1)
dimids=(/xid,yid,tid/)
call check( nf90_get_att(ncid,nf90_global,"stepnum",stepnum) )

count=(/ifpe-ifps+1,jfpe-jfps+1,1/)
counts=(/ifpe-ifps+2,jfpe-jfps+2,1/)
start=(/1,1,tdim/)
starts=(/1,1,tdim/)

#define _SD (ifps:ifpe+1,jfps:jfpe+1)
#define _D  (ifps:ifpe,jfps:jfpe)
call check( nf90_get_var(ncid,lfnid,lfn _SD,start=starts,count=counts) )
call check( nf90_get_var(ncid,tignid,tign _D,start=starts,count=count) )
call check( nf90_get_var(ncid,vfxid,vx _SD,start=starts,count=counts) )
call check( nf90_get_var(ncid,vfyid,vy _SD,start=starts,count=counts) )
call check( nf90_get_var(ncid,grnhftid,grnhfx _D,start=start,count=count) )
call check( nf90_get_var(ncid,grnqftid,grnqfx _D,start=start,count=count) )
call check( nf90_get_var(ncid,fuel_fracid,fuel_frac _D,start=start,count=count) )
call check( nf90_get_var(ncid,fueltimeid,fuel_time _D,start=start,count=count) )
endif


count=(/ifpe-ifps+1,jfpe-jfps+1,1/)
counts=(/ifpe-ifps+2,jfpe-jfps+2,1/)

! fuel data
ifuelread=0
nfuel_cat0=3

! for matlab
open(1,file='model_test_out.txt',form='formatted')
1   format(e25.12e3)
sm=ifpe-ifps+1
sn=jfpe-jfps+1
write(1,1)1.,1.,fdx
write(1,1)1.,1.,fdy

time_start=t0
num_ignitions=2
ignition_start_x(1)=0.5*fdx*(ifde-ifds)
ignition_start_y(1)=0.5*fdy*(jfde-jfds)
ignition_end_x(1)=0.5*fdx*(ifde-ifds)*0.9999999
ignition_end_y(1)=0.5*fdy*(jfde-jfds)*1.0000001
ignition_radius(1) = 0.5*max(5.0,6*max(fdx,fdy))
ignition_time(1)=1
ignition_start_x(2)=1000
ignition_start_y(2)=500
ignition_end_x(2)=1500
ignition_end_y(2)=1500
! at least 6 by 6 cells but no less than 5 m
ignition_radius(2) = 0.5*max(5.0,6*max(fdx,fdy))
ignition_time(2)=2

t0=tvals
time_start=t0
stepnum=stepnum+1
sstep=stepnum

do istep=sstep,msteps
    ifun_start=1
    if(istep.ne.1)ifun_start=3
    do ifun=ifun_start,6
!OMP    PARALLEL DO PRIVATE(ij,ifts,ifte,jfts,jfte)        
        do ij=1,num_tiles
            ifts= if_start(ij)          
            ifte= if_end(ij)
            jfts= jf_start(ij)
            jfte= jf_end(ij)

            call   sfire_model (10*istep+ij,ifun,         &
                need_lfn_update, num_ignitions,              & 
                ifuelread,nfuel_cat0,                   &
                ifds,ifde,jfds,jfde,                    & ! fire domain dims - the whole domain
                ifms,ifme,jfms,jfme,                    & ! fire memory dims - how declared
                ifts,ifte,jfts,jfte,                    & ! fire tile dims  - this thread
                time_start,dt,                          & ! time and increment
                fdx,fdy,                                & ! fire mesh spacing
                2,                                      & ! boundary_guard
                ignition_start_x,ignition_start_y,      &
                ignition_end_x,ignition_end_y,          &  
                ignition_radius,                        &
                ignition_time,                          &
                zsf,                                    & ! terrain height (for gradient)
                vx,vy,                                  & ! input: wind
                lfn,lfn_out,tign,fuel_frac,                     & ! state: level function, ign time, fuel left
                grnhfx,grnqfx,                          & ! output: heat fluxes
                nfuel_cat,                              & ! fuel data per point 
                fuel_time,                              & ! save derived internal data
                dzfsdx,dzfsdy,bbb,betafl,phiwc,r_0,fgip,ischap &
            )
tdim=tdim+1
start=(/1,1,tdim/)
atmp(1)=tvals+dt
call check( &
nf90_put_var(ncid,tvalsid,atmp,start=(/start(3)/),count=(/count(3)/)) &
           )
call check( nf90_put_var(ncid,lfnid,lfn _SD,start=start,count=counts) )
call check( nf90_put_var(ncid,tignid,tign _D,start=start,count=count) )
call check( nf90_put_var(ncid,vfxid,vx _SD,start=start,count=counts) )
call check( nf90_put_var(ncid,vfyid,vy _SD,start=start,count=counts) )
call check( nf90_put_var(ncid,grnhftid,grnhfx _D,start=start,count=count) )
call check( nf90_put_var(ncid,grnqftid,grnqfx _D,start=start,count=count) )
call check( nf90_put_var(ncid,fuel_fracid,fuel_frac _D,start=start,count=count) )
call check( nf90_put_var(ncid,tvalsid,atmp,start=(/tdim/),count=(/1/)) )
call check( nf90_put_var(ncid,fueltimeid,fuel_time _D,start=start,count=count) )
call check( nf90_put_att(ncid,nf90_global,"stepnum",istep) )
call check( nf90_sync(ncid) )
time_start=time_start+dt 
    
!close(1)
call check( nf90_close(ncid) )

        enddo 
        
    enddo !OMP PARALLEL
    
    if(istep.le.10.or.mod(istep,10).eq.0)then
        write(1,1)1.,1.,time_start
        write(1,1)sm,sn,((lfn(i,j),i=ifps,ifpe),j=jfps,jfpe)
        write(1,1)sm,sn,((tign(i,j),i=ifps,ifpe),j=jfps,jfpe)
        write(1,1)sm,sn,((vx(i,j),i=ifps,ifpe),j=jfps,jfpe)
        write(1,1)sm,sn,((vy(i,j),i=ifps,ifpe),j=jfps,jfpe)
        write(1,1)sm,sn,((grnhfx(i,j),i=ifps,ifpe),j=jfps,jfpe) 
    endif
    print *,'test_main: step ',istep,' of ',msteps,' time ',time_start
    time_start=time_start+dt 
enddo

close(1)

end subroutine model_test

!
!******************************
!

subroutine set_tiles(itiles,jtiles,ids,ide,jds,jde,num_tiles,i_start,i_end,j_start,j_end)
!*** set tiles for standalone/testing
implicit none
!*** arguments
integer,intent(in)::itiles,jtiles,ids,ide,jds,jde
integer,intent(out)::num_tiles
integer,intent(out),dimension(itiles*jtiles)::i_start,i_end,j_start,j_end
!*** local
integer::i,j,istep,jstep,ij
num_tiles=itiles*jtiles
istep=(ide-ids+itiles)/itiles
jstep=(jde-jds+jtiles)/jtiles
do i=1,itiles
    do j=1,jtiles
        ij=j+(i-1)*jtiles
        i_start(ij)=min(ide,ids+(i-1)*istep)
        i_end(ij)  =min(ide,ids+(i  )*istep-1)
        j_start(ij)=min(jde,jds+(j-1)*jstep)
        j_end(ij)  =min(jde,jds+(j  )*jstep-1)
    enddo
enddo
call check_tiles(ids,ide,jds,jde,num_tiles,i_start,i_end,j_start,j_end)
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
if (num_tiles.eq.1) then
    if(i_start(1).ne.ips.or.i_end(1).ne.ipe.or.j_start(1).ne.jps.or.j_end(1).ne.jpe)ie=1
else
    do ij=1,num_tiles
        if(i_start(ij).lt.ips.or.i_end(ij).gt.ipe &
        .or.j_start(ij).lt.jps.or.j_end(ij).gt.jpe)ie=ij
    enddo
endif
if(ie.ne.0)then        
    write(msg,*)'bad tile ',ie
    call message(msg)
    write(msg,*)'patch dimensions:',ips,ipe,jps,jpe
    call message(msg)
    do ij=1,num_tiles
        write(msg,*)'tile',ij,i_start(ij),i_end(ij),j_start(ij),j_end(ij)
        call message(msg)
    enddo
    call crash('bad tile bounds')
endif
end subroutine check_tiles


end module module_model_test

!
!******************************
!


program model_test_main
use module_model_test
call  main_sub
end program model_test_main
