! Test the dynamical core of the fire model
! Jan Mandel jmandel@ucar.edu August 2007

!   test driver for core, arguments communicated by files

module module_core
use module_fr_sfire_core
use module_fr_sfire_util
contains

subroutine sfire_core(id, &
    ids,ide,jds,jde,  &
    ims,ime,jms,jme,  &
    its,ite,jts,jte,  &
    ts,dt,dx,dy,fuel_time,   &
    tbound,           &
    lfn_in,lfn_out,tign,fuel_frac,fuel_frac_burnt &
#ifdef SPEED_VARS_ARGS      /* extra arguments for normal_spread */
#include SPEED_VARS_ARGS
#else
#include "fr_sfire_params_args.h"
#endif
                             ) 

!*** purpose
!
! Dynamical core of the fire spread model, insulated from the physics. 
! References the fire grid ONLY. Keep separated from the coupling
! with atm. model and data input.
!

implicit none

!*** arguments

integer, intent(in) :: id,ids,ide,jds,jde,ims,ime,jms,jme,its,ite,jts,jte
real, intent(inout), dimension(ims:ime,jms:jme)::lfn_in,tign
real, intent(out), dimension(ims:ime,jms:jme)::fuel_frac,lfn_out,fuel_frac_burnt
real, intent(in), dimension(ims:ime,jms:jme)::fuel_time
real, intent(in):: ts,dt,dx,dy
real, intent(out)::tbound
#ifdef SPEED_VARS_DECL      /* extra arguments for normal_spread */
#include SPEED_VARS_DECL
#else
#include "fr_sfire_params_decl.h"
#endif


! argument          intent  description                     (unit)  lives at
!
! id                in   unique number for debug prints and dumps
! ifun              in   2=one timestep, 3=compute fuel_frac
! ids,ide,jds,jde   in   mesh domain dimensions (cells)      (1)
! ims,ime,jms,jme   in   mesh aray dimensions                (1)
! its,ite,jts,jte   in   mesh tile dimensions                (1)
! ts                in   the starting time                   (s)
! dt                in   the timestep                        (s)
! dx,dy             in   fire mesh spacings                  (m)
! fuel_time         in   time fuel burns down to 1/e         (s)      cells
! lfn_in,lfn_out    inout  level set function (state)              (1)      nodes
! tign            inout  ignition time (state)               (s)      nodes
! fuel_frac       inout  the fuel fraction at the end        (1)      cells
! fuel_frac_burnt   out  the fuel fraction burnt this call   (1)      cells 

! A speed function is called for the spread rate at selected points.
! fuel_frac is currently ignored on input and computed from tign but it may not
! be in future.

!*** description

! This is a dynamical core of the fire spread model, insulated from the physics. 
! The physics should be done in the pre- and postprocessing, and in 
! the speed function for the fireline propagation (future).
!
! The state of the model is the level function lfn, which determines the fire 
! area, and the ignition time tign, both interpolated from values at nodes.
! The fire area is the level set where lfn <= 0.  The fireline is where lfn=0. 
! The array tign outside of the fire area is not set or referenced. 
! The state should be preserved between the calls, and it can be modified by
! data assimilation. All other quantities are derived from the state in each call.
!
! The level function evolves the fireline with the speed in the normal direction given by
! the spread rate r and the normal component of the wind. The level set method
! takes care of of various special cases automagically, such as ignition of a cell 
! surrounded by cells that all completely burning, and merging of approaching firelines.
!
! The fuel fraction is estimated from the the ignition times assuming 
! exponential decrease since ignition with decrease of fule fraction to 1/e in
! fuel_tim. The ignition times at nodes are interpolated linearly from the 
! evolving level function at the start and at the end, i.e. by assuming that the 
! value of the level function at a point varies linearly with time.
!
! ***NOTE: If a narrow band scheme is used to advance the level function in time
! then the firelines at time_start and time_end must fit within the band 
! at either time, and level function values away from the band should be set
! to some large positive and negative constants to assure that the ignition times
! in the area between the firelines are set reasonably. This is also important 
! when the level function is used for data assimilation. ***
!
! It is the responsibility of the caller to:
!
! before the call
!    - polulate all cells with the proper fuel_time coefficient
!    - interpolate and correct atmospheric winds
!
! after the call
!    - compute the fluxes from the fuel fraction burned
!      and sum up the the fluxes over atmopheric grid cells
!
!*** local
real, dimension(ims:ime,jms:jme)::fuel_frac_start
integer i,j



!*** calls
intrinsic epsilon

!*** executable

    ! compute the fuel fraction at ts 
    call get_fuel_left( &
        ims,ime,jms,jme, &
        its,ite,jts,jte, &
        dx,dy,                        &
        lfn_in,tign,fuel_time,ts,fuel_frac_start) ! 

    ! do one time step
    call prop_ls( id, &                                ! for debug
                ids,ide,jds,jde, &                       ! domain dims
                ims,ime,jms,jme, &                       ! memory dims
                its,ite,jts,jte, &                       ! tile dims
                ts,dt,dx,dy,     & ! scalars in
                tbound,          &                       ! scalars out
                lfn_in,lfn_out,tign         &                     ! arrays inout                   
#ifdef SPEED_VARS_ARGS      /* extra arguments for normal_spread */
#include SPEED_VARS_ARGS
#else
#include "fr_sfire_params_args.h"
#endif
             )

    ! compute the fuel fraction at time_end 
    call get_fuel_left( &
        ims,ime,jms,jme, &
        its,ite,jts,jte, &
        dx,dy,                        &
        lfn_in,tign,fuel_time,ts+dt,fuel_frac) ! replace by time_now later
        
    do j=jts,jte    
        do i=its,ite
            fuel_frac_burnt(i,j)=fuel_frac_start(i,j)-fuel_frac(i,j)
        enddo
    enddo
        

call print_2d_stats(its,ite,jts,jte,ims,ime,jms,jme,fuel_frac_burnt,'core: fuel frac burnt')

end subroutine sfire_core

end module module_core
!
!****************************************
!
program core_test
use module_core
use module_fr_sfire_core
use module_fr_sfire_util
 

implicit none
!*** variables
integer :: m1,m2,n1,n2,i,j,m1m,m2m,n1m,n2m,halo
real,allocatable,dimension(:,:):: &
    lfn,tign,fuel_time,frac_lost,frac_end,vx,vy,r
real:: rm1,rm2,rn1,rn2,time_start,time_end,fire_dx,fire_dy
integer::num_tiles
integer, dimension(100)::i_start,i_end,j_start,j_end
real::normal_spread_c,normal_spread_e
!*** executable

! read the inputs from a file
open(1,form='formatted',file='core_test_in.txt')
read(1,1)rm1,rm2,rn1,rn2,time_start,time_end,fire_dx,fire_dy
print *,rm1,rm2,rn1,rn2,time_start,time_end,fire_dx,fire_dy
m1=nint(rm1)  ! domain dimensions
n1=nint(rn1)
m2=nint(rm2)
n2=nint(rn2)
halo=0
m1m=m1-halo      ! memory dims
m2m=m2+halo
n1m=n1-halo
n2m=n2+halo
print *,'core_test:',m1,m2,n1,n2
allocate (lfn(m1m:m2m,n1m:n2m),tign(m1m:m2m,n1m:n2m), &
  vx(m1m:m2m,n1m:n2m),vy(m1m:m2m,n1m:n2m),            &
  fuel_time(m1m:m2m,n1m:n2m),r(m1m:m2m,n1m:n2m),  &
  frac_lost(m1m:m2m,n1m:n2m),frac_end(m1m:m2m,n1m:n2m))
read(1,1)((lfn(i,j),i=m1,m2),j=n1,n2)
read(1,1)((tign(i,j),i=m1,m2),j=n1,n2)
read(1,1)((fuel_time(i,j),i=m1,m2-1),j=n1,n2-1)
read(1,1)((r(i,j),i=m1,m2),j=n1,n2)
read(1,1)((vx(i,j),i=m1,m2),j=n1,n2)
read(1,1)((vy(i,j),i=m1,m2),j=n1,n2)
close(1)

call set_tiles(1,1,m1,m2-1,n1,n2-1,num_tiles,i_start,i_end,j_start,j_end)

! run the code

!normal_spread_c=0.185060861
!normal_spread_e=1.310758329
normal_spread_c=1
normal_spread_e=1
call sfire_core( &
m1,m2-1,n1,n2-1, &
m1m,m2m,n1m,n2m,          &
m1,m2-1,n1,n2-1,  &
num_tiles,i_start,i_end,j_start,j_end, &
time_start,time_end-time_start,fire_dx,fire_dy,fuel_time,           & 
lfn,tign,frac_end,frac_lost &
#ifdef SPEED_VARS_ARGS      /* extra arguments for normal_spread */
#include SPEED_VARS_ARGS
#endif
)

! write the output to a file

open (1,form='formatted',file='core_test_out.txt')
write(1,1)((lfn(i,j),i=m1,m2),j=n1,n2)
write(1,1)((tign(i,j),i=m1,m2),j=n1,n2)
write(1,1)((frac_lost(i,j),i=m1,m2-1),j=n1,n2-1)
write(1,1)((frac_end(i,j),i=m1,m2-1),j=n1,n2-1)
close(1)

1   format((e25.15,1x)) 

end program core_test
