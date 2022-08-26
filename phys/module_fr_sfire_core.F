!
!*** Jan Mandel August-October 2007 email: jmandel@ucar.edu or Jan.Mandel@gmail.com
!
! With contributions by Minjeong Kim.
#define DEBUG_OUT
#define DEBUG_PRINT
!#define DEBUG_OUT_FUEL_LEFT    

module module_fr_sfire_core

use module_fr_sfire_phys, only: fire_params , fire_ros
use module_fr_sfire_util

implicit none

! The mathematical core of the fire spread model. No physical constants here.
! 
! subroutine sfire_core: only this routine should be called from the outside.
! subroutine fuel_left:  compute remaining fuel from time of ignition.
! subroutine prop_ls: propagation of curve in normal direction.


contains

!
!****************************************
!
    
subroutine init_no_fire(&
    ifds,ifde,jfds,jfde, &
    ifms,ifme,jfms,jfme, &
    ifts,ifte,jfts,jfte, &
    fdx,fdy,time_start,dt,  & ! scalars in
    fuel_frac,fire_area,lfn,tign_in, & !arrays out 
    tign)    ! arrays inout            
implicit none
             
!*** purpose: initialize model to no fire

!*** arguments
integer, intent(in):: ifds,ifde,jfds,jfde   ! fire domain bounds
integer, intent(in):: ifts,ifte,jfts,jfte   ! fire tile bounds
integer, intent(in):: ifms,ifme,jfms,jfme   ! array bounds
real, intent(in) :: fdx,fdy,time_start,dt   ! mesh spacing, time
real, intent(out), dimension (ifms:ifme,jfms:jfme) :: & 
                   fuel_frac,fire_area,lfn,           &  !  model state
                   tign_in       
real, intent(inout), dimension (ifms:ifme,jfms:jfme) :: & 
                   tign       ! model state, on input values tign_in

!*** calls
intrinsic epsilon
                                                
!*** local
integer:: i,j
real:: lfn_init,time_init,time_now
character(len=128):: msg

time_now = time_start+dt
time_init = time_start+2*dt  ! could be time_start+dt if not for rounding errors
lfn_init = 2*max((ifde-ifds+1)*fdx,(jfde-jfds+1)*fdy)      ! more than domain diameter

call print_2d_stats(ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme, tign,'init_no_fire start: tign')

if (fire_perimeter_time > 0. .and. fire_tign_in_time > 0.)then
   call crash('fire_perimeter_time > 0 and fire_tign_in_time > 0')
endif

if(fire_perimeter_time > 0.)then
  
  write(msg,'(a,g14.6)')'init_no_fire: using given ignition time to replay history to fire_perimeter_time=',fire_perimeter_time
!$OMP CRITICAL(SFIRE_CORE_CRIT)
  call message(msg,level=1)
!$OMP END CRITICAL(SFIRE_CORE_CRIT)
  do j=jfts,jfte
     do i=ifts,ifte
        fuel_frac(i,j)=1.          ! fuel at start is 1 by definition
        fire_area(i,j)=0.          ! nothing burning
        lfn(i,j) = tign(i,j)-time_now       ! use specified ignition time as level set function 
     enddo
  enddo
elseif(fire_tign_in_time > 0.)then
  call message('init_no_fire: using ignition from given max fire arrival time',level=1)
  do j=jfts,jfte
     do i=ifts,ifte
        tign_in(i,j)=tign(i,j)
        fuel_frac(i,j)=1.          ! fuel at start is 1 by definition
        fire_area(i,j)=0.          ! nothing burning
        lfn(i,j) = tign_in(i,j)-time_now    ! use specified ignition time as level set function 
     enddo
  enddo
  call print_2d_stats(ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme, tign_in,'init_no_fire: tign_in')
else
  call message('init_no_fire: setting state to no fire',level=1)
  do j=jfts,jfte
     do i=ifts,ifte
        fuel_frac(i,j)=1.          ! fuel at start is 1 by definition
        fire_area(i,j)=0.          ! nothing burning
        tign(i,j) = time_init      ! ignition in future
        lfn(i,j) = lfn_init        ! no fire 
     enddo
  enddo
endif

call check_lfn_tign('init_no_fire',time_now,ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme,lfn, tign)
call print_2d_stats(ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme, tign,'init_no_fire: tign')
call print_2d_stats(ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme, tign,'init_no_fire: lfn')

do j=jfts,jfte
   do i=ifts,ifte
      if(.not.lfn(i,j) > 0.)then
         write(msg,*)'i,j=',i,j,' tign=',tign(i,j),' <= time_now =',time_now
         call message(msg,level=0)
         call crash('init_no_fire: ignition time must be after the end of the first time step')
      endif
   enddo
enddo

end subroutine init_no_fire

!
!******************
!

subroutine ignite_from_tign_in(                                   &
                        ifds,ifde,jfds,jfde,                      & ! fire domain dims - the whole domain
                        ifms,ifme,jfms,jfme,                      &
                        ifts,ifte,jfts,jfte,                      &
                        start_ts,end_ts,                          &
                        tign_in,                                  &
                        lfn,tign,ignited)
implicit none
!*** purpose: ignite from given gridded ignition times
!*** arguments
integer, intent(in):: ifds,ifde,jfds,jfde   ! fire domain bounds
integer, intent(in):: ifts,ifte,jfts,jfte   ! fire tile bounds
integer, intent(in):: ifms,ifme,jfms,jfme   ! array bounds
real, intent(in):: start_ts,end_ts          ! the time step start and end 
real, intent(in), dimension (ifms:ifme,jfms:jfme) :: tign_in     ! the given ignition time 
real, intent(inout), dimension (ifms:ifme,jfms:jfme) :: & 
                   lfn, tign                ! level function, ignition time (state)
integer, intent(out)::ignited
!*** global
! fire_tign_in_time 
!*** local
integer:: i,j
real:: lfn_old

!*** executable
ignited=0
if(.not. start_ts < fire_tign_in_time)return   ! too late, nothing to do

call check_lfn_tign('ignite_from_tign_in start',end_ts,ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme,lfn, tign)

do j=jfts,jfte   
    do i=ifts,ifte
       if (.not. tign_in(i,j)>end_ts )then
          lfn_old = lfn(i,j) 
          tign(i,j)=min(tign(i,j), tign_in(i,j))
          lfn(i,j)=min(lfn(i,j),tign(i,j)-end_ts)
          call check_lfn_tign_ij(i,j,'ignite_from_tign_in',end_ts,lfn(i,j),tign(i,j))
          if(lfn_old > 0 .and. .not. lfn(i,j) > 0.) ignited = ignited + 1
       endif
    enddo
enddo 


call check_lfn_tign('ignite_from_tign_in end',end_ts,ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme,lfn, tign)
call print_2d_stats(ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme, tign,'init_no_fire: tign')
call print_2d_stats(ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme, tign,'init_no_fire: lfn')

end subroutine ignite_from_tign_in

 

subroutine ignite_fire( ifds,ifde,jfds,jfde,                    & ! fire domain dims - the whole domain
                        ifms,ifme,jfms,jfme,                      &
                        ifts,ifte,jfts,jfte,                      &
                        ignition_line,                            &
                        start_ts,end_ts,                    &
                        coord_xf,coord_yf,                &     
                        unit_xf,unit_yf,                  &
                        lfn,tign,ignited)
implicit none

!*** purpose: ignite a circular/line fire 

!*** description
! ignite fire in the region within radius r from the line (sx,sy) to (ex,ey).
! the coordinates of nodes are given as the arrays coord_xf and coord_yf
! r is given in m
! one unit of coord_xf is unit_xf m 
! one unit of coord_yf is unit_yf m 
! so a node (i,j) will be ignited iff for some (x,y) on the line
! || ( (coord_xf(i,j) - x)*unit_xf , (coord_yf(i,j) - y)*unit_yf ) || <= r 


!*** arguments
integer, intent(in):: ifds,ifde,jfds,jfde   ! fire domain bounds
integer, intent(in):: ifts,ifte,jfts,jfte   ! fire tile bounds
integer, intent(in):: ifms,ifme,jfms,jfme   ! array bounds
type(line_type), intent(in):: ignition_line    ! description of the ignition line
real, intent(in):: start_ts,end_ts          ! the time step start and end 
real, dimension(ifms:ifme, jfms:jfme), intent(in):: & 
    coord_xf,coord_yf                       !  node coordinates  
real, intent(in):: unit_xf,unit_yf          !  coordinate units in m
real, intent(inout), dimension (ifms:ifme,jfms:jfme) :: & 
                   lfn, tign                ! level function, ignition time (state)
integer, intent(out):: ignited              ! number of nodes newly ignited
                        
!*** local
integer:: i,j
real::lfn_new,tign_new,time_ign,ax,ay,rels,rele,d
! real:: lfn_new_chk,tign_new_chk,tmperr
real:: sx,sy                    ! start of ignition line, from lower left corner
real:: ex,ey                    ! end of ignition line, or zero
real:: st,et                    ! start and end of time of the ignition line
character(len=128):: msg
real::cx2,cy2,dmax,axmin,axmax,aymin,aymax,dmin
real:: start_x,start_y          ! start of ignition line, from lower left corner
real:: end_x,end_y              ! end of ignition line, or zero
real:: radius                   ! all within the radius of the line will ignite
real:: start_time,end_time      ! the ignition time for the start and the end of the line
real:: ros,tos                  ! ignition rate and time of spread
integer:: msglevel=4,smsg=2     ! print at this level
real:: lfn_min

!*** executable

call check_lfn_tign('ignite_fire start',end_ts,ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme,lfn, tign)

! copy ignition line fields to local variables
start_x    = ignition_line%start_x ! x coordinate of the ignition line start point (m, or long/lat)
start_y    = ignition_line%start_y ! y coordinate of the ignition line start point
end_x      = ignition_line%end_x   ! x coordinate of the ignition line end point
end_y      = ignition_line%end_y   ! y coordinate of the ignition line end point
start_time = ignition_line%start_time ! ignition time for the start point from simulation start (s)
end_time   = ignition_line%end_time! ignition time for the end poin from simulation start (s)
radius     = ignition_line%radius  ! all within this radius ignites immediately
ros        = ignition_line%ros     ! rate of spread


tos        = radius/ros            ! time of spread to the given radius
st         = start_time            ! the start time of ignition considered in this time step
et         = min(end_ts,end_time)  ! the end time of the ignition segment in this time step

ignited = 0
! this should be called whenever (start_ts, end_ts) \subset (start_time, end_time + tos) 
if(start_ts>et+tos .or. end_ts<st)return   ! too late or too early, nothing to do

if(start_time < end_time)then  ! we really want to test start_time .ne. end_time, but avoiding test of floats on equality
        ! segment of nonzero length
        !rels =  (st - start_time) / (end_time - start_time)  ! relative position of st in the segment (start,end)
	!sx = start_x + rels * (end_x - start_x)
	!sy = start_y + rels * (end_y - start_y)
        rels = 0.
        sx = start_x
        sy = start_y
        rele =  (et - start_time) / (end_time - start_time)    ! relative position of et in the segment (start,end)
	ex = start_x + rele * (end_x - start_x)
	ey = start_y + rele * (end_y - start_y)
else
        ! just a point
        rels = 0.
        rele = 1.
	sx = start_x
	sy = start_y
	ex = end_x
	ey = end_y
endif


cx2=unit_xf*unit_xf
cy2=unit_yf*unit_yf

axmin=coord_xf(ifts,jfts)
aymin=coord_yf(ifts,jfts)
axmax=coord_xf(ifte,jfte)
aymax=coord_yf(ifte,jfte)
if(fire_print_msg.ge.smsg)then
!$OMP CRITICAL(SFIRE_CORE_CRIT)
write(msg,'(a,2f11.6,a,2f11.6)')'IGN from ',sx,sy,' to ',ex,ey
call message(msg,level=smsg)
write(msg,'(a,2f10.2,a,2f10.2,a)')'IGN timestep [',start_ts,end_ts,'] in [',start_time,end_time,']'
call message(msg,level=smsg)
write(msg,'(a,2g13.6,a,2g13.6)')'IGN tile coord from  ',axmin,aymin,' to ',axmax,aymax
call message(msg,level=smsg)
!$OMP END CRITICAL(SFIRE_CORE_CRIT)
endif
ignited=0
dmax=0
dmin=huge(dmax)
11      format('IGN ',6(a,g17.7,1x)) 
12      format('IGN ',4(a,2g13.7,1x))
! tmperr=0.
do j=jfts,jfte   
    do i=ifts,ifte
        call check_lfn_tign_ij(i,j,'ignite_fire start',end_ts,lfn(i,j),tign(i,j))
        ax=coord_xf(i,j)
        ay=coord_yf(i,j)

        ! get d= distance from the nearest point on the ignition segment
        ! and time_ign = the ignition time there
        call nearest(d,time_ign,ax,ay,sx,sy,st,ex,ey,et,cx2,cy2)


        ! collect statistics
        dmax=max(d,dmax)
        dmin=min(d,dmin)

        ! tentative new level set function and time of ignition
        ! lfn_new_chk=d - min( radius, ros*(end_ts - time_ign) )  ! lft at end_ts
        ! tign_new_chk = time_ign + d/ros
        if(radius < ros*(end_ts - time_ign))then
            lfn_new = d - radius
            tign_new= time_ign + d/ros
        else
            lfn_new=d-ros*(end_ts - time_ign)
            tign_new=lfn_new/ros + end_ts  ! consistent even with rounding errors
        endif
        ! tmperr=max(tmperr,abs(lfn_new - lfn_new_chk)+abs(tign_new-tign_new_chk))
        !lfn_new=lfn_new_chk
        !tign_new=tign_new_chk
        
        if(fire_print_msg.ge.msglevel)then
!$OMP CRITICAL(SFIRE_CORE_CRIT)
            write(msg,*)'IGN1 i,j=',i,j,' lfn(i,j)=',lfn(i,j),' tign(i,j)=',tign(i,j)
            call message(msg,level=0)
            write(msg,*)'IGN2 i,j=',i,j,' lfn_new= ',lfn_new, ' time_ign= ',time_ign,' d=',d
            call message(msg,level=msglevel)
!$OMP END CRITICAL(SFIRE_CORE_CRIT)
        endif
        if(.not.lfn_new>0.) then
            ignited=ignited+1   ! count
        endif
        if(.not. lfn(i,j)<0. .and. lfn_new < 0.) then
            tign(i,j)=tign_new  ! newly ignited now
            if(fire_print_msg.ge.msglevel)then
!$OMP CRITICAL(SFIRE_CORE_CRIT)
                write(msg,'(a,2i6,a,2g13.6,a,f10.2,a,2f10.2,a)')'IGN ignited cell ',i,j,' at',ax,ay, &
                    ' time',tign(i,j),' in [',start_ts,end_ts,']'
                call message(msg,level=0)
                write(msg,'(a,g10.3,a,f10.2,a,2f10.2,a)')'IGN distance',d,' from ignition line at',time_ign,' in [',st,et,']'
                call message(msg,level=msglevel)
!$OMP END CRITICAL(SFIRE_CORE_CRIT)
            endif
            if(tign(i,j) < start_ts .or. tign(i,j) > end_ts )then
!$OMP CRITICAL(SFIRE_CORE_CRIT)
                write(msg,'(a,2i6,a,f11.6,a,2f11.6,a)')'WARNING ',i,j, &
                ' fixing ignition time ',tign(i,j),' outside of the time step [',start_ts,end_ts,']'
                call message (msg)
!$OMP END CRITICAL(SFIRE_CORE_CRIT)
                tign(i,j) = min(max(tign(i,j),start_ts),end_ts)
            endif
        endif
        lfn(i,j)=min(lfn(i,j),lfn_new)  ! update the level set function
        if(fire_print_msg.ge.msglevel)then
!$OMP CRITICAL(SFIRE_CORE_CRIT)
            write(msg,*)'IGN3 i,j=',i,j,' lfn(i,j)=',lfn(i,j),' tign(i,j)=',tign(i,j)
            call message(msg,level=msglevel)
!$OMP END CRITICAL(SFIRE_CORE_CRIT)
        endif
        call check_lfn_tign_ij(i,j,'ignite_fire end',end_ts,lfn(i,j),tign(i,j))
    enddo
enddo
! write(msg,*)'ignite_fire: max error ',tmperr
! call message(msg)
call check_lfn_tign("ignite_fire end:",end_ts,ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme,lfn, tign)
if(fire_print_msg .ge. smsg)then
lfn_min=huge(lfn_min)
do j=jfts,jfte   
    do i=ifts,ifte
        lfn_min=min(lfn_min,lfn(i,j))
    enddo
enddo
!$OMP CRITICAL(SFIRE_CORE_CRIT)
write(msg,'(a,2g13.2,a,g10.2,a,g10.2)')'IGN units ',unit_xf,unit_yf,' m max dist ',dmax,' min',dmin
call message(msg,level=smsg)
write(msg,'(a,f6.1,a,f8.1,a,i10,a,g10.2)')'IGN radius ',radius,' time of spread',tos, &
   ' ignited nodes',ignited,' lfn min',lfn_min
call message(msg,level=smsg)
!$OMP END CRITICAL(SFIRE_CORE_CRIT)
endif
end subroutine ignite_fire

!
!***
!

subroutine check_lfn_tign(s,time_now,ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme,lfn, tign)
!*** purpose: check consistency of lfn and ignition
implicit none
!*** arguments
character(len=*),intent(in)::s              ! for print
real, intent(in)::time_now                    ! end of timestep = time_now
integer, intent(in):: ifts,ifte,jfts,jfte   ! fire tile bounds
integer, intent(in):: ifms,ifme,jfms,jfme   ! array bounds
real, intent(in), dimension (ifms:ifme,jfms:jfme) :: & 
                   lfn, tign                ! level function, ignition time (state)
!*** local
integer:: i,j
character(len=128)::msg
!***    executable

do j=jfts,jfte   
    do i=ifts,ifte
       call check_lfn_tign_ij(i,j,s,time_now,lfn(i,j),tign(i,j))
    enddo
enddo

! !$OMP CRITICAL(SFIRE_CORE_CRIT)
! 1 format('check_lfn_tign OK',f9.2,1x,a,1x,i5,':',i5,i5,':',i5)
! write(msg,1)time_now,trim(s),ifts,ifte,jfts,jfte
! !$OMP END CRITICAL(SFIRE_CORE_CRIT)
! call message(msg)

end subroutine check_lfn_tign 

subroutine check_lfn_tign_ij(i,j,s,time_now,lfnij,tignij)
!*** purpose: check consistency of lfn and ignition
implicit none
!*** arguments
integer, intent(in)::i,j                    ! indices
character(len=*),intent(in)::s              ! for print
real, intent(in)::time_now                  ! end of timestep = time_now
real, intent(in):: lfnij, tignij                ! level function, ignition time (state)
!*** local
character(len=128):: msg,msg1
!***    executable

if(.not. lfnij<0. .and. tignij<time_now)then
!$OMP CRITICAL(SFIRE_CORE_CRIT)
1  format(a,i5,',',i5,a,g13.6,a,g15.8,a,g15.8)
   msg1 = trim(s)//':check_lfn_tign_ij: inconsistent state'
   if (lfnij > 0.)then
       write(msg,1)'i,j=',i,j,' lfn=',lfnij,'>0 tign=',tignij,'< time_now=',time_now
       call message(msg,level=0)
       call crash(msg1) ! instead of just warning in check_lfn_tign_ij
   else
       write(msg,1)'i,j=',i,j,' lfn=',lfnij,'>=0 tign=',tignij,'< time_now=',time_now
       call message(msg,level=0)
       call warning(msg1,level=0)
   endif
!$OMP END CRITICAL(SFIRE_CORE_CRIT)
endif

end subroutine check_lfn_tign_ij


! called from the inside of a loop, inline if possible
!DEC$ ATTRIBUTES FORCEINLINE
SUBROUTINE nearest(d,t,ax,ay,sx,sy,st,ex,ey,et,cx2,cy2)
        implicit none
!***    arguments
        real, intent(out):: d,t
        real, intent(in):: ax,ay,sx,sy,st,ex,ey,et,cx2,cy2
        ! input:
        ! ax, ay       coordinates of point a
        ! sx,sy,ex,ey  coordinates of endpoints of segment [x,y]
        ! st,et        values at the endpoints x,y
        ! cx2,cy2      x and y unit squared for computing distance 

        ! output
        ! d            the distance of a and the nearest point z on the segment [x,y]
        ! t            linear interpolation from the values st,et to the point z
        !
        ! method: compute d as the distance (ax,ay) from the midpoint (mx,my)
        ! minus a correction (because of rounding errors): |a-c|^2 = |a-m|^2 - |m-c|^2
        ! when |m-c| >= |s-e|/2 the nearest point is one of the endpoints
        ! the computation work also for the case when s=e exactly or approximately 
        !
        !
        !           a    
        !          /| \
        !     s---m-c--e
        !
        ! |m-c| = |a-m| cos (a-m,e-s) 
        !       = |a-m| (a-m).(e-s))/(|a-m|*|e-s|)
!***    local
        real:: mx,my,dam2,dames,am_es,cos2,dmc2,mcrel,mid_t,dif_t,des2,cx,cy
        character(len=128):: msg
        integer::msglevel=4
!***    executable

11      format('IGN ',6(a,g17.7,1x))
12      format('IGN ',4(a,2g13.7,1x))

        ! midpoint m = (mx,my)
        mx = (sx + ex)*0.5
        my = (sy + ey)*0.5
        dam2=(ax-mx)*(ax-mx)*cx2+(ay-my)*(ay-my)*cy2      ! |a-m|^2
        des2 = (ex-sx)*(ex-sx)*cx2+(ey-sy)*(ey-sy)*cy2          ! des2 = |e-s|^2
        dames = dam2*des2
        am_es=(ax-mx)*(ex-sx)*cx2+(ay-my)*(ey-sy)*cy2       ! am_es = (a-m).(e-s)
        if(dames>0)then
            cos2 = (am_es*am_es)/dames                  ! cos2 = cos^2 (a-m,e-s)
        else ! point a already is the midpoint
            cos2 = 0.
        endif
        dmc2 = dam2*cos2                                ! dmc2 = |m-c|^2
        if(4.*dmc2 < des2)then                          ! if |m-c|<=|e-s|/2
            ! d = sqrt(max(dam2 - dmc2,0.))               ! d=|a-m|^2 - |m-c|^2, guard rounding
            mcrel = sign(sqrt(4.*dmc2/des2),am_es)      ! relative distance of c from m 
        elseif(am_es>0)then                             ! if cos > 0, closest is e
            mcrel = 1.0 
        else                                            ! closest is s
            mcrel = -1.0
        endif
	cx = (ex + sx)*0.5 + mcrel*(ex - sx)*0.5     ! interpolate to c by going from m
	cy = (ey + sy)*0.5 + mcrel*(ey - sy)*0.5     ! interpolate to c by going from m
        d=sqrt((ax-cx)*(ax-cx)*cx2+(ay-cy)*(ay-cy)*cy2) ! |a-c|^2
	t = (et + st)*0.5 + mcrel*(et - st)*0.5     ! interpolate to c by going from m
        if(fire_print_msg.ge.msglevel)then
!$OMP CRITICAL(SFIRE_CORE_CRIT)
            write(msg,12)'find nearest to [',ax,ay,'] from [',sx,sy,'] [',ex,ey,']' ! DEB
            call message(msg,level=msglevel)
            write(msg,12)'end times',st,et,' scale squared',cx2,cy2 ! DEB
            call message(msg,level=msglevel)
            write(msg,11)'nearest at mcrel=',mcrel,'from the midpoint, t=',t ! DEB
            call message(msg,level=msglevel)
            write(msg,12)'nearest is [',cx,cy,'] d=',d ! DEB
            call message(msg,level=msglevel)
            write(msg,11)'dam2=',dam2,'des2=',des2,'dames=',dames
            call message(msg,level=msglevel)
            write(msg,11)'am_es=',am_es,'cos2=',cos2,'dmc2=',dmc2 ! DEB
            call message(msg)
!$OMP END CRITICAL(SFIRE_CORE_CRIT)
        endif
END SUBROUTINE nearest


!
!**********************
!            

subroutine fuel_left( &
    ifds,ifde,jfds,jfde, &
    ims,ime,jms,jme, &
    its,ite,jts,jte, &
    ifs,ife,jfs,jfe, &
    lfn, tign, fuel_time, time_now, fuel_frac, fire_area)
implicit none

!*** purpose: determine fraction of fuel remaining
!*** NOTE: because variables are cell centered, need halo/sync width 1 before

!*** Jan Mandel August 2007 email: jmandel@ucar.edu or Jan.Mandel@gmail.com

!*** arguments

integer, intent(in) ::ifds,ifde,jfds,jfde,its,ite,jts,jte,ims,ime &
                      ,jms,jme,ifs,ife,jfs,jfe
real, intent(in), dimension(ims:ime,jms:jme)::lfn,tign,fuel_time
real, intent(in):: time_now
real, intent(out), dimension(ifs:ife,jfs:jfe)::fuel_frac
real, intent(out), dimension(ims:ime,jms:jme):: fire_area

! ims,ime,jms,jme   in   memory dimensions
! its,ite,jts,jte   in   tile dimensions (cells where fuel_frac computed)
! ifs,ife,jfs,jfe   in   fuel_frac memory dimensions
! lfn               in   level function, at nodes at midpoints of cells
! tign              in   ignition time, at nodes at nodes at midpoints of cells
! fuel_time         in   time constant of fuel, per cell
! time_now          in   time now
! fuel_frac         out  fraction of fuel remaining, per cell
! fire_area         out  fraction of cell area on fire

!*** local

integer::i,j,ir,jr,icl,jcl,isubcl,jsubcl,i2,j2,ii,jj,its1,jts1,ite1,jte1
real::fmax,frat,helpsum1,helpsum2,fuel_left_ff,fire_area_ff,rx,ry,tignf(2,2)
real,dimension(3,3)::tff,lff
! help variables instead of arrays fuel_left_f and fire_area_f 
real::lffij,lffi1j,lffij1,lffi1j1,tifij,tifi1j,tifij1,tifi1j1,tx,ty,txx,tyy
! variables for calculation instead of lff(i,j) and tif(i,j)is lffij,tifij etc..#define IFCELLS (ite-its+1)*fuel_left_irl
#define JFCELLS (jte-jts+1)*fuel_left_jrl
character(len=128)::msg
integer::m,omp_get_thread_num
     

call check_mesh_2dim(its-1,ite+1,jts-1,jte+1,ims,ime,jms,jme)
call check_mesh_2dim(its,ite,jts,jte,ifs,ife,jfs,jfe)
call check_lfn_tign('fuel_left start',time_now,its,ite,jts,jte,ims,ime,jms,jme,lfn, tign)

! refinement
ir=fuel_left_irl
jr=fuel_left_jrl

if ((ir.ne.2).or.(jr.ne.2)) then 
   call crash('fuel_left: ir.ne.2 or jr.ne.2 ')
endif

rx=1./ir 
ry=1./jr

! interpolate level set function to finer grid
#ifdef DEBUG_OUT_FUEL_LEFT    
    call write_array_m(1,IFCELLS+1,1,JFCELLS+1,1,IFCELLS+1,1,JFCELLS+1,lff,'lff',0)
    call write_array_m(1,IFCELLS+1,1,JFCELLS+1,1,IFCELLS+1,1,JFCELLS+1,tif,'tif',0)
#endif

!
! example for ir=2:
!
!                     |      coarse cell      |
!      its-1                     its                                   ite             ite+1
! -------X------------|-----.-----X-----.-----|--........----|----------X----------|---------X
!           fine node 1           2           3                   2*(ite-its+1) 
!                fine cell  1           2                             cell 2*(ite-its+1)



!  Loop over cells in Tile 
!  Changes made by Volodymyr Kondratenko 09/24/2009
its1=max(its,ifds+1)
ite1=min(ite,ifde-1)
jts1=max(jts,jfds+1)
jte1=min(jte,jfde-1)

! jm: initialize fuel_frac - we do not compute it near the domain boundary becau!se we cannot interpolate!
do j=jts,jte
   do i=its,ite
      fuel_frac(i,j)=1.
      fire_area(i,j)=0.
   enddo
enddo

call check_lfn_tign('fuel_left',time_now,its,ite,jts,jte,ims,ime,jms,jme,lfn, tign)

do icl=its1,ite1
  do jcl=jts1,jte1
    helpsum1=0
    helpsum2=0

    call tign_lfn_interpolation(time_now,icl,jcl,ims,ime,jms,jme, &
                                  tign,lfn,tff,lff)

    do isubcl=1,ir
      do jsubcl=1,jr 
        if(fuel_left_method.eq.1)then
          call fuel_left_cell_1( fuel_left_ff, fire_area_ff, &
     lff(isubcl,jsubcl),lff(isubcl,jsubcl+1),lff(isubcl+1,jsubcl),lff(isubcl+1,jsubcl+1), &
     tff(isubcl,jsubcl),tff(isubcl,jsubcl+1),tff(isubcl+1,jsubcl),tff(isubcl+1,jsubcl+1), &
     time_now, fuel_time(icl,jcl))
        elseif(fuel_left_method.eq.2)then
          call fuel_left_cell_2( fuel_left_ff, fire_area_ff, &
     lff(isubcl,jsubcl),lff(isubcl,jsubcl+1),lff(isubcl+1,jsubcl),lff(isubcl+1,jsubcl+1), &
     tff(isubcl,jsubcl),tff(isubcl,jsubcl+1),tff(isubcl+1,jsubcl),tff(isubcl+1,jsubcl+1), &
     time_now, fuel_time(icl,jcl)) 
! dont forget to change fire_area_ff here
        else
          call crash('fuel_left: unknown fuel_left_method')
        endif

        ! consistency check
        if(fire_area_ff.lt.-1e-6 .or.  &
          (fire_area_ff.eq.0. .and. fuel_left_ff.lt.1.-1e-6))then
!$OMP CRITICAL(SFIRE_CORE_CRIT)
           write(msg,'(a,2i6,2(a,f11.8))')'fuel_left: at node',i,j, &
              ' of refined mesh fuel burnt',1-fuel_left_ff,' fire area',fire_area_ff
!$OMP END CRITICAL(SFIRE_CORE_CRIT)
           call crash(msg)
        endif

        helpsum1=helpsum1+fuel_left_ff
        helpsum2=helpsum2+fire_area_ff
        ! write(*,*)icl,jcl,isubcl,jsubcl,fuel_left_ff,fire_area_ff
      enddo
    enddo
    ! write(*,*)icl,jcl,helpsum1,helpsum2
    fuel_frac(icl,jcl)=helpsum1 / (ir*jr) ! multiply by weight for averaging over coarse cell
    fire_area(icl,jcl)=helpsum2 / (ir*jr) ! multiply by weight for averaging over coarse cell
  enddo 
enddo



#ifdef DEBUG_OUT_FUEL_LEFT
    call write_array_m(its,ite,jts,jte,ims,ime,jms,jme,fire_area,'fire_area',0)
    call write_array_m(its,ite,jts,jte,ims,ime,jms,jme,fuel_frac,'fuel_frac',0)
    call write_array_m(1,IFCELLS,1,JFCELLS,1,IFCELLS,1,JFCELLS,fuel_left_f,'fuel_left_f',0)
    call write_array_m(1,IFCELLS,1,JFCELLS,1,IFCELLS,1,JFCELLS,fire_area_f,'fire_area_f',0)
#endif

! consistency check after sum
!fmax=0
!do j=jts,jte
!    do i=its,ite        
!       if(fire_area(i,j).eq.0.)then
!           if(fuel_frac(i,j).lt.1.-1e-6)then
!!$OMP CRITICAL(SFIRE_CORE_CRIT)
!               write(msg,'(a,2i6,2(a,f11.8))')'fuel_left: at node',i,j, &
!                   ' fuel burnt',1-fuel_frac(i,j),' but fire area',fire_area(i,j)
!!$OMP END CRITICAL(SFIRE_CORE_CRIT)
!               call crash(msg)
!           endif
!       else
!           frat=(1-fuel_frac(i,j))/fire_area(i,j)
!           fmax=max(fmax,frat)
!       endif
!    enddo
!enddo
!$OMP CRITICAL(SFIRE_CORE_CRIT)
!write(msg,'(a,4i6,a,e15.7)')'fuel_left: tile',its,ite,jts,jte,' max fuel burnt/area',fmax 
!$OMP END CRITICAL(SFIRE_CORE_CRIT)
!call message(msg)

return


end subroutine fuel_left

!
!************************
!


!
!*************************
!Subroutine that is calculating tign and lfn of four endpoints of the subcell
! that is located at isubcl,jsubcl of the cell -(icl,jcl)
!
subroutine tign_lfn_interpolation(time_now,icl,jcl,ims,ime,jms,jme, &
                                  tign,lfn,tff,lff)
real, intent(in):: time_now    ! not ignited nodes will have tign set to >= time_now
integer, intent(in) :: icl,jcl
integer, intent(in) :: ims,ime,jms,jme  ! memory dimensions of all arrays
real, intent(in), dimension(ims:ime,jms:jme)::lfn,tign
real, intent(out),dimension(3,3)::tff,lff

!     |                 |                 |       
!  -(3,1)-------------(3,2)-------------(3,3)
!     |                 |                 |
!     |   (2,1)         |      (2,2)      |       
!     |                 |                 |
!     |                 |                 |
!     |                 |                 |
!     |                 |                 |
!   (2,1)--------node-(icl,jcl)---------(2,3)-----------(icl,jcl+1)-------------|
!     |          sub-node (2,2)           |
!     |                 |                 |
!     |   (1,1)         |      (1,2)      |       each fire mesh cell is decomposed in 4
!     |                 |                 |       tff,lff is computed at the nodes of 
!     |                 |                 |       the subcells, numbered (1,1)...(3,3)
!   (1,1)-------------(1,2)-------------(1,3)--     
!     |                 |                 |       
!

!**********************

      call check_lfn_tign('tign_lfn_interpolation',time_now,icl-1,icl,jcl-1,jcl,ims,ime,jms,jme,lfn, tign)
	
!       Direct calculation tif and lff, avoiding big auxiliary arrays, just for case ir=jr=2
! Checking whether icl or jcl is on the boundary

      call tign_lfn_four_pnts_interp(tign(icl-1,jcl-1),tign(icl-1,jcl),tign(icl,jcl-1),  & 
                                     tign(icl,jcl),lfn(icl-1,jcl-1),lfn(icl-1,jcl),      &
                                     lfn(icl,jcl-1),lfn(icl,jcl),lff(1,1),tff(1,1),time_now)

      call tign_lfn_line_interp(tign(icl-1,jcl),tign(icl,jcl),lfn(icl-1,jcl),lfn(icl,jcl), &
                             lff(1,2),tff(1,2),time_now,icl-1,jcl,icl,jcl)


      call tign_lfn_four_pnts_interp(tign(icl-1,jcl),tign(icl-1,jcl+1),tign(icl,jcl),  & 
                                     tign(icl,jcl+1),lfn(icl-1,jcl),lfn(icl-1,jcl+1),      &
                                     lfn(icl,jcl),lfn(icl,jcl+1),lff(1,3),tff(1,3),time_now)

      call tign_lfn_line_interp(tign(icl,jcl-1),tign(icl,jcl),lfn(icl,jcl-1),lfn(icl,jcl), &
                             lff(2,1),tff(2,1),time_now,icl,jcl-1,icl,jcl)

      lff(2,2)=lfn(icl,jcl)
      tff(2,2)=tign(icl,jcl)

      call tign_lfn_line_interp(tign(icl,jcl+1),tign(icl,jcl),lfn(icl,jcl+1),lfn(icl,jcl), &
                             lff(2,3),tff(2,3),time_now,icl,jcl+1,icl,jcl)

      call tign_lfn_four_pnts_interp(tign(icl,jcl-1),tign(icl,jcl),tign(icl+1,jcl-1),  & 
                                     tign(icl+1,jcl),lfn(icl,jcl-1),lfn(icl,jcl),      &
                                     lfn(icl+1,jcl-1),lfn(icl+1,jcl),lff(3,1),tff(3,1),time_now)

      call tign_lfn_line_interp(tign(icl+1,jcl),tign(icl,jcl),lfn(icl+1,jcl),lfn(icl,jcl), &
                             lff(3,2),tff(3,2),time_now,icl+1, jcl, icl, jcl)

      call tign_lfn_four_pnts_interp(tign(icl,jcl),tign(icl,jcl+1),tign(icl+1,jcl),  & 
                                     tign(icl+1,jcl+1),lfn(icl,jcl),lfn(icl,jcl+1),      &
                                     lfn(icl+1,jcl),lfn(icl+1,jcl+1),lff(3,3),tff(3,3),time_now)

end subroutine tign_lfn_interpolation

!
!************************
!

subroutine tign_lfn_line_interp(tign1,tign2,lfn1,lfn2,lfn_subcl,tign_subcl,time_now,i1,j1,i2,j2)
!***purpose: computes time of ignition of the point(*_subcl) that lies on the line
! between two points, whose lfn and tign is known

!*** arguments
!
real,intent(in) :: tign1,tign2   ! ignition times at the two endpoints
real,intent(in) :: lfn1,lfn2     ! level set function at the two endpoints
real,intent(in) :: time_now      ! tign>=time_now => no fire
real,intent(out) :: lfn_subcl,tign_subcl ! interpolated to midpoint
integer,intent(in):: i1,j1,i2,j2 ! indices of the endpoints of the line, for checking only

! Case 1: both points are on fire -> tign is interpolated linearly
! Case 2: both points are not on fire -> tign_subcl=time_now, not burning
! Case 3: One point is on fire, another is not -> tign_subcl set from
!         the equation  tign - time_now = c*lfn, where the proportionality
!         constant is found from the values of tign and lfn at the burning endpoint.
!         In particular, lfn(x)=0 implies tign(x)=time_now, so the representations of
!         the fireline by lfn and tign are consistent.
!         This is a special case of Case 3 from subroutine tign_lfn_four_pnts_interp.

!*** local
real :: c
character(len=128)::msg

!*** executable

! check consistency of inputs
call check_lfn_tign_ij(i1,j1,'tign_lfn_line_interp',time_now,lfn1,tign1) 
call check_lfn_tign_ij(i2,j2,'tign_lfn_line_interp',time_now,lfn2,tign2) 

lfn_subcl=0.5*(lfn1+lfn2)  ! interpolate lfn linearly

if (.not. lfn_subcl < 0.) then ! never test floats on equality
   tign_subcl=time_now         ! midpoint not on fire => none on fire
elseif ((lfn1 < 0.).AND.(lfn2 < 0.)) then
   tign_subcl=0.5*(tign1+tign2) ! both on fire, interpolate linearly
else  ! one is on fire one is not
   if (lfn1 < 0.) then  ! 1 is on fire, 2 is not
      c = (tign1-time_now)/lfn1
   elseif (lfn2 < 0.) then
      c = (tign2-time_now)/lfn2
   else
      call crash('tign_lfn_line_interp: one of lfn1 or lfn2 should be < 0')
   endif
   if( c < 0.)then
!$OMP CRITICAL(SFIRE_CORE_CRIT)
       write(msg,*)'tign1,tign2,time_now ',tign1,tign2,time_now 
       call message(msg,level=0)
       write(msg,*)'tign1-time_now,lfn1 ',tign1-time_now,lfn1 
       call message(msg,level=0)
       write(msg,*)'tign2-time_now,lfn2 ',tign2-time_now,lfn2 
       call message(msg,level=0)
       write(msg,*)'c ',c 
       call message(msg,level=0)
!$OMP END CRITICAL(SFIRE_CORE_CRIT)
       call crash('tign_lfn_line_interp: bad ignition times, c<0')
   endif
   tign_subcl=c*lfn_subcl+time_now;
endif
end subroutine tign_lfn_line_interp
!
!************************
!

subroutine tign_lfn_four_pnts_interp(tign1,tign2,tign3,tign4, &
  lfn1,lfn2,lfn3,lfn4,lfn_subcl,tign_subcl,time_now)

!***purpose: computes time of ignition of the point(*_subcl) that lies on the middle
! of square with given 4 points on its ends, whose lfn and tign is known
! since lfn is interpolated linearly, lfn_subcl is known 

!***arguments
real,intent(in) :: tign1,tign2,tign3,tign4 ! time of ignition of all corner points
real,intent(in) :: lfn1,lfn2,lfn3,lfn4 ! lfn of central and all corner points
real,intent(in) :: time_now
real,intent(out) :: lfn_subcl,tign_subcl

! Case 1: all 4 points are on fire -> tign is interpolated linearly
!   tign_subcl=0.25*(tign1+tign2+tign3+tign4)
! Case 2: all points are not on fire -> subcl - not burning
! Case 3: some points are on fire, others are not
! Here we assume that when lfn(x)=0 -> tign(x)=time_now (values interpolated to point x)
! For this case, tign of central point is approximated by
! tign~=c*lfn+time_now, which for our case is
! tign_subcl~=c*lfn_subcl+time_now
! where c is computed by least squares from the values at the points that are on fire
!
! sum(c*lfn(Ai)+time_now-tign(Ai))^2 ---> min
! for all lfn(Ai)<0
!
! solution for 'c' would be 
!
!         sum(tign(Ai)-time_now)*lfn(Ai)
!    c=   ------------------------------, both sums are over Ai, s.t lfn(Ai)<0
!         sum(lfn(Ai)*lfn(Ai))
!


real :: a,b,c,scale
character(len=128)::msg

call check_lfn_tign_ij(0,1,'tign_lfn_four_pnts_interp',time_now,lfn1,tign1) 
call check_lfn_tign_ij(0,2,'tign_lfn_four_pnts_interp',time_now,lfn2,tign2) 
call check_lfn_tign_ij(0,3,'tign_lfn_four_pnts_interp',time_now,lfn3,tign3) 
call check_lfn_tign_ij(0,4,'tign_lfn_four_pnts_interp',time_now,lfn4,tign4) 
 
lfn_subcl=0.25*(lfn1+lfn2+lfn3+lfn4)

if(.not. lfn_subcl < 0.) then ! midpoint not on fire, most frequent
      ! Case 2  
      tign_subcl=time_now
elseif((lfn1 < 0.).AND.(lfn2 < 0.).AND.(lfn3 < 0.).AND.(lfn4 < 0.)) then
      ! Case 1
      tign_subcl=0.25*(tign1+tign2+tign3+tign4)  ! all on fire, interpolate
else ! some on fire
      ! Case 3
      ! tign_subcl~=c*lfn_subcl+time_now
      scale = -minval((/lfn1,lfn2,lfn3,lfn4/))
      a=0.
      b=0.
      if (lfn1 < 0.) then
         a=a+(lfn1/scale)**2
         b=b+(tign1-time_now)*(lfn1/scale)
      endif
      if (lfn2 < 0.) then
           a=a+(lfn2/scale)**2
           b=b+(tign2-time_now)*(lfn2/scale)
      endif
      if (lfn3 < 0.) then
           a=a+(lfn3/scale)**2
           b=b+(tign3-time_now)*(lfn3/scale)
      endif
      if (lfn4 < 0.) then
           a=a+(lfn4/scale)**2
           b=b+(tign4-time_now)*(lfn4/scale)
      endif 
      if (.not. a>0. .or. .not. scale>0.)then
!$OMP CRITICAL(SFIRE_CORE_CRIT)
          write(msg,*)'lfn1,2,3,4,_subcl,ssq = ',lfn1,lfn2,lfn3,lfn3,lfn_subcl,a,' scale=',scale
!$OMP END CRITICAL(SFIRE_CORE_CRIT)
          call message(msg,level=0)
          call crash('tign_lfn_four_pnts_interp: internal: not a>0, none of lfn < 0 ?')
      endif
      if( b < 0.)call crash('tign_lfn_four_pnts_interp: bad ignition times, b<0') ! can have 0 because of rounding
      c=b/a
      tign_subcl=c*(lfn_subcl/scale)+time_now;
endif

end subroutine tign_lfn_four_pnts_interp


!
!************************
!


subroutine fuel_left_cell_1( fuel_frac_left, fire_frac_area, &
    lfn00,lfn01,lfn10,lfn11, &
    tign00,tign01,tign10,tign11,&
    time_now, fuel_time_cell)
!*** purpose: compute the fuel fraction left in one cell
implicit none
!*** arguments
real, intent(out):: fuel_frac_left, fire_frac_area ! 
real, intent(in)::lfn00,lfn01,lfn10,lfn11    ! level set function at 4 corners of the cell
real, intent(in)::tign00,tign01,tign10,tign11! ignition time at the  4 corners of the cell
real, intent(in)::time_now                   ! the time now
real, intent(in)::fuel_time_cell            ! time to burns off to 1/e
!*** Description
! The area burning is given by the condition L <= 0, where the function P is
! interpolated from lfn(i,j)
!
! The time since ignition is the function T, interpolated in from tign(i,j)-time_now.
! The values of tign(i,j) where lfn(i,j)>=0 are ignored, tign(i,j)=0 is taken 
! when lfn(i,j)=0.
!
! The function computes an approxmation  of the integral
!
!
!                                  /\
!                                  |              
! fuel_frac_left  =      1   -     | 1 -  exp(-T(x,y)/fuel_time_cell)) dxdy
!                                  |            
!                                 \/
!                                0<x<1
!                                0<y<1
!                             L(x,y)<=0
!
! When the cell is not burning at all (all lfn>=0), then fuel_frac(i,j)=1.
! Because of symmetries, the result should not depend on the mesh spacing dx dy
! so dx=1 and dy=1 assumed.
!
! Example:
!
!        lfn<0         lfn>0
!      (0,1) -----O--(1,1)            O = points on the fireline, T=tnow
!            |      \ |               A = the burning area for computing
!            |       \|                        fuel_frac(i,j)
!            |   A    O 
!            |        |
!            |        |
!       (0,0)---------(1,0)
!       lfn<0          lfn<0
!
! Approximations allowed: 
! The fireline can be approximated by straight line(s).
! When all cell is burning, approximation by 1 point Gaussian quadrature is OK.
! 
! Requirements:
! 1. The output should be a continuous function of the arrays lfn and
!  tign whenever lfn(i,j)=0 implies tign(i,j)=tnow.  
! 2. The output should be invariant to the symmetries of the input in each cell.
! 3. Arbitrary combinations of the signs of lfn(i,j) should work.
! 4. The result should be at least 1st order accurate in the sense that it is
!    exact if the time from ignition is a linear function.
!
! If time from ignition is approximated by polynomial in the burnt
! region of the cell, this is integral of polynomial times exponential
! over a polygon, which can be computed exactly.
!
! Requirement 4 is particularly important when there is a significant decrease
! of the fuel fraction behind the fireline on the mesh scale, because the
! rate of fuel decrease right behind the fireline is much larger 
! (exponential...). This will happen when
!
! change of time from ignition within one mesh cell / fuel_time_cell is not << 1
!
! This is the same as
!
!               mesh cell size
!  X =    -------------------------      is not << 1
!       fireline speed * fuel_time_cell
!         
!
! When X is large then the fuel burnt in one timestep in the cell is
! approximately proportional to length of  fireline in that cell.
!
! When X is small then the fuel burnt in one timestep in the cell is
! approximately proportional to the area of the burning region.
!

!*** calls
intrinsic tiny

!*** local
real::ps,aps,area,ta,out
real::t00,t01,t10,t11
real,parameter::safe=tiny(aps)
character(len=128)::msg

! the following algorithm is a very crude approximation

! minus time since ignition, 0 if no ignition yet
! it is possible to have 0 in fire region when ignitin time falls in 
! inside the time step because lfn is updated at the beginning of the time step

t00=tign00-time_now
if(lfn00>0. .or. t00>0.)t00=0.
t01=tign01-time_now
if(lfn01>0. .or. t01>0.)t01=0.
t10=tign10-time_now
if(lfn10>0. .or. t10>0.)t10=0.
t11=tign11-time_now
if(lfn11>0. .or. t11>0.)t11=0.

! approximate burning area, between 0 and 1   
ps = lfn00+lfn01+lfn10+lfn11   
aps = abs(lfn00)+abs(lfn01)+abs(lfn10)+abs(lfn11)
aps=max(aps,safe)
area =(-ps/aps+1.)/2.
area = max(area,0.) ! make sure area is between 0 and 1
area = min(area,1.)
    
! average negative time since ignition
ta=0.25*(t00+t01+t10+t11)

! exp decay in the burning area
out=1.
!if(area>0.)out=1. - area*(1. - exp(ta/fuel_time_cell))
if(area>0)out=area*exp(ta/fuel_time_cell) + (1. - area)

if(out>1.)then
!$OMP CRITICAL(SFIRE_CORE_CRIT)
    write(msg,*)'out=',out,'>1 area=',area,' ta=',ta
    call message(msg)
    write(msg,*)'tign=',tign00,tign01,tign10,tign11,' time_now=',time_now
!$OMP END CRITICAL(SFIRE_CORE_CRIT)
    call message(msg)
    !call message('WARNING: fuel_left_cell_1: fuel fraction > 1')
    call crash('fuel_left_cell_1: fuel fraction > 1')
endif

!out = max(out,0.) ! make sure out is between 0 and 1
!out = min(out,1.)

fuel_frac_left = out
fire_frac_area = area

end subroutine fuel_left_cell_1

!
!****************************************
!
! function calculation fuel_frac made by Volodymyr Kondratenko on the base of
! the Matlab code made by Jan Mandel and the fortran port by Minjeong Kim

subroutine fuel_left_cell_2( fuel_frac_left, fire_frac_area, &
                             lfn00,lfn01,lfn10,lfn11, &
                             tign00,tign01,tign10,tign11,&
                             time_now, fuel_time_cell)
!*** purpose: compute the fuel fraction left in one cell
implicit none
!*** arguments
real, intent(out):: fuel_frac_left, fire_frac_area ! 
real, intent(in)::lfn00,lfn01,lfn10,lfn11    ! level set function at 4 corners of the cell
real, intent(in)::tign00,tign01,tign10,tign11! ignition time at the  4 corners of the cell
real, intent(in)::time_now                   ! the time now
real, intent(in)::fuel_time_cell             ! burns time to 1/e
!*** Description
! The burning area is given by the condition L <= 0, where the function L is
! interpolated from lfn(i,j)
!
! The time since ignition is the function T, interpolated in from tign(i,j)-time_now.
! The values of tign(i,j) where lfn(i,j)>=0 are ignored, tign(i,j)=0 is taken 
! when lfn(i,j)=0.
!
! The function computes an approxmation  of the integral
!
!
!                                  /\
!                                  |              
! fuel_frac_left  =      1   -     | 1 -  exp(-T(x,y)/fuel_time_cell)) dxdy
!                                  |            
!                                 \/
!                                0<x<1
!                                0<y<1
!                             L(x,y)<=0
!
! When the cell is not burning at all (all lfn>=0), then fuel_frac(i,j)=1.
! Because of symmetries, the result should not depend on the mesh spacing dx dy
! so dx=1 and dy=1 assumed.
!
! Example:
!
!        lfn<0         lfn>0
!      (0,1) -----O--(1,1)            O = points on the fireline, T=tnow
!            |      \ |               A = the burning area for computing
!            |       \|                        fuel_frac(i,j)
!            |   A    O 
!            |        |
!            |        |
!       (0,0)---------(1,0)
!       lfn<0          lfn<0
!
! Approximations allowed: 
! The fireline can be approximated by straight line(s).
! When all cell is burning, approximation by 1 point Gaussian quadrature is OK.
! 
! Requirements:
! 1. The output should be a continuous function of the arrays lfn and
!  tign whenever lfn(i,j)=0 implies tign(i,j)=tnow.  
! 2. The output should be invariant to the symmetries of the input in each cell.
! 3. Arbitrary combinations of the signs of lfn(i,j) should work.
! 4. The result should be at least 1st order accurate in the sense that it is
!    exact if the time from ignition is a linear function.
!
! If time from ignition is approximated by polynomial in the burnt
! region of the cell, this is integral of polynomial times exponential
! over a polygon, which can be computed exactly.
!
! Requirement 4 is particularly important when there is a significant decrease
! of the fuel fraction behind the fireline on the mesh scale, because the
! rate of fuel decrease right behind the fireline is much larger 
! (exponential...). This will happen when
! 
! change of time from ignition within one mesh cell * fuel speed is not << 1
!
! This is the same as
!
!         mesh cell size*fuel_speed 
!         -------------------------      is not << 1
!             fireline speed
!         
!
! When X is large then the fuel burnt in one timestep in the cell is
! approximately proportional to length of  fireline in that cell.
!
! When X is small then the fuel burnt in one timestep in the cell is
! approximately proportional to the area of the burning region.

!*** calls
intrinsic tiny
#define DREAL real(kind=8)
!*** local
DREAL::ps,aps,area,ta,out
DREAL::t00,t01,t10,t11
DREAL,parameter::safe=tiny(aps)
DREAL::dx,dy ! mesh sizes
integer::i,j,k

DREAL,dimension(3)::u

DREAL::tweight,tdist
integer::kk,ll,ss
DREAL::rnorm
DREAL,dimension(8,2)::xylist,xytlist
DREAL,dimension(8)::tlist,llist,xt
DREAL,dimension(5)::xx,yy
DREAL,dimension(5)::lfn,tign

integer:: npoint
DREAL::tt,x0,y0,xts,xte,yts,yte,xt1,xt2
DREAL::lfn0,lfn1,dist,nr,s,errQ,ae,ce,ceae,a0,a1,a2,d,cet
DREAL::s1,s2,s3
DREAL::upper,lower,ah,ch,aa,cc,aupp,cupp,alow,clow
DREAL,dimension(2,2)::mQ
DREAL,dimension(2)::ut
character(len=128)::msg

!calls
intrinsic epsilon

DREAL, parameter:: zero=0.,one=1.,eps=epsilon(zero)

!!!! For finite differences by VK
DREAL::tign_middle,dt_dx,dt_dy,lfn_middle,a,b,c
DREAL:: alg_err

!*** executable

call crash('fuel_left_method=2 not working, please use fuel_left_method=1')

! check consistency
call check_lfn_tign_ij(0,0,'fuel_left_cell_2',time_now,lfn00,tign00) 
call check_lfn_tign_ij(0,1,'fuel_left_cell_2',time_now,lfn01,tign01) 
call check_lfn_tign_ij(1,0,'fuel_left_cell_2',time_now,lfn10,tign10) 
call check_lfn_tign_ij(1,1,'fuel_left_cell_2',time_now,lfn11,tign11) 

alg_err=0
dx=1
dy=1
t00=time_now-tign00
if(lfn00>=0. .or. t00<0.)t00=0.
t01=time_now-tign01
if(lfn01>=0. .or. t01<0.)t01=0.
t10=time_now-tign10
if(lfn10>=0. .or. t10<0.)t10=0.
t11=time_now-tign11
if(lfn11>=0. .or. t11<0.)t11=0.

! approximate burning area, between 0 and 1  
! was taken from fuel_left_cell_1 made by Jan, will need to be changed 
ps = lfn00+lfn01+lfn10+lfn11
aps = abs(lfn00)+abs(lfn01)+abs(lfn10)+abs(lfn11)
aps=max(aps,safe)
area =(-ps/aps+1.)/2.
area = max(area,zero) ! make sure area is between 0 and 1
area = min(area,one)

!*** case0 Do nothing
if ( lfn00>=0 .and. lfn10>=0 .and. lfn01>=0 .and. lfn11>=0 ) then
    out = 1.0 !  fuel_left, no burning
    area= 0.
!*** case4 all four coners are burning
else if (lfn00<=0 .and. lfn10<=0 .and. lfn01<=0 .and. lfn11<=0) then
! All burning
! T=u(1)*x+u(2)*y+u(3)
! t(0,0)=tign(1,1)
! t(0,fd(2))=tign(1,2)
! t(fd(1),0)=tign(2,1)
! t(fd(1),fd(2))=tign(2,2)
! t(g/2,h/2)=sum(tign(i,i))/4
! dt/dx=(1/2h)*(t10-t00+t11-t01)
! dt/dy=(1/2h)*(t01-t00+t11-t10)
! approximate T(x,y)=u(1)*x+u(2)*y+u(3) Using finite differences
! t(x,y)=t(h/2,h/2)+(x-h/2)*dt/dx+(y-h/2)*dt/dy
! u(1)=dt/dx
! u(2)=dt/dy
! u(3)=t(h/2,h/2)-h/2(dt/dx+dt/dy)

 tign_middle=(t00+t01+t10+t11)/4

 ! since mesh_size is 1 we replace fd(1) and fd(2) by 1
 dt_dx=(t10-t00+t11-t01)/2
 dt_dy=(t01-t00+t11-t10)/2

 u(1)=dt_dx
 u(2)=dt_dy
 u(3)=tign_middle-(dt_dx+dt_dy)/2

! integrate
u(1)=-u(1)/fuel_time_cell
u(2)=-u(2)/fuel_time_cell
u(3)=-u(3)/fuel_time_cell
    s1=u(1)
    s2=u(2)            
    out=exp(u(3))*intexp(s1)*intexp(s2)
    area=1
    if ( out<0 .or. out>1.0 ) then
        call message('WARNING: fuel_left_cell: case all burning: out should be between 0 and 1')
    end if
!*** case 1,2,3- other cases
!*** part of cell is burning 
else
    ! set xx, yy for the coner points
    ! move these values out of i and j loop to speed up
    ! comments for xx, yy - make center [0,0], cyclic, counterclockwise
    ! comments for lfn,tign - cyclic, counterclockwise
    xx(1) = -0.5
    xx(2) =  0.5
    xx(3) =  0.5
    xx(4) = -0.5
    xx(5) = -0.5
    yy(1) = -0.5
    yy(2) = -0.5
    yy(3) =  0.5
    yy(4) =  0.5
    yy(5) = -0.5     
    lfn(1)=lfn00
    lfn(2)=lfn10
    lfn(3)=lfn11
    lfn(4)=lfn01
    lfn(5)=lfn00
    tign(1)=t00
    tign(2)=t10
    tign(3)=t11
    tign(4)=t01
    tign(5)=t00
    npoint = 0 ! number of points in polygon

    do k=1,4
        lfn0=lfn(k  )
        lfn1=lfn(k+1)
        if ( lfn0 <= 0.0 ) then
            npoint = npoint + 1
            xylist(npoint,1)=xx(k) ! add corner to list
            xylist(npoint,2)=yy(k)
            tlist(npoint)=-tign(k) ! time since ignition
            llist(npoint)=lfn0
        end if
        if ( lfn0*lfn1 < 0 ) then
            npoint = npoint + 1
! coordinates of intersection of the fire line with segment k k+1
! lfn(t)=lfn0+t*(lfn1-lfn0)=0            
            tt=lfn0/(lfn0-lfn1)
            x0=xx(k)+( xx(k+1)-xx(k) )*tt
            y0=yy(k)+( yy(k+1)-yy(k) )*tt
            xylist(npoint,1)=x0
            xylist(npoint,2)=y0
            tlist(npoint)=0 ! now at ignition
            llist(npoint)=0 ! on fireline
        end if
    end do

    ! make the list circular and trim to size
    tlist(npoint+1)=tlist(1)
    llist(npoint+1)=llist(1)   
    xylist(npoint+1,1)=xylist(1,1)
    xylist(npoint+1,2)=xylist(1,2)

! approximate L(x,y)=u(1)*x+u(2)*y+u(3)
 lfn_middle=(lfn00+lfn01+lfn10+lfn11)/4
 dt_dx=(lfn10-lfn00+lfn11-lfn01)/2
 dt_dy=(lfn01-lfn00+lfn11-lfn10)/2
 u(1)=dt_dx
 u(2)=dt_dy
 u(3)=lfn_middle-(dt_dx+dt_dy)/2
! finding the coefficient c, reminder we work over one subcell only
! T(x,y)=c*L(x,y)+time_now
    a=0
    b=0

    if (lfn00 <= 0) then
                        a=a+lfn00*lfn00
                        if (t00 < 0) then
                        call crash('fuel_burnt_fd: tign(i1) should be less then time_now')
                        else
                        b=b+t00*lfn00
                        end if
    end if


    if (lfn01 <= 0) then
                        a=a+lfn01*lfn01
                        if (t01< 0) then
                        call crash('fuel_burnt_fd: tign(i1) should be less then time_now')
                        else
                        b=b+t01*lfn01
                        end if
    end if


    if (lfn10<=0) then
                        a=a+lfn10*lfn10
                        if (t10<0) then
                        call crash('fuel_burnt_fd: tign(i1) should be less then time_now')
                        else
                        b=b+t10*lfn10
                        end if
  end if

    if (lfn11<=0) then
                        a=a+lfn11*lfn11
                        if (t11<0) then
                        call crash('fuel_burnt_fd: tign(i1) should be less then time_now')
                        else
                        b=b+t11*lfn11
                        end if
    end if


                     if (a==0) then
                        call crash('fuel_burnt_fd: if c is on fire then one of cells should be on fire')
                     end if
                        c=b/a
                        u(1)=u(1)*c
                        u(2)=u(2)*c
                        u(3)=u(3)*c

    ! rotate to gradient on x only
    nr = sqrt(u(1)**2+u(2)**2)
    if(.not.nr.gt.eps)then
        out=1.
        goto 900
    endif
    c = u(1)/nr
    s = u(2)/nr
! rotation such that Q*u(1:2)-[something;0]
    mQ(1,1)=c
    mQ(1,2)=s
    mQ(2,1)=-s
    mQ(2,2)=c            
    ! mat vec multiplication
    call matvec(mQ,2,2,u,3,ut,2,2,2)            
    errQ = ut(2) ! should be zero            
    ae = -ut(1)/fuel_time_cell
    ce = -u(3)/fuel_time_cell ! -T(xt,yt)/fuel_time=ae*xt+ce      
    cet=ce!keep ce for later
    call matmatp(xylist,8,2,mQ,2,2,xytlist,8,2,npoint+1,2,2)            
    call sortxt( xytlist, 8,2, xt,8,npoint ) !sort ascending in x           
    out=0.0
    aupp=0.0
    cupp=0.0
    alow=0.0
    clow=0.0
    do k=1,npoint-1
        xt1=xt(k)
        xt2=xt(k+1)
        upper=0
        lower=0
        ah=0
        ch=0
        if ( xt2-xt1 > eps*100 ) then
! slice of nonzero width                
! find slice height as h=ah*x+ch
            do ss=1,npoint ! pass counterclockwise
                xts=xytlist(ss,1) ! start point of the line
                yts=xytlist(ss,2)
                xte=xytlist(ss+1,1) ! end point of the line
                yte=xytlist(ss+1,2)
                  
                if ( (xts>xt1 .and. xte>xt1) .or. &
                     (xts<xt2 .and. xte<xt2) ) then
                    aa = 0 ! do nothing
                    cc = 0
                else ! line y=a*x+c through (xts,yts), (xte,yte)
                    aa = (yts-yte)/(xts-xte)
                    cc = (xts*yte-xte*yts)/(xts-xte)                    
                    if (xte<xts) then ! upper boundary
                        aupp = aa
                        cupp = cc
                        ah=ah+aa
                        ch=ch+cc
                        upper=upper+1
                    else              ! lower boundary
                        alow = aa
                        clow = cc
                        lower=lower+1
                    end if
                end if!(xts>xt1 .and. xte>xt1)              
            end do ! ss
            ce=cet !use stored ce
!shift small amounts exp(-**) to avoid negative fuel burnt
            if (ae*xt1+ce > 0 ) then
              ce=ce-(ae*xt1+ce)!
            end if
            if (ae*xt2+ce > 0) then
            ce=ce-(ae*xt2+ce)
            end if

            ah = aupp-alow
            ch = cupp-clow  
            ! integrate (ah*x+ch)*(1-exp(ae*x+ce) from xt1 to xt2
            ! numerically sound for ae->0, ae -> infty
            ! this can be important for different model scales
            ! esp. if someone runs the model in single precision!!
            ! s1=int((ah*x+ch),x,xt1,xt2)
            s1 = (xt2-xt1)*((1./2.)*ah*(xt2+xt1)+ch)            
            ! s2=int((ch)*(-exp(ae*x+ce)),x,xt1,xt2)
            ceae=ce/ae;
            s2 = -ch*exp(ae*(xt1+ceae))*(xt2-xt1)*intexp(ae*(xt2-xt1))                
            ! s3=int((ah*x)*(-exp(ae*x+ce)),x,xt1,xt2)
            ! s3=int((ah*x)*(-exp(ae*(x+ceae))),x,xt1,xt2)
            ! expand in Taylor series around ae=0
            ! collect(expand(taylor(int(x*(-exp(ae*(x+ceae))),x,xt1,xt2)*ae^2,ae,4)/ae^2),ae)
            ! =(1/8*xt1^4+1/3*xt1^3*ceae+1/4*xt1^2*ceae^2-1/8*xt2^4-1/3*xt2^3*ceae-1/4*xt2^2*ceae^2)*ae^2
            !     + (-1/3*xt2^3-1/2*xt2^2*ceae+1/3*xt1^3+1/2*xt1^2*ceae)*ae 
            !     + 1/2*xt1^2-1/2*xt2^2
            !
            ! coefficient at ae^2 in the expansion, after some algebra            
            a2=(xt1-xt2)*((1./4.)*(xt1+xt2)*ceae**2+(1./3.)* &
               (xt1**2+xt1*xt2+xt2**2)*ceae+(1./8.)* &
               (xt1**3+xt1*(xt2**2)+xt1**2*xt2+xt2**3))               
            d=(ae**4)*a2
            
            if (abs(d)>eps) then
            ! since ae*xt1+ce<=0 ae*xt2+ce<=0 all fine for large ae
            ! for ae, ce -> 0 rounding error approx eps/ae^2
                s3=( exp(ae*(xt1+ceae))*(ae*xt1-1)-&
                     exp(ae*(xt2+ceae))*(ae*xt2-1) )/(ae**2)
                
            !we do not worry about rounding as xt1 -> xt2, then s3 -> 0
            else
                ! coefficient at ae^1 in the expansion
                a1=(xt1-xt2)*((1./2.)*ceae*(xt1+xt2)+(1./3.)*&
                   (xt1**2+xt1*xt2+xt2**2))
                ! coefficient at ae^0 in the expansion for ae->0
                a0=(1./2.)*(xt1-xt2)*(xt1+xt2)
                s3=a0+a1*ae+a2*ae**2; ! approximate the integral
                            end if

            s3=ah*s3                                                
            out=out+s1+s2+s3
            if(out<0. .or. out>1.) then                                  
             write(msg,'(a,g14.4,a)')'WARNING::fuel_fraction ',out,' should be between 0 and 1'
            end if!print
            if(out.ne.out.or..not.out.le.huge(out).or..not.out.ge.-huge(out)) call crash('fuel_fraction out is not a valid number')
                
        end if
    end do ! k     

            out=1-out !fuel_left
end if ! if case0, elseif case4 ,else case123

900 continue 
fuel_frac_left = out
fire_frac_area= area
if(isnotfinite(fuel_frac_left)) call crash('fuel_frac_left is not a valid number')
if(isnotfinite(fire_frac_area)) call crash('fire_frac_area is not a valid number')
end subroutine fuel_left_cell_2






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function intexp(ab)
implicit none
DREAL::ab
!calls
intrinsic epsilon

real, parameter:: zero=0.,one=1.,eps=epsilon(zero)

!eps = 2.2204*(10.0**(-8))!from matlab
if ( eps < abs(ab)**3/6. ) then
    intexp=(exp(ab)-1)/ab
  else
    intexp=1+ab/2.
end if
end function intexp
!
!****************************************
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         
!**************************************** 
!       


!
!****************************************


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine sortxt(xytlist,nrow,ncolumn,xt,nxt,nvec)
implicit none
integer::nrow,ncolumn,nxt,nvec
DREAL,dimension(nrow,ncolumn)::xytlist
DREAL,dimension(nxt)::xt

integer::i,j
DREAL::temp

do i=1,nvec
  xt(i)=xytlist(i,1)
end do

do i=1,nvec-1
  do j=i+1,nvec
    if ( xt(i) > xt(j) ) then
         temp = xt(i)
         xt(i)=xt(j)
         xt(j)=temp
    end if
  end do
end do

end subroutine sortxt
!
!****************************************
!
subroutine matvec(A,m,n,V,nv,out,nout,nrow,ncolumn)
implicit none
integer::m,n,nv,nout,nrow,ncolumn
DREAL,dimension(m,n)::A   ! allocated m by n 
DREAL,dimension(nv)::V    ! allocated nv
DREAL,dimension(nout)::out! allocated nout 

integer::i,j

do i=1,nrow
  out(i)=0.0
  do j=1,ncolumn
    out(i)=out(i)+A(i,j)*V(j)
  end do
end do
end subroutine
!
!****************************************
!
subroutine matmatp(A,mA,nA,B,mB,nB,C,mC,nC,nrow,ncolumn,nP)
implicit none
integer::mA,nA,mB,nB,mC,nC,nrow,ncolumn,nP
DREAL,dimension(mA,nA)::A   ! allocated m by n 
DREAL,dimension(mB,nB)::B   ! allocated m by n 
DREAL,dimension(mC,nC)::C   ! allocated m by n 
integer::i,j,k
do i=1,nrow  
  do j=1,ncolumn
    C(i,j)=0.0
  do k=1,nP
    C(i,j)=C(i,j)+A(i,k)*B(j,k) ! B'
  end do
end do
end do
end subroutine

!
!****************************************
!

subroutine prop_ls( id, &                                ! for debug
                ipart,  &                                ! 1:compute lfn_out 2:tign & speedcheck
                ids,ide,jds,jde, &                       ! domain dims
                ims,ime,jms,jme, &                       ! memory dims
                ips,ipe,jps,jpe, &                ! patch - nodes owned by this process 
                its,ite,jts,jte, &                       ! tile dims
                ts,dt,dx,dy,     &                       ! scalars in
                tbound,          &                       ! scalars out
                lfn_in,lfn_out,tign,ros,  &              ! arrays inout          
                fp               &
                   )
implicit none

!*** purpose: advance level function in time

! Jan Mandel August 2007 - February 2008

!*** description
!
! Propagation of closed curve by a level function method. The level function
! lfn is defined by its values at the nodes of a rectangular grid. 
! The area where lfn < 0 is inside the curve. The curve is 
! described implicitly by lfn=0. Points where the curve intersects gridlines
! can be found by linear interpolation from nodes.
!
! The level function is advanced from time ts to time ts + dt. 
!
! The level function should be initialized to (an approximation of) the signed
! distance from the curve. If the initial curve is a circle, the initial level
! function is simply the distance from the center minus the radius.
! 
! The curve moves outside with speed given by function speed_func.
!   
! Method: Godunov/ENO method for the normal motion. The timestep is checked for
! CFL condition. For a straight segment in a constant field and locally linear
! level function, the method reduces to the exact normal motion. The advantage of 
! the level set method is that it treats automatically special cases such as
! the curve approaching itself and merging components of the area inside the curve.
!
! Based on S. Osher and R. Fedkiw, Level set methods and dynamic implicit surfaces,
! Springer, 2003, Sec. 6.4, as implemented in toolboxLS for Matlab by 
! I. Mitchell, A toolbox of Level Set Methods (Version 1.1), TR-2007-11,
! Dept. Computer Science, University of British Columbia, 2007
! http://www.cs.ubc.ca/\~mitchell/Toolbo\LS
! 
  
!*** arguments 

! id                in    unique identification for prints and dumps
! ipart             in    1:compute lfn_out 2:tign & speedcheck 
!                         *** must call halo on lfn_out between part 1 and part 2***
! ids,ide,jds,jde   in    domain dimensions
! ims,ime,jms,jme   in    memory dimensions
! its,ite,jts,jte   in    tile dimensions
! ts                in    start time
! dt                in    time step
! dx,dy             in    grid spacing
! tbound            out   bound on stable time step from CFL condition, if tbound>=dt then OK
! lfn_in,lfn_out    inout,out the level set function at nodes
! tign              inout the ignition time at nodes

! The dimensions are cell-based, the nodal value is associated with the south-west corner.
! The whole computation is on domain indices ids:ide+1,jds:jde+1.
!
! The region where new lfn and tign are computed is the tile its:ite,jts:jte 
! except when the tile is at domain upper boundary, an extra band of points is added:
! if ite=ide then region goes up to ite+1, if jte=jde then region goes up to jte+1.

! The time step requires values from 2 rows of nodes beyond the region except when at the 
! domain boundary one-sided derivatives are used. This is implemented by extending the input
! beyond the domain boundary so sufficient memory bounds must be allocated. 
! The update on all tiles can be done in parallel. To avoid the race condition (different regions
! of the same array updated by different threads), the in and out versions of the
! arrays lft and tign are distinct. If the time step dt is larger
! that the returned tbound, the routine should be called again with timestep td<=tbound, and then
! having distinct inputs and outputs comes handy.

!*** calls
!
! tend_ls
!

integer,intent(in)::id,ipart,ims,ime,jms,jme,ids,ide,jds,jde,its,ite,jts,jte,ips,ipe,jps,jpe 
real,dimension(ims:ime,jms:jme),intent(inout)::lfn_in,tign
real,dimension(ims:ime,jms:jme),intent(out)::lfn_out,ros
real,intent(in)::dx,dy,ts,dt
real,intent(out)::tbound
type(fire_params),intent(in)::fp

!*** local 
! arrays
#define IMTS its-1
#define IMTE ite+1
#define JMTS jts-1
#define JMTE jte+1
real,dimension(IMTS:IMTE,JMTS:JMTE):: tend, lfn1 ! region-sized with halo
! scalars
real::grad2,rr,tbound2,a,a1 ! a=0 euler, a=0.5 heun

real::gradx,grady,aspeed,err,aerr,time_future,time_now
real::tmp,t0,t1,t2
integer::ihs,ihe,jhs,jhe
integer::ihs2,ihe2,jhs2,jhe2
integer::itso,iteo,jtso,jteo
integer::i,j,its1,ite1,jts1,jte1,k,kk,id1
character(len=128)::msg
integer::nfirenodes,nfireline,ierrx
real::sum_err,min_err,max_err,sum_aerr,min_aerr,max_aerr   

! constants
integer,parameter :: mstep=1000, printl=1
real, parameter:: zero=0.,one=1.,eps=epsilon(zero),tol=100*eps, &
    safe=2.,rmin=safe*tiny(zero),rmax=huge(zero)/safe

! f90 intrinsic function

intrinsic max,min,sqrt,nint,epsilon,tiny,huge
  
!*** executable

!$OMP CRITICAL(SFIRE_CORE_CRIT)
write(msg,'(a,i6,a,i2,4(a,i5))')'prop_ls:',id,' part',ipart,' tile',its,':',ite,',',jts,':',jte
!$OMP END CRITICAL(SFIRE_CORE_CRIT)
call message(msg)

if (ipart==1)then
!   lfn_in -> lfn_out

    call check_lfn_tign('prop_ls start',ts,its,ite,jts,jte,ims,ime,jms,jme,lfn_in,tign)

    a=fire_back_weight ! from module_fr_sfire_util
    a1=1. - a
    
    ! tend = F(lfn)

    ihs2=max(its-2,ids)   ! need lfn two beyond the tile but not outside the domain 
    ihe2=min(ite+2,ide)
    jhs2=max(jts-2,jds) 
    jhe2=min(jte+2,jde)

    ihs=max(its-1,ids)   ! compute tend one beyond the tile but not outside the domain 
    ihe=min(ite+1,ide)
    jhs=max(jts-1,jds) 
    jhe=min(jte+1,jde)


#ifdef DEBUG_OUT    
    call write_array_m(ihs,ihe,jhs,jhe,ims,ime,jms,jme,lfn_in,'lfn_in',id)
#endif

    ! check array dimensions
    call check_mesh_2dim(ihs2,ihe2,jhs2,jhe2,ims,ime,jms,jme)
    call print_2d_stats(ihs2,ihe2,jhs2,jhe2,ims,ime,jms,jme,lfn_in,'prop_ls: lfn in')
    
    ! NOTE: tend_ls will extrapolate to one node strip at domain boundaries
    ! so that it can compute gradient at domain boundaries.
    ! To avoid copying, lfn_in is declared inout.
    ! At tile boundaries that are not domain boundaries values of lfn_in two nodes
    ! outside of the tile are needed.
    id1 = id  ! for debug prints
    if(id1.ne.0)id1=id1+1000
    call  tend_ls( id1, &
    ims,ime,jms,jme, &                       ! memory dims for lfn_in
    IMTS,IMTE,JMTS,JMTE, &                   ! memory dims for tend 
    ids,ide,jds,jde, &                       ! domain dims - where lfn exists
    ips,ipe,jps,jpe, &                       ! patch - nodes owned by this process 
    ihs,ihe,jhs,jhe, &                       ! where tend computed
    ims,ime,jms,jme, &                       ! memory dims for ros 
    its,ite,jts,jte, &                       ! tile dims - where to set ros
    ts,dt,dx,dy,      &                      ! scalars in
    lfn_in, &                                ! arrays in
    tbound, &                                ! scalars out 
    tend, ros, &                              ! arrays out        
    fp         &                             ! params
)

#ifdef DEBUG_OUT    
    call write_array_m(ihs,ihe,jhs,jhe,IMTS,IMTE,JMTS,JMTE,tend,'tend1',id)
#endif

    ! Euler method, the half-step, same region as ted
    do j=jhs,jhe
        do i=ihs,ihe
            lfn1(i,j) = lfn_in(i,j) + dt*tend(i,j)
        enddo
    enddo
    
    call print_2d_stats(ihs,ihe,jhs,jhe,IMTS,IMTE,JMTS,JMTE,lfn1,'prop_ls: lfn1')
    ! tend = F(lfn1) on the tile (not beyond)

    if(id1.ne.0)id1=id1+1000
    call  tend_ls( id1,&
    IMTS,IMTE,JMTS,JMTE, &                   ! memory dims for lfn
    IMTS,IMTE,JMTS,JMTE, &                   ! memory dims for tend 
    ids,ide,jds,jde,     &                   ! domain dims - where lfn exists
    ips,ipe,jps,jpe, &                       ! patch - nodes owned by this process 
    its,ite,jts,jte, &                       ! tile dims - where is tend computed
    ims,ime,jms,jme, &                       ! memory dims for ros 
    its,ite,jts,jte, &                       ! tile dims - where is ros set
    ts+dt,dt,dx,dy,      &                   ! scalars in
    lfn1, &                                  ! arrays in
    tbound2, &                               ! scalars out 
    tend,ros, &                               ! arrays out        
    fp &
)

#ifdef DEBUG_OUT    
    call write_array_m(its,ite,jts,jte,IMTS,IMTE,JMTS,JMTE,tend,'tend2',id)
#endif

    call print_2d_stats(its,ite,jts,jte,ims,ime,jms,jme, ros,'prop_ls: ros')
    call print_2d_stats(its,ite,jts,jte,IMTS,IMTE,JMTS,JMTE,tend,'prop_ls: tend2')
        
    tbound=min(tbound,tbound2)

!$OMP CRITICAL(SFIRE_CORE_CRIT)
    write(msg,'(a,f10.2,4(a,f7.2))')'prop_ls: time',ts,' dt=',dt,' bound',min(tbound,999.99), &
        ' dx=',dx,' dy=',dy
!$OMP END CRITICAL(SFIRE_CORE_CRIT)
    call message(msg)
    if(dt>tbound)then
!$OMP CRITICAL(SFIRE_CORE_CRIT)
        write(msg,'(2(a,f10.2))')'prop_ls: WARNING: time step ',dt, &
        ' > bound =',tbound
!$OMP END CRITICAL(SFIRE_CORE_CRIT)
        call message(msg)
    endif
    
    ! combine lfn1 and lfn_in + dt*tend -> lfn_out
    
    do j=jts,jte
        do i=its,ite
            lfn_out(i,j) = a1*lfn1(i,j) + a*(lfn_in(i,j) + dt*tend(i,j))
            lfn_out(i,j) = min(lfn_out(i,j),lfn_in(i,j)) ! fire area can only increase
        enddo
    enddo     

elseif(ipart==2)then

    call continue_at_boundary(1,1,zero, &   !ext 
    ims,ime,jms,jme, &                ! memory dims
    ids,ide,jds,jde, &                    ! domain - nodes where lfn defined
    ips,ipe,jps,jpe, &                    ! patch - nodes owned by this process 
    its,ite,jts,jte, &                ! tile - nodes update by this thread
    itso,iteo,jtso,jteo, &                ! where set now
    lfn_out )                                  ! array

    ! copy to use consistently old values
    do j=jts,jte
        do i=its,ite
           lfn1(i,j)=lfn_out(i,j)
        enddo
    enddo

    !fix up lfn_out - no new local minima
    do j=jts,jte
        do i=its,ite
            t0=min(lfn_in(i+1,j),lfn_in(i-1,j),lfn_in(i,j+1),lfn_in(i,j-1))
            t1=min(lfn1(i+1,j),lfn1(i-1,j),lfn1(i,j+1),lfn1(i,j-1))
            if(.not.t0>lfn_in(i,j) .and. t1>lfn_out(i,j))then
                if(fire_print_msg>2)then
!$OMP CRITICAL(SFIRE_CORE_CRIT)
                    write(msg,'(a,2i6)')'prop_ls: new local min',i,j
                    call message(msg)
                    write(msg,'((a,g13.5,1x))') &
                       'prop_ls: new local min',lfn_out(i,j),'fixed to',t1,'incr',t1-lfn_in(i,j)
                    call message(msg)
!$OMP END CRITICAL(SFIRE_CORE_CRIT)
                endif
                lfn_out(i,j)=t1
            endif
        enddo
    enddo

    ! compute ignition time by interpolation
    ! the node was not burning at start but it is burning at end
    ! interpolate from the level functions at start and at end
    ! lfn_in   is the level set function value at time ts
    ! lfn_out  is the level set function value at time ts+dt
    ! 0        is the level set function value at time tign(i,j)
    ! thus assuming the level function is approximately linear =>
    ! tign(i,j)= ts + ((ts + td) - ts) * lfn_in / (lfn_in - lfn_out)
    !        =        ts + dt * lfn_in  / (lfn_in - lfn_out)
    !        = (ts + dt) + dt * lfn_out / (lfn_in - lfn_out)
    ! the second form is better because we want make sure tign < ts+dt (used only for lfn_out < 0)

    time_now=ts+dt
    time_future=ts+2*dt
    do j=jts,jte
        do i=its,ite
            call check_lfn_tign_ij(i,j,'prop_ls before',ts,lfn_in(i,j),tign(i,j))
            if (lfn_out(i,j) < 0. ) then  ! now on fire
               if( .not. lfn_in(i,j) < 0)  then ! was not before, set ignition time by interpolation
                  tign(i,j) = time_now + dt * lfn_out(i,j) / (lfn_in(i,j) - lfn_out(i,j))
                  ! tmp = max(tmp,abs(ts + dt * lfn_in(i,j) / (lfn_in(i,j) - lfn_out(i,j)) - tign(i,j)))
               endif  !  was already burning, leave ignition time unchanged
            else ! set the ignition time outside of burning region
               tign(i,j)=time_future
            endif
            call check_lfn_tign_ij(i,j,'prop_ls after',time_now,lfn_out(i,j),tign(i,j))
            ! lfn_out(i,j) = max(lfn_out(i,j),-100.0)  ! stability fix
        enddo
    enddo
    ! write(msg,*)'prop_ls: tign err',tmp
    ! call message(msg)
    
    ! check local speed error and stats 
    ! may not work correctly in parallel
    ! init stats
    nfirenodes=0
    nfireline=0
    sum_err=0.
    min_err=rmax
    max_err=rmin     
    sum_aerr=0.
    min_aerr=rmax
    max_aerr=rmin    
    its1=its+1
    jts1=jts+1
    ite1=ite-1
    jte1=jte-1
    ! loop over right inside of the domain
    ! cannot use values outside of the domain, would have to reflect and that
    ! would change lfn_out
    ! cannot use values outside of tile, not synchronized yet
    ! so in parallel mode the statistics is not accurate
    do j=jts1,jte1
        do i=its1,ite1
            if(lfn_out(i,j)>0.0)then   ! a point out of burning region
                if(lfn_out(i+1,j)<=0.or.lfn_out(i,j+1)<=0.or. & ! neighbor in burning region
                   lfn_out(i-1,j)<=0.or.lfn_out(i,j-1)<=0)then ! point next to fireline
                   gradx=(lfn_out(i+1,j)-lfn_out(i-1,j))/(2.0*dx) ! central differences
                   grady=(lfn_out(i,j+1)-lfn_out(i,j-1))/(2.0*dy)
                   grad2=sqrt(gradx*gradx+grady*grady)
                   aspeed = (lfn_in(i,j)-lfn_out(i,j))/(dt*max(grad2,rmin))                   
                   rr = speed_func(gradx,grady,dx,dy,i,j,fp,ierrx,msg)
                   err=aspeed-rr
                   sum_err=sum_err+err
                   min_err=min(min_err,err)
                   max_err=max(max_err,err)     
                   aerr=abs(err)
                   sum_aerr=sum_aerr+aerr
                   min_aerr=min(min_aerr,aerr)
                   max_aerr=max(max_aerr,aerr)
                   nfireline=nfireline+1
                endif
            else
                nfirenodes=nfirenodes+1
            endif
        enddo
    enddo
!$OMP CRITICAL(SFIRE_CORE_CRIT)
    write(msg,'(2(a,i6,f8.4))')'prop_ls: nodes burning',nfirenodes, &
      (100.*nfirenodes)/((ite1-its1+1)*(jte1-jts1+1)),'% next to fireline',nfireline
!$OMP END CRITICAL(SFIRE_CORE_CRIT)
   call message(msg)
!    if(nfireline>0)then
!        call print_stat_line('speed error',its1,ite1,jts1,jte1,min_err,max_err,sum_err/nfireline)
!        call print_stat_line('abs(speed error)',its1,ite1,jts1,jte1,min_aerr,max_aerr,sum_aerr/nfireline)
!    endif

    ! check if the fire did not get to the domain boundary
    do k=-1,1,2
        ! do kk=1,(boundary_guard*(ide-ids+1))/100  ! in %
        do kk=1,boundary_guard   ! measured in cells
            i=ids+k*kk
            if(i.ge.its.and.i.le.ite)then
                do j=jts,jte
                    if(lfn_out(i,j)<=0.)goto 9
                enddo
            endif
    enddo
        ! do kk=1,(boundary_guard*(jde-jds+1))/100
        do kk=1,boundary_guard    ! measured in cells
            j=jds+k*kk
            if(j.ge.jts.and.j.le.jte)then
                do i=its,ite
                    if(lfn_out(i,j)<=0.)goto 9
                enddo
            endif
        enddo
    enddo
    goto 10
9   continue
!$OMP CRITICAL(SFIRE_CORE_CRIT)
    write(msg,'(a,i2,a,2i8)')'prop_ls: fire',boundary_guard, &
        ' cells from domain boundary at node ',i,j
!$OMP END CRITICAL(SFIRE_CORE_CRIT)
    call message(msg)     
    call crash('prop_ls: increase the fire region')
10  continue

    call print_2d_stats(its,ite,jts,jte,ims,ime,jms,jme, lfn_out,'prop_ls: lfn out')
    call print_2d_stats(its,ite,jts,jte,ims,ime,jms,jme, tign,'prop_ls: tign out')

    call check_lfn_tign('prop_ls end',ts+dt,its,ite,jts,jte,ims,ime,jms,jme,lfn_out, tign)

else
    call crash('prop_ls: ipart must be 1 or 2')
endif
end subroutine prop_ls

!
!*****************************
!

subroutine tend_ls( id, &
    lims,lime,ljms,ljme, &                   ! memory dims for lfn
    tims,time,tjms,tjme, &                   ! memory dims for tend 
    ids,ide,jds,jde, &                       ! domain - nodes where lfn defined
    ips,ipe,jps,jpe, &                       ! patch - nodes owned by this process 
    ints,inte,jnts,jnte, &                   ! region - nodes where tend computed
    ims,ime,jms,jme, &                       ! memory dims for ros 
    its,ite,jts,jte, &                       ! tile dims - where is ros set
    t,dt,dx,dy,      &                       ! scalars in
    lfn, &                                   ! arrays in
    tbound, &                                ! scalars out 
    tend, ros,  &                              ! arrays out
    fp &
)

implicit none
! purpose
! compute the right hand side of the level set equation

!*** arguments
integer,intent(in)::id,lims,lime,ljms,ljme,tims,time,tjms,tjme
integer,intent(in)::ims,ime,jms,jme,its,ite,jts,jte
integer, intent(in)::ids,ide,jds,jde,ints,inte,jnts,jnte,ips,ipe,jps,jpe 
real,intent(in)::t                                     ! time
real,intent(in)::dt,dx,dy                                 ! mesh step
real,dimension(lims:lime,ljms:ljme),intent(inout)::lfn ! level set function
real,intent(out)::tbound                               ! max allowed time step
real,dimension(tims:time,tjms:tjme),intent(out)::tend  ! tendency (rhs of the level set pde)
real,dimension(ims:ime,jms:jme),intent(out)::ros  ! rate of spread 
type(fire_params),intent(in)::fp

!*** local 
real:: te,diffLx,diffLy,diffRx,diffRy, & 
   diffCx,diffCy,diff2x,diff2y,grad,rr, &
   ros_back,ros_wind,ros_slope,advx,advy,scale,nvx,nvy,speed
integer::i,j,itso,iteo,jtso,jteo,ierrx,nerr
integer::upwind_err=0
character(len=128)msg,msg2

! constants
real, parameter:: eps=epsilon(0.0)
!intrinsic epsilon
real, parameter:: zero=0.,one=1.,tol=100*eps, &
    safe=2.,rmin=safe*tiny(zero),rmax=huge(zero)/safe


! f90 intrinsic function

intrinsic max,min,sqrt,nint,tiny,huge

!*** executable
    
    ! check array dimensions
    call check_mesh_2dim(ints-1,inte+1,jnts-1,jnte+1,lims,lime,ljms,ljme)
    call check_mesh_2dim(ints,inte,jnts,jnte,tims,time,tjms,tjme)
    
    call continue_at_boundary(1,1,fire_lfn_ext_up, &   !extend by extrapolation but never down 
    lims,lime,ljms,ljme, &                ! memory dims
    ids,ide,jds,jde, &                    ! domain - nodes where lfn defined
    ips,ipe,jps,jpe, &                    ! patch - nodes owned by this process 
    ints,inte,jnts,jnte, &                ! tile - nodes update by this thread
    itso,iteo,jtso,jteo, &                ! where set now
    lfn)                                  ! array

    call print_2d_stats(itso,iteo,jtso,jteo,lims,lime,ljms,ljme, &
                   lfn,'tend_ls: lfn cont')

#ifdef DEBUG_OUT
    call write_array_m(ints-1,inte+1,jnts-1,jnte+1,lims,lime,ljms,ljme,lfn,'tend_lfn_in',id)
#endif
    
    nerr=0
    tbound=0    
    do j=jnts,jnte
        do i=ints,inte
            ! one sided differences
            diffRx = (lfn(i+1,j)-lfn(i,j))/dx
            diffLx = (lfn(i,j)-lfn(i-1,j))/dx
            diffRy = (lfn(i,j+1)-lfn(i,j))/dy
            diffLy = (lfn(i,j)-lfn(i,j-1))/dy
            diffCx = diffLx+diffRx   ! 	TWICE CENTRAL DIFFERENCE
            diffCy = diffLy+diffRy
    
            !upwinding - select right or left derivative
            select case(sfire_upwinding)
            case(0)  ! none
                grad=sqrt(diffCx**2 + diffCy**2)
            case(1) ! standard
                diff2x=select_upwind(diffLx,diffRx)
                diff2y=select_upwind(diffLy,diffRy)
                grad=sqrt(diff2x*diff2x + diff2y*diff2y)
            case(2) ! godunov per osher/fedkiw
                diff2x=select_godunov(diffLx,diffRx)
                diff2y=select_godunov(diffLy,diffRy)
                grad=sqrt(diff2x*diff2x + diff2y*diff2y)
            case(3) ! eno
                diff2x=select_eno(diffLx,diffRx)
                diff2y=select_eno(diffLy,diffRy)
                grad=sqrt(diff2x*diff2x + diff2y*diff2y)
            case(4) ! Sethian - twice stronger pushdown of bumps
                grad=sqrt(max(diffLx,0.)**2+min(diffRx,0.)**2   &
                        + max(diffLy,0.)**2+min(diffRy,0.)**2)
            case default
                grad=0.
                upwind_err=1
            end select
  
            ! normal direction, from central differences
            scale=sqrt(diffCx*diffCx+diffCy*diffCy+eps) 
            nvx=diffCx/scale
            nvy=diffCy/scale
                      
            ! wind speed in direction of spread
            ! speed =  vx(i,j)*nvx + vy(i,j)*nvy
        
    
            ! get rate of spread from wind speed and slope

            call fire_ros(ros_back,ros_wind,ros_slope, &
            nvx,nvy,i,j,fp,ierrx,msg2)
            nerr = nerr + ierrx

            rr=ros_back + ros_wind + ros_slope
            if(fire_grows_only.gt.0)rr=max(rr,0.)

            ! set ros for output
            if(i.ge.its.and.i.le.ite.and.j.ge.jts.and.j.le.jte)ros(i,j)=rr

            ! get rate of spread
            te = -rr*grad   ! normal term 

            ! cfl condition
            if (grad > 0.) then
                 tbound = max(tbound,rr*(abs(diff2x)/dx+abs(diff2y)/dy)/grad)
            endif

            ! add numerical viscosity
            te=te + fire_viscosity*abs(rr)*((diffRx-diffLx)+(diffRy-diffLy))

            tend(i,j)=te
        enddo
    enddo        

    if(upwind_err>0)call crash('sfire_upwinding value not supported, 3 recommended')

    if(nerr>0)then
!$OMP CRITICAL(SFIRE_CORE_CRIT)
        write(msg,'(a,i6,1x,a)')'tend_ls:',nerr,'messages in rate of spread computations. Last message:'
!$OMP END CRITICAL(SFIRE_CORE_CRIT)
        call warning(msg)
        call warning(msg2)
    endif

    call print_2d_stats(its,ite,jts,jte,ims,ime,jms,jme, ros,'tend_ls: ros')
    call print_2d_stats(ints,inte,jnts,jnte,tims,time,tjms,tjme, &
                   tend,'tend_ls: tend out')

    ! the final CFL bound
    tbound = 1/(tbound+tol)

end subroutine tend_ls

!
!**************************
!

real function select_upwind(diffLx,diffRx)
implicit none
real, intent(in):: diffLx, diffRx
real diff2x

! upwind differences, L or R if bith same sign, otherwise zero    

diff2x=0
if (diffLx>0.and.diffRx>0.)diff2x=diffLx
if (diffLx<0.and.diffRx<0.)diff2x=diffRx

select_upwind=diff2x
end function select_upwind


!
!**************************
!


real function select_godunov(diffLx,diffRx)
implicit none
real, intent(in):: diffLx, diffRx
real diff2x,diffCx

! Godunov scheme: upwind differences, L or R or none    
! always test on > or < never = , much faster because of IEEE
! central diff >= 0 => take left diff if >0, ortherwise 0
! central diff <= 0 => take right diff if <0, ortherwise 0

diff2x=0
diffCx=diffRx+diffLx
if (diffLx>0.and..not.diffCx<0)diff2x=diffLx
if (diffRx<0.and.     diffCx<0)diff2x=diffRx

select_godunov=diff2x
end function select_godunov

!
!**************************
!

real function select_eno(diffLx,diffRx)
implicit none
real, intent(in):: diffLx, diffRx
real diff2x

! 1st order ENO scheme

! diffLx <= 0 and diffRx <= 0: diffRx
! diffLx >= 0 and diffRx >= 0: diffLx
! diffLx >= 0 and diffRx <= 0:
!        diffRx + diffLx <= 0: diffRx
!        diffRx + diffLx >  0: diffLx


if    (.not.diffLx>0 .and. .not.diffRx>0)then
    diff2x=diffRx
elseif(.not.diffLx<0 .and. .not.diffRx<0)then
    diff2x=diffLx
elseif(.not.diffLx<0 .and. .not.diffRx>0)then
    if(.not. abs(diffRx) < abs(diffLx))then
        diff2x=diffRx
    else
        diff2x=diffLx
    endif
else
    diff2x=0.
endif

select_eno=diff2x
end function select_eno
      
!
!**************************
!

real function speed_func(diffCx,diffCy,dx,dy,i,j,fp,ierrx,msg)
!*** purpose
!    the level set method speed function
implicit none
!*** arguments
real, intent(in)::diffCx,diffCy  ! x and y coordinates of the direction of propagation
real, intent(in)::dx,dy  ! x and y coordinates of the direction of propagation
integer, intent(in)::i,j         ! indices of the node to compute the speed at
type(fire_params),intent(in)::fp
integer, intent(out)::ierrx
character(len=*), intent(out)::msg
!*** local
real::scale,nvx,nvy,r
real::ros_back , ros_wind , ros_slope
real, parameter:: eps=epsilon(0.0)
!*** executable
            ! normal direction, from central differences
            scale=sqrt(diffCx*diffCx+diffCy*diffCy+eps) 
            nvx=diffCx/scale
            nvy=diffCy/scale
                      
            ! get rate of spread from wind speed and slope

            call fire_ros(ros_back,ros_wind,ros_slope, &
            nvx,nvy,i,j,fp,ierrx,msg)

            r=ros_back + ros_wind + ros_slope
            if(fire_grows_only.gt.0)r=max(r,0.)
            speed_func=r

end function speed_func

end module module_fr_sfire_core
