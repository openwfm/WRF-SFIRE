program burn_test
use module_fr_sfire_core
implicit none
integer :: ids,ide,jds,jde,ims,ime,jms,jme
parameter(ids=1,ide=1,jds=1,jde=1,ims=1,ime=2,jms=1,jme=2)
real,  dimension(ims:ime,jms:jme)::lfn,tign
real,  dimension(ims:ime,jms:jme)::w
real  tnow
real,  dimension(ims:ime,jms:jme)::fuel_frac
integer ::k,i,j
parameter(k=3) ! number of test cases 
real,  dimension(ims:ime,jms:jme,k)::lfn_t,tign_t
real,  dimension(ims:ime,jms:jme,k)::frac_t
real,  dimension(ims:ime,jms:jme,k)::w_t
real,  dimension(k)::tnow_t
character(len=20),dimension(k)::desc

data desc(1)/'all burning'/
data ((lfn_t(i,j,1),i=1,2),j=1,2)/-1.0,-2.0,-3.0,-4.0/
data ((tign_t(i,j,1),i=1,2),j=1,2)/1.0,2.0,3.0,4.0/
data frac_t(1,1,1)/0.142829357916163/ !(exp(-1.0)+exp(-2.0)+exp(-3.0)+exp(-4.0))/4
data w_t(1,1,1)/1.0/
data tnow_t(1)/5.0/

data desc(2)/'nothing burning'/
data ((lfn_t(i,j,2),i=1,2),j=1,2)/1.0,2.0,3.0,4.0/
data ((tign_t(i,j,2),i=1,2),j=1,2)/1.0,2.0,3.0,4.0/
data w_t(1,1,2)/2.0/
data frac_t(1,1,2)/1.0/
data tnow_t(2)/1.0/

data desc(3)/'fireline on diagonal'/
data ((lfn_t(i,j,3),i=1,2),j=1,2)/-1.0,0.0,0.0,1.0/
data ((tign_t(i,j,3),i=1,2),j=1,2)/1.0,0.0,0.0,0.0/
data w_t(1,1,3)/2.0/
data frac_t(1,1,3)/0.894646573528574/  ! 0.5*(exp(-1)+2)/3+0.5
data tnow_t(3)/2.0/
! keep adding more tests here
! keep adding more tests here

do i=1,k
    lfn=lfn_t(:,:,i)
    tign=tign_t(:,:,i)
    w=w_t(1,1,i)
    tnow=tnow_t(i)
    call fuel_left(1,ids,ide,jds,jde,ims,ime,jms,jme,lfn,tign,w,tnow,fuel_frac)
    write(*,1) i, desc(i), fuel_frac(1,1), frac_t(1,1,i), (fuel_frac(1,1)-frac_t(1,1,i)) 
1   format( i3,1x, a,' fuel_frac', f8.5,' should be', f8.5,' error', g11.3)
enddo

end program burn_test
