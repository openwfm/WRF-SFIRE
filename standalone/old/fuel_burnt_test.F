program fuel_burnt_test
use module_fr_sfire_core
use module_fr_sfire_phys
implicit none

! ims,ime,jms,jme   in   memory dimensions
! its,ite,jts,jte   in   tile dimensions (cells where fuel_frac computed)
! ifs,ife,jfs,jfe   in   fuel_frac memory dimensions
! ir,jr,            in   refinement - quadrature cells per fire cell, even, at least 2
! lfn               in   level function, at nodes at midpoints of cells
! tign              in   ignition time, at nodes at nodes at midpoints of cells
! fuel_time         in   time constant of fuel, per cell
! time_now          in   time now
! fuel_frac         out  fraction of fuel remaining, per cell

!*** local
integer:: ids,ide,jds,jde,ims,ime,jms,jme,ir,jr,k
parameter(ims=1,ime=2,jms=1,jme=2,ir=2,jr=2)
real:: time_now,off,fuel_left_f
real,dimension(ims:ime,jms:jme)::lff,tif,lfn0,fuel_time
integer::i,iu

lfn0(1,1)=0.419279 
lfn0(1,2)=-1.9766893E-02
lfn0(2,1)=6.310914  
lfn0(2,2)=5.983462
time_now=2.0
fuel_time(1,1)= 8.23529
iu=1
open(iu,file='tmp.txt')
do k=1,23

  off = (k-1)*0.05
  lff(1,1)=lfn0(1,1)-off
  lff(1,2)=lfn0(1,2)-off
  lff(2,1)=lfn0(2,1)-off
  lff(2,2)=lfn0(2,2)-off

  tif(1,1)=2.000000
  tif(1,2)=1.565244
  tif(2,1)=2.000000
  tif(2,2)=2.000000

  fuel_left_f=fuel_left_cell_1( &
                 lff(1,1),lff(1,2),lff(2,1),lff(2,2),&
                 tif(1,1),tif(1,2),tif(2,1),tif(2,2),&
                 time_now, fuel_time(1,1)) 
                                              ! premultiply by weight
  fuel_left_f=1-fuel_left_f                                            
  print *,'fuel_left_f=',fuel_left_f 
  write(iu,*) fuel_left_f 
end do
close(iu)   
end program 

