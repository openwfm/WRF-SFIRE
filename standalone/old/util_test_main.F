module util_test
contains
subroutine util_test_sub
use module_fr_sfire_util
implicit none
call test(sum_2d_cells,'sum_2d_cells',0,0, &
    0,1,0,1,-1,2,0,2,      &
    1,1,1,1,0,3,0,3)
call testn2n(2,2,1,1,1,1,     &
    1,2,1,2,1,2,1,2,      &
    1,3,1,3,1,3,1,3)
call testn2n(2,2,1,1,1,1,     &
    1,2,1,3,1,2,1,3,      &
    1,3,1,5,1,3,1,5)
call testn2n(2,3,1,1,2,2,     &
    1,2,1,3,1,5,1,4,      &
    1,3,1,10,1,3,1,10)
call testn2n(2,3,1,1,2,2,     &
    1,2,1,3,-1,5,1,4,      &
    1,3,2,6,1,3,1,8)
!call test(interpolate_2d_cells2cells,'interpolate_2d_cells2cells',0,0,  &
!    1,2,1,2,1,2,1,2,      &
!    1,4,1,4,1,4,1,4)    
!call test(interpolate_2d_cells2cells,'interpolate_2d_cells2cells',0,0,  &
!    1,2,1,2,1,2,1,2,      &
!    1,6,1,6,1,6,1,6)    
!call test(interpolate_2d_cells2nodes,'interpolate_2d_cells2nodes',0,1,  &
!    1,2,1,2,1,2,1,2,      &
!    1,4,1,4,1,5,1,5)    
!call test(interpolate_2d_cells2nodes,'interpolate_2d_cells2nodes',0,1,  &
!    1,2,1,2,1,2,1,2,      &
!    1,6,1,6,1,7,1,7)    
end subroutine util_test_sub


subroutine test(sub,s,n1,n2,              &
       ids1,ide1,jds1,jde1,ims1,ime1,jms1,jme1, &
       ids2,ide2,jds2,jde2,ims2,ime2,jms2,jme2)
use module_fr_sfire_util
implicit none

external sub
character(len=*)s
integer:: n1,n2 ! 0 if cell based, 1 if node based
integer, intent(in)::ids1,ide1,jds1,jde1,ims1,ime1,jms1,jme1
real::v1(ims1:ime1,jms1:jme1)
integer, intent(in)::ids2,ide2,jds2,jde2,ims2,ime2,jms2,jme2
real::v2(ims2:ime2,jms2:jme2)
integer:: i,j,k

write(*,'(a,a)')'test ',s
write(*,2)'mesh 1 size',ids1,ide1,jds1,jde1
k=0
do i=ids1,ide1+n1
    do j=jds1,jde1+n1
       k=k+1
       v1(i,j)=k
    enddo
    write(*,1)(v1(i,j),j=jds1,jde1+n1)
enddo
1 format(20f7.3)
2 format(a,i4,':',i4,' by ',i4,':',i4)
call sub(      &
       ims1,ime1,jms1,jme1,ids1,ide1,jds1,jde1,v1, &
       ims2,ime2,jms2,jme2,ids2,ide2,jds2,jde2,v2)
write(*,2)'mesh 2, interpolated, size',ids2,ide2,jds2,jde2
do i=ids2,ide2+n2
    write(*,1)(v2(i,j),j=jds2,jde2+n2)
enddo
end subroutine test

subroutine testn2n(ir,jr,ip1,jp1,ip2,jp2,       &
       ids1,ide1,jds1,jde1,ims1,ime1,jms1,jme1, &
       ids2,ide2,jds2,jde2,ims2,ime2,jms2,jme2)
use module_fr_sfire_util
implicit none

integer, intent(in):: ir,jr ! refinement ratios
integer, intent(in):: ip1,jp1,ip2,jp2 ! offset of start of 2 in 1
integer, intent(in)::ids1,ide1,jds1,jde1,ims1,ime1,jms1,jme1
real::v1(ims1:ime1,jms1:jme1)
integer, intent(in)::ids2,ide2,jds2,jde2,ims2,ime2,jms2,jme2
real::v2(ims2:ime2,jms2:jme2)
integer:: i,j,k

write(*,'(a)')'test interpolate_2d_nodes2nodes'
write(*,2)'mesh 1 size',ids1,ide1,jds1,jde1
k=0
do i=ids1,ide1
    do j=jds1,jde1
       k=k+1
       v1(i,j)=k
    enddo
    write(*,1)(v1(i,j),j=jds1,jde1)
enddo
v2=-1
1 format(20f7.3)
2 format(a,i4,':',i4,' by ',i4,':',i4)

call interpolate_2d_nodes2nodes(  &
    ims1,ime1,jms1,jme1, & ! array coarse grid
    ids1,ide1,jds1,jde1, & ! dimensions fine grid
    ims2,ime2,jms2,jme2, & ! array coarse grid
    ids2,ide2,jds2,jde2, & ! dimensions coarse grid
    ir,jr,               & ! refinement ration
    ip1,jp1,ip2,jp2,     & ! 
    v1,                  & ! in coarse grid  
    v2  )                  ! out fine grid

write(*,2)'mesh 2, interpolated, size',ids2,ide2,jds2,jde2
do i=ids2,ide2
    write(*,1)(v2(i,j),j=jds2,jde2)
enddo
end subroutine testn2n

end module util_test

program main
use util_test
call util_test_sub
end program main
