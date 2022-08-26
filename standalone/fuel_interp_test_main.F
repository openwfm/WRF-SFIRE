	program test
	use module_fr_sfire_core
        implicit none
real :: time_now
integer :: ims,ime,jms,jme,icl,jcl,i,j,num_tests,flag
real,dimension(3,3)::tff,lff,lfn,tign,real_tff,real_lff
real,dimension(3,3)::tff0,lff0,lfn0,tign0,real_tff0,real_lff0
real,dimension(3,3)::tff1,lff1,lfn1,tign1,real_tff1,real_lff1
real,dimension(3,3)::tff2,lff2,lfn2,tign2,real_tff2,real_lff2
real,dimension(3,3)::tff3,lff3,lfn3,tign3,real_tff3,real_lff3
real,dimension(3,3)::tff8,lff8,lfn8,tign8
real,dimension(3,5)::tff9,lff9,lfn9,tign9
real, dimension(5,5)::d
real,dimension(2,2)::result
real::test_err,error
real::lfn00,lfn01,lfn10,lfn11,tign00,tign01,tign10,tign11,fuel_time_cell
real:: result1,Result2,Result3,Result4,Result5,result6
real:: res,res1,res2,res3,res4!
real::fire_area
time_now=3.
101 format(a/(3f7.3))
ims=1
ime=3
jms=1
jme=3
icl=2
jcl=2
error=0.0001
num_tests=0
fire_area=0;
!!!!! Test #1
!lfn - plane Z=1.5-X;
! tign -- everywhere where lfn>=0, tign=3;
! everywhere where lfn<0, tign=Z+3;
data lfn/0.5, -0.5, -1.5, &
         0.5, -0.5,-1.5, &
         0.5, -0.5,-1.5/

data tign/3.0, 2.5, 1.5, &
          3.0, 2.5, 1.5, &
          3.0, 2.5, 1.5 /

data real_tff/3.0, 2.5, 2.0, &
          3.0, 2.5, 2.0, &
          3.0, 2.5, 2.0 /

data real_lff/0.0, -0.5, -1.0, &
          0.0, -0.5, -1.0, &
          0.0, -0.5, -1.0 /

        call tign_lfn_interpolation(time_now,icl,jcl,ims,ime,jms,jme, &
                                    tign,lfn,tff,lff)
test_err=1
	test_err=maxval(abs(real_tff-tff))
if (test_err>error) then
	print*,'Test #1 is not working, error ',test_err
        print*,'is greater than 0.0001'
else
num_tests=num_tests+1;
endif
!!!! Test #2
!lfn - plane Z=4-X-Y;
! tign -- everywhere where lfn>=0, tign=3;
! everywhere where lfn<0, tign=Z+3;
data lfn1/2.0, 1.0, 0.0, &
         1.0, 0.0,-1.0, &
         0.0, -1.0,-2.0/

data tign1/3.0, 3.0, 3.0, &
          3.0, 3.0, 2.0, &
          3.0, 2.0, 1.0 /

data real_tff1/3.0, 3.0, 3.0, &
          3.0, 3.0, 2.5, &
          3.0, 2.5, 2.0 /

data real_lff1/1.0, 0.5, 0.0, &
          0.5, 0.0, -0.5, &
          0.0, -0.5, -1.0 /

        call tign_lfn_interpolation(time_now,icl,jcl,ims,ime,jms,jme, &
                                    tign1,lfn1,tff1,lff1)

test_err=1
	test_err=maxval(abs(real_tff1-tff1))

if (test_err>error) then
	print*,'Test #2 is not working, error ',test_err
        print*,'is greater than 0.0001'
else
num_tests=num_tests+1;
endif
!!!! Test #3

!=======
!lfn - plane Z=3-X-Y;
! tign -- everywhere where lfn>=0, tign=3;
! everywhere where lfn<0, tign=Z+3;
data lfn2/1.0, 0.0, -1.0, &
         0.0, -1.0,-2.0, &
         -1.0, -2.0,-3.0/

data tign2/3.0, 3.0, 2.0, &
          3.0, 2.0, 1.0, &
          2.0, 1.0, 0.0 /

data real_tff2/3.0, 2.5, 2.0, &
          2.5, 2.0, 1.5, &
          2.0, 1.5, 1.0 /

data real_lff2/0.0, -0.5, -1.0, &
          -0.5, -1.0, -1.5, &
          -1.0, -1.5, -2.0 /



        call tign_lfn_interpolation(time_now,icl,jcl,ims,ime,jms,jme, &
                                    tign2,lfn2,tff2,lff2)

test_err=1
test_err=maxval(abs(real_tff2-tff2))
if (test_err>error) then
	print*,'Test #3 is not working, error ',test_err
        print*,'is greater than 0.0001'
else
num_tests=num_tests+1;
endif
!!!! Test #4

!=======
!lfn - plane Z=5-X-Y;
! tign -- everywhere where lfn>=0, tign=3;
! everywhere where lfn<0, tign=Z+3;
data lfn3/3.0, 2.0, 1.0, &
         2.0, 1.0, 0.0, &
         1.0, 0.0,-1.0/

data tign3/3.0, 3.0, 3.0, &
          3.0, 3.0, 3.0, &
          3.0, 3.0, 2.0 /

data real_tff3/3.0, 3.0, 3.0, &
          3.0, 3.0, 3.0, &
          3.0, 3.0, 3.0 /

data real_lff3/2.0, 1.5, 1.0, &
          1.5, 1.0, 0.5, &
          1.0, 0.5, 0.0 /



        call tign_lfn_interpolation(time_now,icl,jcl,ims,ime,jms,jme, &
                                    tign3,lfn3,tff3,lff3)
test_err=1
test_err=maxval(abs(real_tff3-tff3))
if (test_err>error) then
	print*,'Test #4 is not working, error ',test_err
        print*,'is greater than 0.0001'
else
num_tests=num_tests+1;
endif
print*,''
!!!Tests for fuel_left calculation!!!
!!!Test #5, result in Matlab u=(-2;0;1) f =0.91691!
res=0.9169
lfn00=-1.
lfn01=-1.
lfn10=1.
lfn11=1.
tign00=1.
tign01=1.
tign10=2.
tign11=2.
time_now=2
fuel_time_cell= 8.235294 ;

call fuel_left_cell_3(  &
    result1,fire_area,  &
    lfn00,lfn01,lfn10,lfn11, &
    tign00,tign01,tign10,tign11,&
    time_now, fuel_time_cell)

test_err=1
test_err=res-result1

if (test_err>error) then
	print*,'Test #5 is not working, error ',test_err
        print*,'is greater than 0.0001'
else
num_tests=num_tests+1;
endif



!Test#6, result in Matlab u=(1;0;1) f =0.1660
res=0.166
lfn00=-2.
lfn01=-2.
lfn10=-1.
lfn11=-1.
tign00=1.
tign01=1.
tign10=2.
tign11=2.
time_now=3
fuel_time_cell= 8.235294 ;

call fuel_left_cell_3(  &
    result1,fire_area,  & 
    lfn00,lfn01,lfn10,lfn11, &
    tign00,tign01,tign10,tign11,&
    time_now, fuel_time_cell)

test_err=1
test_err=res-result1

if (test_err>error) then
	print*,'Test #6 is not working, error ',test_err
        print*,'is greater than 0.0001'
else
num_tests=num_tests+1;
endif

!Test7, result in Matlab u=-(2;0;4) f =0.4538
res=0.4538
lfn00=-2.
lfn01=-2.
lfn10=-4.
lfn11=-4.
tign00=6.
tign01=6.
tign10=4.
tign11=4.
time_now=10
fuel_time_cell= 8.235294 ;

call fuel_left_cell_3(  &
    result1,fire_area,  &
    lfn00,lfn01,lfn10,lfn11, &
    tign00,tign01,tign10,tign11,&
    time_now, fuel_time_cell)

test_err=1
test_err=res-result1
if (test_err>error) then
	print*,'Test #7 is not working, error ',test_err
        print*,'is greater than 0.0001'
else
num_tests=num_tests+1;
endif


!Test#8, result in Matlab u=(-4;-4;6) f =0.5336
res=0.5336
lfn00=-3.
lfn01=-1.
lfn10=-1.
lfn11=1.
tign00=4.
tign01=8.
tign10=8.
tign11=10.
time_now=10
fuel_time_cell= 8.235294 ;

call fuel_left_cell_3(  &
    result1,fire_area,  &
    lfn00,lfn01,lfn10,lfn11, &
    tign00,tign01,tign10,tign11,&
    time_now, fuel_time_cell)

test_err=1
test_err=res-result1
if (test_err>error) then
	print*,'Test #8 is not working, error ',test_err
        print*,'is greater than 0.0001'
else
num_tests=num_tests+1;
endif

!Test9, result in Matlab u=(-2;-2;1) f =0.9692
res=0.9692
lfn00=-1.
lfn01=1.
lfn10=1.
lfn11=3.
tign00=1.
tign01=2.
tign10=2.
tign11=2.
time_now=2
fuel_time_cell= 8.235294 ;

call fuel_left_cell_3(  &
    result1,fire_area,  &
    lfn00,lfn01,lfn10,lfn11, &
    tign00,tign01,tign10,tign11,&
    time_now, fuel_time_cell)

test_err=1
test_err=res-result1
if (test_err>error) then
	print*,'Test #9 is not working, error ',test_err
        print*,'is greater than 0.0001'
else
num_tests=num_tests+1;
endif

print*,""

print*,"Main Test (#10), we calculate the fuel_left of the center cell"
print*,"(we have 9 cells initially)"
print*,"The initial values are"
print*,"lfn is described by the plane z=y-2; tign is described by the plane z=y"


time_now=2.

data lfn8/-1., 0., 1., &
         -1., 0.,1., &
         -1., 0.,1./

data tign8/1.0, 2., 2., &
          1.0, 2., 2., &
          1.0, 2., 2. /
res1=2.9752314E-02
res2=2.9752314E-02
res3=1.000
res4=1.000
print*,"time_now",time_now
write(*,101)'tign',tign8
write(*,101)'lfn',lfn8
print*,"icl,jcl should be equal 2 here"
write(*,*)'icl,jcl',icl,jcl
print*,"calling tign_lfn_interpolation, that calculates"
        call tign_lfn_interpolation(time_now,icl,jcl,ims,ime,jms,jme, &
                                    tign8,lfn8,tff8,lff8)

print*,"as result tff,lff are"
write(*,101)'tff',tff8
write(*,101)'lff',lff8

d=1e20
d(1:5:2,1:5:2)=tign8
d(2:4,2:4)=tff8
write(*,*)'tign and tff, here you can check whether they were  & 
           calculated properly'
write(*,'(5f7.3)')d
write(*,*)' '
d(1:5:2,1:5:2)=lfn8
d(2:4,2:4)=lff8
write(*,*)'lfn and lff, here you can check whether they were calculated properly'
write(*,'(5f7.3)')d
write(*,*)' '


fuel_time_cell= 8.235294 ;
print*,"Calculation of fuel_frac over 4 subcells"

do i=1,2  ! Loop over 4 subcells
   do j=1,2
print*,""
print*,"meaning of tff and lff in the subcell",i,j
write(*,*)'lff',lff8(i,j),lff8(i,j+1),lff8(i+1,j),lff8(i+1,j+1)
write(*,*)'tff',tff8(i,j),tff8(i,j+1),tff8(i+1,j),tff8(i+1,j+1)
call fuel_left_cell_3(  &
    result(i,j),fire_area,  &
    lff8(i,j),lff8(i,j+1),lff8(i+1,j),lff8(i+1,j+1), &
    tff8(i,j),tff8(i,j+1),tff8(i+1,j),tff8(i+1,j+1),&
    time_now, fuel_time_cell)
write(*,*)'result after step i and j',result(i,j),i,j
end do
end do

if (res1-result(1,1)>error.or.res2-result(1,2)>error.or.res3-result(2,1)>error.or.res4-result(2,2)>error) then
	print*,'Test #10 is not working, one of errors '
        print*,'For cell (1,1)',res1-result(1,1)
        print*,'For cell (1,2)',res2-result(1,2)
        print*,'For cell (2,1)',res3-result(2,1)
        print*,'For cell (2,2)',res4-result(2,2)
        print*,'is greater than 0.0001',error
else
num_tests=num_tests+1;
endif

print*,""
print*,""
!!!! Test #2

write(*,*)'MAIN TEST #2 (#11)'
print*,"Check for consistency, for 3 inner cells with the same results"
time_now=2.
ims=1
ime=5
jms=1
jme=5
res1=2.9752314E-02
res2=1.000
res3=2.9752314E-02
res4=1.000
flag=1
data lfn9/-1., -1.,-1., &
         -1., -1.,0., &
         0., 0.,0., &
         0., 1., 1., &
         1.,1.,1./

data tign9/1.0, 1., 1., &
          1.0, 1., 2., &
          2.0, 2., 2.,  &
          2.,2.,2., &
          2.,2.,2./

do icl=2,4
        write(*,*)"icl,jcl",icl,jcl
       call tign_lfn_interpolation(time_now,icl,jcl,ims,ime,jms,jme, &
                                    tign9,lfn9,tff9,lff9)

fuel_time_cell= 8.235294 ;
        write(*,*)"Calculation of fuel_frac over 4 subcells"
result6=0;
do i=1,2
   do j=1,2


call fuel_left_cell_3(  &
    result(i,j),fire_area,  &
    lff9(i,j),lff9(i,j+1),lff9(i+1,j),lff9(i+1,j+1), &
    tff9(i,j),tff9(i,j+1),tff9(i+1,j),tff9(i+1,j+1),&
    time_now, fuel_time_cell)
write(*,*)'result after step i and j',result(i,j),i,j

end do
end do

if (res1-result(1,1)>error.or.res2-result(1,2)>error.or.res3-result(2,1)>error.or.res4-result(2,2)>error) then
	print*,'Test #10 is not working, one of errors '
        print*,'For cell (1,1)',res1-result(1,1)
        print*,'For cell (1,2)',res2-result(1,2)
        print*,'For cell (2,1)',res3-result(2,1)
        print*,'For cell (2,2)',res4-result(2,2)
        print*,'is greater than 0.0001',error
        flag=0 
end if

end do
if (flag.eq.1) then 
    num_tests=num_tests+1;
endif
print*,"If the results for each correspondent subcell are equal in all cells"
print*,"then the code is consistent"

print*,num_tests,'of 11 tests worked succesfully'


 end program test
