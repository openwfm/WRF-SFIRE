!   test driver for subroutine prop_ls
 
    program prop_test
    use module_fr_sfire_core
    use module_fr_sfire_util
    
    implicit none
    integer::m1,m2,n1,n2,i,j,m1m,m2m,n1m,n2m,halo
    real::rm,rn,ts,te,dx,dy
    integer::num_tiles
    integer, dimension(100)::i_start,i_end,j_start,j_end
    real,allocatable,dimension(:,:)::lfn,tc,vx,vy,r
    real::normal_spread_c,normal_spread_e
    
    !external spread_rate
    open(1,form='formatted',file='prop_test_in.txt')
    read(1,1)rm,rn
    m1=1
    m2=nint(rm)
    n1=1
    n2=nint(rn)
    halo=1
    m1m=m1-halo      ! memory dims
    m2m=m2+halo
    n1m=n1-halo
    n2m=n2+halo
    

    print *,'prop_test_m:',m1,m2,n1,n2 
    allocate (lfn(m1m:m2m,n1m:n2m),&
              tc(m1m:m2m,n1m:n2m),&
              r(m1m:m2m,n1m:n2m),&
              vx(m1m:m2m,n1m:n2m),&
              vy(m1m:m2m,n1m:n2m))
    read(1,1)ts,te
    read(1,1)normal_spread_c,normal_spread_e
    print *,'prop_test_main: ts=',ts,' te=',te,' c=',normal_spread_c,' e=',normal_spread_e
    read(1,1)((r(i,j),i=m1,m2),j=n1,n2)
    !print *,'prop_test_m: read r'
    read(1,1)dx,dy
    read(1,1)((lfn(i,j),i=m1,m2),j=n1,n2)
    !print *,'prop_test_m: read lfn'
    read(1,1)((vx(i,j),i=m1,m2),j=n1,n2)
    !print *,'prop_test_m: read vx'
    read(1,1)((vy(i,j),i=m1,m2),j=n1,n2)
    !print *,'prop_test_m: read vy'
    !print *,'prop_test_m: calling prop_test'
    close(1)
    call set_tiles(1,1,m1,m2-1,n1,n2-1,num_tiles,i_start,i_end,j_start,j_end)
    call prop_ls( &
        m1m,m2m,n1m,n2m,&
        m1,m2-1,n1,n2-1, &
        num_tiles,i_start,i_end,j_start,j_end,&
        ts,te-ts,dx,dy,lfn,tc &
#       ifdef SPEED_VARS_ARGS  /* extra arguments for the speed function */
#           include SPEED_VARS_ARGS
#       endif        
    )
    open (2,form='formatted',file='prop_test_out.txt')
    write(2,1)((lfn(i,j),i=m1,m2),j=n1,n2)
    close(2)

1   format((e25.15)) 

    end program prop_test
