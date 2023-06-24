    ! fake physics module - simple spread formula
    
    module module_fr_sfire_phys
    
    contains
    ! the speed function
    
    subroutine normal_spread(  &
        ims,ime,jms,jme, &  ! memory bounds
        its,ite,jts,jte, &  ! tile bounds
        nvx,nvy,         &  ! in: directions to compute the speed in
        speed            &  ! out: the spread speed
#ifdef SPEED_VARS_ARGS      /* extra arguments for normal_spread */
#include SPEED_VARS_ARGS
#endif
    )
    use module_fr_sfire_util
    implicit none
!   arguments
    integer, intent(in)::its,ite,jts,jte,ims,ime,jms,jme
    real,intent(in),dimension(ims:ime,jms:jme)::nvx,nvy
    real,intent(out),dimension(ims:ime,jms:jme)::speed
#ifdef SPEED_VARS_DECL      /* extra arguments for normal_spread */
#include SPEED_VARS_DECL
#endif
    real::e,c
    ! calls
    intrinsic nint
    !  internal
    integer i,j
 
 !  compute speed at all nodes
 !
 !  pull coeffiecients
 !
    e=normal_spread_e
    c=normal_spread_c
    write(*,*)'c=',normal_spread_c,' e=',normal_spread_e
    do j=jts,jte+1
        do i=its,ite+1
            speed(i,j) = r(i,j) + max(vx(i,j)*nvx(i,j) + vy(i,j)*nvy(i,j),0.0)**e * c
        enddo
    enddo
    
    end subroutine normal_spread

    end module module_fr_sfire_phys
