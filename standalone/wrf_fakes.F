! wrf_fakes.F
! things that would be normally done in WRF

module module_model_constants
! The model may not define any physical constants, all must be done here.
! copied what needed from WRFV3/share/module_model_constants.F
   REAL   , PARAMETER :: r_d          = 287.
   REAL   , PARAMETER :: cp           = 7.*r_d/2.   ! specific heat of the atmosphere
   REAL   , PARAMETER :: xlv          = 2.5E6       ! latent heat
   REAL   , PARAMETER ::  pi2=2.*3.1415926          ! 2*pi
   REAL   , PARAMETER :: reradius     = 1./6370.0e03 ! 1/earth radius
   REAL   , PARAMETER :: g            = 9.81        ! gravity acceleration
end module module_model_constants

module module_state_description
   integer, parameter:: num_tracer=1
   integer ::             p_tr17_1=1
   integer ::             p_tr17_2=1
   integer ::             p_tr17_3=1
   integer ::             p_tr17_4=1
   integer ::             p_tr17_5=1
   integer ::             p_tr17_6=1
   integer ::             p_tr17_7=1
   integer ::             p_tr17_8=1
end module module_state_description

!
!***
!

module module_wrf_error
implicit none
contains

! mock-up of various wrf utility functions

subroutine wrf_error_fatal(s)
!*** purpose: abort with a message
implicit none
character(len=*), intent(in)::s
call latch ! food for debugger so you can say "stop at latch" and not worry about module name
write(6,*)s
call abort()
end subroutine wrf_error_fatal

!
!***
!

subroutine wrf_debug(level,s)
!*** purpose: print a message

implicit none
character(len=*), intent(in)::s
integer, intent(in):: level
write(6,*)s
end subroutine wrf_debug

!
!***
!

subroutine wrf_message(s)
character(len=*), intent(in)::s
integer i
do i=len(s),2,-1
    select case(s(i:i))
    case(' ')
    case default
        goto 1
    end select
enddo
1 write(6,'(a)')s(1:i)
end subroutine wrf_message

end module module_wrf_error

!
!***
!

! just for testing
subroutine latch
end

subroutine wrf_error_fatal(s)
! ESMF has this subroutine outside the module for some reason
use module_wrf_error, ONLY : wrf_error_fatal_mod => wrf_error_fatal
implicit none
character(len=*),intent(in)::s
call wrf_error_fatal_mod(s)
end subroutine wrf_error_fatal

!
!*** various stubs, mostly doing nothing
!

module module_dm
use module_wrf_error
implicit none
contains
subroutine wrf_get_nproc (nprocs)
integer nprocs
nprocs=1
end subroutine wrf_get_nproc

!
!***
!

subroutine wrf_get_myproc( myproc )
integer myproc
myproc=1
end subroutine wrf_get_myproc

!
!***
!

subroutine wrf_dm_maxval_integer( val, idex, jdex )
integer::val
integer::idex,jdex
call wrf_error_fatal('wrf_dm_maxval_integer:not implemented')
end subroutine wrf_dm_maxval_integer

subroutine wrf_dm_sum_reals(t,s)
real, intent(in)::t(:)
real, intent(out)::s(:)
s=t
end subroutine wrf_dm_sum_reals

end module module_dm

!
!*** external stubs
!


SUBROUTINE wrf_dm_bcast_real( BUF, N1 )
   IMPLICIT NONE
   INTEGER n1
   REAL  buf(*)
   RETURN
END SUBROUTINE wrf_dm_bcast_real

!
!***
!

LOGICAL FUNCTION wrf_dm_on_monitor()
  wrf_dm_on_monitor = .true.
END FUNCTION wrf_dm_on_monitor

!
!***
!

SUBROUTINE wrf_dm_bcast_integer( BUF, N1 )
   IMPLICIT NONE
   INTEGER n1
   INTEGER  buf(*)
   RETURN
END SUBROUTINE wrf_dm_bcast_integer


!WRF:DRIVER_LAYER:UTIL
!

MODULE module_timing
   use module_wrf_error 

   INTEGER, PARAMETER, PRIVATE :: cnmax = 30
   INTEGER, PRIVATE, DIMENSION(cnmax) :: count_int1 , count_rate_int1 , count_max_int1
   INTEGER, PRIVATE, DIMENSION(cnmax) :: count_int2 , count_rate_int2 , count_max_int2
   INTEGER, PRIVATE :: cn = 0 
   REAL, PRIVATE    :: elapsed_seconds , elapsed_seconds_total = 0
   REAL, PRIVATE    :: cpu_1 , cpu_2 , cpu_seconds , cpu_seconds_total = 0

CONTAINS

   SUBROUTINE init_module_timing
      cn = 0
   END SUBROUTINE init_module_timing


   SUBROUTINE start_timing

      IMPLICIT NONE

      cn = cn + 1
      IF ( cn .gt. cnmax ) THEN
        CALL  wrf_error_fatal( 'module_timing: clock nesting error (too many nests)' )
        RETURN
      ENDIF
      CALL SYSTEM_CLOCK ( count_int1(cn) , count_rate_int1(cn) , count_max_int1(cn) )
!     CALL CPU_TIME ( cpu_1 )

   END SUBROUTINE start_timing


   SUBROUTINE end_timing ( string )
   
      IMPLICIT NONE

      CHARACTER *(*) :: string

      IF ( cn .lt. 1 ) THEN
        CALL  wrf_error_fatal( 'module_timing: clock nesting error, cn<1' ) 
      ELSE IF ( cn .gt. cnmax ) THEN
        CALL  wrf_error_fatal( 'module_timing: clock nesting error, cn>cnmax' ) 
      ENDIF

      CALL SYSTEM_CLOCK ( count_int2(cn) , count_rate_int2(cn) , count_max_int2(cn) )
!     CALL CPU_TIME ( cpu_2 )

      IF ( count_int2(cn) < count_int1(cn) ) THEN
         count_int2(cn) = count_int2(cn) + count_max_int2(cn)
      ENDIF

      count_int2(cn) = count_int2(cn) - count_int1(cn)
      elapsed_seconds = REAL(count_int2(cn)) / REAL(count_rate_int2(cn))
      elapsed_seconds_total = elapsed_seconds_total + elapsed_seconds

      WRITE(6,'(A,A,A,F10.5,A)') 'Timing for ',TRIM(string),': ',elapsed_seconds,' elapsed seconds.'
#if defined(DM_PARALLEL) && ! defined(STUBMPI)
      WRITE(0,'(A,A,A,F10.5,A)') 'Timing for ',TRIM(string),': ',elapsed_seconds,' elapsed seconds.'
#endif

!     cpu_seconds = cpu_2 - cpu_1
!     cpu_seconds_total = cpu_seconds_total + cpu_seconds
!     PRINT '(A,A,A,F10.5,A)' ,'Timing for ',TRIM(string),': ',cpu_seconds,' cpu seconds.'

      cn = cn - 1

   END SUBROUTINE end_timing

END MODULE module_timing

