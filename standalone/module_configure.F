module module_configure
use module_wrf_error, only: wrf_error_fatal

implicit none

private
public read_namelist,grid_config_rec_type,nmlfile


! type for all namelist fields, domain 1 only
integer,parameter::max_domains=10,max_moves=1,max_eta=501,max_bogus=5
type grid_config_rec_type
#include "namelist_defines2.inc"
end type

! namelist file name
character(len=*),parameter::nmlfile='namelist.input'

contains

subroutine read_namelist(config_flags)
implicit none

!*** arguments
type(grid_config_rec_type), intent(out)::config_flags

!*** local
character(len=32)::grp
integer::iou=1
#include "namelist_defines.inc"
#include "namelist_statements.inc"

integer::placeholder
namelist /standalone/ placeholder

!*** executable

#include "namelist_defaults.inc"

write(*,'(2a)')'reading flags from file ',trim(nmlfile)

open(iou,file=nmlfile,err=999)

grp='time_control'
read(iou,time_control,err=998)
rewind(iou)

grp='domains'
read(iou,domains,err=998)
rewind(iou)

grp='fire'
read(iou,fire,err=998)
rewind(iou)

grp='standalone'
read(iou,standalone,err=997,end=997)
close(iou)

997 continue
 
#define DEST_RECORD config_flags %
#define SOURCE_RECORD 
#define SOURCE_REC_DEX (1)
#include "config_assigns.inc"

return

998 continue
call wrf_error_fatal('error reading namelist group '//grp)

999 continue
call wrf_error_fatal('error opening namelist '//trim(nmlfile))

end subroutine read_namelist

end module module_configure
