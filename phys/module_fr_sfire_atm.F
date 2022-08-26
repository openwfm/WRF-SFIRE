
!*** Jan Mandel August 2007 - February 2010
!*** email: Jan.Mandel@gmail.com

! Routines dealing with the atmosphere


module module_fr_sfire_atm

use module_model_constants, only: cp,xlv,g
use module_fr_sfire_phys, only: fire_wind_height,have_fuel_cats,fcz0,fcwh,nfuelcats,no_fuel_cat,no_fuel_cat2,windrf,itree,fueldepthm,mfuelcats
use module_fr_sfire_util
USE module_dm , only : wrf_dm_sum_reals


use module_fr_sfire_phys, only:  mfuelcats, nfuelcats    ! for emissions
USE module_state_description, only: num_tracer
use module_fr_sfire_phys, only: fuel_name

#ifdef WRF_CHEM
USE module_state_description, only: num_chem
USE module_configure, only: &
p_co,  &
p_ch4, &
p_h2, &
p_no,  &
p_no2,  &
p_so2,  &
p_nh3,  &
p_p25,  &
p_p25i, &
p_p25j, &
p_oc1,  &
p_oc2,  &
p_bc1,  &
p_bc2,  &
p_ald,  &
p_csl,  &
p_eth,  &
p_hc3,  &
p_hc5,  &
p_hcho,  &
p_iso,  &
p_ket,  &
p_mgly,  &
p_ol2,  &
p_olt, &
p_oli, &
p_ora2,&
p_tol, &
p_xyl, &
p_bigalk, &
p_bigene, &
p_c10h16, &
p_c2h4, &
p_c2h5oh, &
p_c2h6, &
p_c3h6, &
p_c3h8, &
p_ch3cooh, &
p_ch3oh, &
p_cres, &
p_glyald, &
! p_hyac, &
p_isopr, &
p_macr, &
p_mek, &
p_mvk, &
p_smoke, &   ! tracer smoke exists only with CHEM
p_sulf, &
p_dms, &
p_msa, &
p_dust_1, &
p_dust_2, &
p_dust_3, &
p_dust_4, &
p_dust_5, &
p_seas_1, &
p_seas_2, &
p_seas_3, &
p_seas_4, &
p_p10 
#endif

USE module_state_description, only: &
p_tr17_1, &
p_tr17_2, &
p_tr17_3, &
p_tr17_4, &
p_tr17_5, &
p_tr17_6, &
p_tr17_7, &
p_tr17_8

implicit none

#ifndef WRF_CHEM
integer, parameter, private:: num_chem=0
#endif

logical, save :: have_wind_log_interpolation = .false.  ! status

! emission tables
REAL, dimension(mfuelcats), save:: &
co=0.,  &
ch4=0., &
h2=0., &
no=0.,  &
no2=0.,  &
so2=0.,  &
nh3=0.,  &
oc1=0.,  &
oc2=0.,  &
bc1=0.,  &
bc2=0.,  &
ald=0.,  &
csl=0.,  &
eth=0.,  &
hc3=0.,  &
p25=0., &
p25i=0., &
p25j=0., &
hc5=0.,  &
hcho=0.,  &
iso=0.,  &
ket=0.,  &
mgly=0.,  &
ol2=0., &
olt=0., &
oli=0., &
ora2=0.,&
tol=0., &
xyl=0., &
bigalk=0., &
bigene=0., &
c10h16=0., &
c2h4=0., &
c2h5oh=0., &
c2h6=0., &
c3h6=0., &
c3h8=0., &
ch3cooh=0., &
ch3oh=0., &
cres=0., &
glyald=0., &
! hyac=0., &
isopr=0., &
macr=0., &
mek=0., &
mvk=0., &
smoke=0., &
sulf=0., &
dms=0., &
msa=0., &
dust_1=0., &
dust_2=0., &
dust_3=0., &
dust_4=0., &
dust_5=0., &
seas_1=0., &
seas_2=0., &
seas_3=0., &
seas_4=0., &
p10=0., &
tr17_1=0., &
tr17_2=0., &
tr17_3=0., &
tr17_4=0., &
tr17_5=0., &
tr17_6=0., &
tr17_7=0., &
tr17_8=0.

real, parameter:: & ! reciprocal molecular weights  (mol/g)
             imw_co    = 1./28.010,&
             imw_ch4   = 1./16.04,&
             imw_h2    = 1./2.016,&
             imw_no    = 1./30.006,&
             imw_no2   = 1./46.006,& 
             imw_so2   = 1./64.066,&
             imw_nh3   = 1./17.031    

real, parameter:: mw_air=28.97   ! molecular weight of air g/mol

integer, parameter:: max_tracer=10, max_chem=64

! cumulative emissions - should be stored in the state, will be wrong on restarts
real, save, dimension(max_chem)::c_chem
real, save, dimension(mfuelcats)::c_fuel
real, save, dimension(max_tracer)::c_tracer

logical, save:: emis_read = .false.
integer, save:: msglevel =1,printsums=0 ! when to print sums
integer, parameter:: line=5 ! number of species, emissions, etc. per line

contains

! chem arrays are chem tracer
! indices p_species are generated in inc/scalar_indices.inc and included in frame/module_configure.F

subroutine read_emissions_table(chem_opt,tracer_opt)
   implicit none
   integer, intent(in)::chem_opt,tracer_opt
   logical, external:: wrf_dm_on_monitor
   external::wrf_dm_bcast_integer , wrf_dm_bcast_real
   integer, dimension(10)::compatible_chem_opt
   integer:: iounit,ierr,i
   character(len=128)::msg
   namelist/emissions/ compatible_chem_opt, printsums, &
co,  &
ch4, &
h2, &
no,  &
no2,  &
so2,  &
nh3,  &
p25, &
p25i, &
p25j, &
oc1,  &
oc2,  &
bc1,  &
bc2,  &
ald,  &
csl,  &
eth,  &
hc3,  &
hc5,  &
hcho,  &
iso,  &
ket,  &
mgly,  &
ol2, &
olt, &
oli, &
ora2,&
tol, &
xyl, &
bigalk, &
bigene, &
c10h16, &
c2h4, &
c2h5oh, &
c2h6, &
c3h6, &
c3h8, &
ch3cooh, &
ch3oh, &
cres, &
glyald, &
! hyac, &
isopr, &
macr, &
mek, &
mvk, &
smoke, &
sulf, &
dms, &
msa, &
dust_1, &
dust_2, &
dust_3, &
dust_4, &
dust_5, &
seas_1, &
seas_2, &
seas_3, &
seas_4, &
p10, &
tr17_1, &
tr17_2, &
tr17_3, &
tr17_4, &
tr17_5, &
tr17_6, &
tr17_7, &
tr17_8

!$  if (OMP_GET_THREAD_NUM() .ne. 0)then
!$     call crash('read_emissions_table: must be called from master thread')
!$  endif

IF ( wrf_dm_on_monitor() ) THEN
    ! we are the master task

    iounit=open_text_file('namelist.fire_emissions','read')
    compatible_chem_opt=0
    read(iounit,emissions,iostat=ierr)
    if(ierr.ne.0)call crash('read_emissions_table: error reading namelist emissions in file namelist.fire_emissions')
    CLOSE(iounit)
    write(msg,'(a,i3,a)')'reading emissions table for',nfuelcats,' fuel categories'
    call message(msg,level=0)
    if (.not.any(compatible_chem_opt.eq.chem_opt))then
        write(msg,'(a,i4,a)')'read_emissions_table: chem_opt=',chem_opt,' not between given compatible_chem_opt in namelist.fire_emissions'
        call  message(msg,level=0)
        write(msg,'(a,10i4)')'compatible_chem_opt=', compatible_chem_opt
        call  message(msg,level=0)
        call crash('chem_opt in namelist.input is not consistent with namelist.fire_emissions')
    endif
ENDIF
call wrf_dm_bcast_integer(printsums, 1)
call wrf_dm_bcast_real(co,  nfuelcats)
call wrf_dm_bcast_real(ch4,  nfuelcats)
call wrf_dm_bcast_real(h2,  nfuelcats)
call wrf_dm_bcast_real(no,  nfuelcats)
call wrf_dm_bcast_real(no2,  nfuelcats)
call wrf_dm_bcast_real(so2,  nfuelcats)
call wrf_dm_bcast_real(nh3,  nfuelcats)
call wrf_dm_bcast_real(p25,  nfuelcats)
call wrf_dm_bcast_real(p25i,  nfuelcats)
call wrf_dm_bcast_real(p25j,  nfuelcats)
call wrf_dm_bcast_real(oc1,  nfuelcats)
call wrf_dm_bcast_real(oc2,  nfuelcats)
call wrf_dm_bcast_real(bc1,  nfuelcats)
call wrf_dm_bcast_real(bc2,  nfuelcats)
call wrf_dm_bcast_real(ald,  nfuelcats)
call wrf_dm_bcast_real(csl,  nfuelcats)
call wrf_dm_bcast_real(eth,  nfuelcats)
call wrf_dm_bcast_real(hc3,  nfuelcats)
call wrf_dm_bcast_real(hc5,  nfuelcats)
call wrf_dm_bcast_real(hcho,  nfuelcats)
call wrf_dm_bcast_real(iso,  nfuelcats)
call wrf_dm_bcast_real(ket,  nfuelcats)
call wrf_dm_bcast_real(mgly,  nfuelcats)
call wrf_dm_bcast_real(ol2, nfuelcats)
call wrf_dm_bcast_real(olt, nfuelcats)
call wrf_dm_bcast_real(oli, nfuelcats)
call wrf_dm_bcast_real(ora2,nfuelcats)
call wrf_dm_bcast_real(tol, nfuelcats)
call wrf_dm_bcast_real(xyl, nfuelcats)
call wrf_dm_bcast_real(bigalk, nfuelcats)
call wrf_dm_bcast_real(bigene, nfuelcats)
call wrf_dm_bcast_real(c10h16, nfuelcats)
call wrf_dm_bcast_real(c2h4, nfuelcats)
call wrf_dm_bcast_real(c2h5oh, nfuelcats)
call wrf_dm_bcast_real(c2h6, nfuelcats)
call wrf_dm_bcast_real(c3h6, nfuelcats)
call wrf_dm_bcast_real(c3h8, nfuelcats)
call wrf_dm_bcast_real(ch3cooh, nfuelcats)
call wrf_dm_bcast_real(ch3oh, nfuelcats)
call wrf_dm_bcast_real(cres, nfuelcats)
call wrf_dm_bcast_real(glyald, nfuelcats)
! call wrf_dm_bcast_real(hyac, nfuelcats)
call wrf_dm_bcast_real(isopr, nfuelcats)
call wrf_dm_bcast_real(macr, nfuelcats)
call wrf_dm_bcast_real(mek, nfuelcats)
call wrf_dm_bcast_real(mvk, nfuelcats)
call wrf_dm_bcast_real(smoke, nfuelcats)
call wrf_dm_bcast_real(sulf, nfuelcats)
call wrf_dm_bcast_real(dms, nfuelcats)
call wrf_dm_bcast_real(msa, nfuelcats)
call wrf_dm_bcast_real(dust_1, nfuelcats)
call wrf_dm_bcast_real(dust_2, nfuelcats)
call wrf_dm_bcast_real(dust_3, nfuelcats)
call wrf_dm_bcast_real(dust_4, nfuelcats)
call wrf_dm_bcast_real(dust_5, nfuelcats)
call wrf_dm_bcast_real(seas_1, nfuelcats)
call wrf_dm_bcast_real(seas_2, nfuelcats)
call wrf_dm_bcast_real(seas_3, nfuelcats)
call wrf_dm_bcast_real(seas_4, nfuelcats)
call wrf_dm_bcast_real(p10, nfuelcats) 
call wrf_dm_bcast_real(tr17_1, nfuelcats)
call wrf_dm_bcast_real(tr17_2, nfuelcats)
call wrf_dm_bcast_real(tr17_3, nfuelcats)
call wrf_dm_bcast_real(tr17_4, nfuelcats)
call wrf_dm_bcast_real(tr17_5, nfuelcats)
call wrf_dm_bcast_real(tr17_6, nfuelcats)
call wrf_dm_bcast_real(tr17_7, nfuelcats)
call wrf_dm_bcast_real(tr17_8, nfuelcats)

if(fire_print_msg .ge. msglevel .and.printsums .gt. 0)then
   ! should be stored in the registry future because of restarts
   write(msg,'(3(a,i3,1x))')'c_chem size',num_chem,'c_tracer size',num_tracer,'c_fuel size',nfuelcats
   call message(msg,level=2)
   if(num_chem>max_chem)then
      call crash('read_emissions_table: increase max_chem')
   endif
   c_chem=0.  ! cumulative burnt
   if(num_tracer>max_tracer)then
      call crash('read_emissions_table: increase max_tracer')
   endif
   c_tracer=0.  ! cumulative burnt
   c_fuel=0.  ! total per timestep, rate burnt, cumulative burnt
endif

emis_read=.true.

end subroutine read_emissions_table

subroutine add_fire_emissions(chem_opt,tracer_opt,dt,dx,dy,      &
       ifms,ifme,jfms,jfme,    &
       ifts,ifte,jtfs,jfte,    &
       ids,ide,kds,kde,jds,jde,          &
       ims,ime,kms,kme,jms,jme,          &
       its,ite,kts,kte,jts,jte,          &
       rho,dz8w,                         & ! input on atmosphere mesh 
       fgip, fuel_frac_burnt, nfuel_cat, & ! input on fire mesh
       chem,tracer)                        ! output


implicit none

!*** purpose
! average fire emissions from fire mesh to coarser atmosphereic mesh and add to chemistry arrays

!*** arguments
! the dimensions are in cells, not nodes!

! input
integer, intent(in)::chem_opt,tracer_opt
real, intent(in):: dt,dx,dy                                           ! time step & mesh spacing
integer, intent(in)::its,ite,kts,kte,jts,jte,ims,ime,kms,kme,jms,jme &! atm grid dims
       ,ids,ide,kds,kde,jds,jde          
integer, intent(in)::ifts,ifte,jtfs,jfte,ifms,ifme,jfms,jfme          ! fire grid dims
real, intent(in)::rho(ims:ime,kms:kme,jms:jme),  &                ! air density  kg/m^3
                  dz8w(ims:ime,kms:kme,jms:jme)                       ! layer height
real, intent(in), dimension(ifms:ifme,jfms:jfme):: fgip, &            ! initial fuel load kg/m^2 
                                            fuel_frac_burnt, &        ! fuel fraction burned this step kg/kg
                                            nfuel_cat                 ! fuel category (Anderson= 1 to 13)
! update
real, intent(inout)::chem(ims:ime,kms:kme,jms:jme,num_chem),tracer(ims:ime,kms:kme,jms:jme,num_tracer)

!*** local
integer:: i,i_f,j,j_f,ir,jr,isz1,isz2,jsz1,jsz2,ioff,joff,ibase,jbase,cat,k1,areaw,m,k,k_p,errors
real::fuel_burnt,vol,air,conv, avgw,emis 
character(len=128)msg

real, dimension(mfuelcats)::s_fuel,t_fuel,r_fuel  ! total per timestep, rate burnt
#ifdef WRF_CHEM
integer, parameter:: chem_np=58
#else
integer, parameter:: chem_np=0
#endif
real, dimension(num_chem) ::s_chem,t_chem,r_chem  ! total per timestep, rate burnt
real, dimension(num_chem) ::a_chem,g_chem  ! concentration in ground level 1 
integer:: chem_pointers(chem_np)
character(len=8)::chem_names(chem_np)
integer, parameter:: tracer_np_max=8
integer:: tracer_pointers(tracer_np_max)

real, dimension(num_tracer) ::s_tracer,t_tracer,r_tracer  ! total per timestep, rate burnt
character(len=8)::tracer_names(tracer_np_max)
integer::tracer_np  ! actual number of tracer pointers

!*** executable

!check mesh dimensions and domain dimensions
call check_mesh_2dim(its,ite,jts,jte,ims,ime,jms,jme)
call check_mesh_2dim(ifts,ifte,jtfs,jfte,ifms,ifme,jfms,jfme)

if(.not.emis_read)call crash('add_fire_emissions: read_emissions_table must be called first') 

#ifdef WRF_CHEM
chem_pointers= (/ &
p_co,  &
p_ch4, &
p_h2, &
p_no,  &
p_no2,  &
p_so2,  &
p_nh3,  &
p_p25, &
p_p25i, &
p_p25j, &
p_oc1,  &
p_oc2,  &
p_bc1,  &
p_bc2,  &
p_ald,  &
p_csl,  &
p_eth,  &
p_hc3,  &
p_hc5,  &
p_hcho,  &
p_iso,  &
p_ket, &
p_mgly,  &
p_ol2,  &
p_olt, &
p_oli, &
p_ora2,&
p_tol, &
p_xyl, &
p_bigalk, &
p_bigene, &
p_c10h16, &
p_c2h4, &
p_c2h5oh, &
p_c2h6, &
p_c3h6, &
p_c3h8, &
p_ch3cooh, &
p_ch3oh, &
p_cres, &
p_glyald, &
! p_hyac, &
p_isopr, &
p_macr, &
p_mek, &
p_mvk, &
! p_smoke, &
p_sulf, &
p_dms, &
p_msa, &
p_dust_1, &
p_dust_2, &
p_dust_3, &
p_dust_4, &
p_dust_5, &
p_seas_1, &
p_seas_2, &
p_seas_3, &
p_seas_4, &
p_p10 /)

chem_names= (/ &
'co      ',  &
'ch4     ', &
'h2      ', &
'no      ',  &
'no2     ',  &
'so2     ',  &
'nh3     ',  &
'p25     ', &
'p25i    ', &
'p25j    ', &
'oc1     ',  &
'oc2     ',  &
'bc1     ',  &
'bc2     ',  &
'ald     ',  &
'csl     ',  &
'eth     ',  &
'hc3     ',  &
'hc5     ',  &
'hcho    ',  &
'iso     ',  &
'ket     ', &
'mgly    ',  &
'ol2     ',  &
'olt     ', &
'oli     ', &
'ora2    ',&
'tol     ', &
'xyl     ', &
'bigalk  ', &
'bigene  ', &
'c10h16  ', &
'c2h4    ', &
'c2h5oh  ', &
'c2h6    ', &
'c3h6    ', &
'c3h8    ', &
'ch3cooh ', &
'ch3oh   ', &
'cres    ', &
'glyald  ', &
! 'hyac     ', &
'isopr   ', &
'macr    ', &
'mek     ', &
'mvk     ', &
! 'smoke   ', &
'sulf    ', &
'dms     ', &
'msa     ', &
'dust_1  ', &
'dust_2  ', &
'dust_3  ', &
'dust_4  ', &
'dust_5  ', &
'seas_1  ', &
'seas_2  ', &
'seas_3  ', &
'seas_4  ', &
'p10     ' /)

#endif

if (tracer_opt == 0) then
    call message('add_fire_emissions: no tracers')
    tracer_np = 0
#ifdef WRF_CHEM
elseif (tracer_opt == 1) then
    tracer_pointers(1:1)= (/ p_smoke  /)
    tracer_names(1:1)= (/ 'smoke  ' /)
    tracer_np = 1
#endif
elseif (tracer_opt == 2) then
    tracer_pointers(1:8)= (/ &
       p_tr17_1, &
       p_tr17_2, &
       p_tr17_3, &
       p_tr17_4, &
       p_tr17_5, &
       p_tr17_6, &
       p_tr17_7, &
       p_tr17_8  /)
    tracer_names(1:8)= (/ &
       'tr17_1  ', &
       'tr17_2  ', &
       'tr17_3  ', &
       'tr17_4  ', &
       'tr17_5  ', &
       'tr17_6  ', &
       'tr17_7  ', &
       'tr17_8  ' /)
    tracer_np = 8
else
    call crash('add_fire_emissions: tracer_opt value not supported')
endif


if(fire_print_msg .ge. msglevel)then
!$OMP CRITICAL(SFIRE_ATM_CRIT)
900 format(6(a,"=",i3,1x))
    write(msg,900)'add_fire_emissions: chem_opt',chem_opt,'tracer_opt',tracer_opt, &
        'num_chem',num_chem,'num_tracer',num_tracer, &
#ifdef WRF_CHEM
        'max chem',maxval(chem_pointers), &
#endif
        'max tracer',maxval(tracer_pointers)
    call message(msg,level=msglevel)
    call print_dims_4d('tracer',tracer)
    call print_dims_4d('chem',chem)
!$OMP END CRITICAL(SFIRE_ATM_CRIT)
endif

#ifdef WRF_CHEM
    call check_pointers('chem',chem,chem_np,chem_names,chem_pointers)
#endif
if(tracer_opt>0)then
    call check_pointers('tracer',tracer,tracer_np,tracer_names,tracer_pointers)
endif

! compute mesh sizes
isz1 = ite-its+1
jsz1 = jte-jts+1
isz2 = ifte-ifts+1
jsz2 = jfte-jtfs+1


! check mesh sizes
if(isz1.le.0.or.jsz1.le.0.or.isz2.le.0.or.jsz2.le.0)then
    call message('all mesh sizes must be positive',level=0)
    goto 9
endif

! compute mesh ratios
ir=isz2/isz1
jr=jsz2/jsz1

if(isz2.ne.isz1*ir .or. jsz2.ne.jsz1*jr)then
    call message('input mesh size must be multiple of output mesh size',level=0)
    goto 9
endif

avgw = 1.0/(ir*jr)     ! averaging weight = 1/number of fire cells per atm cell

! initialize emissions statistics per timestep
#ifdef WRF_CHEM
t_chem = 0.
#endif
t_fuel = 0.
do i=1,num_tracer
   t_tracer(i)=0.
enddo

! conversion fuel_burnt kg/m^2 ->  chem_X 1e6*mol/mol
! emis_X(g/kg)*fuel_burnt(kg/m^2)/mw_X(g/mol) = emissions (mol/m^2)
! rho(kg/m^3)*dz8w(m))/mw_air(28.97e-3 kg/mol) = dry air in the 1st layer (mol/m^2)

k1 = kts 
!$OMP CRITICAL(SFIRE_ATM_CRIT)
write(msg,'(a,i3)')'Fire tracers into atmosphere level',k1
!$OMP END CRITICAL(SFIRE_ATM_CRIT)
call message(msg,level=msglevel)

#ifdef WRF_CHEM
if(fire_print_msg .ge. msglevel .and.printsums .gt. 0)then
    ! sum ground concentrations and check for nans
    a_chem=0.0
    g_chem=0.0
    errors=0
    do j=jts,jte
       do k=1,chem_np
          k_p=chem_pointers(k)
          do i=its,ite
             if(chem(i,k1,j,k_p ).ne.chem(i,k1,j,k_p ))errors=errors+1
             a_chem(k_p)=a_chem(k_p)+chem(i,k1,j,k_p )
          enddo
       enddo
    enddo
    if(errors>0)call crash('NaN before chem update')
    call wrf_dm_sum_reals(a_chem,g_chem)
    call message('Layer1 raw sums before adding fire emissions',level=msglevel)
    do i=1,chem_np,line
        m=min(i+line-1,chem_np)
        write(msg,80)'Emissions ',(trim(chem_names(j)),j=i,m)
        call message(msg,level=msglevel)
        write(msg,81)'Layer1 beg',(g_chem(chem_pointers(j)),j=i,m)
        call message(msg,level=msglevel)
        call message(' ',level=msglevel)
    enddo
endif
#endif

call message("Inserting chemical species",level=msglevel)
do j=max(jds+1,jts),min(jte,jde-1)          ! safe distance from domain boundary
    jbase=jtfs+jr*(j-jts)                   ! indexing
    do i=max(ids+1,its),min(ite,ide-1)
       ibase=ifts+ir*(i-its)                ! indexing
       !air = 1e6*mw_air/rho(i,kds,j)    ! 1e6*mw_air/air density 
       !vol = avgw/dz8w(i,kds,j)             ! averaging volume factor / 1st layer depth

       do joff=0,jr-1
           j_f=joff+jbase
           do ioff=0,ir-1
               i_f=ioff+ibase


               !*** fire cell (i_f,j_f) contributes to atmosphere cell (i,j) at ground level

               fuel_burnt = fgip(i_f,j_f) * fuel_frac_burnt(i_f,j_f) ! kg/m^2 
               cat = nfuel_cat(i_f, j_f)                       ! usually 1 to 13

               if(cat.lt.no_fuel_cat)t_fuel(cat)=t_fuel(cat) + fuel_burnt


               !*** chem compounds emissions given in g/kg

               ! fuel_burnt kg/m^2 *  table g/kg ->  ppmv = 1e6*mol/mol in 1st layer
!AK rho is in kg/m3 so the conversion factor must be scaled by a 1000 to match
!emissions in grams
!                conv = avgw*1e6*mw_air/(rho(i,k1,j)*dz8w(i,k1,j))
                 conv = avgw*1e3*mw_air/(rho(i,k1,j)*dz8w(i,k1,j))
#ifdef WRF_CHEM

               emis=co  (cat)*fuel_burnt                                    ! emission from fire cell in g/m^2
               t_chem(p_co) = t_chem(p_co) + emis                           ! add to total
               ! if(isnan(chem(i,k1,j,p_co  )))call crash('NaN before')
               chem(i,k1,j,p_co  )=chem(i,k1,j,p_co  ) + emis*conv*imw_co   ! add to chem
               ! if(isnan(chem(i,k1,j,p_co  )))call crash('NaN after')

               emis=ch4  (cat)*fuel_burnt                                    ! emission from fire cell in g/m^2
               t_chem(p_ch4) = t_chem(p_ch4) + emis                           ! add to total
               chem(i,k1,j,p_ch4  )=chem(i,k1,j,p_ch4  ) + emis*conv*imw_ch4   ! add to chem

               emis=h2  (cat)*fuel_burnt                                    ! emission from fire cell in g/m^2
               t_chem(p_h2) = t_chem(p_h2) + emis                           ! add to total
               chem(i,k1,j,p_h2  )=chem(i,k1,j,p_h2  ) + emis*conv*imw_h2   ! add to chem

               emis=no  (cat)*fuel_burnt                                    ! emission from fire cell in g/m^2
               t_chem(p_no) = t_chem(p_no) + emis                           ! add to total
               chem(i,k1,j,p_no  )=chem(i,k1,j,p_no  ) + emis*conv*imw_no   ! add to chem

               emis=no2 (cat)*fuel_burnt                                    ! emission from fire cell in g/m^2
               t_chem(p_no2) = t_chem(p_no2) + emis                         ! add to total
               chem(i,k1,j,p_no2 )=chem(i,k1,j,p_no2 ) + emis*conv*imw_no2  ! add to chem

               emis=so2 (cat)*fuel_burnt                                    ! emission from fire cell in g/m^2
               t_chem(p_so2) = t_chem(p_so2) + emis                         ! ads to total
               chem(i,k1,j,p_so2 )=chem(i,k1,j,p_so2 ) + emis*conv*imw_so2  ! add to chem

               emis=nh3 (cat)*fuel_burnt                                    ! emission from fire cell in g/m^2
               t_chem(p_nh3) = t_chem(p_nh3) + emis                         ! add to total
               chem(i,k1,j,p_nh3 )=chem(i,k1,j,p_nh3 ) + emis*conv*imw_nh3  ! add to chem


               !*** other emissions already given in mol/kg
               ! fuel_burnt kg/m^2 *  table mol/kg -> ppmv = 1e6*mol/mol in 1st layer dry air in 1st layer  
                
               ! same conversion factor but we will not divide by the molecular weight of the compound
               ! conv = avgw*1e6*mw_air/(rho(i,kds,j)*dz8w(i,kds,j)) 

               emis=ald (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_ald) = t_chem(p_ald) + emis                         ! add to total
               chem(i,k1,j,p_ald )=chem(i,k1,j,p_ald )+emis*conv

               emis=csl (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_csl) = t_chem(p_csl) + emis                         ! add to total
               chem(i,k1,j,p_csl )=chem(i,k1,j,p_csl )+emis*conv

               emis=eth (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_eth) = t_chem(p_eth) + emis                         ! add to total
               chem(i,k1,j,p_eth )=chem(i,k1,j,p_eth )+emis*conv

               emis=hc3 (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_hc3) = t_chem(p_hc3) + emis                         ! add to total
               chem(i,k1,j,p_hc3 )=chem(i,k1,j,p_hc3 )+emis*conv

               emis=hc5 (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_hc5) = t_chem(p_hc5) + emis                         ! add to total
               chem(i,k1,j,p_hc5 )=chem(i,k1,j,p_hc5 )+emis*conv

               emis=hcho (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_hcho) = t_chem(p_hcho) + emis                         ! add to total
               chem(i,k1,j,p_hcho)=chem(i,k1,j,p_hcho)+emis*conv

               emis=iso (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_iso) = t_chem(p_iso) + emis                         ! add to total
               chem(i,k1,j,p_iso )=chem(i,k1,j,p_iso )+emis*conv

               emis=ket (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_ket) = t_chem(p_ket) + emis                         ! add to total
               chem(i,k1,j,p_ket )=chem(i,k1,j,p_ket )+emis*conv

               emis=mgly (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_mgly) = t_chem(p_mgly) + emis                         ! add to total
               chem(i,k1,j,p_mgly)=chem(i,k1,j,p_mgly)+emis*conv

               emis=ol2 (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_ol2) = t_chem(p_ol2) + emis                         ! add to total
               chem(i,k1,j,p_ol2 )=chem(i,k1,j,p_ol2 )+emis*conv

               emis=olt (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_olt) = t_chem(p_olt) + emis                         ! add to total
               chem(i,k1,j,p_olt )=chem(i,k1,j,p_olt )+emis*conv

               emis=oli (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_oli) = t_chem(p_oli) + emis                         ! add to total
               chem(i,k1,j,p_oli )=chem(i,k1,j,p_oli )+emis*conv

               emis=ora2 (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_ora2) = t_chem(p_ora2) + emis                         ! add to total
               chem(i,k1,j,p_ora2)=chem(i,k1,j,p_ora2)+emis*conv

               emis=tol (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_tol) = t_chem(p_tol) + emis                         ! add to total
               chem(i,k1,j,p_tol )=chem(i,k1,j,p_tol )+emis*conv

               emis=xyl (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_xyl) = t_chem(p_xyl) + emis                         ! add to total
               chem(i,k1,j,p_xyl )=chem(i,k1,j,p_xyl )+emis*conv

               emis=bigalk (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_bigalk) = t_chem(p_bigalk) + emis                         ! add to total
               chem(i,k1,j,p_bigalk )=chem(i,k1,j,p_bigalk )+emis*conv

               emis=bigene (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_bigene) = t_chem(p_bigene) + emis                         ! add to total
               chem(i,k1,j,p_bigene )=chem(i,k1,j,p_bigene )+emis*conv

               emis=c10h16 (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_c10h16) = t_chem(p_c10h16) + emis                         ! add to total
               chem(i,k1,j,p_c10h16 )=chem(i,k1,j,p_c10h16 )+emis*conv

               emis=c2h4 (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_c2h4) = t_chem(p_c2h4) + emis                         ! add to total
               chem(i,k1,j,p_c2h4 )=chem(i,k1,j,p_c2h4 )+emis*conv

               emis=c2h5oh (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_c2h5oh) = t_chem(p_c2h5oh) + emis                         ! add to total
               chem(i,k1,j,p_c2h5oh )=chem(i,k1,j,p_c2h5oh )+emis*conv

               emis=c2h6 (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_c2h6) = t_chem(p_c2h6) + emis                         ! add to total
               chem(i,k1,j,p_c2h6 )=chem(i,k1,j,p_c2h6 )+emis*conv

               emis=c3h6 (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_c3h6) = t_chem(p_c3h6) + emis                         ! add to total
               chem(i,k1,j,p_c3h6 )=chem(i,k1,j,p_c3h6 )+emis*conv

               emis=c3h8 (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_c3h8) = t_chem(p_c3h8) + emis                         ! add to total
               chem(i,k1,j,p_c3h8 )=chem(i,k1,j,p_c3h8 )+emis*conv

               emis=ch3cooh (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_ch3cooh) = t_chem(p_ch3cooh) + emis                         ! add to total
               chem(i,k1,j,p_ch3cooh )=chem(i,k1,j,p_ch3cooh )+emis*conv

               emis=ch3oh (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_ch3oh) = t_chem(p_ch3oh) + emis                         ! add to total
               chem(i,k1,j,p_ch3oh )=chem(i,k1,j,p_ch3oh )+emis*conv

               emis=cres (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_cres) = t_chem(p_cres) + emis                         ! add to total
               chem(i,k1,j,p_cres )=chem(i,k1,j,p_cres )+emis*conv

               emis=glyald (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_glyald) = t_chem(p_glyald) + emis                         ! add to total
               chem(i,k1,j,p_glyald )=chem(i,k1,j,p_glyald )+emis*conv

               emis=isopr (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_isopr) = t_chem(p_isopr) + emis                         ! add to total
               chem(i,k1,j,p_isopr )=chem(i,k1,j,p_isopr )+emis*conv

               emis=macr (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_macr) = t_chem(p_macr) + emis                         ! add to total
               chem(i,k1,j,p_macr )=chem(i,k1,j,p_macr )+emis*conv

               emis=mek (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_mek) = t_chem(p_mek) + emis                         ! add to total
               chem(i,k1,j,p_mek )=chem(i,k1,j,p_mek )+emis*conv

               emis=mvk (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_mvk) = t_chem(p_mvk) + emis                         ! add to total
               chem(i,k1,j,p_mvk )=chem(i,k1,j,p_mvk )+emis*conv

               ! aerosols
               ! fuel_burnt kg/m^2 *  table g/kg -> ug/kg dry air in 1st layer  

               ! see also chem/emissions_driver.F line 515

               conv = avgw*1e6/(rho(i,k1,j)*dz8w(i,k1,j))

               emis=p25 (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_p25) = t_chem(p_p25) + emis                         ! add to total
               chem(i,k1,j,p_p25 )=chem(i,k1,j,p_p25 )+emis*conv

               emis=p25i (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_p25i) = t_chem(p_p25i) + emis                         ! add to total
               chem(i,k1,j,p_p25i )=chem(i,k1,j,p_p25i )+emis*conv

               emis=p25j (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_p25j) = t_chem(p_p25j) + emis                         ! add to total
               chem(i,k1,j,p_p25j )=chem(i,k1,j,p_p25j )+emis*conv

               emis=oc1  (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_oc1 ) = t_chem(p_oc1 ) + emis                         ! add to total
               chem(i,k1,j,p_oc1  )=chem(i,k1,j,p_oc1  )+emis*conv

               emis=oc2  (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_oc2 ) = t_chem(p_oc2 ) + emis                         ! add to total
               chem(i,k1,j,p_oc2  )=chem(i,k1,j,p_oc2  )+emis*conv

               emis=bc1  (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_bc1 ) = t_chem(p_bc1 ) + emis                         ! add to total
               chem(i,k1,j,p_bc1  )=chem(i,k1,j,p_bc1  )+emis*conv

               emis=bc2  (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_bc2 ) = t_chem(p_bc2 ) + emis                         ! add to total
               chem(i,k1,j,p_bc2  )=chem(i,k1,j,p_bc2  )+emis*conv

               emis=sulf  (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_sulf ) = t_chem(p_sulf ) + emis                         ! add to total
               chem(i,k1,j,p_sulf  )=chem(i,k1,j,p_sulf  )+emis*conv

               emis=dms  (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_dms ) = t_chem(p_dms ) + emis                         ! add to total
               chem(i,k1,j,p_dms  )=chem(i,k1,j,p_dms  )+emis*conv

               emis=msa  (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_msa ) = t_chem(p_msa ) + emis                         ! add to total
               chem(i,k1,j,p_msa  )=chem(i,k1,j,p_msa  )+emis*conv

               emis=dust_1  (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_dust_1 ) = t_chem(p_dust_1 ) + emis                         ! add to total
               chem(i,k1,j,p_dust_1  )=chem(i,k1,j,p_dust_1  )+emis*conv

               emis=dust_2  (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_dust_2 ) = t_chem(p_dust_2 ) + emis                         ! add to total
               chem(i,k1,j,p_dust_2  )=chem(i,k1,j,p_dust_2  )+emis*conv

               emis=dust_3  (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_dust_3 ) = t_chem(p_dust_3 ) + emis                         ! add to total
               chem(i,k1,j,p_dust_3  )=chem(i,k1,j,p_dust_3  )+emis*conv

               emis=dust_4 (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_dust_4 ) = t_chem(p_dust_4 ) + emis                         ! add to total
               chem(i,k1,j,p_dust_4  )=chem(i,k1,j,p_dust_4  )+emis*conv

               emis=dust_5  (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_dust_5 ) = t_chem(p_dust_5 ) + emis                         ! add to total
               chem(i,k1,j,p_dust_5  )=chem(i,k1,j,p_dust_5  )+emis*conv

               emis=seas_1  (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_seas_1 ) = t_chem(p_seas_1 ) + emis                         ! add to total
               chem(i,k1,j,p_seas_1  )=chem(i,k1,j,p_seas_1  )+emis*conv

               emis=seas_2  (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_seas_2 ) = t_chem(p_seas_2 ) + emis                         ! add to total
               chem(i,k1,j,p_seas_2  )=chem(i,k1,j,p_seas_2  )+emis*conv

               emis=seas_3  (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_seas_3 ) = t_chem(p_seas_3 ) + emis                         ! add to total
               chem(i,k1,j,p_seas_3  )=chem(i,k1,j,p_seas_3  )+emis*conv

               emis=seas_4  (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_seas_4 ) = t_chem(p_seas_4 ) + emis                         ! add to total
               chem(i,k1,j,p_seas_4  )=chem(i,k1,j,p_seas_4  )+emis*conv

               emis=p10  (cat)*fuel_burnt                                    ! emission from fire cell
               t_chem(p_p10 ) = t_chem(p_p10 ) + emis                         ! add to total
               chem(i,k1,j,p_p10  )=chem(i,k1,j,p_p10  )+emis*conv

#endif
               if (tracer_opt > 0) then

                     ! treat tracers exactly the same as aerosols, emissions g/kg fuel burned, tracer concentration ug/kg dry air
                     conv = avgw*1e6/(rho(i,k1,j)*dz8w(i,k1,j))
               endif

#ifdef WRF_CHEM
               if (tracer_opt == 1) then
                     emis=smoke (cat)*fuel_burnt                                    ! emission from fire cell
                     t_tracer(p_smoke) = t_tracer(p_smoke) + emis                         ! add to total
                     tracer(i,k1,j,p_smoke )=tracer(i,k1,j,p_smoke )+emis*conv
               endif
#endif
               if (tracer_opt == 2) then
                     emis=tr17_1 (cat)*fuel_burnt                                    ! emission from fire cell
                     t_tracer(p_tr17_1) = t_tracer(p_tr17_1) + emis                         ! add to total
                     tracer(i,k1,j,p_tr17_1 )=tracer(i,k1,j,p_tr17_1 )+emis*conv

                     emis=tr17_2 (cat)*fuel_burnt                                    ! emission from fire cell
                     t_tracer(p_tr17_2) = t_tracer(p_tr17_2) + emis                         ! add to total
                     tracer(i,k1,j,p_tr17_2 )=tracer(i,k1,j,p_tr17_2 )+emis*conv

                     emis=tr17_3 (cat)*fuel_burnt                                    ! emission from fire cell
                     t_tracer(p_tr17_3) = t_tracer(p_tr17_3) + emis                         ! add to total
                     tracer(i,k1,j,p_tr17_3 )=tracer(i,k1,j,p_tr17_3 )+emis*conv

                     emis=tr17_4 (cat)*fuel_burnt                                    ! emission from fire cell
                     t_tracer(p_tr17_4) = t_tracer(p_tr17_4) + emis                         ! add to total
                     tracer(i,k1,j,p_tr17_4 )=tracer(i,k1,j,p_tr17_4 )+emis*conv

                     emis=tr17_5 (cat)*fuel_burnt                                    ! emission from fire cell
                     t_tracer(p_tr17_5) = t_tracer(p_tr17_5) + emis                         ! add to total
                     tracer(i,k1,j,p_tr17_5 )=tracer(i,k1,j,p_tr17_5 )+emis*conv

                     emis=tr17_6 (cat)*fuel_burnt                                    ! emission from fire cell
                     t_tracer(p_tr17_6) = t_tracer(p_tr17_6) + emis                         ! add to total
                     tracer(i,k1,j,p_tr17_6 )=tracer(i,k1,j,p_tr17_6 )+emis*conv

                     emis=tr17_7 (cat)*fuel_burnt                                    ! emission from fire cell
                     t_tracer(p_tr17_7) = t_tracer(p_tr17_7) + emis                         ! add to total
                     tracer(i,k1,j,p_tr17_7 )=tracer(i,k1,j,p_tr17_7 )+emis*conv

                     emis=tr17_8 (cat)*fuel_burnt                                    ! emission from fire cell
                     t_tracer(p_tr17_8) = t_tracer(p_tr17_8) + emis                         ! add to total
                     tracer(i,k1,j,p_tr17_8 )=tracer(i,k1,j,p_tr17_8 )+emis*conv
                endif
           enddo
       enddo
    enddo
enddo

#ifdef WRF_CHEM
if(fire_print_msg .ge. msglevel .and.printsums .gt. 0)then
    ! sum ground concentrations and check for nans
    a_chem=0.0
    g_chem=0.0
    errors=0
    do j=jts,jte
       do k=1,chem_np
          k_p=chem_pointers(k)
          do i=its,ite
             if(chem(i,k1,j,k_p ).ne.chem(i,k1,j,k_p ))errors=errors+1
             a_chem(k_p)=a_chem(k_p)+chem(i,k1,j,k_p )
          enddo
       enddo
    enddo
    if(errors>0)call crash('NaN after chem update')
    call wrf_dm_sum_reals(a_chem,g_chem)
    call message('Layer1 raw sums after adding fire emissions',level=msglevel)
    do i=1,chem_np,line
        m=min(i+line-1,chem_np)
        write(msg,80)'Emissions ',(trim(chem_names(j)),j=i,m)
        call message(msg,level=msglevel)
        write(msg,81)'Layer1 end',(g_chem(chem_pointers(j)),j=i,m)
        call message(msg,level=msglevel)
        call message(' ',level=msglevel)
    enddo
endif
#endif

if(fire_print_msg .ge. msglevel .and.printsums .gt. 0)then
    call message("Computing totals over all processes",level=msglevel)
    ! sum over processes and add to cumulative sums

    ! fuel burned
    call wrf_dm_sum_reals(t_fuel,s_fuel)
    ! scale
    s_fuel = s_fuel*dx*dy 
    ! get rates
    r_fuel = s_fuel/dt 
    ! add to cumulative totals
    if(mfuelcats < nfuelcats)call crash('add_fire_emissions: bad size c_fuel')
    do i=1,nfuelcats
        c_fuel(i) = c_fuel(i) + s_fuel(i)
    enddo


#ifdef WRF_CHEM
    call wrf_dm_sum_reals(a_chem,g_chem)
    ! chem
    call wrf_dm_sum_reals(t_chem,s_chem)
    s_chem = s_chem*dx*dy 
    r_chem = s_chem/dt 
    if(max_chem<num_chem)call crash('add_fire_emissions: bad size c_chem')
    do i=1,num_chem
        c_chem(i) = c_chem(i) + s_chem(i)
    enddo

    call message('Total emissions in g or mol per the table',level=msglevel)
    do i=1,chem_np,line
        m=min(i+line-1,chem_np)
        write(msg,80)'Emissions ',(trim(chem_names(j)),j=i,m)
        call message(msg,level=msglevel)
        write(msg,81)'Cumulative',(c_chem(chem_pointers(j)),j=i,m)
        call message(msg,level=msglevel)
        write(msg,81)'Rate per s',(r_chem(chem_pointers(j)),j=i,m)
        call message(msg,level=msglevel)
        call message(' ',level=msglevel)
    enddo
#endif

    if(num_tracer >0)then
        call message("Computing totals of tracers over all processes",level=msglevel)
        ! tracer
        call wrf_dm_sum_reals(t_tracer,s_tracer)
        s_tracer = s_tracer*dx*dy 
        r_tracer = s_tracer/dt 
        if(max_tracer<num_tracer)call crash('add_fire_emissions: bad size c_tracer')
        do i=1,num_tracer
            c_tracer(i) = c_tracer(i) + s_tracer(i)
        enddo
    
        call message('Total emissions in g',level=msglevel)
        do i=1,tracer_np,line
            m=min(i+line-1,tracer_np)
            write(msg,80)'Emissions ',(trim(tracer_names(j)),j=i,m)
            call message(msg,level=msglevel)
            write(msg,81)'Cumulative',(c_tracer(tracer_pointers(j)),j=i,m)
            call message(msg,level=msglevel)
            write(msg,81)'Rate per s',(r_tracer(tracer_pointers(j)),j=i,m)
            call message(msg,level=msglevel)
            call message(' ',level=msglevel)
        enddo
    endif


    do cat=1,nfuelcats
        if(c_fuel(cat) > 0.)then 
            write(msg,83)fuel_name(cat),' burned',c_fuel(cat),'kg, rate',r_fuel(cat),' kg/s'
            call message(msg,level=msglevel)
        endif
    enddo
    write(msg,83)'Total fuel',' burned',sum(c_fuel),'kg, rate',sum(r_fuel),' kg/s'
    call message(msg,level=msglevel)

endif
80 format(a,8a11)
81 format(a,8e11.3)
83 format(a30,a,g14.4,a,g14.4,a,a)

return

9 continue
!$OMP CRITICAL(SFIRE_ATM_CRIT)
write(msg,91)ifts,ifte,jtfs,jfte,ifms,ifme,jfms,jfme
call message(msg,level=0)
write(msg,91)its,ite,jts,jte,ims,ime,jms,jme
call message(msg,level=0)
write(msg,92)'input  mesh size:',isz2,jsz2
call message(msg,level=0)
91 format('dimensions: ',8i8)
write(msg,92)'output mesh size:',isz1,jsz1
call message(msg,level=0)
92 format(a,2i8)
!$OMP END CRITICAL(SFIRE_ATM_CRIT)
call crash('add_fire_emissions: bad mesh sizes')

end subroutine add_fire_emissions

!
!***
!

subroutine check_pointers(array_name,array,np,pointer_names,pointers)
implicit none

!*** arguments
character(len=*),intent(in)::array_name
real, dimension(:,:,:,:),intent(in)::array
integer, intent(in)::np
character(len=*), dimension(:),intent(in)::pointer_names
integer, dimension(:),intent(in)::pointers

!*** local
integer::i,m,j,ps,pe
character(len=256)::msg

!** executable
!$OMP CRITICAL(SFIRE_ATM_CRIT)
do i=1,np,line
   m=min(i+line-1,np)
   write(msg,'(a7,8(1x,a8))')array_name,(trim(pointer_names(j)),j=i,m)
   call message(msg,level=msglevel)
   write(msg,'(a,8i9)')     'Pointer',(pointers(j),j=i,m)
   call message(msg,level=msglevel)
   call message(' ')
enddo

ps = lbound(array,4)
pe = ubound(array,4)
do i=1,np
    write(msg,995)array_name,'species',trim(pointer_names(i)),'pointer',pointers(i)
995 format(a,1x,a,1x,a,1x,a,i3)
    call message(msg,level=2)
    if(pointers(i) < ps .or. pointers(i) > pe) then 
        call print_dims_4d(array_name,array)
        write(msg,995)array_name,'species',trim(pointer_names(i)),'pointer',pointers(i)
        call message(msg,level=-1)
        call crash('add_fire_emissions: pointer out of bounds')
    endif
enddo
!$OMP END CRITICAL(SFIRE_ATM_CRIT)
end subroutine check_pointers

!
!***
!

subroutine print_dims_4d(array_name,array)
! print line with 4D array name and dimensions

!*** arguments
character(len=*), intent(in) ::array_name
real, dimension(:,:,:,:), intent(in) ::array

!*** local
character(len=256)::msg

!*** executable
993 format(3a,4(1x,i3,':',i3))
write(msg,993)'array ',array_name,' has dimensions ',&
         lbound(array,1),ubound(array,1), &
         lbound(array,2),ubound(array,2), &
         lbound(array,3),ubound(array,3), &
         lbound(array,4),ubound(array,4)
call message(msg,level=-1)
end subroutine print_dims_4d

!
!***
!


SUBROUTINE fire_tendency( &
    ids,ide, kds,kde, jds,jde,   & ! dimensions
    ims,ime, kms,kme, jms,jme,   &
    its,ite, kts,kte, jts,jte,   &
    grnhfx,grnqfx,canhfx,canqfx, & ! heat fluxes summed up to  atm grid 
    alfg,alfc,z1can,             & ! coeffients, properties, geometry 
    zs,z_at_w,dz8w,mu,rho,       &
    rthfrten,rqvfrten)             ! theta and Qv tendencies 

! This routine is atmospheric physics 
! it does NOT go into module_fr_sfire_phys because it is not fire physics

! taken from the code by Ned Patton, only order of arguments change to the convention here
! --- this routine takes fire generated heat and moisture fluxes and
!     calculates their influence on the theta and water vapor 
! --- note that these tendencies are valid at the Arakawa-A location

   IMPLICIT NONE

! --- incoming variables

   INTEGER , INTENT(in) :: ids,ide, kds,kde, jds,jde, &
                           ims,ime, kms,kme, jms,jme, &
                           its,ite, kts,kte, jts,jte

   REAL, INTENT(in), DIMENSION( ims:ime,jms:jme ) :: grnhfx,grnqfx  ! W/m^2
   REAL, INTENT(in), DIMENSION( ims:ime,jms:jme ) :: canhfx,canqfx  ! W/m^2
   REAL, INTENT(in), DIMENSION( ims:ime,jms:jme ) :: zs  ! topography (m abv sealvl)
   REAL, INTENT(in), DIMENSION( ims:ime,jms:jme ) :: mu  ! dry air mass (Pa)

   REAL, INTENT(in), DIMENSION( ims:ime,kms:kme,jms:jme ) :: z_at_w ! m abv sealvl
   REAL, INTENT(in), DIMENSION( ims:ime,kms:kme,jms:jme ) :: dz8w   ! dz across w-lvl
   REAL, INTENT(in), DIMENSION( ims:ime,kms:kme,jms:jme ) :: rho    ! density

   REAL, INTENT(in) :: alfg ! extinction depth ground fire heat (m)
   REAL, INTENT(in) :: alfc ! extinction depth crown  fire heat (m)
   REAL, INTENT(in) :: z1can    ! height of crown fire heat release (m)

! --- outgoing variables

   REAL, INTENT(out), DIMENSION( ims:ime,kms:kme,jms:jme ) ::   &
       rthfrten, & ! theta tendency from fire (in mass units)
       rqvfrten    ! Qv tendency from fire (in mass units)
! --- local variables

   INTEGER :: i,j,k
   INTEGER :: i_st,i_en, j_st,j_en, k_st,k_en

   REAL :: cp_i
   REAL :: rho_i
   REAL :: xlv_i
   REAL :: z_w
   REAL :: fact_g, fact_c
   REAL :: alfg_i, alfc_i

   REAL, DIMENSION( its:ite,kts:kte,jts:jte ) :: hfx,qfx
   
!!   character(len=128)::msg


        do j=jts,jte
            do k=kts,min(kte+1,kde)
               do i=its,ite
                   rthfrten(i,k,j)=0.
                   rqvfrten(i,k,j)=0.
               enddo
            enddo
        enddo


! --- set some local constants
   

   cp_i = 1./cp     ! inverse of specific heat
   xlv_i = 1./xlv   ! inverse of latent heat
   alfg_i = 1./alfg
   alfc_i = 1./alfc

!!write(msg,'(8e11.3)')cp,cp_i,xlv,xlv_i,alfg,alfc,z1can
!!call message(msg)

   call print_2d_stats(its,ite,jts,jte,ims,ime,jms,jme,grnhfx,'fire_tendency:grnhfx')
   call print_2d_stats(its,ite,jts,jte,ims,ime,jms,jme,grnqfx,'fire_tendency:grnqfx')

! --- set loop indicies : note that 

   i_st = MAX(its,ids+1)
   i_en = MIN(ite,ide-1)
   k_st = kts
   k_en = MIN(kte,kde-1)
   j_st = MAX(jts,jds+1)
   j_en = MIN(jte,jde-1)

! --- distribute fluxes

   DO j = j_st,j_en
      DO k = k_st,k_en
         DO i = i_st,i_en

            ! --- set z (in meters above ground)

            z_w = z_at_w(i,k,j) - zs(i,j) ! should be zero when k=k_st

            ! --- heat flux

            fact_g = cp_i * EXP( - alfg_i * z_w )
            IF ( z_w < z1can ) THEN
               fact_c = cp_i
            ELSE
               fact_c = cp_i * EXP( - alfc_i * (z_w - z1can) )
            END IF
            hfx(i,k,j) = fact_g * grnhfx(i,j) + fact_c * canhfx(i,j) 

!!            write(msg,2)i,k,j,z_w,grnhfx(i,j),hfx(i,k,j)
!!2           format('hfx:',3i4,6e11.3)
!!            call message(msg)

            ! --- vapor flux

            fact_g = xlv_i * EXP( - alfg_i * z_w )
            IF (z_w < z1can) THEN
               fact_c = xlv_i
            ELSE
               fact_c = xlv_i * EXP( - alfc_i * (z_w - z1can) )
            END IF
            qfx(i,k,j) = fact_g * grnqfx(i,j) + fact_c * canqfx(i,j) 
            
!!            if(hfx(i,k,j).ne.0. .or. qfx(i,k,j) .ne. 0.)then
!!                write(msg,1)i,k,j,hfx(i,k,j),qfx(i,k,j)
!!1               format('tend:',3i6,2e11.3)
!!                call message(msg)
!            endif

         END DO
      END DO
   END DO

! --- add flux divergence to tendencies
!
!   multiply by dry air mass (mu) to eliminate the need to 
!   call sr. calculate_phy_tend (in dyn_em/module_em.F)

   ! print *,'fire_tendency:',i_st,i_en,j_st,j_en
   DO j = j_st,j_en
      DO k = k_st,k_en-1
         DO i = i_st,i_en

            rho_i = 1./rho(i,k,j)

            rthfrten(i,k,j) = - mu(i,j) * rho_i * (hfx(i,k+1,j)-hfx(i,k,j)) / dz8w(i,k,j)
            rqvfrten(i,k,j) = - mu(i,j) * rho_i * (qfx(i,k+1,j)-qfx(i,k,j)) / dz8w(i,k,j)
 
            ! print *,i,j,k,rthfrten(i,k,j)

         END DO
      END DO
   END DO

   call print_3d_stats(its,ite,kts,kte,jts,jte,ims,ime,kms,kme,jms,jme,rthfrten,'fire_tendency:rthfrten')
   call print_3d_stats(its,ite,kts,kte,jts,jte,ims,ime,kms,kme,jms,jme,rqvfrten,'fire_tendency:rqvfrten')

   RETURN

END SUBROUTINE fire_tendency

!
!***
!
subroutine interpolate_atm2fire(id,               & ! for debug output, <= 0 no output
    ids,ide, kds,kde, jds,jde,                    & ! atm grid dimensions
    ims,ime, kms,kme, jms,jme,                    &
    ips,ipe,jps,jpe,                              &
    its,ite,jts,jte,                              &
    ifds, ifde, jfds, jfde,                       & ! fire grid dimensions
    ifms, ifme, jfms, jfme,                       &
    ifps, ifpe, jfps, jfpe,                       & ! fire patch bounds
    ifts,ifte,jfts,jfte,                          &
    ir,jr,                                        & ! atm/fire grid ratio
    u_frame, v_frame,                             & ! velocity frame correction
    u,v,                                          & ! atm grid arrays in
    ph,phb,                                       &
    z0,zs,                                        &
    uah,vah,                                      &
    uf,vf)                                          ! fire grid arrays out
    
implicit none
! Jan Mandel, October 2010
!*** purpose: interpolate winds and height

!*** arguments
integer, intent(in)::id                           
integer, intent(in)::                             &
    ids,ide, kds,kde, jds,jde,                    & ! atm domain bounds
    ims,ime, kms,kme, jms,jme,                    & ! atm memory bounds 
    ips,ipe,jps,jpe,                              &
    its,ite,jts,jte,                              & ! atm tile bounds
    ifds, ifde, jfds, jfde,                       & ! fire domain bounds
    ifms, ifme, jfms, jfme,                       & ! fire memory bounds
    ifps, ifpe, jfps, jfpe,                       & ! fire patch bounds
    ifts,ifte,jfts,jfte,                          & ! fire tile bounds
    ir,jr                                         ! atm/fire grid refinement ratio
real, intent(in):: u_frame, v_frame                 ! velocity frame correction
real,intent(in),dimension(ims:ime,kms:kme,jms:jme)::&
    u,v,                                          & ! atm wind velocity, staggered  
    ph,phb                                          ! geopotential
real,intent(in),dimension(ims:ime,jms:jme)::&
    z0,                                           & ! roughness height
    zs                                              ! terrain height
real,intent(out),dimension(ims:ime,jms:jme)::&
    uah,                                           & ! atm wind at fire_wind_height, diagnostics
    vah                                              ! 
real,intent(out), dimension(ifms:ifme,jfms:jfme)::&
    uf,vf                                           ! wind velocity fire grid nodes 
    
    
!*** local
character(len=256)::msg
#define TDIMS its-2,ite+2,jts-2,jte+2
real, dimension(its-2:ite+2,jts-2:jte+2):: ua,va   ! atm winds, interpolated over height, still staggered grid
real, dimension(its-2:ite+2,kds:kde,jts-2:jte+2):: altw,altub,altvb,hgtu,hgtv ! altitudes
integer:: i,j,k,ifts1,ifte1,jfts1,jfte1,ite1,jte1
integer::itst,itet,jtst,jtet,itsu,iteu,jtsu,jteu,itsv,itev,jtsv,jtev
integer::kdmax,its1,jts1,ips1,jps1
integer::itsou,iteou,jtsou,jteou,itsov,iteov,jtsov,jteov
real:: ground,loght,loglast,logz0,logfwh,ht,zr
real::r_nan
integer::i_nan
equivalence (i_nan,r_nan)

!*** executable

! debug init local arrays
i_nan=2147483647
ua=r_nan
va=r_nan
altw=r_nan
altub=r_nan
hgtu=r_nan
hgtv=r_nan


if(kds.ne.1)then
!$OMP CRITICAL(SFIRE_DRIVER_CRIT)
  write(msg,*)'WARNING: bottom index kds=',kds,' should be 1?'
  call message(msg)
!$OMP END CRITICAL(SFIRE_DRIVER_CRIT)
endif

!                            ^ j
!        ------------        |
!        |          |         ----> i
!        u    p     |
!        |          |    nodes in cell (i,j)
!        ------v-----    view from top     
!
!             v(ide,jde+1)
!            -------x------        
!            |            |      
!            |            |      
! u(ide,jde) x      x     x u(ide+1,jde) 
!            | p(ide,hde) |   
!            |            |   p=ph,phb,z0,...
!            -------x------        
!              v(ide,jde)
!
! staggered values set u(ids:ide+1,jds:jde) v(ids:ide,jds:jde+1)
! p=ph+phb set at (ids:ide,jds:jde) 
! location of u(i,j) needs p(i-1,j) and p(i,j)
! location of v(i,j) needs p(i,j-1) and p(i,j)
! *** NOTE need HALO in ph, phb
! so we can compute only u(ids+1:ide,jds:jde) v(ids:ide,jds+1,jde)
! unless we extend p at the boundary
! but because we care about the fire way in the inside it does not matter
! if the fire gets close to domain boundary the simulation is over anyway

    ite1=snode(ite,ide,1)
    jte1=snode(jte,jde,1)
    ! do this in any case to check for nans
    call print_3d_stats(its,ite1,kds,kde,jts,jte,ims,ime,kms,kme,jms,jme,u,'wind U in')
    call print_3d_stats(its,ite,kds,kde,jts,jte1,ims,ime,kms,kme,jms,jme,v,'wind V in')

    if(fire_print_msg.gt.0)then
!$OMP CRITICAL(SFIRE_DRIVER_CRIT)
       write(msg,'(a,f7.2,a)')'interpolate_atm2fire: log-interpolation of wind to',fire_wind_height,' m'
       call message(msg)
!$OMP END CRITICAL(SFIRE_DRIVER_CRIT)
    endif

! indexing
       
! for w 
    itst=ifval(ids.eq.its,its,its-1)
    itet=ifval(ide.eq.ite,ite,ite+1)
    jtst=ifval(jds.ge.jts,jts,jts-1)
    jtet=ifval(jde.eq.jte,jte,jte+1)

! for u
    itsu=ifval(ids.eq.its,its+1,its)  ! staggered direction
    iteu=ifval(ide.eq.ite,ite,ite+1)  ! staggered direction
    jtsu=ifval(jds.ge.jts,jts,jts-1)
    jteu=ifval(jde.eq.jte,jte,jte+1)

! for v
    jtsv=ifval(jds.eq.jts,jts+1,jts)  ! staggered direction
    jtev=ifval(jde.eq.jte,jte,jte+1)  ! staggered direction
    itsv=ifval(ids.ge.its,its,its-1)
    itev=ifval(ide.eq.ite,ite,ite+1)


if(fire_print_msg.ge.1)then
!$OMP CRITICAL(SFIRE_DRIVER_CRIT)
  write(msg,7001)'atm input  ','tile',its,ite,jts,jte
  call message(msg)
  write(msg,7001)'altw       ','tile',itst,itet,jtst,jtet
  call message(msg)
!$OMP END CRITICAL(SFIRE_DRIVER_CRIT)
endif
7001 format(a,' dimensions ',a4,':',i6,' to ',i6,' by ',i6,' to ',i6)

!**********************************************************
!*                                                        *
!*           find the altitude of the w points            *
!*                                                        *
!**********************************************************
!! in future, replace by z8w & test if the same

    kdmax=kde-1   ! max layer to interpolate from, can be less

    do j = jtst,jtet
      do k=kds,kdmax+1
        do i = itst,itet       
          altw(i,k,j) = (ph(i,k,j)+phb(i,k,j))/g             ! altitude of the bottom w-point 
        enddo
      enddo
    enddo

! values at u points
if(fire_print_msg.ge.1)then
!$OMP CRITICAL(SFIRE_DRIVER_CRIT)
  write(msg,7001)'u interp at','tile',itsu,iteu,jtsu,jteu
  call message(msg)
!$OMP END CRITICAL(SFIRE_DRIVER_CRIT)
endif

!**********************************************************
!*                                                        *
!*  interpolate the altitude from w points to u points    *
!*                                                        *
!**********************************************************

    do j = jtsu,jteu          
      do k=kds,kdmax+1
        do i = itsu,iteu       
          altub(i,k,j)= 0.5*(altw(i-1,k,j)+altw(i,k,j))      ! altitude of the bottom point under u-point
        enddo
      enddo
      do k=kds,kdmax
        do i = itsu,iteu       
          hgtu(i,k,j) =  0.5*(altub(i,k,j)+altub(i,k+1,j)) - altub(i,kds,j)  ! height of the u-point above the ground
        enddo
      enddo
    enddo

! values at v points
if(fire_print_msg.ge.1)then
!$OMP CRITICAL(SFIRE_DRIVER_CRIT)
  write(msg,7001)'v interp at','tile',itsv,itev,jtsv,jtev
  call message(msg)
!$OMP END CRITICAL(SFIRE_DRIVER_CRIT)
endif

!**********************************************************
!*                                                        *
!*  interpolate the altitude from w points to v points    *
!*                                                        *
!**********************************************************

    do j = jtsv,jtev          
      do k=kds,kdmax+1
        do i = itsv,itev       
          altvb(i,k,j)= 0.5*(altw(i,k,j-1)+altw(i,k,j))      ! altitude of the bottom point under v-point
        enddo
      enddo
      do k=kds,kdmax
        do i = itsv,itev       
          hgtv(i,k,j) =  0.5*(altvb(i,k,j)+altvb(i,k+1,j)) - altvb(i,kds,j)  ! height of the v-point above the ground
        enddo
      enddo
    enddo

#ifdef DEBUG_OUT
        call write_array_m3(itsu,iteu,kds,kdmax,jtsu,jteu,its-2,ite+2,kds,kde,jts-2,jte+2,altub,'altub',id)
        call write_array_m3(itsv,itev,kds,kdmax,jtsv,jtev,its-2,ite+2,kds,kde,jts-2,jte+2,altvb,'altvb',id)
        call write_array_m3(itsu,iteu,kds,kdmax,jtsu,jteu,its-2,ite+2,kds,kde,jts-2,jte+2,hgtu,'hgtu',id)
        call write_array_m3(itsv,itev,kds,kdmax,jtsv,jtev,its-2,ite+2,kds,kde,jts-2,jte+2,hgtv,'hgtv',id)
#endif

    logfwh = log(fire_wind_height)

!**********************************************************
!*                                                        *
!*  interpolate u vertically to fire_wind_height          *
!*                                                        *
!**********************************************************

    ! interpolate u, staggered in X

    do j = jtsu,jteu            ! compute on domain by 1 smaller, otherwise z0 and ph not available
      do i = itsu,iteu        ! compute with halo 2
        zr = 0.5*(z0(i,j)+z0(i-1,j)) ! interpolated roughness length under this u point
        if(fire_wind_height > zr)then       !  
          do k=kds,kdmax
            ht = hgtu(i,k,j)      ! height of this u point above the ground
            if( .not. ht < fire_wind_height) then ! found layer k this point is in
              loght = log(ht)
              if(k.eq.kds)then               ! first layer, log linear interpolation from 0 at zr
                logz0 = log(zr)
                ua(i,j)= u(i,k,j)*(logfwh-logz0)/(loght-logz0)
              else                           ! log linear interpolation
                loglast=log(hgtu(i,k-1,j))
                ua(i,j)= u(i,k-1,j) + (u(i,k,j) - u(i,k-1,j)) * ( logfwh - loglast) / (loght - loglast)
              endif
              goto 10
            endif
            if(k.eq.kdmax)then                 ! last layer, still not high enough
              ua(i,j)=u(i,k,j) 
            endif
          enddo
10        continue
        else  ! roughness higher than the fire wind height
          ua(i,j)=0.
        endif
      enddo
    enddo


!**********************************************************
!*                                                        *
!*  interpolate v vertically to fire_wind_height          *
!*                                                        *
!**********************************************************

    ! interpolate v, staggered in Y

    do j = jtsv,jtev
      do i = itsv,itev
        zr = 0.5*(z0(i,j-1)+z0(i,j)) ! roughness length under this v point
        if(fire_wind_height > zr)then       !  
          do k=kds,kdmax
            ht = hgtv(i,k,j)      ! height of this u point above the ground
            if( .not. ht < fire_wind_height) then ! found layer k this point is in
              loght = log(ht)
              if(k.eq.kds)then               ! first layer, log linear interpolation from 0 at zr
                logz0 = log(zr)
                va(i,j)= v(i,k,j)*(logfwh-logz0)/(loght-logz0)
              else                           ! log linear interpolation
                loglast=log(hgtv(i,k-1,j))
                va(i,j)= v(i,k-1,j) + (v(i,k,j) - v(i,k-1,j)) * ( logfwh - loglast) / (loght - loglast)
              endif
              goto 11
            endif
            if(k.eq.kdmax)then                 ! last layer, still not high enough
              va(i,j)=v(i,k,j) 
            endif
          enddo
11        continue
        else  ! roughness higher than the fire wind height
          va(i,j)=0.
        endif
      enddo
    enddo

#ifdef DEBUG_OUT
!   store the output for diagnostics
    do j = jts,jte1
      do i = its,ite1
        uah(i,j)=ua(i,j)
        vah(i,j)=va(i,j)
      enddo
    enddo

    call write_array_m(its,ite1,jts,jte,ims,ime,jms,jme,uah,'uah_n',id) ! no reflection
    call write_array_m(its,ite,jts,jte1,ims,ime,jms,jme,vah,'vah_n',id)
#endif

!**********************************************************
!*                                                        *
!*  interpolate ua,va vertically to the fire mesh         *
!*                                                        *
!**********************************************************


    ips1 = ifval(ips.eq.ids,ips+1,ips)
    call continue_at_boundary(1,1,0., & ! x direction 
       TDIMS,                  &! memory dims atm grid tile
       ids+1,ide,jds,jde, &     ! domain dims - where u defined
       ips1,ipe,jps,jpe, &     ! patch dims 
       itsu,iteu,jtsu,jteu, & ! tile dims - in nonextended direction one beyond if at patch boundary but not domain
       itsou,iteou,jtsou,jteou, & ! out, where set
       ua)                           ! array

    jps1 = ifval(jps.eq.jds,jps+1,jps)
    call continue_at_boundary(1,1,0., & ! y direction 
       TDIMS,                  & ! memory dims atm grid tile
       ids,ide,jds+1,jde, &      ! domain dims - where v wind defined
       ips,ipe,jps1,jpe, &        ! patch dims 
       itsv,itev,jtsv,jtev, & ! tile dims
       itsov,iteov,jtsov,jteov, & ! where set
       va)                           ! array

!   store the output for diagnostics
    do j = jts,jte1
      do i = its,ite1
        uah(i,j)=ua(i,j)
        vah(i,j)=va(i,j)
      enddo
    enddo

#ifdef DEBUG_OUT
        call write_array_m(itsou,iteou,jtsou,jteou,TDIMS,ua,'ua',id)
        call write_array_m(itsov,iteov,jtsov,jteov,TDIMS,va,'va',id)
#endif

!$OMP CRITICAL(SFIRE_DRIVER_CRIT)
    ! don't have all values valid, don't check
     write(msg,12)'atm mesh wind U at',fire_wind_height,' m'
     call print_2d_stats(itsou,iteou,jtsou,jteou,TDIMS,ua,msg)
     write(msg,12)'atm mesh wind V at',fire_wind_height,' m'
     call print_2d_stats(itsov,iteov,jtsov,jteov,TDIMS,va,msg)
12  format(a,f6.2,a)
    call print_2d_stats(its,ite1,jts,jte,ims,ime,jms,jme,uah,'UAH')
    call print_2d_stats(its,ite,jts,jte1,ims,ime,jms,jme,vah,'VAH')
    !call write_array_m(its,ite1,jts,jte,ims,ime,jms,jme,uah,'uah',id)
    !call write_array_m(its,ite,jts,jte1,ims,ime,jms,jme,vah,'vah',id)
!$OMP END CRITICAL(SFIRE_DRIVER_CRIT)

!      ---------------
!     | F | F | F | F |   Example of atmospheric and fire grid with
!     |-------|-------|   ir=jr=4.
!     | F | F | F | F |   Winds are given at the midpoints of the sides of the atmosphere grid,
!     ua------z-------|   interpolated to midpoints of the cells of the fine fire grid F.
!     | F | F | F | F |   This is (1,1) cell of atmosphere grid, and [*] is the (1,1) cell of the fire grid.
!     |---------------|   ua(1,1) <--> uf(0.5,2.5)
!     | * | F | F | F |   va(1,1) <--> vf(2.5,0.5)
!      -------va------    za(1,1) <--> zf(2.5,2.5)
!
!   ^ x2
!   |  --------va(1,2)---------
!   | |            |           |   Example of atmospheric and fire grid with
!   | |            |           |   ir=jr=1.
!   | |          za,zf         |   Winds are given at the midpoints of the sides of the atmosphere grid,
!   | ua(1,1)----uf,vf-----ua(2,1) interpolated to midpoints of the cells of the (the same) fire grid 
!   | |           (1,1)        |   ua(1,1) <--> uf(0.5,1) 
!   | |            |           |   va(1,1) <--> vf(1,0.5) 
!   | |            |           |   za(1,1) <--> zf(1,1)
!   |  --------va(1,1)---------
!   |--------------------> x1 
!
! Meshes are aligned by the lower left cell of the domain. Then in the above figure
! u = node with the ua component of the wind at (ids,jds), midpoint of side
! v = node with the va component of the wind at (ids,jds), midpoint of side
! * = fire grid node at (ifds,jfds)
! z = node with height, midpoint of cell
! 
! ua(ids,jds)=uf(ifds-0.5,jfds+jr*0.5-0.5)         = uf(ifds-0.5,jfds+(jr-1)*0.5)
! va(ids,jds)=vf(ifds+ir*0.5-0.5,jfds-0.5)         = vf(ifds+(ir-1)*0.5,jfds-0.5)
! za(ids,jds)=zf(ifds+ir*0.5-0.5,jfds+jr*0.5-0.5)  = zf(ifds+(ir-1)*0.5,jfds+(jr-1)*0.5)
    
    ! ifts1=max(snode(ifts,ifps,-1),ifds) ! go 1 beyond patch boundary but not at domain boundary
    ! ifte1=min(snode(ifte,ifpe,+1),ifde)
    ! jfts1=max(snode(jfts,jfps,-1),jfds)
    ! jfte1=min(snode(jfte,jfpe,+1),jfde)

    call interpolate_2d(  &
        TDIMS,                  & ! memory dims atm grid tile
        itsou,iteou,jtsou,jteou,& ! where set
        ifms,ifme,jfms,jfme,    & ! array dims fire grid
        ifts,ifte,jfts,jfte,& ! dimensions on the fire grid to interpolate to
        ir,jr,                  & ! refinement ratio
        real(ids),real(jds),ifds-0.5,jfds+(jr-1)*0.5, & ! line up by lower left corner of domain
        ua,                     & ! in atm grid     
        uf)                      ! out fire grid

    call interpolate_2d(  &
        TDIMS,                  & ! memory dims atm grid tile
        itsov,iteov,jtsov,jteov,& ! where set 
        ifms,ifme,jfms,jfme,    & ! array dims fire grid
        ifts,ifte,jfts,jfte,& ! dimensions on the fire grid to interpolate to
        ir,jr,                  & ! refinement ratio
        real(ids),real(jds),ifds+(ir-1)*0.5,jfds-0.5, & ! line up by lower left corner of domain
        va,                     & ! in atm grid     
        vf)                      ! out fire grid

call print_2d_stats_vec(ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme,uf,vf,'fire wind (m/s)')
! call print_2d_stats_vec(ifts1,ifte1,jfts1,jfte1,ifms,ifme,jfms,jfme,uf,vf,'fire wind extended')
#ifdef DEBUG_OUT
        call write_array_m(ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme,uf,'uf1',id)
        call write_array_m(ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme,vf,'vf1',id)
        ! call write_array_m(ifts1,ifte1,jfts1,jfte1,ifms,ifme,jfms,jfme,uf,'uf1',id)
        ! call write_array_m(ifts1,ifte1,jfts1,jfte1,ifms,ifme,jfms,jfme,vf,'vf1',id)
#endif

return

end subroutine interpolate_atm2fire

subroutine apply_windrf(                            &
            ifms,ifme,jfms,jfme,                    &
            ifts,ifte,jfts,jfte,                    &
            nfuel_cat,uf,vf)
integer::                                           &
    ifms, ifme, jfms, jfme,                       & ! fire memory bounds
    ifts, ifte, jfts, jfte                          ! fire tile bounds
real,intent(in),dimension(ifms:ifme,jfms:jfme)::nfuel_cat
real,intent(inout),dimension(ifms:ifme,jfms:jfme)::uf,vf
!*** local
    integer::i,j,k
!*** executable

    do j = jfts,jfte
      do i = ifts,ifte
          k=int( nfuel_cat(i,j) )
          if(k.lt.no_fuel_cat)then 
              uf(i,j)=uf(i,j)*windrf(k)
              vf(i,j)=vf(i,j)*windrf(k)
          else
              uf(i,j)=0.
              vf(i,j)=0.
          endif
      enddo
    enddo

end subroutine apply_windrf

!
!***
!

subroutine setup_wind_log_interpolation(          &
            ids,ide,  jds,jde,                    & ! atm grid dimensions
            ims,ime,  jms,jme,                    &
            ips,ipe,jps,jpe,                              &
            its,ite,jts,jte,                              &
            ifds, ifde, jfds, jfde,                       & ! fire grid dimensions
            ifms, ifme, jfms, jfme,                       &
            ifts, ifte, jfts, jfte,                       &
            ir,jr,                                        & ! atm/fire grid ratio
            z0,                                           & ! atm grid arrays in
            nfuel_cat,                                    & ! fire array in
            fz0,fwh,wz0)                                    ! fire arrays out
!*** arguments 
integer, intent(in)::                    &
    ids,ide, jds,jde,                    & ! atm domain bounds
    ims,ime, jms,jme,                    & ! atm memory bounds 
    ips,ipe,jps,jpe,                     &
    its,ite,jts,jte,                     & ! atm tile bounds
    ifds, ifde, jfds, jfde,                       & ! fire domain bounds
    ifms, ifme, jfms, jfme,                       & ! fire memory bounds
    ifts, ifte, jfts, jfte,                       & ! fire tile bounds
    ir,jr                                         ! atm/fire grid refinement ratio
real,intent(in),dimension(ims:ime,jms:jme)::z0    ! landuse roughness length
real,intent(in),dimension(ifms:ifme,jfms:jfme)::nfuel_cat   ! fuel category
real,intent(out),dimension(ifms:ifme,jfms:jfme)::&
    fz0,                                          & ! roughness height (fire)
    fwh,                                          & ! height to read the wind from
    wz0                                             ! WRF roughness height

!*** local
integer::i,j,ii,jj,k,id=0
character(len=128)::msg
real::r
!*** executable
 
     if(.not.have_fuel_cats)call crash('setup_wind_log_interpolation: fuel categories not yet set')

     select case(fire_wind_log_interp)

     case(1)
       call message('fire_wind_log_interp=1: log interpolation on fire mesh, roughness and wind height from fuel categories')
       do j=jfts,jfte
          do i=ifts,ifte
            k = int(nfuel_cat(i,j))
            if(k.ge.no_fuel_cat.and.k.le.no_fuel_cat2)then ! no fuel
                fz0(i,j)=-1.
                fwh(i,j)=-1.
            elseif(k < 1 .or. k  > nfuelcats) then
!$OMP CRITICAL(SFIRE_ATM_CRIT)
                 write(msg,*)'i,j,nfuel_cat,nfuelcats=',i,j,k,nfuelcats
!$OMP END CRITICAL(SFIRE_ATM_CRIT)
                 call message(msg)
                 call crash('setup_wind_log_interpolation: fuel category out of bounds')
            else
                fz0(i,j)=fcz0(k)
                fwh(i,j)=fcwh(k)
            endif
          enddo
       enddo
       do j=jts,jte
         do i=its,ite
            do jj=(j-1)*jr+1,(j-1)*jr+jr
               do ii=(i-1)*ir+1,(i-1)*ir+ir
                  wz0(ii,jj)=z0(i,j)
               enddo
            enddo
         enddo
       enddo
     case(2)
        call message('fire_wind_log_interp=2: log interpolation on fire mesh' // &
          'piecewise constant roughness from landuse, constant fire_wind_height')
        do j=jts,jte
           do i=its,ite
              do jj=(j-1)*jr+1,(j-1)*jr+jr
                 do ii=(i-1)*ir+1,(i-1)*ir+ir
                     fz0(ii,jj)=z0(i,j)
                     wz0(ii,jj)=z0(i,j)
                 enddo
              enddo
           enddo
        enddo
       do j=jfts,jfte
          do i=ifts,ifte
            k = int(nfuel_cat(i,j))
            if(k.lt.no_fuel_cat)then ! no fuel, interpolation does not matter
                fwh(i,j)=fcwh(k)
            else
                fz0(i,j)=-1.
                fwh(i,j)=-1.
            endif
          enddo
       enddo

     case(3)
       call message('fire_wind_log_interp=3: log interpolation on fire mesh,' // &
           ' interpolated roughness from landuse, constant fire_wind_height')
       call interpolate_z2fire(id,1,                & ! for debug output, <= 0 no output
            ids,ide,  jds,jde,                    & ! atm grid dimensions
            ims,ime,  jms,jme,                    &
            ips,ipe,jps,jpe,                              &
            its,ite,jts,jte,                              &
            ifds, ifde, jfds, jfde,                       & ! fire grid dimensions
            ifms, ifme, jfms, jfme,                       &
            ifts,ifte,jfts,jfte,                          &
            ir,jr,                                        & ! atm/fire grid ratio
            z0,                                       & ! atm grid arrays in
            fz0)                                      ! fire grid arrays out
       do j=jfts,jfte
          do i=ifts,ifte
            k = int(nfuel_cat(i,j))
            if(k.ne.no_fuel_cat)then ! no fuel, interpolation does not matter
                fwh(i,j)=fcwh(k)
            else
                fz0(i,j)=-1.
                fwh(i,j)=-1.
            endif
            wz0(i,j) = 0      ! Not using WRF roughness length here so set to 0.
          enddo
       enddo

     case(4)
        call message('fire_wind_log_interp=4: log interpolation on atmospheric' // &
           ' mesh, roughness from landuse, constant fire_wind_height')
        return

     case default
        !$OMP CRITICAL(SFIRE_ATM_CRIT)
        write(msg,*)'setup_wind_log_interpolation: invalid fire_wind_log_interp=',fire_wind_log_interp
        !$OMP END CRITICAL(SFIRE_ATM_CRIT)
        call crash('msg') 

     end select 

     select case(fire_use_windrf)

     case(0)  
       call message('setup_wind_log_interpolation: not using wind reduction factors')

     case(1)
       call message('setup_wind_log_interpolation: multiplying wind by reduction factors')

     case(2)
       call message('setup_wind_log_interpolation: resetting wind interpolation height from wind reduction factors')
       do j=jfts,jfte
          do i=ifts,ifte
            k = int(nfuel_cat(i,j))
            if(k.ne.no_fuel_cat)then 
                fwh(i,j) = fz0(i,j) ** (1.-windrf(k)) * fire_wind_height ** windrf(k) ! GMD paper eq. (26)

                if (.not. fz0(i,j) > 0. .or. .not. fwh(i,j) > fz0(i,j))then
!$OMP CRITICAL(SFIRE_ATM_CRIT)
                  write(msg,*)'category ',k,'windrf=',windrf(k),' fire_wind_height=',fire_wind_height
                  call message(msg,level=-1)
                  write(msg,*)'i=',i,' j=',j,' fz0(i,j)=',fz0(i,j),' fwh(i,j)=',fwh(i,j)
                  call message(msg,level=-1)
!$OMP END CRITICAL(SFIRE_ATM_CRIT)
                  call crash('setup_wind_log_interpolation: must have fwh > fz0 > 0')
                endif

            endif
          enddo
       enddo

     case(3)
      if(fire_wind_log_interp.eq.2.or.fire_wind_log_interp.eq.3)then
       call message('setup_wind_log_interpolation: adjusting wind interpolation height for LANDUSE roughness height')
       do j=jfts,jfte
          do i=ifts,ifte
            k = int(nfuel_cat(i,j))
            if(k.lt.no_fuel_cat)then 
                r = log(fcwh(k)/fcz0(k))/log(fire_wind_height/fcz0(k))! GMD paper eq. (25)
                fwh(i,j) = fz0(i,j) ** (1.-r) * fire_wind_height ** r ! GMD paper eq. (26)

                if (.not. fz0(i,j) > 0. .or. .not. fwh(i,j) > fz0(i,j))then
!$OMP CRITICAL(SFIRE_ATM_CRIT)
                  write(msg,*)'category ',k, 'roughness ',fcz0(k),' midflame height ',fcwh(k),' fire_wind_height=',fire_wind_height
                  call message(msg,level=-1)
                  write(msg,*)'computed wind reduction factor ',r
                  call message(msg,level=-1)
                  write(msg,*)'i=',i,' j=',j,' fz0(i,j)=',fz0(i,j),' fwh(i,j)=',fwh(i,j)
                  call message(msg,level=-1)
!$OMP END CRITICAL(SFIRE_ATM_CRIT)
                  call crash('setup_wind_log_interpolation: must have fwh > fz0 > 0')
                endif

            endif
          enddo
       enddo
      else
       call message('setup_wind_log_interpolation: not using wind reduction factors')
      endif

     case(4)
      call message('setup_wind_log_interpolation: resetting wind interpolation height from wind reduction factors (DVM NEW)')
      do j=jfts,jfte
         do i=ifts,ifte
           k = int(nfuel_cat(i,j))
           if(k.ne.no_fuel_cat)then 
               fwh(i,j) = fz0(i,j) ** (1.-windrf(k)) * fire_wind_height ** windrf(k) ! GMD paper eq. (26)

               if (.not. fz0(i,j) > 0. .or. .not. fwh(i,j) > fz0(i,j))then
!$OMP CRITICAL(SFIRE_ATM_CRIT)
                 write(msg,*)'category ',k,'windrf=',windrf(k),' fire_wind_height=',fire_wind_height
                 call message(msg,level=-1)
                 write(msg,*)'i=',i,' j=',j,' fz0(i,j)=',fz0(i,j),' fwh(i,j)=',fwh(i,j)
                 call message(msg,level=-1)
!$OMP END CRITICAL(SFIRE_ATM_CRIT)
                 call crash('setup_wind_log_interpolation: must have fwh > fz0 > 0')
               endif

           endif
         enddo
      enddo

     case default
        !$OMP CRITICAL(SFIRE_ATM_CRIT)
        write(msg,*)'setup_wind_log_interpolation: invalid fire_use_windrf=',fire_use_windrf
        !$OMP END CRITICAL(SFIRE_ATM_CRIT)
        call crash('msg') 

     end select 

!      consistency check
       do j=jfts,jfte
          do i=ifts,ifte
            k = int(nfuel_cat(i,j))
            if(k.lt.no_fuel_cat)then 
              if (.not. fz0(i,j) > 0. .or. .not. fwh(i,j) > fz0(i,j))then
!$OMP CRITICAL(SFIRE_ATM_CRIT)
                 write(msg,*)'i=',i,' j=',j,' fz0(i,j)=',fz0(i,j),' fwh(i,j)=',fwh(i,j)
                 call message(msg,level=-1)
!$OMP END CRITICAL(SFIRE_ATM_CRIT)
                 call crash('setup_wind_log_interpolation: must have fwh > fz0 > 0')
              endif
            else
              if(.not.fwh(i,j)<0.)then
!$OMP CRITICAL(SFIRE_ATM_CRIT)
                 write(msg,*)'i=',i,' j=',j,' fz0(i,j)=',fz0(i,j),' fwh(i,j)=',fwh(i,j)
                 call message(msg,level=-1)
!$OMP END CRITICAL(SFIRE_ATM_CRIT)
                 call crash('setup_wind_log_interpolation: no fuel must be signalled by fwh<0')
              endif
            endif
          enddo
       enddo

     have_wind_log_interpolation = .true.

     call print_2d_stats(ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme,fz0,'setup_wind_log:fz0')        
     call print_2d_stats(ifts,ifte,jfts,jfte,ifms,ifme,jfms,jfme,fz0,'setup_wind_log:fwh')        

end subroutine setup_wind_log_interpolation

!
!***
!

subroutine interpolate_wind2fire_height(id,       & ! to identify debugging prints and files if needed
    ids,ide, kds,kde, jds,jde,                    & ! atm grid dimensions
    ims,ime, kms,kme, jms,jme,                    &
    ips,ipe,jps,jpe,                              &
    its,ite,jts,jte,                              &
    ifds, ifde, jfds, jfde,                       & ! fire grid dimensions
    ifms, ifme, jfms, jfme,                       &
    ifps, ifpe, jfps, jfpe,                       & ! fire patch bounds
    ifts,ifte,jfts,jfte,                          &
    ir,jr,                                        & ! atm/fire grid ratio
    u_frame, v_frame,                             & ! velocity frame correction
    u,v,ph,phb,                                   & ! input atmospheric arrays
    fz0,fwh,                                      & ! input fire arrays
    uf,vf)                                          ! output fire arrays

  implicit none
!*** arguments
integer, intent(in):: id,                         & ! debug identification
    ids,ide, kds,kde, jds,jde,                    & ! atm domain bounds
    ims,ime, kms,kme, jms,jme,                    & ! atm memory bounds 
    ips,ipe,jps,jpe,                              &
    its,ite,jts,jte,                              & ! atm tile bounds
    ifds, ifde, jfds, jfde,                       & ! fire domain bounds
    ifms, ifme, jfms, jfme,                       & ! fire memory bounds
    ifps, ifpe, jfps, jfpe,                       & ! fire patch bounds
    ifts,ifte,jfts,jfte,                          & ! fire tile bounds
    ir,jr                                         ! atm/fire grid refinement ratio
real, intent(in):: u_frame, v_frame                 ! velocity frame correction
real,intent(in),dimension(ims:ime,kms:kme,jms:jme)::&
    u,v,                                          & ! atm wind velocity, staggered  
    ph,phb                                          ! geopotential
real,intent(in),dimension(ifms:ifme,jfms:jfme)::&
    fz0,                                          & ! roughness height
    fwh                                             ! height to read the wind from
real,intent(out),dimension(ifms:ifme,jfms:jfme)::&
    uf,                                           & ! atm wind at fire_wind_height, diagnostics
    vf                                              ! 

!*** local
  integer:: i,j,k,jcb,jcm,icb,icm,kdmax,kmin,kmax
  integer::itst,itet,jtst,jtet
  integer::iftst,iftet,jftst,jftet
  real:: wjcb,wjcm,wicb,wicm,ht,i_g,loght,zr,ht_last,logwh,wh,loght_last,uk,vk,uk1,vk1,z0,logz0
  real, dimension (its-1:ite+1,kds:kde,jts-1:jte+1)::z,zw
  character(len=128)::msg

!*** executable

  if(.not. have_wind_log_interpolation) call crash('interpolate_wind2fire_height: wind_log_interpolation must be set up first')

  ! print *,'interpolate_wind2fire_height start, id=',id

  kdmax=kde-1                                       ! max layer to use

! find the altitude of atm cell centers            

! index bounds for cell centers - need to go one beyond if at end of tile but not domain  
  itst=ifval(ids.eq.its,its,its-1)
  itet=ifval(ide.eq.ite,ite,ite+1)
  jtst=ifval(jds.ge.jts,jts,jts-1)
  jtet=ifval(jde.eq.jte,jte,jte+1)
  
  ! print *,'its, ite, jts, jte    =',its, ite, jts, jte
  ! print *,'itst, itet, jtst, jtet=',itst, itet, jtst, jtet

! get altitudes
  i_g = 1./g
  do j = jtst,jtet
    do k=kds,kdmax+1
      do i = itst,itet       
        zw(i,k,j) = (ph(i,k,j)+phb(i,k,j))*i_g      ! altitude of the bottom w-point 
      enddo
    enddo
    do k=kds,kdmax
      do i = itst,itet       
        z(i,k,j) = (zw(i,k,j)+zw(i,k+1,j))*0.5 - zw(i,kds,j) ! height of the cell center above the bottom
      enddo
    enddo
  enddo

! index bounds for fire mesh cell centers
! to prevent setting values from uninitialized memory
  iftst=ifval(ifds.eq.ifts,ifts+ir/2,ifts)
  iftet=ifval(ifde.eq.ifte,ifte-ir/2,ifte)
  jftst=ifval(jfds.ge.jfts,jfts+jr/2,jfts)
  jftet=ifval(jfde.eq.jfte,jfte-jr/2,jfte)

  ! print *,'iftst, iftet, jftst, jftet=',iftst, iftet, jftst, jftet

! zero out first, to prevent unitialized values on strips along domain boundaries
! it would be faster but longer code to do just cleanup loop on the strips
   do j = jfts,jfte
     do i = ifts,ifte
       uf(i,j)=0.
       vf(i,j)=0.
     enddo
  enddo

! vertical and horizontal interpolation

  kmin=kde  ! init stats
  kmax=kds

  loop_j: do j = jftst,jftet
    call coarse(j,jr,-2,jcb,wjcb)    ! get interpolation coefficients from the boundary 
    call coarse(j,jr,ir,jcm,wjcm)    ! get interpolation coefficients from the midpoint 
    loop_i: do i = iftst,iftet
      call coarse(i,ir,-2,icb,wicb)    ! get interpolation coefficients from the boundary 
      call coarse(i,ir,ir,icm,wicm)    ! get interpolation coefficients from the midpoint 
      z0 = fz0(i,j)                    ! roughness length 
      wh = fwh(i,j)                    ! wind height

      if( wh > z0 .and. z0 > 0)then

      ht_last=z0                     ! initialize starting height of this layer
      loop_k: do k=kds,kdmax         ! search for layer k such that ht(k-1)<=wh<ht(k), ht(0)=z0
        ! interpolate height from atmospheric cell midpoints
        ht=interpolate_h(its-1,ite+1,kds,kde,jts-1,jte+1,icm,k,jcm,wicm,wjcm,z)
        if(.not. ht < wh) exit loop_k ! found layer k this point is in
        ht_last = ht
      enddo loop_k

      if(k .gt. kdmax) then
          goto 91  ! run out of vertical levels, this must be wrong
      endif
      kmin=min(k,kmin)
      kmax=max(k,kmax)
      ! found layer k, ht_last < wh <= ht 
      logz0 = log(z0)
      logwh= log(wh)
      loght_last = log(ht_last)
      loght = log(ht)

      ! interpolate u at level k from staggered coarse grid: boundary in i, midpoints in j
      uk=interpolate_h(ims,ime,kms,kme,jms,jme,icb,k,jcm,wicb,wjcm,u) - u_frame
      ! interpolate v at level k from staggered coarse grid: midpoints in i, boundary in j
      vk=interpolate_h(ims,ime,kms,kme,jms,jme,icm,k,jcb,wicm,wjcb,v) - v_frame

      if(k.gt.kds)then           ! interpolate u,v horizontaly at the previous layer,  k-1
        ! interpolate u at level k-1 from staggered coarse grid: boundary in i, midpoints in j
        uk1=interpolate_h(ims,ime,kms,kme,jms,jme,icb,k-1,jcm,wicb,wjcm,u)
        ! interpolate v at level k-1 from staggered coarse grid: midpoints in i, boundary in j
        vk1=interpolate_h(ims,ime,kms,kme,jms,jme,icm,k-1,jcb,wicm,wjcb,v)
      else ! there is no previous layer, wind at roughness is assumed 0
        uk1=0. 
        vk1=0.
      endif

      ! log interpolate the wind to wh
      uf(i,j)= uk1 + (uk - uk1) * ( logwh - loght_last) / (loght - loght_last)
      vf(i,j)= vk1 + (vk - vk1) * ( logwh - loght_last) / (loght - loght_last)

    else
      ! below roughness or roughness zero, no wind 
      uf(i,j) = 0.
      vf(i,j) = 0.
    endif

    enddo loop_i
  enddo loop_j

  ! print *,'interpolate_wind2fire_height complete, id=',id

!$OMP CRITICAL(SFIRE_ATM_CRIT)
  write(msg,*)'wind interpolated from layers',kmin,' to ',kmax
  call message(msg)
!$OMP END CRITICAL(SFIRE_ATM_CRIT)

  return

  91 call crash('interpolate_wind2fire_height: fire wind height too large, increase kdmax or atm height')
  92 continue
!$OMP CRITICAL(SFIRE_ATM_CRIT)
  write(msg,*)'fz0(',i,j,')=',fz0(i,j),'fwh(',i,j,')=',fwh(i,j)
  call message(msg)
!$OMP END CRITICAL(SFIRE_ATM_CRIT)
  call crash('interpolate_wind2fire_height: must have fire wind height > roughness height > 0')

contains

real function interpolate_h(ims,ime,kms,kme,jms,jme,ic,kc,jc,wic,wjc,a)
!*** purpose: bilinear interpolation from a(ic:ic+1,k,jc:jc+1) with weights wicm wjcm
  implicit none
!*** arguments
  integer, intent(in)::ims,ime,jms,kms,kme,jme,ic,kc,jc
  real, intent(in)::wic,wjc,a(ims:ime,kms:kme,jms:jme)
!*** executable
  interpolate_h = &
    a(ic,kc,jc)    *wic     *wjc      + &
    a(ic,kc,jc+1)  *wic     *(1.-wjc) + &
    a(ic+1,kc,jc)  *(1.-wic)*wjc      + &
    a(ic+1,kc,jc+1)*(1.-wic)*(1.-wjc)
end function interpolate_h


subroutine coarse(ix,nr,ia,ic,w)
!*** find coarse mesh index and interpolation weight
!*** arguments
    implicit none
    integer, intent(in)::ix,nr,ia
    integer, intent(out)::ic
    real, intent(out)::w
!*** description
! given fine mesh with nr cells for each coarse mesh cell and such that
!   coarse mesh node 1 is aligned with the fine mesh at (na+1)/2
! for fine mesh node ix find its coarse mesh coordinate c and return 
!   ic=floor(c), the index of the nearest coarse mesh node below 
!   w =1-(c-ic), the interpolation weight from coarse mesh node ic to fine mesh node ix 
!
! Intended use:
! fine mesh nodes are always at the middle of their cells
!
! the alignment when the coarse nodes are on the boundary of coarse cells:
! |---1---|---2---|.......|--nr---|   fine mesh
! 1-------------------------------2   coarse mesh
! ia = -2  because coarse node 1 is aligned with the fine mesh at -1/2 = (-2 + 1)/2
!
! the alignment when the coarse node is at the midpoint of coarse  cell:
! |---1---|---2---|---3---|---4---|   fine mesh, here nr=4
! |---------------1---------------|   coarse mesh
! ia = nr because coarse node 1 is aligned with the fine mesh at (nr + 1)/2
! here,  (4 + 1)/2 = 2.5
!

!*** local
    real:: c,a

!*** executable

  a = (ia + 1)*0.5    ! the location on the fine grid where coarse node 1 is aligned
  c = 1 + (ix - a)/nr ! coarse mesh coordinate of ix
  ic= floor(c)        ! nearest coarse node to the left
  w = (1 + ic) - c    ! interpolation weight, 1-(c-ic)

end subroutine coarse

end subroutine interpolate_wind2fire_height


subroutine find_trees_fmesh(id,       & ! to identify debugging prints and files if needed
  ids,ide, kds,kde, jds,jde,                    & ! atm grid dimensions
  ims,ime, kms,kme, jms,jme,                    &
  ips,ipe,jps,jpe,                              &
  its,ite,jts,jte,                              &
  ifds, ifde, jfds, jfde,                       & ! fire grid dimensions
  ifms, ifme, jfms, jfme,                       &
  ifps, ifpe, jfps, jfpe,                       & ! fire patch bounds
  ifts,ifte,jfts,jfte,                          &
  ir,jr,                                        & ! atm/fire grid ratio
  can_top,nfuel_cat)                              ! canopy top and fuel cat (input)


  implicit none

  !*** arguments
  integer, intent(in):: id,                     & ! debug identification
  ids,ide, kds,kde, jds,jde,                    & ! atm domain bounds
  ims,ime, kms,kme, jms,jme,                    & ! atm memory bounds 
  ips,ipe,jps,jpe,                              &
  its,ite,jts,jte,                              & ! atm tile bounds
  ifds, ifde, jfds, jfde,                       & ! fire domain bounds
  ifms, ifme, jfms, jfme,                       & ! fire memory bounds
  ifps, ifpe, jfps, jfpe,                       & ! fire patch bounds
  ifts,ifte,jfts,jfte,                          & ! fire tile bounds
  ir,jr                                         

  real,intent(in),dimension(ifms:ifme,jfms:jfme)::nfuel_cat
  real,intent(inout),dimension(ifms:ifme,jfms:jfme)::can_top


  !*** local variables
  
  integer :: i, j, iveg, tree_type
  real, dimension(8) :: z_max
  

  !Canopy height lookup table

  z_max   = (/10.,20.,15.,24.,22.5,10.,2.2,0.73/)


  !Loop through each i j point within our patch and relate our nfuel category # to a 
  !canopy height, which we specify in our fire namelist file. 

  do j = jfts,jfte
    do i = ifts,ifte

      iveg = nfuel_cat(i,j)
      tree_type = itree(iveg)

      !If the nfuel category is associated with a fuel type that is not associated with a
      !tree (designated as a '0' within the fire.namelist file), set the canopy top height 
      !as 0.

      if(tree_type /= 0) then
        
        !write(*,*) z_max(tree_type)
        can_top(i,j) = z_max(tree_type)

      else 

        can_top(i,j) = 0 

      end if 

    enddo
  enddo


end subroutine find_trees_fmesh


subroutine massman_fwh(id,       & ! to identify debugging prints and files if needed
    ids,ide, kds,kde, jds,jde,                    & ! atm grid dimensions
    ims,ime, kms,kme, jms,jme,                    &
    ips,ipe,jps,jpe,                              &
    its,ite,jts,jte,                              &
    ifds, ifde, jfds, jfde,                       & ! fire grid dimensions
    ifms, ifme, jfms, jfme,                       &
    ifps, ifpe, jfps, jfpe,                       & ! fire patch bounds
    ifts,ifte,jfts,jfte,                          &
    ir,jr,                                        & ! atm/fire grid ratio
    cuf,cvf,                                      & ! canopy top winds (input)
    nfuel_cat,                                    & ! fire fuel category (input)
    can_top,                                      & ! canopy top and FWH heights (input)
    uf,vf)                                          ! fwh winds (input and output)

    implicit none

    !*** arguments
    integer, intent(in):: id,                     & ! debug identification
    ids,ide, kds,kde, jds,jde,                    & ! atm domain bounds
    ims,ime, kms,kme, jms,jme,                    & ! atm memory bounds 
    ips,ipe,jps,jpe,                              &
    its,ite,jts,jte,                              & ! atm tile bounds
    ifds, ifde, jfds, jfde,                       & ! fire domain bounds
    ifms, ifme, jfms, jfme,                       & ! fire memory bounds
    ifps, ifpe, jfps, jfpe,                       & ! fire patch bounds
    ifts,ifte,jfts,jfte,                          & ! fire tile bounds
    ir,jr                                         

    real,intent(in),dimension(ifms:ifme,jfms:jfme)::&
    cuf,cvf,                                       & ! canopy winds on fire mesh
    nfuel_cat,                                     & ! fuel categories on fire mesh
    can_top                                          ! cantop top and fire wind heights on fire mesh       
    
    real,intent(inout),dimension(ifms:ifme,jfms:jfme)::uf,vf

    !*** local variables

    real :: xi_og, k, c1, c3, c2, mid_fh
    real :: ha, z_og, drag_index, cantop_frac, C_surf, N, cd, zeta_sum, u_scale,z_top
    real, dimension(8) :: xi_max, sigma_u, sigma_l, PAI, big_cd, z_max, d1
    real, dimension(20) :: z_agl, xi, fa, zeta, Ub, Ut, u
    integer :: iveg, i, j, z, znum


    !---------------------------------------------------------------------------------------!
    !   Fortran code for the Massman et al. 2017 parameterzation for within-canopy winds.   !
    !   This parameterization was specifically greared toward downscaling winds above the   !
    !   canopy displacement height to the fire flame height. However, this parameterization !
    !   could be applied for other applications involving canopy winds and forested areas.  !
    !                                                                                       !
    !   Written by DVM on 5/3/2019                                                          !
    !   W.J. Massman, J.M. Forthofer, and M.A. Finney, 2017 : An improved canopy wind       !
    !   model for predicting wind. Can. J. For. Res. 47: 594–603.                           !                     
    !---------------------------------------------------------------------------------------!


    !iveg = 5               !Hardcoded value for trouble shooting


    !Model constants used by the Massman canopy parameterization 

    xi_og = 0.0025          !"Ground" roughness length (in xi coordinates) [unitless]
    k = 0.40                !Von Karman constant
    c1 = 0.38               !used to compute ratio beteen u* and u at canopy top height [Wang 2012]
    c3 = 15                 !used to compute ratio beteen u* and u at canopy top height [Wang 2012]
    c2 = c1+k/log(xi_og)    !used to compute ratio beteen u* and u at canopy top height [Wang 2012]

    znum = 20               !Number of vertical levels for our canopy model? (20 is the default)


    !Computations for an arbituary function for 0 < xi < 1 that describes the 'desired mathematical 
    !shape of the plant surface distribution'. The assigned values below corresponds to:
    ! 
    !        Aspen, Scots Pine, Jack Pine, Loblolly Pine, Hardwood, Spruce, Corn and Rice.

    xi_max  = (/0.36,0.60,0.58,0.60,0.84,0.60,0.94,0.62/)
    sigma_u = (/0.60,0.30,0.20,0.10,0.13,0.38,0.03,0.50/)
    sigma_l = (/0.20,0.10,0.20,0.27,0.30,0.16,0.60,0.45/)
    PAI     = (/3.28,2.41,2.14,3.78,4.93,5.73,2.94,3.10/)
    big_cd  = (/0.20,0.25,0.27,0.20,0.09,0.20,0.30,0.30/)
    z_max   = (/10.,20.,15.,24.,22.5,10.,2.2,0.73/)
    d1      = (/0.,0.,-0.7,0.,2.,0.,0.,0./)


    !Loop through each i j point within our patch
    do j = jfts,jfte
      do i = ifts,ifte


        !What is the height of the tree canopy for the gridcell we are working on?
        !If z_top equals 0, we are not working with a gridcell with trees, thus 
        !we don't want to modify our u and v winds on the fire mesh so cycle to 
        !the next iteration.  
        
        z_top = can_top(i,j)

        if(z_top == 0) then 
          cycle
        end if 
        

        !Determine tree type based on whats prescribed in the namelist.fire file. 
        
        iveg = itree(nfuel_cat(i,j))


        !Compute the heights for our vertical coordinate based on the vegetation 
        !height (z_max) that was assigned. Heights are in mAGL. In addition, 
        !compute a unitless height coordinate (xi).

        do z = 1, znum
          if(z>1) then
              z_agl(z) = z_agl(z-1) + z_top/(znum-1)
          else
              z_agl(z) = 0
          end if 
        end do 

        xi = z_agl/z_top
      

        !Let compute fa, from equation 2 of the Massman paper.

        do z = 1, znum
          if(xi(z)>=xi_max(iveg)) then
              fa(z) = exp((-(xi(z)-xi_max(iveg))**2)/(sigma_u(iveg)**2))
          else
              fa(z) = exp((-(xi_max(iveg)-xi(z))**2)/(sigma_l(iveg)**2))
          end if   
        end do 
      
      
        !Equation 1 of Massman et al. [2017] for 'nondimensional plant surface distribution
      
        ha = sum((fa/sum(fa))*PAI(iveg))
      
      
        !"Ground" roughness length (in height coordinates) [mAGL], according to Massman et al. [2017]
      
        z_og = xi_og*z_top  
      
      
        !Drag area index, which simplifies down to the drag coefficient (Cd) x the plant area index 
        !(PAI) Massman et al. [2017]. This is needed to compute the ratio between u* and u 
        !(cantop_frac, eqn. 10), N (variables from Massman, 1997). 
      
        drag_index = big_cd(iveg)*PAI(iveg)
      
        cantop_frac = c1 - c2*exp(-c3*drag_index)
      
        C_surf = (cantop_frac**2)*2
      
        N = drag_index/C_surf
      
      
        !cd represents the 'bulk drag coefficient' at a given height. Since d1 and p1 are set zero
        !, equation 3 in Massman et al. 2017 simplifies down to cd(xi) = Cd... easy peasy!
      
        cd = big_cd(iveg)
      
      
        !Lastly, we need to compute the cumulative drag area as zeta which is a function of xi
      
        zeta_sum = 0
      
        do z = 1, znum
          zeta(z) = cd*(fa(z)/sum(fa))*PAI(iveg)+zeta_sum
          zeta_sum = zeta(z)
        end do 
      
      
        !Canopy wind and stress profiles (computed from eqn 6 and 7 Massman et al [2017])
      
        do z = 1, znum
          if(xi(z)>=xi_og) then
            Ub(z) = log(xi(z)/xi_og)/log(1/xi_og)
          else
            Ub(z) = 0
          end if   
        end do 
      
        Ut  = cosh((N*zeta)/drag_index)/cosh(N)
      
      
        !Compute the non-dimensional wind profile (u(z)/u(h)) in eqn. 11 Massman et al [2017])
      
        u = Ut*Ub      


        !What is the mid-flame height? We assume that the flame height is fuel depthx2,
        !thus the mid-flame height is set as the fuel depth. See Andrews et al. 2012.
        !Fuel depth is a namelist.fire parameter, which can be obtained by the fuel
        !category of the gridcell of interest.

        mid_fh = fueldepthm(nfuel_cat(i,j))


        !Finally, find the nearest two vertical levels to our mid-flame height, and interpolate our
        !wind to the mid-flame height using a simple interpolation scheme. 

        do z = 1, znum
          if(z_agl(z)>=mid_fh) then
              u_scale = (((z_agl(z-1)-mid_fh)*(u(z-1)-u(z)))/(z_agl(z)-z_agl(z-1)))+u(z-1)
              EXIT
          end if 
        end do 
      
        !Scale our u and v winds at the mid-flame height using our non-dimensional wind 
        !profile, which was interpolated to the mid-flame height. 

        uf(i,j) = cuf(i,j)*u_scale
        vf(i,j) = cvf(i,j)*u_scale


      enddo
    enddo


end subroutine massman_fwh


end module module_fr_sfire_atm
