# fefines FC, FCFLAGS, FFREE
include make.inc
DEFS =
WRF = ../wrfv2_fire
FPFLAGS= -I$(NETCDF)/include -I.
LDFLAGS = -L$(NETCDF)/lib -lnetcdff -lnetcdf
#LDFLAGS = `$(NETCDF)/bin/nc-config --flibs`
FFLAGS=$(FPFLAGS) $(FCFLAGS) $(FFREE)
INC = $(WRF)/inc


MAKEFILE=Makefile

%.o: $(WRF)/phys/%.F 
	$(FP) $(FFLAGS) $(DEFS) -E $(WRF)/phys/$*.F | sed 's/^#.*//' > $*.f90
	$(FC) -c $(FFLAGS) $*.f90

%.o: %.F 
	$(FP) $(FFLAGS) $(DEFS) -E $*.F | sed 's/^#.*//' > $*.f90
	$(FC) -c $(FFLAGS) $*.f90

ATM_OBJ = \
    module_domain.o \
    wrf_netcdf.o \
    module_configure.o \
    module_fr_sfire_model.o \
    module_fr_sfire_core.o \
    module_fr_sfire_phys.o \
    module_fr_sfire_atm.o \
    module_fr_sfire_util.o  \
    wrf_fakes.o


FIRE_OBJ = \
    module_domain.o \
    wrf_netcdf.o \
    module_configure.o \
    module_fr_sfire_model.o \
    module_fr_sfire_core.o \
    module_fr_sfire_phys.o \
    module_fr_sfire_atm.o \
    module_fr_sfire_util.o  \
	module_fr_sfire_driver.o \
    wrf_fakes.o

INIT_OBJ = \
    wrf_netcdf.o \
    module_configure.o \
    module_fr_sfire_util.o  \
    wrf_fakes.o

FUELS_OBJ = \
    fuels_main.o \
    module_fr_sfire_phys.o \
    module_fr_sfire_util.o  \
    wrf_fakes.o
    

NAMELIST= \
    namelist_statements.inc \
    namelist_defaults.inc \
    namelist_defines.inc \
    namelist_defines2.inc \
    config_assigns.inc

FUEL_INTERP_OBJ = \
    fuel_interp_test_main.o \
    module_fr_sfire_core.o \
    module_fr_sfire_phys.o \
    module_fr_sfire_util.o  \
    wrf_fakes.o

MOISTURE_TEST_OBJ = \
    moisture_test_main.o \
    module_fr_sfire_phys.o \
    module_fr_sfire_util.o  \
    wrf_fakes.o
    
MOISTURE_OBJ = \
    moisture_main.o \
    module_fr_sfire_phys.o \
    module_fr_sfire_util.o  \
    wrf_fakes.o

.PHONY: all


#.PHONY: all

default: fire

moisture: moisture.exe
moisture.exe: $(MAKEFILE) libesmf_time.a $(MOISTURE_OBJ)
	$(FC) $(FCFLAGS) -o moisture.exe $(MOISTURE_OBJ) -L. -lesmf_time $(LDFLAGS)
moisture_test: moisture_test.exe
moisture_test.exe: $(MAKEFILE) $(MOISTURE_TEST_OBJ)
	$(FC) $(FCFLAGS) -o moisture_test.exe $(MOISTURE_TEST_OBJ)


fuel: fuel_interp_test.exe
fuel_interp_test: fuel_interp_test.exe
 
fuel_interp_test.exe: $(MAKEFILE) $(FUEL_INTERP_OBJ)
	$(FC) $(FCFLAGS) -o fuel_interp_test.exe $(FUEL_INTERP_OBJ)

all: fire init 

atm: atm.exe
	cd ../wrfv2_fire/test/em_fire; rm -f atm.exe; ln -s ../../../standalone/atm.exe atm.exe

fire: fire.exe
	cd ../wrfv2_fire/test/em_fire; rm -f fire.exe; ln -s ../../../standalone/fire.exe fire.exe

init: init.exe
	cd ../wrfv2_fire/test/em_fire; rm -f init.exe; ln -s ../../../standalone/init.exe init.exe

atm.exe: $(MAKEFILE) atm.o $(FIRE_OBJ)
	$(FC) $(FCFLAGS) -o atm.exe fire.o $(FIRE_OBJ) $(LDFLAGS)

fire.exe: $(MAKEFILE) fire.o $(FIRE_OBJ)
	$(FC) $(FCFLAGS) -o fire.exe fire.o $(FIRE_OBJ) $(LDFLAGS)

init.exe: $(MAKEFILE) init.o $(INIT_OBJ)
	$(FC) $(FCFLAGS) -o init.exe init.o $(INIT_OBJ) $(LDFLAGS)

fuels: fuels.exe

fuels.exe: $(MAKEFILE) $(FUELS_OBJ)
	$(FC) $(FCFLAGS) $(LDFLAGS) -o fuels.exe $(FUELS_OBJ) $(LDFLAGS)
 
# for the order of generating modules

moisture_test_main.o: $(MAKEFILE) \
    module_fr_sfire_phys.o  \
    module_fr_sfire_util.o 

moisture_main.o: $(MAKEFILE) \
    module_fr_sfire_core.o  \
    module_fr_sfire_util.o 

fire.o: $(MAKEFILE) \
    module_fr_sfire_model.o \
    module_fr_sfire_core.o  \
	module_fr_sfire_driver.o \
    module_domain.o \
    module_configure.o \
    wrf_fakes.o \
    wrf_netcdf.o

init.o: $(MAKEFILE) \
    module_configure.o \
    wrf_netcdf.o \
    module_fr_sfire_util.o

fuel_interp_test_main.o: $(MAKEFILE) \
    module_fr_sfire_phys.o \
    module_fr_sfire_util.o 

model_fuel_interp_test.o: $(MAKEFILE) \
    module_fr_sfire_core.o  \
    module_fr_sfire_util.o \
    wrf_namelist.o \
    wrf_netcdf.o

module_fr_sfire_model.o: $(MAKEFILE) \
    module_fr_sfire_core.o  \
    module_fr_sfire_util.o 

module_fr_sfire_core.o: $(MAKEFILE) \
    module_fr_sfire_phys.o \
    module_fr_sfire_util.o 

module_fr_sfire_atm.o: \
    module_fr_sfire_phys.o \
    module_fr_sfire_util.o
      
module_fr_sfire_phys.o: $(MAKEFILE) \
    module_fr_sfire_util.o 

module_fr_sfire_util.o: $(MAKEFILE) \
    wrf_fakes.o

wrf_netcdf.o: $(MAKEFILE) \
    module_fr_sfire_util.o

module_fr_sfire_driver.o: $(MAKEFILE) \
	sfire_id.inc \
	module_domain.o \
	module_configure.o \
    module_fr_sfire_phys.o \
	module_fr_sfire_model.o \
	module_fr_sfire_atm.o \
	module_fr_sfire_util.o

module_configure.o: $(MAKEFILE) $(NAMELIST) \
    wrf_fakes.o

sfire_id.inc: __always
	chmod +x commit_hash
	./commit_hash > sfire_id.inc
	cat sfire_id.inc

fuels_main.o: $(MAKEFILE) \
    module_fr_sfire_phys.o \
    module_fr_sfire_util.o  \
    wrf_fakes.o

__always:


# build rules for ESMF time support
ESMF_OBJS=ESMF_Alarm.o ESMF_BaseTime.o ESMF_Time.o \
        ESMF_Base.o ESMF_Calendar.o ESMF_Fraction.o   \
        ESMF_TimeInterval.o ESMF_Stubs.o ESMF_Clock.o Meat.o \
		ESMF_AlarmClock.o ESMF_Mod.o

libesmf_time.a: $(ESMF_OBJS)
	$(RM) $@
	$(AR) ru $@ $^

ESMF_Alarm.o : ESMF_BaseTime.o ESMF_Time.o ESMF_TimeInterval.o
ESMF_BaseTime.o : ESMF_Base.o
ESMF_Time.o : ESMF_BaseTime.o ESMF_TimeInterval.o ESMF_Calendar.o \
              ESMF_Stubs.o
ESMF_Base.o :
ESMF_Calendar.o : ESMF_BaseTime.o
ESMF_Fraction.o : ESMF_BaseTime.o
ESMF_TimeInterval.o : ESMF_BaseTime.o ESMF_Calendar.o ESMF_Fraction.o
Meat.o : ESMF_Alarm.o ESMF_BaseTime.o ESMF_Clock.o ESMF_Time.o \
        ESMF_Base.o ESMF_Calendar.o ESMF_Fraction.o    \
        ESMF_TimeInterval.o
ESMF_Mod.o : ESMF_Alarm.o ESMF_BaseTime.o ESMF_Clock.o ESMF_Time.o \
        ESMF_Base.o ESMF_Calendar.o ESMF_Fraction.o    \
        ESMF_TimeInterval.o Meat.o ESMF_Stubs.o ESMF_AlarmClock.o
ESMF_AlarmClock.o : ESMF_Alarm.o ESMF_Clock.o
ESMF_Clock.o : ESMF_BaseTime.o ESMF_Time.o ESMF_TimeInterval.o
$(ESMF_OBJS) : %.o : $(WRF)/external/esmf_time_f90/%.F90
	$(FC) -c $< -I$(WRF)/external/esmf_time_f90 -I. -DMAX_HISTORY=25 -o $@

INCFILES = $(INC)/config_assigns.inc $(INC)/namelist_defaults.inc $(INC)/namelist_defines.inc $(INC)/namelist_defines2.inc $(INC)/namelist_statements.inc
sync_inc: $(INCFILES)
	cp $(INCFILES) .

clean:
	rm -f *.mod *.o *.f90 *.s *.exe *.a
