This directory contains SFIRE (the fire component of wrf-fire) test driver.
It builds the fire executable from the files sources WRF-SFIRE, files generated
in a WRF-SFIRE build, and the files in this directory.

How to build standalone:

1. build WRF-SFIRE as usual: cd ..; ./configure; compile em_fire (not required for make fire_ros.exe)
2. go back here: cd standalone
3. Select compiler: ln -s make.inc.ifort make.inc or one of the others or make your own
4. make  or make <nameofexecutable>.exe

All standalone executables should be run in the simulation directory where wrf.exe runs.
They use the inputs available to wrf.exe such as namelist.input and namelist.fire, etc.

fire.exe is the complete fire model running with atmosphere from WRF-SFIRE output. To use:

1. Create atmospheric forcing: cd test/em_fire/hill; ideal.exe; wrf.exe
2. link atmospheric forcing: ln -s the_wrfout_just_created fire_input.nc
3. Run the standalone in the simulation directoryt: <WRF SFIRE ROOT>/standalone/fire.exe

fire_ros.exe calls the fuels and rate of spread (ROS) subsystem. To use:
1. Run the standalone in the simulation directory: <WRF SFIRE ROOT>/standalone/fire_ros.ex
   This will create Matlab file fuels.m. Follow the directions to use it in Matlab at
   https://wiki.openwfm.org/wiki/How_to_diagnose_fuel_properties_in_WRF-SFIRE
    

Files:

fire.F             the main fire program 
module_domain.F    fake WRF declaration of grid
module_configure.F fake WRF declaration of namelist 
wrf_fakes.F        fake WRF subroutines
wrf_netcdf.F       read and write files in WRF compatible format
util_netcdf.F      convenience stubs to netcdf API   

The following main programs almost certainly do not work:

init.F                     generate a basic fire_input.nc 
moisture_main.F            run the fuel  moisture model
moisture_test_main.F       another one
fuels_main.F               not sure
atm.F                      not sure
fuel_interp_test_main.F    not sure

See also old guide for the original 2012 pre-WRF-SFIRE version at 
http://wiki.openwfm.org/wiki/How_to_run_the_standalone_fire_model_in_WRF-Fire

If you still get errors, it means that the standalone fell behind the fire
code (again) and some hand updating is needed.

Jan Mandel
June 18, 2010
Updated December 22, 2012
Updated August 1, 2022
