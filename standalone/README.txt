This directory contains SFIRE (the fire component of wrf-fire) test drivers.
It builds executables from the files sources in WRF-SFIRE, files generated
in a WRF-SFIRE build, and the files in this directory.

A. fire_ros.exe

fire_ros.exe calls the fuels and rate of spread (ROS) subsystem. 
It reads namelist.fire if present in the current directory.
It does not read the namelist.input so WRF does not need to be built first. 

1. Select compiler: ln -s make.inc.ifort make.inc or one of the others or make your own
2. Build the executables: make fire_ros 
3. Run: cd ../test/em_fire/balbi or cd ../test/em_fire/hill;  ./fire_ros.exe
4. The first time, fire_ros.exe will create a file namelist_standalone.output,
   which has the same format as the input and can be used as a template.
5. cp namelist_standalone.output namelist_standalone.input
6. Edit namelist_standalone.input to change the inputs to fire_ros_balbi as desired, run again, repeat.  
7. If call_write_fuels_m=T in in input, fire_ros.exe will create file fuels.m. Follow the directions in
   https://wiki.openwfm.org/wiki/How_to_diagnose_fuel_properties_in_WRF-SFIRE
   to examine fuels.m and make graphics.

B. fire.exe 
fire.exe is the complete fire model running with atmosphere from WRF-SFIRE output.

1. build WRF-SFIRE as usual, options do not matter: cd ..; ./configure; compile em_fire
2. Select compiler: ln -s make.inc.ifort make.inc or one of the others or make your own
3. go back here: cd standalone
4. build the executable: make fire 
5. Create the atmospheric forcing: cd test/em_fire/hill; ideal.exe; wrf.exe
6. link the atmospheric forcing: ln -s <the_wrfout_just_createdi> fire_input.nc
7. Run the standalone: ../../../standalone/fire.exe


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
