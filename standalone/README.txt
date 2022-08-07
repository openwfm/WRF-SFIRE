This directory contains SFIRE (the fire component of wrf-fire) test driver.
It builds the fire executable from the files sources WRF-SFIRE, files generated
in a WRF-SFIRE build, and the files in this directory.

How to build and test:

1. build WRF-SFIRE as usual: cd ..; ./configure; compile em_fire
2. back here: cd standalone
3. Select compiler: ln -s make.inc.ifort make.inc or one of the others or make your own
4. Build fire.exe: make 
5. Create atmospheric forcing: cd test/em_fire/hill; ideal.exe; wrf.exe
6. link atmospheric forcing: ln -s the_wrfout_just_created fire_input.nc
7. Fun the standalone: ../fire.exe

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
