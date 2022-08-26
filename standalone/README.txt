This directory contains SFIRE (the fire component of wrf-fire) test driver.
It links the fire model from the files wrfv2_fire/phys/module_fr_sfire_*.F with the 
files in this directory:

model_test_main.F	the main program that calls the model
wrf_fakes.F             stubs that perform various wrf functions

To build the standalone driver, create your own make.inc, or link an existing one, and type make. 
This will create file model_test_main.exe, then just execute this file.

For more information see http://www.openwfm.org/wiki/How_to_run_the_standalone_fire_model_in_WRF-Fire

Keeping the standalone up to date (for programmers)

If you get errors about missing configure_flags%something it means that the registry 
has changed and the include files need to be updated as follows:

 - build wrf in the wrf2_fire directory
 - in this directory: make sync_inc
 - commit the resulting include files

If you still get errors, it means that the standalone fell behind (again) and some hand updating is
needed. Let me know.

Spread rate calculation interface:

The spread rate is computed in module_fr_sfire_phys.F (maybe it should be renamed to 
module_fr_sfire_spread.F?). This module defines derived type fire_params which contains only
pointers.  The driver declares an object type fire_params and assigns the pointers to 
parameter arrays. This object is passed down the call chain to the spread rate calculation.
This way, when the parameters passed to the spread rate calculation change, no changes
are required in the code between the driver and the spread rate calculation.

Jan Mandel
June 18, 2010
Updated December 22, 2012
