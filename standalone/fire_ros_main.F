program fire_ros_main

use module_fr_sfire_phys

! create file fuels.m read and visualize the results as described in 
! https://wiki.openwfm.org/wiki/How_to_diagnose_fuel_properties_in_WRF-SFIRE

call init_fuel_cats(init_fuel_moisture=.false.) 

end program fire_ros_main
