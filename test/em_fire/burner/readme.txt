This is an idealized simulation with a stationary heat source AKA burner.
The burner.m script generates input_tign_g text file, that defines the 
heart source. This file is used by ideal.exe which processes the data and 
put them into the wrfinput_d01 file. The burner script creates a linear burner
placed along Y axis (N-S) direction, 25m thick and 750m long located 400m from the
western domain boundary, with the heat relaese of 1800KW/m2
The simulation is set for 20 minutes, in 25m resolution domain covering an area of 3000x6000m.
240x120x120 grid point in x,y, and z direction. Fuel type is set to category 3 (tall grass). 

The follwing setup is included in the burner.m:
burner_size_m     = [750,25]    % burner is 25m thick and 750m long 
domain_size_m     = [3000,6000] % that is the size of the domain in m
mesh_step_m       = [25,25]     % that is the fire mesh resolution
burner_dist_m     = 400         % that is the distance fron the weatern domain edge where the burner will be placed
burner_start_s    = 10          % that is the start time when the burner starts emitting heat
burner_end_s      = 1e6         % that is the end time when the heat release ends

The test burn setup uses a isothermal atmosphere with the temperature of 20C (273.15K), RH=33% and a 1/7 power law vertical wind profile:
U(z) = 5 (z/2)^1/7 m/s, surface pressure is set to 10^5 Pa, surface temperature is set to 20C as well. 
The surface properties are defined as custom category 28, in LANDUSE.TBL (USGS Winter) with the roughnes length z0=0.03m. 

The burner is turned on by:
   fire_update_fuel_frac =2, in namelist.input
   weight and combustion constant (cmbcnst) in namelist.fire determines the heat flux.

How to define the heat flux:
The heat flux is controlled by the weighting parameter that determines the slope of the mass lose curve of the fuel, which is used to compute the fuel_time:
  fuel_time=weight/0.85 
fuel_time is used to compute the burn rate as: burn_rate=1/fuel_time
the heat flux is computed as: burn_rate * cmbcnst, so:
  heat_flux = 0.85/weight * cmbcnst, 
to get a custom heat flux one has to modify the weighting parameter according to the formula below:
   weight=0.85/desired_heat_flux * cmbcnst
for instance:
   for weight = 7.4
   heat_flux = 0.85/7.4 * 17.433e+06 = 2e+06W/m2 = 2000KW/m2 

How to prepare wrfinput file with the predefined heat source:

1. Modify burner.m according to your setup, remeber to set the domain_size_m and mesh_step (fire grid resolution).
2. start matlab in this directory, or cd here and run startup
3. in matlab, run burner to create file input_tign_g 
4. run ./ideal.exe , this should generate wrfinput_d01 with TIGN_G mtrix defining the burner position and timming
5. run ./wrf.exe
