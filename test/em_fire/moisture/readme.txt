To set initial  moisture:

1. set in namelist_input

 fmoist_run = .false.,  
 fmoist_interp = .true., 

2. run ideal.exe or real.exe

3. start matlab in the fire_em directory, or cd to fire_em and run startup
to set the path

4. prepare in matlab the desired values  array x size (n1,n2,5) and run
ncreplace('wrfinput_d01','FMC_GC',x)

5. run ./wrf.exe 

To use constant moisture in some moisture class: set in namelist.fire
drying_lag and wetting lag to a large number, such as 1e10

The proportions of moisture classes in fuel categories are in the arrays
fmc_gw01,..., fmc_gw05 in namelist.fire


