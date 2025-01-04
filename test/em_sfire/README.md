# Test cases

This is the folder with the WRF-SFIRE-specific idealized test cases. It has been created to avoid conflicts with the legacy (NCAR's) fire implementation using ifire = 2 option.

The subdirectories contain all the necessary files (besides the executables that should be built when you compile wrf) to run idealized WRF-SFIRE simulations. To link the executables to your run, perform the following command after you have changed the directory into a test simulation:

```shell
ln -s ../../main/*.exe .
```

That will link the executables to your directory. Next you need to copy in your job script, run ideal.exe (assuming your job script does not do that already), and then run wrf.exe. 

## Test case description
### fireflux_small
This is a test case built upon the FireFlux I experiment (Clements, C. B., S. Zhong, S. Goodrick, J. Li, X. Bian, B.E. Potter, W. E. Heilman, J.J. Charney, R. Perna, M. Jang, D. Lee, M.Patel, S. Street and G. Aumann, 2007: Observing the Dynamics of Wildland Grass Fires: FireFlux- A Field Validation Experiment). This test case serves as an example of how custom fuel, topography, and land use data can be used in idealized fire cases.

### rain
This test case includes a simulation of a convective cloud that generates precipitation. This is the test case illustrating the use of the fuel moisture model and shows how the fuel moisture content responds to changing meteorological conditions.

### simpleflat
This is a simple test case illustrating an evolution of the parabolic fire front in response to a uniform westerly wind. 

### simplehill
This is a simple test case illustrating the effect of topography on the fire rate of spread. Both direct and indirect effects are included here. The effects of the topography on the wind speed, and the impact of slope on the ROS. This test case is analogous to simple flat, but has a hill extending from north to south. The terrain in this case is specified directly in the namelist.input.

