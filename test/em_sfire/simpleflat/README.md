# simpleflat idealized case

This is a simple test case that demonstrates the evolution of a grassfire burning in tall grass, classified as fuel category 3, in an environment characterized by a constant westerly wind flow. The wind profile specified in the `input_sounding` file is logarithmic, with a roughness length of 0.036 m and a wind speed of 5.5 m/s at a reference height of 6.1 m. This profile adheres to the LOG profile outlined by Kochanski et al. (2013) (https://doi.org/10.1002/jgrd.50436).
 
The simulation uses the surface layer model initialization activated using the `sfc_full_init =.true.` flag, which allows the use of the Revised MM5 Monin-Obukhov surface layer scheme (`sf_sfclay_physics = 1`), and the 5-layer thermal diffusion land surface model (`sf_surface_physics = 1`).
The critical variables describing the land surface type, surface temperature, and mean soil temperature are specified using the `sfc_lu_index`, `sfc_tsk`, and `sfc_tmn` variables respectively.

The simulation employs vertical hyperbolic grid stretching, activated by the flags `stretch_grd = .true.` and `stretch_hyp = .true.` The arrangement of the near-surface model layers is managed using the `z_grd_scale` parameter. A higher value of `z_grd_scale` results in a shallower first model layer and increases the number of model layers near the surface.

The simulation begins with two walking ignitions traveling north-south, perpendicular to the wind, starting from a common point located 1000 meters from the western edge of the domain, and lasts for 20 minutes. The simulation is run in Large Eddy Simulation (LES) mode using a mesh size of 82 x 82 x 41 (X x Y x Z) with a horizontal resolution of 60 m. Open boundary conditions are applied at the eastern and western boundaries, while periodic boundary conditions are used at the northern and southern boundaries.
