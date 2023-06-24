Test example for assimilation of NDWI

1. run ideal.exe
2. in matlab, cd up directory (em_fire), startup, cd ndwi
3. still in matlab, set_wrfinput, note the value of fmc_g at the end
3. look at output - the values of fmc_g should be the same
4. in matlab, confirm fmc_g=ncread('wrfout_d01_0001-01-01_00:00:00','FMG_C') 
for fmc_g(:,:,2) (fmc_g is set only after the initial wrfout is made) 