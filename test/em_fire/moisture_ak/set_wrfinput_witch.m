function set_wrfinput_witch
% first run ideal.exe to create wrfinput_d01
% then run this script to modify some variables in wrfinput_d01
ddd=pwd;
cd ..
startup
cd(ddd)
load wrfinput_witch.mat
ncreplace('wrfinput_d01','NFUEL_CAT',nfuel_cat)
ncreplace('wrfinput_d01','FMC_GC',fmc_gc)
end