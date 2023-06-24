#!/bin/csh -f
echo usage: $0 wrf.exe or something like $0 mpirun -np 5 wrf.exe for a parallel run 
setenv CMD "$*"
setenv DIFFWRF ../../../external/io_netcdf/diffwrf
setenv NLISTV namelist.input.var
setenv DIR _regtest_restart_tmp.dir
setenv WRFRST1 wrfrst_d01_0001-01-01_00:01:00 
setenv WRFRST5 wrfrst_d01_0001-01-01_00:02:00 

echo complete run 

cat $NLISTV | sed 's/RESTART/\.false\./' | sed 's/START_M/00/g' > namelist.input

rm -rf wrfrst_* wrfout_* rsl.* > /dev/null

./ideal.exe

$CMD
if ($status) then
	exit(1)
endif

rm -rf $DIR > /dev/null ; mkdir $DIR 
mv wrfrst_* wrfout_* rsl.* $DIR

echo restart run 

cat $NLISTV | sed 's/RESTART/\.true\./' | sed 's/START_M/01/g' > namelist.input
cp $DIR/$WRFRST1 .
$CMD
if ($status) then
	exit(1)
endif

$DIFFWRF $WRFRST5 $DIR/$WRFRST5 | grep -v '^ *TS_'





