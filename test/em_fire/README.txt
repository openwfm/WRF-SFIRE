ifire=1:
All files to run fire test cases are in this directory em_fire.
The files necessary for running the hill_simple case are linked directly 
in this directory.

Subdirectories for the other test cases and the links necessary for
running the test cases there are created by the script create_links.sh,
which is executed when running ./compile em_fire in the top level WRF 
directory, The subdirectories are destroyed when running ./clean -a.

If you want to make your own test case, you will simply need to 
create the files namelist.input and input_sounding that will
correspond to your case (best by modifying a copy from another case).

ifire=2:
The fire test problems are in subdirectories. At the moment, these include

rain_1
two_fires_1
fireflux

If you want to make your own test case subdirectory, all you need to do 
is to create soft links for ideal.exe and wrf.exe pointing the the parent 
directory, and to create the files namelist.input and input_sounding 
(best by modifying a copy from another subdirectory).

Do not just copy one of the existing subdirectories, the soft links might not 
be copied properly.
