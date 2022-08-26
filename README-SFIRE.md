### WRF-SFIRE ###
This repository is https://github.com/wrf-sfire/WRF-SFIRE.git master branch. 
A mirror with a graphical log is at https://repo.or.cz/git-browser/by-commit.html?r=WRF-SFIRE.git

## Branches ##
* master, develop - the merged code.
* WRF-track/master, WRF-track/develop, etc. - tracking the official WRF repository including release tags
* filtered/* - rebased commits from the original wrf-fire repository, directory wrfv2_fire mapped to root:
_  filtered/balbi
-  filtered/devel
-  filtered/dvm_branch
-  filtered/master
-  filtered/submitted-to-3.3
-  filtered/wrf-fire-branch-for-fmc-merge
-  fuel-moisture-model
* wrf-fire-track/* - the original wrf-fire
-  wrf-fire-track/master
-  wrf-fire-track/dvm_branch
  
## How to upgrade WRF version ##
* git checkout WRF/master
* git pull git@github.com:wrf-model/WRF.git master
* git checkout develop
* git merge WRF/master
* test in case something got broken, when happy:
* git checkout master
* git merge develop

## How to run ##
* In namelist input, ifire=1 is WRF-SFIRE, ifire=2  the fire code in WRF we put there in 2012 (see branch filtered/submitted-to-3.3) with changes at NCAR. Most namelist flags are the same but each version has some of its own. 
* Test problems available in test/em_fire:
- rain ifire=2
- hill ifire=1
_ two_fires ifire=2
* The problems with ifire=2 should give the same result in branch master (the merged code) and in branch WRF/develop (the WRF distribution)
* It seems to run serial and MPI and so far passed some limited testing on cheyenne.

## Notes ##
Standalone is included but have not started on it yet.
Branches balbi and dvm_branch were carried over but they are not merged into master or tested yet because they were not merged into master in wrf-fire. Also I do not know how to test dvm_branch. 
I can't do real problems yet, that requires pushing data through current WPS, which Adam and Angel are working on.
