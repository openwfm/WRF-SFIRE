# Test cases

These directories contain all the necessary files (besides the executables that should be built when you compile wrf) to run idealized WRF-SFIRE simulations. To link the executables to your run, perform the following command after you have changed directory into a test simulation:

```shell
ln -s ../../main/*.exe .
```

That will link the executables to your directory. Next you need to copy in your job script, run ideal.exe (assuming your job script does not do that already), and then run wrf.exe. 
