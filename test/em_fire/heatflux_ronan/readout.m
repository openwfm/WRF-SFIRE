function t=readout(name)
t0=ncread('wrfout_d01_0001-01-01_00:00:00',name);
t1=ncread('wrfout_d01_0001-01-01_00:10:00',name);
t2=ncread('wrfout_d01_0001-01-01_00:20:00',name);
t3=ncread('wrfout_d01_0001-01-01_00:30:00',name);
t=cat(ndims(t0),t0,t1,t2,t3);
end