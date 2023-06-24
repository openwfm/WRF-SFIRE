function run_moisture_test(varargin)
if nargin>=1,
    ncfile=varargin{1};
else
    ncfile=[];
end
if isempty(ncfile),
    ncfile='moisture100.nc';
end
if nargin>=2,
    hr=varargin{2};
else
    hr=[0,realmax];
end
ncload(ncfile)
Times=ncread(ncfile,'Times');
[d,hours]=make_moisture_input(Times,t2,psfc,q2,rain);
setenv('DYLD_LIBRARY_PATH','/opt/local/lib:/opt/local/lib/gcc46')
system('./moisture_test.exe');
load moisture_output.txt
out.t2   =moisture_output(:,3);
out.t2_mid=midpoints(out.t2);
out.psfc =moisture_output(:,4);
out.q2   =moisture_output(:,5);
% rh_fire is computed from midpoint values in the fortran code
out.rain=moisture_output(:,6);
out.rh_fire=moisture_output(:,7);
out.fmc_equi =moisture_output(:,8);
out.tlag     =moisture_output(:,9);
out.fmc_gc   =moisture_output(:,10);
out.id='from moisture_test.exe';
out.id='';
out.hours=hours;
out.hours_mid=midpoints(hours);
out.r=find(hours>=hr(1) & hours<=hr(2)); % range to display
orig.t2=t2(~d);
orig.t2_mid=midpoints(orig.t2);
orig.psfc=psfc(~d);
orig.q2=q2(~d);
orig.rain=rain(~d);
orig.rh_fire=rh_fire(~d);
orig.fmc_equi=fmc_equi(1,~d)';
orig.fmc_gc=fmc_gc(1,~d)';
orig.id=['from ',ncfile];
orig.hours=hours;
orig.hours_mid=midpoints(hours);
orig.r=out.r;
figure(1)
plot(hours,orig.fmc_gc,'r--',hours,out.fmc_gc,'k-')
xlabel hours
ylabel kg/kg
title('Fuel moisture')
h=legend(orig.id,out.id);
set(h,'Interpreter','none')
figure(2)
plot(hours,orig.fmc_equi,'r--',hours,out.fmc_equi,'k-')
xlabel hours
set(ylabel('kg/kg'),'Interpreter','none')
title('Equilibrium moisture')
h=legend(orig.id,out.id);
set(h,'Interpreter','none')
plot_moisture(3,orig)
plot_moisture(4,out)
plot_all(5,orig)
plot_all(6,out)
end

function plot_moisture(f,s)
figure(f);
[d,w]=equilibrium_moisture(s.rh_fire,s.t2_mid);
plot(s.hours_mid(s.r),d(s.r),'b--',...
     s.hours_mid(s.r),w(s.r),'b-.',...
     s.hours_mid(s.r),s.fmc_equi(s.r),'r--',...
     s.hours(s.r),s.fmc_gc(s.r),'k-')
xlabel hours
set(ylabel('Fuel moisture contents (kg / kg)'),'Interpreter','none')
h=title(['Fuel moisture ',s.id]);
set(h,'Interpreter','none')
h=legend('Drying','Wetting','Equilibrium','Actual');
set(h,'Interpreter','none')
setvmax(0.3)
end

function plot_all(f,s)
figure(f)
[d,w]=equilibrium_moisture(s.rh_fire,s.t2_mid);
plot(s.hours(s.r),s.t2(s.r)/max(s.t2),...
    s.hours(s.r),s.psfc(s.r)/max(s.psfc),...
    s.hours(s.r),s.q2(s.r)/max(s.q2),...
    s.hours(s.r),s.rh_fire(s.r),...
    s.hours_mid(s.r),s.fmc_equi(s.r),...
    s.hours_mid(s.r),d(s.r),'g-.',...
    s.hours_mid(s.r),w(s.r),'b-.',...
    s.hours(s.r),s.fmc_gc(s.r))
h=title(['All variables ',s.id]);
set(h,'Interpreter','none');
xlabel hours
h=legend('Scaled T2',...
    'Scaled PSFC',...
    'Scaled Q2',...
    'RH',...
    'FMC_EQUI',...
    'Drying equilibrium',...
    'Wetting equilibrium',...
    'FMC_GC')
set(h,'Interpreter','none');
setvmax(1)
end

function setvmax(v)
a=axis;
a(3:4)=[0,v];
axis(a)
end

function T=midpoints(T2)
T=[T2(1);0.5*(T2(1:end-1)+T2(2:end))];
end