% to create conus.kml:
% download http://firemapper.sc.egov.usda.gov/data_viirs/kml/conus_hist/conus_20120914.kmz
% and gunzip 
% 
% to create w.mat:
% run Adam's simulation, currently results in
% /home/akochans/NASA_WSU/wrf-fire/WRFV3/test/em_barker_moist/wrfoutputfiles_live_0.25
% then in Matlab
% f='wrfout_d05_2012-09-15_00:00:00'; 
% t=nc2struct(f,{'Times'},{});  n=size(t.times,2);  w=nc2struct(f,{'TIGN_G','FXLONG','FXLAT','UNIT_FXLAT','UNIT_FXLONG'},{},n);
% save ~/w.mat w    
%
% to create c.mat
% c=nc2struct(f,{'NFUEL_CAT'},{},1);
% save ~/c.mat c
%
% to create s.mat:
% s=read_wrfout_sel({'wrfout_d05_2012-09-09_00:00:00','wrfout_d05_2012-09-12_00:00:00','wrfout_d05_2012-09-15_00:00:00'},{'FGRNHFX'}); 
% save ~/s.mat s 
% 

% ****** REQUIRES Matlab 2013a - will not run in earlier versions *******

v=read_fire_kml('conus_viirs.kml');

% v=read_fire_kml('conus_modis.kml');
load w
load s
load c
fuels

% establish boundaries

min_lat = min(w.fxlat(:))
max_lat = max(w.fxlat(:))
min_lon = min(w.fxlong(:))
max_lon = max(w.fxlong(:))

% select fire detection within the domain and time
bi=(v.lon > min_lon & v.lon < max_lon & v.lat > min_lat & v.lat < max_lat);
u = unique(v.tim(bi));
tim_ref = u(1);
bi = bi & v.tim == tim_ref;
lon=v.lon(bi);
lat=v.lat(bi);
res=v.res(bi);
tim=v.tim(bi);

% plot satellite detection points

% plot3(lon,lat,tim,'ko'),
% m=500; n=500;
[m,n]=size(w.fxlong);
for j=1:length(fuel), W(j)=fuel(j).weight;end
for j=length(fuel)+1:max(c.nfuel_cat(:)),W(j)=NaN;end
T = zeros(m,n);
for j=1:n, for i=1:m
        T(i,j)=W(c.nfuel_cat(i,j));
end,end

% rebase time on the largest tign_g = the time of the last frame

last_time=datenum(char(w.times)');
max_tign_g=max(w.tign_g(:));

tim = tim - last_time;
tim_ref = tim_ref - last_time;
tign = (w.tign_g - max_tign_g)/(24*60*60);
tign(tign==0)=NaN;

while 1
tscale = 1000; % multiply fuel.weight by this to get detection time scale 
a=input('enter [tscale dir add1p add1m add2p add2m] (1,rad,s/m,s/m)\n');
tscale=a(1)
if tscale<=0, break, end
dir=a(2)
add1p=a(3)
add1m=a(4)
add2p=a(5)
add2m=a(6)

v1=[cos(dir),sin(dir)]; % main direction
v2=[cos(dir+pi/2),sin(dir+pi/2)]; % secondary direction

% find ignition point
[i_ign,j_ign]=find(w.tign_g == min(w.tign_g(:)));
    
% vector field (x,y) - (x_ign,y_ign) 
VDx=(w.fxlong-w.fxlong(i_ign,j_ign))*w.unit_fxlong;
VDy=(w.fxlat-w.fxlat(i_ign,j_ign))*w.unit_fxlat;
    
p1 = VDx*v1(1)+VDy*v1(2); % projections on main direction
p2 = VDx*v2(1)+VDy*v2(2); % projections on secondary direction
[theta,rho]=cart2pol(p1,p2); % in polar coordinates, theta between [-pi,pi]
thetas=pi*[-3/2,-1,-1/2,0,1/2,1,3/2];
deltas=[add2p add1m add2m add1p add2p add1m add2m];
delta = interp1(thetas,deltas,theta,'pchip').*rho;

tign_mod = tign + delta;


pmap = p_map(tscale*T/(24*60*60),tign_mod-tim_ref);

[mm,nn]=size(w.fxlong);
mi=1:ceil(mm/m):mm;
ni=1:ceil(nn/n):nn;
mesh_fxlong=w.fxlong(mi,ni);
mesh_fxlat=w.fxlat(mi,ni);
mesh_fgrnhfx=s.fgrnhfx(mi,ni,:);
mesh_pmap = pmap(mi,ni); 
mesh_lon = mesh_fxlong(mi,ni);
mesh_lat = mesh_fxlat(mi,ni);
mesh_tign= tign(mi,ni); 

% draw the fireline
figure(1)
hold off, clf
%contour(mesh_fxlong,mesh_fxlat,mesh_tign,[tim_ref,tim_ref]);

% add the detection probability map
%hold on
h=pcolor(mesh_fxlong,mesh_fxlat,pmap);
set(h,'EdgeAlpha',0,'FaceAlpha',0.5);
colorbar
cc=colormap; cc(1,:)=1; colormap(cc);
hold off

% plot detection squares

C=0.5*ones(1,length(res));
rlon=0.5*res/w.unit_fxlong;
rlat=0.5*res/w.unit_fxlat;
X=[lon-rlon,lon+rlon,lon+rlon,lon-rlon]';
Y=[lat-rlat,lat-rlat,lat+rlat,lat+rlat]';
hold on
hh=fill(X,Y,C);
hold off

grid,drawnow

% hold on, (mesh_lon,mesh_lat,mesh_tim2),grid on
title('Barker Canyon fire VIIRS fire detection')
ylabel('latitude')
xlabel('longitude')

end