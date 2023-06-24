function prop_test
% a jm/mkim version of propagate
% follows Osher-Fedkiw and Mitchell toolbox

% defaults
% grid
%xs=15; % domain size
xs=10;
ys=10;
%m=121 % domain mesh
m=101
n=101
time0=0
time1=10
r=0
speed=2  % wind speed
wind='c'
alpha=0.75*pi
dt=0.2;
plotSteps=50
normal_spread_c=1  % speed = r + wind^e * c
normal_spread_e=1

example='p'

switch example
    case {'circlespin','c'}
    wind='c'
    ignition='c'
    case {'gostraight','g'}
    wind='s'
    ignition='c'
    case {'physical','p'}
    wind='s';
    ignition='c';
    alpha=pi;
    r=0.02;
    m=401
    n=401
    %m=151
    %n=151
    speed=10;
    xs=6*(m-1);
    ys=6*(n-1);
    time1=200;
    dt=1;
    plotSteps=400
    normal_spread_c=0.185060861
    normal_spread_e=1.310758329
    %normal_spread_c=1
    %normal_spread_e=1
    case {'windmill','w'}
    time1=3;
    r=0.1;
    wind='c'
    ignition='w'    
    case {'spinstardemo','s'},
    % same as spinStarDemo('low',1)
    xs=2; % domain size
    ys=2;
    m=101 % domain mesh
    n=101
    % normal propagation speed
    r=0.2;
    % max wind speed
    w=1;
    wind='c'
    ignition='s'
    time0=0
    time1=1
    case {'o','only growth'}
    speed=0
    wind='s'
    ignition='c'
    r=0.1
end

phi=zeros(m,n); % allocate level function

dx=xs/(m-1);
dy=ys/(n-1);
x=dx*[0:m-1];
y=dy*[0:n-1];

switch wind
    case {'straight','s'}
        vx=speed*cos(alpha)*ones(m,n);
        vy=speed*sin(alpha)*ones(m,n);
    case {'circular','c'}
        [xx,yy]=ndgrid(x-xs/2,y-ys/2);
        rotation=-0.75 * pi;
        rotation=0.3 * pi;
        vx=yy*rotation;
        vy=-xx*rotation;
    otherwise
        error('no wind')
end

switch ignition
    case {'star','s'}
        points = 7;
        shift = 2.5;
        scale = 0.20;
        [xx,yy]=ndgrid(x-xs/2,y-ys/2);
        [ theta, rad ] = cart2pol(xx, yy);
        phi = rad - scale * (cos(points * theta) + shift);
    case {'circle','c'}
        c=[0.8*xs,0.5*ys]; % center
        d=3*max(dx,dy);      % radius
        for i=1:m
            for j=1:n
%               phi(i,j)=sign(norm([i*dx,j*dy]-c)-d);
                phi(i,j)=(norm([i*dx,j*dy]-c)-d);
            end
        end
    case {'widmill','w'}
        phi=ones(m,n);
        phi(floor(m*0.20):ceil(m*0.45),floor(n*0.49):ceil(n*0.51))=-1;
        phi(floor(m*0.55):ceil(m*0.80),floor(n*0.49):ceil(n*0.51))=-1;
        phi(floor(m*0.49):ceil(m*0.51),floor(n*0.20):ceil(n*0.45))=-1;
        phi(floor(m*0.49):ceil(m*0.51),floor(n*0.55):ceil(n*0.80))=-1;
    otherwise
        error('no ignition')
end


%r=ones(m,n).*r + 0.04*sqrt(vx.*vx+vy.*vy); % shrinkage correction
r=ones(m,n).*r;  % make sure it is array

data2=phi;
data3=phi;
tNow=time0;
dofort=1;
domat=0;
for i=1:plotSteps
    tNext=tNow+dt
%    % call the fortran implementation
    if dofort,
        [data3]=prop_test_f(data3,tNow,tNext,...
            normal_spread_c,normal_spread_e,vx,vy,r,dx,dy);
    end
    %call own imp
    %[ux,uy]=get_advection(phi,r,vx,vy,dx,dy);
    %vis_wind(ux,uy,dx,dy)
    if domat,
        [dummy,data2]=prop_ls_cir(data2,tNow,tNext,vx,vy,r,dx,dy,@spread_rate);
    end
    %[tNow,data2]=prop_ls(data2,tNow,tNext,ux,uy,0,dx,dy);
    if dofort & domat, 
        err_fort=norm(data3(:)-data2(:))
        if err_fort>1e-10,
            i,warning('large difference between matlab and fortran')
        end
    end
    % display
    if dofort,
        vis(data3,vx,vy,dx,dy,tNow);
    else
        vis(data2,vx,vy,dx,dy,tNow);
    end
    tNow = tNext;
end
%err_fort=norm(data3(:)-data2(:))
% call the toolbox routines 
% for that need to be in toolboxls/Examples/OsherFedkiw
if 0
data = normal_advection(phi,time0,time1,vx,vy,r,dx,dy,'low');
fprintf('\n');
err_toolbox=norm(data(:)-data2(:))
end 

end

%-------------------------------------------

function vis(u,vx,vy,dx,dy,tNow)
[m,n]=size(u);
x=[0:m-1]*dx;
y=[0:n-1]*dy;
figure(1);clf
hold off
vis_type= '2d';
drawn=false;
switch vis_type
    case '3d'
    h=surf(y,x,u);
    set(h,'edgecolor','none')
    alpha(0.3)
    hold on
    contour3(y,x,u,[0 0],'k')
    hold off
    drawn=true;
    case '2d'
    [m,n]=size(u);
    contour(y,x,u,[0 0],'k')
    hold on
    vis_wind(vx,vy,dx,dy)
    axis equal
    drawn=true;
    otherwise
        disp('no visualization selected, set vis=type to 2d or 3d')
end
if drawn,
        title(sprintf('t=%g',tNow));
    drawnow;
    pause(0.3); % for ctrl-c
end
end

function [phi_out]=prop_test_f(phi,ts,te,c,e,vx,vy,r,dx,dy)
fmt='%25.15e \n';
if any(any(isnan(phi))), error('phi NaN'), end
[m,n]=size(phi);
fid = fopen('prop_test_in.txt','wt');
fprintf(fid,fmt,m,n,ts,te,c,e,r,dx,dy);
fprintf(fid,fmt,phi,vx,vy);
fclose(fid);
! prop_test_prog.exe
fid = fopen('prop_test_out.txt','rt');
if fid < 0,
    error('cannot open prop_test_out.txt')
end
[in,c]=fscanf(fid,'%g');
fclose(fid);
if c~=m*n,
    fprintf('read %g terms should have %g\n',c,m*n+1)
    error('bad number of terms in output file')
end
%t=in(1);
phi_out=reshape(in(1:c),size(phi));
if any(any(isnan(phi_out))), error('phi_out NaN'), end
end

function speed=spread_rate(r,vx,vy,nvx,nvy,scale);
    speed = r + max(vx.*nvx + vy.*nvy,0);
end

