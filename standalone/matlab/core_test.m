function core_test(exe)


% defaults
if ~exist('exe','var'), exe='core_test_prog.exe';end
exe
% grid
xs=10; % domain size
ys=10;
%ys=12;
%m=111 % domain mesh
m=101
n=101
time0=0
time1=10
r=0
example='w'

switch example
    case {'circlespin','c'}
    wind='c'
    ignition='c'
    case {'gostraight','g'}
    wind='s'
    ignition='c'
    case {'windmill','w'}
    wind='c'
    ignition='w'
    time1=5;
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
end

phi=zeros(m,n); % allocate level function

dx=xs/m;
dy=ys/n;
x=dx*[0:m-1];
y=dy*[0:n-1];

switch wind
    case {'straight','s'}
        alpha=-1.3*pi
        speed=2
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
        c=[0.3*xs,0.5*ys]; % center
        d=0.5;      % radius
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

plotSteps = 10;               % How many intermediate plots to produce?
t0 = time0;                      % Start time.

% Period at which intermediate plots should be produced.
tPlot = (time1 - time0) / (plotSteps - 1);

%r=ones(m,n).*r + 0.04*sqrt(vx.*vx+vy.*vy); % shrinkage correction
r=ones(m,n).*r;  % make sure it is array
tign=-1e10*ones(m,n);
fuel_time=20*ones(m-1,n-1);
tign(phi<=0)=time0;
frac_end=ones(m-1,n-1);
frac_lost=zeros(m-1,n-1);

tNow=time0;
vis(phi,frac_lost,vx,vy,dx,dy,tNow);

for i=1:plotSteps
    tNext=min(time1, tNow + tPlot);

    % write the input file for the fortran code
    fmt='%25.15e \n';
    fid = fopen('core_test_in.txt','wt');
    if(fid<0), error('cannot open file core_test_in.txt'),end
    fprintf(fid,fmt,1,m,1,n,tNow,tNext,dx,dy);
    fprintf(fid,fmt,phi,tign,fuel_time,r,vx,vy);
    fclose(fid);

    %   call the fortran implementation
    eval(['! ',exe])

    % read the output file from the fortran code
    fid = fopen('core_test_out.txt','rt');
    if(fid<0), error('cannot open file core_test_out.txt'),end
    [in,c]=fscanf(fid,'%g');
    fclose(fid);

    % parse the file read into a single array
    count=2*m*n+2*(m-1)*(n-1);
    if(c~=count), error(sprintf('read %g terms should be %g\n',c,count)),end
    last=0;
    terms=m*n;phi=reshape(in(1+last:terms+last),size(phi));last=terms+last;
    terms=m*n;tign=reshape(in(1+last:terms+last),size(tign));last=terms+last;
    tign_d=tign;tign_d(tign<0)=NaN;
    terms=(m-1)*(n-1);frac_lost=reshape(in(1+last:terms+last),size(frac_lost));last=terms+last;
    terms=(m-1)*(n-1);frac_end=reshape(in(1+last:terms+last),size(frac_end));last=terms+last;
    if(last~=count),error('bad count'),end
    frac_tend=frac_lost/(tNext-tNow);

    % display
    clf;figure(1)
    vis(phi,-frac_tend,vx,vy,dx,dy,tNow)
    %subplot(1,3,1);vis(phi,frac_tend,vx,vy,dx,dy,tNow)
    %subplot(1,3,2);mesh(tign_d);title('tign')
    %subplot(1,3,3);mesh(phi);title('phi')
    drawnow
    pause(0.3); % for ctrl-c
    tNow = tNext;
end

end

%-------------------------------------------

function vis(u,f,vx,vy,dx,dy,tNow)
[m,n]=size(u);
x=[0:m-1]*dx;
y=[0:n-1]*dy;
hold off
vis_type= '3d';
drawn=false;
switch vis_type
    case '3d'
    figure(1);clf
    h=surf(y,x,u);
    set(h,'edgecolor','none')
    alpha(0.3)
    hold on
    contour3(y,x,u,[0 0],'k')
    drawn=true;
    figure(2);clf
    xh=[1/2:m-3/2]*dx;
    yh=[1/2:n-3/2]*dy;
    h=pcolor(xh,yh,-f');
    set(h,'edgecolor','none')
    colorbar

    case '2d'
    [m,n]=size(u);
    axis equal
    xh=[1/2:m-3/2]*dx;
    yh=[1/2:n-3/2]*dy;
    h=pcolor(xh,yh,-f');
    set(h,'edgecolor','none')
    colorbar
    hold on
    contour(y,x,u',[0 0],'k')
    s=20;
    ix=1:ceil(m/s):m;
    iy=1:ceil(n/s):n;
    [xx,yy]=ndgrid(x,y);
    quiver(xx(ix,iy),yy(ix,iy),vx(ix,iy),vy(ix,iy))
    drawn=true;
    otherwise
        disp('no visualization selected, set vis=type to 2d or 3d')
end
hold off
if drawn,
    title(sprintf('t=%g',tNow));
end
end

