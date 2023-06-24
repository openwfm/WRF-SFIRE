function [t,phi]=prop_ls(phi0,time0,time1,vvx,vvy,r,dx,dy)
% this version always propagates in the normal direction
% the normal direction is approximated by gradient of level function
% wind is projected into the rate of spread
% but this does not work well for the windmill..

% a jm/mkim version of propagate
% follows Osher-Fedkiw and Mitchell toolbox

% phi       level function out
% phi0      level function in
% time0     starting time
% time1     end time
% vx        component x of velocity field
% vy        component y of velocity field
% r         propagation speed in normal direction
% dx        mesh step in x direction
% dy        mesh step in y direction

% initialize
phi=phi0; % allocate level function
[m,n]=size(phi);
t=time0;  % current time
tol=300*eps;
istep=0;
msteps=1000;
diffLx=zeros(m,n);
diffLy=zeros(m,n);
diffRx=zeros(m,n);
diffRy=zeros(m,n);

% time loop
while t<time1-tol & istep < msteps
    istep=istep+1;
    % one-sided differences
    [diffLx,diffRx,diffLy,diffRy,diffCx,diffCy]=get_diff(phi,dx,dy);     
    % Godunov scheme for normal motion
    flowLx=(diffLx>=0 & diffCx>=0);
    flowRx=(diffRx<=0 & diffCx<0);
    diff2x=diffLx.*flowLx + diffRx.*flowRx; %
    flowLy=(diffLy>=0 & diffCy>=0);
    flowRy=(diffRy<=0 & diffCy<0);
    diff2y=diffLy.*flowLy + diffRy.*flowRy; %
    grad=sqrt(diff2x.*diff2x + diff2y.*diff2y);
    nz=find(grad);
    %tbound_n(1)=max(abs(r.*sqrt(diff2x(nz)))./grad(nz))/dx; % time step bnd
    %tbound_n(2)=max(abs(r.*sqrt(diff2y(nz)))./grad(nz))/dy;
    %tbound_np=r.*(abs(diff2x)./dx+abs(diff2y)./dy); % pointwise time step bnd
    %tbound_n=max(tbound_np(nz)./abs(grad(nz))); % worst case
    %tend_n =-r.*grad;   % the result, tendency
    tbound_n=0; tend_n=0;
    % replace wind by normal advection
    
    % scale [diffCx,diffCy] to normal vector
    scale=sqrt(realmin+diffCx.*diffCx+diffCy.*diffCy);
    nvx = diffCx ./ scale;
    nvy = diffCy ./ scale;

    % given spread rate plus size of velocity field projected on the normal direction
    spread_rate=ros(r,vvx,vvy,nvx,nvy);
    % empirical correction
    spread_rate=max(0,spread_rate - 0.5*max(dx,dy)*scale);

    % get the advection size in the normal direction: 
    % project [vx,vy] on the normal vector
    vx = nvx .* spread_rate;
    vy = nvy .* spread_rate;
    % standard upwinding for advection
    tend_a=-(diffLx.*max(vx,0)+diffRx.*min(vx,0)+...
        diffLy.*max(vy,0)+diffRy.*min(vy,0));
    tbound_a=max(abs(vx(:)))/dx+max(abs(vy(:)))/dy;
    % complete right-hand side
    tend=tend_n + tend_a;
    % complete time step bound
    tbound = 1/(tbound_n+tbound_a+realmin); 
    % decide on timestep
    dt=min(time1-t,0.5*tbound);
    % trailing edge correction - do not allow fireline to go backwards
    %tt=max(dx,dy);
    %ins=find(phi<=0);
    %tend(ins)=min(tend(ins),-0.5);
    % tend=min(tend,0);
    % advance
    phi=phi+dt*tend;
    t=t+dt;
end
fprintf('prop_ls: %g steps from %g to %g\n',istep,time0,t)
end % of the function prop_ls

function spread_rate=ros(r,vx,vy,nvx,nvy)
% get the spread rate from wind vector and
% fireline normal vector
spread_rate=max(r+vx.*nvx+vy.*nvy,0);
end
