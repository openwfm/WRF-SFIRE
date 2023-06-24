function model_test
%a=load('model_test_out.txt');
f=fopen('model_test_out.txt','r');
if f<0, error('cannot open input file'), end
dx=next_mat(f);
dy=next_mat(f);
eof=0;
i=0;
t=next_mat(f);
while ~isempty(t),
    i=i+1;
    lfn=next_mat(f);
    if(i>1),
        fprintf('step %i time %g ',i,t),
        fprintf('lfn change %g\n',...
            max(max(abs(lfn-lfn_last))));
    end
    lfn_last=lfn;
    tign=next_mat(f);
    tign(lfn>0)=NaN;
    vx=next_mat(f);
    vy=next_mat(f);
    flux=next_mat(f);
    vis(lfn,flux,vx,vy,dx,dy,t)
    drawnow
    pause(0.3)
    t=next_mat(f);
end
fclose(f);
    
end

function vis(u,f,vx,vy,dx,dy,tNow)
% u = level function
% f = reaction intensity or flux
clf
[m,n]=size(u);
x=[0:m-1]*dx;
y=[0:n-1]*dy;
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
    drawn=true;
    case '2d'
    xh=[1/2:m-1/2]*dx;
    yh=[1/2:n-1/2]*dy;
    h=pcolor(xh,yh,f');
    % shading('interp')
    set(h,'edgecolor','none')
    colorbar
    hold on
    contour(y,x,u',[0 0],'k')
    s=20;
    ix=1:ceil(m/s):m;
    iy=1:ceil(n/s):n;
    [xx,yy]=ndgrid(x,y);
    quiver(xx(ix,iy),yy(ix,iy),vx(ix,iy),vy(ix,iy))
    axis image
    drawn=true;
    otherwise
        disp('no visualization selected, set vis=type to 2d or 3d')
end
hold off
if drawn,
    title(sprintf('t=%g',tNow));
end
end
