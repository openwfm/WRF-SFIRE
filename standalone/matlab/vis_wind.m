function vis_wind(vx,vy,dx,dy)
    [m,n]=size(vx);
    x=[0:m-1]*dx;
    y=[0:n-1]*dy;
    s=20;
    ix=1:ceil(m/s):m;
    iy=1:ceil(n/s):n;
    [xx,yy]=ndgrid(x,y);
    quiver(yy(ix,iy),xx(ix,iy),vy(ix,iy),vx(ix,iy))
end