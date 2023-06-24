function fuel_fraction=fuel_burnt(lfn,tign,tnow,fd,fuel_time)
% lfn       - level set function, size (2,2)
% tign      - time of ignition, size (2,2)
% tnow      - time now
% fd        - mesh cell size (2)
% fuel_time - time constant of fuel
%
% output: approximation of
%                              /\
%                     1        |                   T(x,y)
%  fuel_burnt  =  ----------   |  ( 1  -  exp( - ----------- ) ) dxdy
%                 fd(1)*fd(2)  |                  fuel_time
%                             \/
%                        (x,y) in R; 
%                         L(x,y)<=0
%
% where R is the rectangle [0,fd(1)] * [0,fd(2)]
%       T is given by tign and L is given by lfn at the 4 vertices of R
%
% note that fuel_left = 1 - fuel_burnt
%
% Requirements:
%       exact in the case when T and L are linear 
%       varies continuously with input
%       value of tign-tnow where lfn>0 ignored
%       assume T=0 when lfn=0
figs=1:4;
for i=figs,figure(i),clf,hold off,end

if all(lfn(:)>=0),
    % nothing burning, nothing to do - most common case, put it first
    fuel_fraction=0;
elseif all(lfn(:)<=0),
    % all burning
    % T=u(1)*x+u(2)*y+u(3)
    % set up least squares(A*u-v)*inv(C)*(A*u-v)->min
    A=[0,      0,    1;
       fd(1),  0,    1;
       0,    fd(2),  1;
       fd(1),fd(2),  1];
    v=tnow-tign(:);    % time from ignition
    rw=2*ones(4,1);
    u = lscov(A,v,rw);  % solve least squares to get coeffs of T
    residual=norm(A*u-v); % zero if T is linear
    % integrate
    uu = -u / fuel_time;
    fuel_fraction = 1 - exp(uu(3)) * intexp(uu(1)*fd(1)) * intexp(uu(2)*fd(2));
    if(fuel_fraction<0 | fuel_fraction>1),
        warning('fuel_fraction should be between 0 and 1')
    end
else
    % part of cell is burning - the interesting case
    % set up a list of points with the corresponding values of T and L
    xylist=zeros(8,2);    % allocate list of points
    Tlist=zeros(8,1);
    Llist=zeros(8,1);
    xx=[-fd(1) fd(1) fd(1) -fd(1) -fd(1)]/2; % make center (0,0)
    yy=[-fd(2) -fd(2) fd(2) fd(2) -fd(2)]/2; % cyclic, counterclockwise
    ii=[1 2 2 1 1]; % indices of corners, cyclic, counterclockwise
    jj=[1 1 2 2 1];
    points=0;
    for k=1:4
        lfn0=lfn(ii(k),jj(k));
        lfn1=lfn(ii(k+1),jj(k+1));
        if(lfn0<=0),
            points=points+1;
            xylist(points,:)=[xx(k),yy(k)]; % add corner to list
            Tlist(points)=tnow-tign(ii(k),jj(k)); % time since ignition
            Llist(points)=lfn0;
        end
        if(lfn0*lfn1<0),
            points=points+1;
            % coordinates of intersection of fire line with segment k k+1
            %lfn(t)=lfn0 + t*(lfn1-lfn0)=0
            t=lfn0/(lfn0-lfn1);
            x0=xx(k)+(xx(k+1)-xx(k))*t;
            y0=yy(k)+(yy(k+1)-yy(k))*t;
            xylist(points,:)=[x0,y0];
            Tlist(points)=0; % now at ignition
            Llist(points)=0; % at fireline
        end
    end
    % make the lists circular and trim to size
    Tlist(points+1)=Tlist(1);
    Tlist=Tlist(1:points+1);
    Llist(points+1)=Llist(1);
    Llist=Llist(1:points+1);
    xylist(points+1,:)=xylist(1,:);
    xylist=xylist(1:points+1,:);
    for k=1:5,lfnk(k)=lfn(ii(k),jj(k));end
    figure(1)
    plot3(xx,yy,lfnk,'k')
    hold on
    plot3(xylist(:,1),xylist(:,2),Tlist,'m--o')
    plot3(xylist(:,1),xylist(:,2),Llist,'g--o')
    plot(xx,yy,'b')
    plot(xylist(:,1),xylist(:,2),'-.r+')
    legend('lfn on whole cell','tign-tnow on burned area',...
        'lfn on burned area', 'cell boundary','burned area boundary')
    patch(xylist(:,1),xylist(:,2),zeros(points+1,1),250,'FaceAlpha',0.3)
    patch(xylist(:,1),xylist(:,2),Tlist,100,'FaceAlpha',0.3)
    hold off,grid on,drawnow,pause(0.5)
    % set up least squares 
    A=[xylist(1:points,1:2),ones(points,1)];
    v=Tlist(1:points);
    for i=1:points
        for j=1:points
            if(j~=i),
                dist(j)=norm(xylist(i,:)-xylist(j,:));
            else
                dist(j)=max(fd);  % large
            end
        end
        rw(i)=1+min(dist);  % weight = 1+min dist from other nodes
    end
    u = lscov(A,v,rw);  % solve least squares to get coeffs of T
    residual=norm(A*u-v); % should be zero if T and L are linear
    nr=sqrt(u(1)*u(1)+u(2)*u(2));
    c=u(1)/nr;
    s=u(2)/nr;
    Q=[c,  s; -s, c];  % rotation such that Q*u(1:2)=[something;0]
    ut=Q*u(1:2);
    errQ=ut(2);  % should be zero
    ae=-ut(1)/fuel_time;
    ce=-u(3)/fuel_time;     %  -T(xt,yt)/fuel_time=ae*xt+ce
    xytlist=xylist*Q'; % rotate the points in the list
    xt=sort(xytlist(1:points,1)); % sort ascending in x
    fuel_fraction=0;
    for k=1:points-1
        % integrate the vertical slice from xt1 to xt2
        figure(2)
        plot(xytlist(:,1),xytlist(:,2),'-o'),grid on,hold on
        xt1=xt(k);
        xt2=xt(k+1);
        plot([xt1,xt1],sum(fd)*[-1,1]/3,'y')
        plot([xt2,xt2],sum(fd)*[-1,1]/3,'y')
        if(xt2-xt1>100*eps*max(fd)), % slice of nonzero width
            % find slice height as h=ah*x+ch
            upper=0;lower=0;
            ah=0;ch=0;
            for s=1:points % pass counterclockwise
                xts=xytlist(s,1); % start point of the line
                yts=xytlist(s,2);
                xte=xytlist(s+1,1); % end point of the line 
                yte=xytlist(s+1,2);
                if (xts>xt1 & xte > xt1) | (xts<xt2 & xte < xt2)
                    plot([xts,xte],[yts,yte],'--k')
                    % that is not the one
                else % line y=a*x+c through (xts,yts), (xte,yte)
                    a=(yts-yte)/(xts-xte);
                    c=(xts*yte-xte*yts)/(xts-xte);
                    if xte<xts % upper boundary
                        aupp=a;
                        cupp=c;
                        plot([xts,xte],[yts,yte],'g')
                        ah=ah+a;
                        ch=ch+c;
                        upper=upper+1;
                    else % lower boundary
                        alow=a;
                        clow=c;
                        plot([xts,xte],[yts,yte],'m')
                        lower=lower+1;
                    end
                end
            end
            if(lower~=1|upper~=1),
                error('slice does not have one upper and one lower line')
            end
            ah=aupp-alow;
            ch=cupp-clow;
            % debug only
            patch([xt1,xt2,xt2,xt1],...
                [alow*[xt1,xt2]+clow,aupp*[xt2,xt1]+cupp,],k*10)
            axis equal,drawnow, pause(0.5)
            figure(3)
            x=[xt1:(xt2-xt1)/100:xt2];
            plot(x,1-exp(ae*x+ce),x,ah*x+ch,x,...
                (ah*x+ch).*(1-exp(ae*x+ce)),[xt1,xt2],[0,0],'+k')
            xlabel x,legend('burned frac','slice height','height*burned','dividers')
            hold on,grid on,drawnow,pause(0.5)
            % integrate (ah*x+ch)*(1-exp(ae*x+ce) from xt1 to xt2
            % numerically sound for ae->0, ae -> infty
            % this can be important for different model scales
            % esp. if someone runs the model in single precision!!
            % s1=int((ah*x+ch),x,xt1,xt2)
            s1=(xt2-xt1)*((1/2)*ah*(xt2+xt1)+ch);
            % s2=int((ch)*(-exp(ae*x+ce)),x,xt1,xt2)
            ceae=ce/ae;
            s2=-ch*exp(ae*(xt1+ceae))*(xt2-xt1)*intexp(ae*(xt2-xt1));
            % s3=int((ah*x)*(-exp(ae*x+ce)),x,xt1,xt2)
            % s3=int((ah*x)*(-exp(ae*(x+ceae))),x,xt1,xt2)
            % expand in Taylor series around ae=0
            % collect(expand(taylor(int(x*(-exp(ae*(x+ceae))),x,xt1,xt2)*ae^2,ae,4)/ae^2),ae)
            % =(1/8*xt1^4+1/3*xt1^3*ceae+1/4*xt1^2*ceae^2-1/8*xt2^4-1/3*xt2^3*ceae-1/4*xt2^2*ceae^2)*ae^2
            %     + (-1/3*xt2^3-1/2*xt2^2*ceae+1/3*xt1^3+1/2*xt1^2*ceae)*ae 
            %     + 1/2*xt1^2-1/2*xt2^2
            %
            % coefficient at ae^2 in the expansion, after some algebra
            a2=(xt1-xt2)*((1/4)*(xt1+xt2)*ceae^2+(1/3)*(xt1^2+xt1*xt2+xt2^2)*ceae+(1/8)*(xt1^3+xt1*xt2^2+xt1^2*xt2+xt2^3));
            d=ae^4*a2;
            if abs(d)>eps
                % since ae*xt1+ce<=0 ae*xt2+ce<=0 all fine for large ae
                % for ae, ce -> 0 rounding error approx eps/ae^2
                s3=(exp(ae*(xt1+ceae))*(ae*xt1-1)-exp(ae*(xt2+ceae))*(ae*xt2-1))/ae^2;
                % we do not worry about rounding as xt1 -> xt2, then s3 -> 0
            else
                % coefficient at ae^1 in the expansion
                a1=(xt1-xt2)*((1/2)*ceae*(xt1+xt2)+(1/3)*(xt1^2+xt1*xt2+xt2^2));
                % coefficient at ae^0 in the expansion for ae->0
                a0=(1/2)*(xt1-xt2)*(xt1+xt2);
                s3=a0+a1*ae+a2*ae^2; % approximate the integral
            end
            s3=ah*s3;
            fuel_fraction=fuel_fraction+s1+s2+s3;
            if(fuel_fraction<0 | fuel_fraction>(fd(1)*fd(2))),
                fuel_fraction,(fd(1)*fd(2)),s1,s2,s3
                warning('fuel_fraction should be between 0 and 1')
            end
        end
    end
    fuel_fraction=fuel_fraction/(fd(1)*fd(2));
end
for i=figs,figure(i),hold off,end
end % function

function s=intexp(ab)
% function s=intexp(ab)
% s = (1/b)*int(exp(a*x),x,0,b)  ab = a*b
%    = 1 if a==0
if eps < abs(ab)^3/6,
    s = (exp(ab)-1)/(ab);  % rounding error approx eps/(a*b) for small a*b
else
    s = 1 + ab/2+ab^2/6;   % last term (a*b)^2/6 for small a*b
end
end
