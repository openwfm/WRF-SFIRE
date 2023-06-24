function [diffLx,diffRx,diffLy,diffRy,diffCx,diffCy]=get_diff(phi,dx,dy)
% get one-sided differences of phi
[m,n]=size(phi);
diffLx(2:m,:)  =(phi(2:m,:)-phi(1:m-1,:))/dx    ;
diffLx(1  ,:)  =diffLx(2,:);
diffRx(1:m-1,:)=diffLx(2:m,:);    
diffRx(m    ,:)=diffLx(m,:);
diffLy(:,2:n)  =(phi(:,2:n)-phi(:,1:n-1))/dy ;   
diffLy(:  ,1)  =diffLy(:,2);
diffRy(:,1:n-1)=diffLy(:,2:n);    
diffRy(:  ,n)  =diffLy(:,n);
diffCx=0.5*(diffLx+diffRx);
diffCy=0.5*(diffLy+diffRy);
end