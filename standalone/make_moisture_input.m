function [d,hours]=make_moisture_input(times,t2,psfc,q2,rainc)
% usage: ncload moisture100.nc; 
% make_moisture_input(times,t2,psfc,q2,rainc)
%
% input: arrays from reading WRF nc files
% times: (:,k)  character string for tiime at step k
% t2            surface temperature (K)
% psfc          surface pressure (Pa)
% q2            surface air moisture contents (kg/kg)
% rainc         accumulated rain (mm)
% i1,i2         indices of a point to take the values from, defaut 1,1
%
% output:
% d             list of time steps to omit
% hours         time from beginning of sim

if ~exist('i2','var'),i2=1;end
if ~exist('i1','var'),i1=1;end

steps=size(times,2);

% convert times

ctimes=char(times');
hours=zeros(steps,1);
for i=1:steps
    hours(i)=24*datenum(ctimes(i,:));
end
hours=hours-hours(1);

% make sure time increases and delete extra steps if not

fmt='step %i %s from start %i hours\n';
m=-1;
d=false(steps,1);
for i=2:steps
    m=max(m,hours(i-1));
    if m>=hours(i),
        if(m==hours(i-1)),
            fprintf(fmt,i-1,ctimes(i-1,:),hours(i-1))
        end
        fprintf(fmt,i,ctimes(i,:),hours(i))
        warning(['Time does not increase, deleting step ',num2str(i)])
        d(i)=true;
    elseif i<steps-1 && hours(i)-hours(i-1) ~= hours(i+1)-hours(i),
        fprintf(fmt,i,ctimes(i,:),hours(i))
        fprintf(fmt,i+1,ctimes(i+1,:),hours(i+1))
        fprintf(fmt,i+2,ctimes(i+2,:),hours(i+2))
        warning('Time step is not uniform')
    end
end
if(any(d)),
    hours(d)=[];
    ctimes(d,:)=[];
    t2(d)=[];
    q2(d)=[];
    rainc(d)=[];
    psfc(d)=[];
end

% create the input file for moisture_test.exe

mm=[hours,t2,psfc,q2,rainc];
save('moisture_input.txt','mm','-ascii')

end