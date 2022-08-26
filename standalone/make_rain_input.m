function make_rain_input(n,mmph,hours)
% create moisture_input.txt for moisture_test.exe 
% usage: make_rain_input(n,mmph)
% input: 
%    n        number of time steps
%    mmph     rain intensity mm/h
%    hours    length of simulation in hours
% example
% make_rain_input(n,100,14) should give 1.723 for 14-hour wetting time lag
%
% note: the last line in the output from moisture_test.exe
%       should not depend on n

mm=[hours*[0:n]'/n,ones(n+1,1)*[300,1e5,0.01],mmph*hours*[0:n]'/n];
save('moisture_input.txt','mm','-ascii')

end