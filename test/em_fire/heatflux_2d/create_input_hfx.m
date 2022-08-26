wd=pwd;                                          % remember work directory
cd ..
startup                                          % set up search path
cd(wd)                                           % return
n1 = 420; n2=420;                                % grid dimensions
width=50;                                        % bump width in grid cells
[i,j]=ndgrid(1:n1,1:n2);                         % create index arrays 
hfx=1e6*exp(-((i-n1/2).^2 + 2*((j-n2/3).^2))/width^2);   % elliptic gaussian bump 
write_array_2d('input_hfx',hfx)                  % write for ideal.exe
mesh(hfx)                                        