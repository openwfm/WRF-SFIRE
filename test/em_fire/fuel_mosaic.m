
% Generate a mosaic fuel file in ./mosaic/input_fc to test
% non-uniform fuel behavior.  This script constructs an array
% containing tiles of variable fuel types.  See below for 
% parameters.

% size of fire domain (from namelist.input)
nx=200;
ny=320;

% size of fuel tiles
tx=20;
ty=20;

% How to fill in fuel categories in tiles.
% This is a 2d array of any size 
% containing fuel types.  This array will be
% repeated if it is smaller than ceil(nx/tx), ceil(ny/ty).
% For example:
%
%  fuels = [ 1 2
%            3 4 ];
% 
% Then the output array will look like:
%
% [ 1 2 1 2 1 2 1 2 ...
%   3 4 3 4 3 4 3 4 ...
%   1 2 1 2 1 2 1 2 ...
%   3 4 3 4 3 4 3 4 ... 
%   ...                 ]
% 
% If empty, then use random categories.
fuels=[];


% fuel categories to use if randomizing
random_fuels=1:13;

% output file
filename='mosaic/input_fc';

%%

% get tile array size
mx=ceil(nx/tx);
my=ceil(ny/ty);

% construct fuel array
if isempty(fuels)
    fuels=ceil(length(random_fuels)*rand( mx,my ));
    fuels=random_fuels(fuels);
end

% extend fuels array if necessary
fuels=repmat(fuels,ceil(nx/mx),ceil(ny/my));
fuels=fuels(1:mx,1:my);

% create output array
fc=zeros(nx,ny);

% fill in fuel categories
for j=1:my
    js=(j-1) * ty + 1;
    je=(j)   * ty;
    for i=1:mx
        is=(i-1) * tx + 1;
        ie=(i)   * tx;
        fc(is:ie,js:je)=fuels(i,j);
    end
end

pcolor(fc); shading interp; colorbar;

write_array_2d(filename,fc);
%system(sprintf('mv fc.txt %s',filename));

