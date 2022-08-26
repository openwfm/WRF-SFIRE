function tign=create_tign
% nodify to your liking
[ii,jj]=meshgrid(1:420);
tign=sqrt(2*(ii-200).^2+(jj-171).^2)+1;
% mesh(tign)
write_array_2d('input_tign_g',tign)
write_array_2d('input_tign_in',tign)
end
function write_array_2d(filename,a)
% write_array_2d(filename,a)
% Purpose: write 2d matrix as input to WRF-Fire
%
% Arguments
% filename  string, the name of the file
% a         2d matrix, the array to be written
%
% Example: write_array_2d('input_ht',ht)
% 
% See also: 
%    read_array_2d  read the file back by a=read_array_2d(filename)
%    image_array_2d visualize the array by image_array_2d(a)
[m,n]=size(a);
h=fopen(filename,'w');
fprintf(h,'%i\n',m,n);
fprintf(h,'%.7g\n',a');
fclose(h);
end
