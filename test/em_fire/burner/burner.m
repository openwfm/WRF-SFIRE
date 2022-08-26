function input_tign_g
% create file input_tign_g for reading by ideal.exe

% edit these numbers as needed
burner_size_m     = [750,25]   %
domain_size_m     = [3000,6000]
mesh_step_m       = [25,25]
burner_dist_m     = 400
burner_start_s    = 10
burner_end_s      = 1e6

% do not edit under this line

domain_size_cells = div(domain_size_m,mesh_step_m)
burner_size_cells = div(burner_size_m,mesh_step_m)

% place burner burner_dist_m from start in direction 1 and symmetrically in direction 2
burner_dist_cells = [(domain_size_cells(1)-burner_size_cells(1))/2,...
    div(burner_dist_m,mesh_step_m(2))]
burner_start_cells = burner_dist_cells+1;
burner_end_cells   = burner_dist_cells+burner_size_cells;

burner_mask = zeros(domain_size_cells);
burner_mask(burner_start_cells(1):burner_end_cells(1),...
            burner_start_cells(2):burner_end_cells(2)) = 1;
spy(burner_mask)
fprintf('burner cells %i from %i\n',nnz(burner_mask),prod(domain_size_cells))

burner_tign = zeros(domain_size_cells);
burner_tign(burner_mask(:)==0) = burner_end_s;
burner_tign(burner_mask(:) >0) = burner_start_s;

file = 'input_tign_g';
write_array_2d(file,burner_tign')
fprintf('Burner ignition written to file %s for ideal.exe\n',file)
% ncreplace('wrfinput_d01','TIGN_G',burner_tign)

end

function size_cells=div(size_m,step_m)
size_cells=size_m./step_m;
if any((size_cells - round(size_cells))>10*eps),
    warning('given size in m should be multiple of step_size_m')
end
size_cells = round(size_cells)
end    
