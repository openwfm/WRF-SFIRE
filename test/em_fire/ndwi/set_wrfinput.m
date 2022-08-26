disp('Set wrfinput for NDWI assimilation testing')
format compact

f='wrfinput_d01'

ndwi = 0.5
live_fmc_prior=0.1

w=ncread(f,'NDWI');
w=ndwi * ones(size(w));
ncreplace(f,'NDWI',w)

% set FMC_GC(:,:,5)=0.1
u=ncread(f,'FMC_GC');
u(:,:,5)=live_fmc_prior;
ncreplace(f,'FMC_GC',u);

% in namelist.fire we have set
fmc_gw05 = 0.5 ,      % fuel is consists of 0.5 live fuel
fmc_gl_stdev = 0.1,   % live fuel moisture prior has error +- 0.1
% ! observation:  FMC = 0.1 + 0.6 * NDWI +- 0.2
fmc_gl_ndwi_0 = 0.1,
fmc_gl_ndwi_rate = 0.6,
fmc_gl_ndwi_stdev = 0.2,

% so we should see
live_fmc_obs = fmc_gl_ndwi_0 + fmc_gl_ndwi_rate * ndwi
w_prior = 1/(fmc_gl_stdev^2)
w_obs   = 1/(fmc_gl_ndwi_stdev^2)
fmc_live = (live_fmc_prior*w_prior + live_fmc_obs*w_obs)/(w_prior+w_obs)
fmc_g = fmc_live * fmc_gw05    % contribution to fmc_g on fire mesh