clear;
clc;
%%

warning('off');
format long e;

disp([datestr(now,'HH:MM:SS'), ' programme started'])

grid_num = 7;

% model_src_dir = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/RADIOCARBON/demo_harvard_forest/src_clm_cen';
model_src_dir = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/TRENDY_CLM5/src_clm_cen';
cd(model_src_dir)

% mac
data_path = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/';

%% simulations from CLM5
cesm2_case_name = 'sasu_f05_g16_checked_step4';
start_year = 661;
end_year = 680;
model_name = 'cesm2_clm5_cen_vr_v2';

%% soil depth information
% width between two interfaces
dz = [2.000000000000000E-002, 4.000000000000000E-002, 6.000000000000000E-002, ...
    8.000000000000000E-002, 0.120000000000000, 0.160000000000000, ...
    0.200000000000000, 0.240000000000000, 0.280000000000000, ...
    0.320000000000000, 0.360000000000000, 0.400000000000000, ...
    0.440000000000000, 0.540000000000000, 0.640000000000000, ...
    0.740000000000000, 0.840000000000000, 0.940000000000000, ...
    1.04000000000000, 1.14000000000000, 2.39000000000000, ...
    4.67553390593274, 7.63519052838329, 11.1400000000000, ...
    15.1154248593737]';

% depth of the interface
zisoi = [2.000000000000000E-002, 6.000000000000000E-002, ...
    0.120000000000000, 0.200000000000000, 0.320000000000000, ...
    0.480000000000000, 0.680000000000000, 0.920000000000000, ...
    1.20000000000000, 1.52000000000000, 1.88000000000000, ...
    2.28000000000000, 2.72000000000000, 3.26000000000000, ...
    3.90000000000000, 4.64000000000000, 5.48000000000000, ...
    6.42000000000000, 7.46000000000000, 8.60000000000000, ...
    10.9900000000000, 15.6655339059327, 23.3007244343160, ...
    34.4407244343160, 49.5561492936897]';

% depth of the node
zsoi = [1.000000000000000E-002, 4.000000000000000E-002, 9.000000000000000E-002, ...
    0.160000000000000, 0.260000000000000, 0.400000000000000, ...
    0.580000000000000, 0.800000000000000, 1.06000000000000, ...
    1.36000000000000, 1.70000000000000, 2.08000000000000, ...
    2.50000000000000, 2.99000000000000, 3.58000000000000, ...
    4.27000000000000, 5.06000000000000, 5.95000000000000, ...
    6.94000000000000, 8.03000000000000, 9.79500000000000, ...
    13.3277669529664, 19.4831291701244, 28.8707244343160, ...
    41.9984368640029]';

% depth between two node
dz_node = zsoi - [0; zsoi(1:end-1)];

% cesm2 resolution
cesm2_resolution_lat = 180/384;
cesm2_resolution_lon = 360/576;
lon_grid = (-180 + cesm2_resolution_lon/2 : cesm2_resolution_lon : 180 - cesm2_resolution_lon/2)';
lat_grid = (90 - cesm2_resolution_lat/2 : -cesm2_resolution_lat : -90 + cesm2_resolution_lat/2)';

month_num = 12;
soil_cpool_num = 7;
soil_decom_num = 20;

%% Simulation input
year_start = 1900;
year_end = 2010;

soc_time_series = (year_start:(year_end-year_start)/((year_end-year_start+1)*12-1):year_end)';

% load simulation from CESM2 (steady state)
load([data_path, 'radiocarbon/demo_harvard_forest/input_data/forcing_steady_state_harvard_forest.mat']);

% load simulation from CESM2 (transient)
load([data_path, 'radiocarbon/demo_harvard_forest/input_data/forcing_transient_harvard_forest.mat']);

obs_tot_soc = reshape(frocing_transient.soc_stock, [length(soc_time_series), 1]);


%% load DA results
simu2_num = 50000;
npara = 21;

load([data_path, 'TRENDY_CLM5/output_data/da_reconstruct/clm5_HF_point_reconstruct.mat'])
soc_mod_trace = reshape(da_summary.soc_mod_trace(1, :, :), [simu2_num, length(soc_time_series)]);
posterior_para = reshape(da_summary.parameters_keep2, [npara, simu2_num]);

valid_num = length(find(isnan(soc_mod_trace(:, 1)) == 0));

soc_mod_mean = mean(soc_mod_trace(round(valid_num/2):valid_num, :), 1);

scatter(obs_tot_soc, soc_mod_mean, 'filled')
hold on
plot([min(obs_tot_soc), max(obs_tot_soc)], [min(obs_tot_soc), max(obs_tot_soc)])
