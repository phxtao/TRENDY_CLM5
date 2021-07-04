close all;
clear;
clc;
%%
year_start = 1700;
year_end = 2019;
seconds_per_month = [31, 30, 31, 28, 31, 30, 31, 31, 30, 31, 30, 31]*24*3600;

input_data_path = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/TRENDY_CLM5/input_data/';
output_data_path = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/TRENDY_CLM5/output_data/';
fig_data_path = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/TRENDY_CLM5/figures/';

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

% grid info
lon_size = 288;
lat_size = 192;
depth_size = 20;
time_size = 3840;
time_series = (year_start:((year_end - year_start)/(time_size - 1)):year_end)';

% lon and lat maps
resolution_lon = 360/lon_size;
resolution_lat = 180/lat_size;
lon_grid = [(0 + resolution_lon/2 : resolution_lon : 180 - resolution_lon/2), (-180 + resolution_lon/2 : resolution_lon : 0 - resolution_lon/2)]';
lat_grid = (-90 + resolution_lat/2: resolution_lat : 90 - resolution_lat/2)';

lon_map = nan(length(lon_grid), length(lat_grid));
lat_map = nan(length(lon_grid), length(lat_grid));

for ilon = 1:length(lon_grid)
    lon_map(ilon, :) = ilon;
    lat_map(ilon, :) = 1:length(lat_grid);
end


% land mask
world_land_mask = ncread([input_data_path, 'TRENDY2020_S3_CO2ClimateLUC_Matrix.clm2.h0.TOTSOMC_1m.170001-201912.nc'], 'landmask');
world_land_mask(:, (lat_grid < -56 | lat_grid > 80)) = 0;


% gpp info unit: gc/m2/s
gpp_data_origin = ncread([input_data_path, 'TRENDY2020_S3_CO2ClimateLUC_Matrix.clm2.h0.GPP.170001-201912.nc'], 'GPP');
gpp_data_1700 = gpp_data_origin(:, :, 1:12);

for imonth = 1:12
    gpp_data_1700(:, :, imonth) = gpp_data_1700(:, :, imonth)*seconds_per_month(imonth); % unit gc/m2
end

gpp_data_1700 = sum(gpp_data_1700, 3);
gpp_data_1700(gpp_data_1700 == 0) = nan;


%% select random grids
rng(8)

sample_size = 10;
grid_interval = 5;

select_sample_index = nan(sample_size, 3);

valid_grid_index = find(gpp_data_1700 > 100 & gpp_data_1700 < quantile(gpp_data_1700, 0.95, 'all') & world_land_mask == 1);

for isample = 1:sample_size
    sample_middle = randsample(valid_grid_index, 1, false);
    sample_middle_lon = lon_map(sample_middle);
    sample_middle_lat = lat_map(sample_middle);
    
    select_sample_index(isample, :) = [sample_middle, sample_middle_lon, sample_middle_lat];
    
    world_land_mask(lon_map >= sample_middle_lon - grid_interval & lon_map <= sample_middle_lon + grid_interval & ...
        lat_map >= sample_middle_lat - grid_interval & lat_map <= sample_middle_lat + grid_interval ) = nan;
    
    valid_grid_index = find(gpp_data_1700 > 50 & gpp_data_1700 < quantile(gpp_data_1700, 0.99, 'all') & world_land_mask == 1);
end


site_label = cell(sample_size, 1);
for isite = 1:sample_size
    site_label{isite} = [num2str(isite), ' Lon: ', num2str(round(lon_grid(select_sample_index(isite, 2)), 2)), ' Lat: ', num2str(round(lat_grid(select_sample_index(isite, 3)), 2)),];
end


save([output_data_path, 'sample_grid_info.mat'], 'select_sample_index');

% distribution
load coastlines
plot(coastlon, coastlat, 'k', 'LineWidth', 0.5);
xlim([-180 180]);
ylim([-60 80]);
hold on
scatter(lon_grid(select_sample_index(:, 2)), lat_grid(select_sample_index(:, 3)), 600, '.r')
hold on
text(lon_grid(select_sample_index(:, 2)), lat_grid(select_sample_index(:, 3)), site_label, 70, 'color', 'blue')
title('', 'FontSize', 20)
set(gca, 'xcolor', 'none', 'ycolor', 'none')

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 12 5];

print([fig_data_path, 'sample_grid_dist.tif'], '-dtiffn', '-r300')
close

