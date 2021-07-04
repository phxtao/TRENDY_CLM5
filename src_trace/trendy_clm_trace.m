close all;
clear;
clc;
%% load nc data
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

resolution_lon = 360/lon_size;
resolution_lat = 180/lat_size;
lon_grid = [(0 + resolution_lon/2 : resolution_lon : 180 - resolution_lon/2), (-180 + resolution_lon/2 : resolution_lon : 0 - resolution_lon/2)]';
lat_grid = (-90 + resolution_lat/2: resolution_lat : 90 - resolution_lat/2)';

% land mask
grid_area = ncread([input_data_path, 'TRENDY2020_S3_CO2ClimateLUC_Matrix.clm2.h0.TOTSOMC_1m.170001-201912.nc'], 'area')*1000^2; % convert unit from km^2 to m^2
grid_land_fraction = ncread([input_data_path, 'TRENDY2020_S3_CO2ClimateLUC_Matrix.clm2.h0.TOTSOMC_1m.170001-201912.nc'], 'landfrac');
world_land_mask = ncread([input_data_path, 'TRENDY2020_S3_CO2ClimateLUC_Matrix.clm2.h0.TOTSOMC_1m.170001-201912.nc'], 'landmask');
grid_area_real = grid_area.*grid_land_fraction;


%% processing scalars (4dim)
var_list = {'T_SCALAR', 'W_SCALAR'};
plot_unit = 'unitless';
for ivar = 1:length(var_list)
    var_name = var_list{ivar};
    
    nc_var_info = ncinfo([input_data_path, 'TRENDY2020_S3_CO2ClimateLUC_Matrix.clm2.h0.', var_name, '.170001-201912.nc'], var_name);
    var_dimension = nc_var_info.Dimensions;
    var_unit = struct2table(nc_var_info.Attributes);
    var_unit = var_unit.Value{strcmp(var_unit.Name, 'units')};
    
    var_lat_time_series = nan(lat_size, time_size, depth_size);
    var_globe_time_series = nan(time_size, depth_size);
    
    for ilayer = 1 : depth_size
        disp(['processing var ', var_name, ' layer ', num2str(ilayer)]);
        var_data_origin = ncread([input_data_path, 'TRENDY2020_S3_CO2ClimateLUC_Matrix.clm2.h0.', var_name, '.170001-201912.nc'], var_name, [1, 1, ilayer, 1], [lon_size, lat_size, 1, time_size]);
        var_data_origin = reshape(var_data_origin, [lon_size, lat_size, time_size]).*dz(ilayer)/sum(dz(1:depth_size)); % weighting
        
        % calculate lat mean
        var_lat_time_series(:, :, ilayer) = reshape(mean(var_data_origin, 1, 'omitnan'), [lat_size, time_size]);
        var_globe_time_series(:, ilayer) = reshape(mean(var_data_origin, [1, 2], 'omitnan'), [time_size, 1]);
    end
    
    save([output_data_path, 'lat_time_series_', var_name, '.mat'], 'var_lat_time_series')
    save([output_data_path, 'var_globe_time_series_', var_name, '.mat'], 'var_globe_time_series')
    
    % calculate depth weighted mean
    var_lat_time_series_tot = sum(var_lat_time_series, 3, 'omitnan');
    var_globe_time_series_tot = sum(var_globe_time_series, 2, 'omitnan');
    
    var_lat_time_series_diff = var_lat_time_series_tot - repmat(var_lat_time_series_tot(:, 1), [1, time_size]);
    var_globe_time_series_diff = var_globe_time_series_tot - var_globe_time_series_tot(1);
    
    fun_plot(var_name, plot_unit, fig_data_path, time_series, year_start, year_end, lat_grid, var_globe_time_series_diff, var_lat_time_series_diff)
end


%% processing scalars (3dim)

var_list = {'FPI', 'O_SCALAR'};
plot_unit = 'unitless';
for ivar = 1:length(var_list)
    var_name = var_list{ivar};
    
    disp(['processing var ', var_name]);
    var_data_origin = ncread([input_data_path, 'TRENDY2020_S3_CO2ClimateLUC_Matrix.clm2.h0.', var_name, '.170001-201912.nc'], var_name);
    var_data_origin = reshape(var_data_origin, [lon_size, lat_size, time_size]);
    
    var_lat_time_series = reshape(mean(var_data_origin, 1, 'omitnan'), [lat_size, time_size]);
    var_globe_time_series = reshape(sum(var_data_origin, [1, 2], 'omitnan'), [time_size, 1]);
    
    save([output_data_path, 'lat_time_series_', var_name, '.mat'], 'var_lat_time_series')
    save([output_data_path, 'var_globe_time_series_', var_name, '.mat'], 'var_globe_time_series')
    
    var_lat_time_series_tot = var_lat_time_series;
    var_globe_time_series_tot = var_globe_time_series;
    
    var_lat_time_series_diff = var_lat_time_series_tot - repmat(var_lat_time_series_tot(:, 1), [1, time_size]);
    var_globe_time_series_diff = var_globe_time_series_tot - var_globe_time_series_tot(1);
    
    fun_plot(var_name, plot_unit, fig_data_path, time_series, year_start, year_end, lat_grid, var_globe_time_series_diff, var_lat_time_series_diff)
end


%% processing flux variables (carbon related total)
var_list = {'GPP', 'NPP', 'NEP', 'HR'};
plot_unit = 'Pg C/month';
for ivar = 1:length(var_list)
    var_name = var_list{ivar};
    
    disp(['processing var ', var_name]);
    var_data_origin = ncread([input_data_path, 'TRENDY2020_S3_CO2ClimateLUC_Matrix.clm2.h0.', var_name, '.170001-201912.nc'], var_name);
    var_data_origin = reshape(var_data_origin, [lon_size, lat_size, time_size]);
    
    var_lat_time_series = reshape(sum(var_data_origin.*repmat(grid_area_real, [1, 1, time_size]).*repmat(reshape(repmat(seconds_per_month, [1, time_size/12]), [1, 1, time_size]), [lon_size, lat_size, 1]), 1, 'omitnan'), [lat_size, time_size])/10^15; % unit: Pg
    var_globe_time_series = reshape(sum(var_data_origin.*repmat(grid_area_real, [1, 1, time_size]).*repmat(reshape(repmat(seconds_per_month, [1, time_size/12]), [1, 1, time_size]), [lon_size, lat_size, 1]), [1, 2], 'omitnan'), [time_size, 1])/10^15;
    
    save([output_data_path, 'lat_time_series_', var_name, '.mat'], 'var_lat_time_series')
    save([output_data_path, 'var_globe_time_series_', var_name, '.mat'], 'var_globe_time_series')
    
    var_lat_time_series_tot = var_lat_time_series;
    var_globe_time_series_tot = var_globe_time_series;
    
    var_lat_time_series_diff = var_lat_time_series_tot - repmat(var_lat_time_series_tot(:, 1), [1, time_size]);
    var_globe_time_series_diff = var_globe_time_series_tot - var_globe_time_series_tot(1);
    
    fun_plot(var_name, plot_unit, fig_data_path, time_series, year_start, year_end, lat_grid, var_globe_time_series_diff, var_lat_time_series_diff)
end


%% processing variables (carbon related total)
var_list = {'CWDC', 'LITR1C', 'LITR2C', 'LITR3C', 'SOIL1C', 'SOIL2C', 'SOIL3C', 'TOTSOMC', 'TOTLITC'};
plot_unit = 'Pg C';
for ivar = 1:length(var_list)
    var_name = var_list{ivar};
    
    
    disp(['processing var ', var_name]);
    var_data_origin = ncread([input_data_path, 'TRENDY2020_S3_CO2ClimateLUC_Matrix.clm2.h0.', var_name, '.170001-201912.nc'], var_name);
    var_data_origin = reshape(var_data_origin, [lon_size, lat_size, time_size]);
    
    var_lat_time_series = reshape(sum(var_data_origin.*repmat(grid_area_real, [1, 1, time_size]), 1, 'omitnan'), [lat_size, time_size])/10^15; % unit: Pg
    var_globe_time_series = reshape(sum(var_data_origin.*repmat(grid_area_real, [1, 1, time_size]), [1, 2], 'omitnan'), [time_size, 1])/10^15;
    
    save([output_data_path, 'lat_time_series_', var_name, '.mat'], 'var_lat_time_series')
    save([output_data_path, 'var_globe_time_series_', var_name, '.mat'], 'var_globe_time_series')
    
    var_lat_time_series_tot = var_lat_time_series;
    var_globe_time_series_tot = var_globe_time_series;
    
    var_lat_time_series_diff = var_lat_time_series_tot - repmat(var_lat_time_series_tot(:, 1), [1, time_size]);
    var_globe_time_series_diff = var_globe_time_series_tot - var_globe_time_series_tot(1);
    
    fun_plot(var_name, plot_unit, fig_data_path, time_series, year_start, year_end, lat_grid, var_globe_time_series_diff, var_lat_time_series_diff)
end

%% processing variables (carbon related vertical)
var_list = {'CWDC_vr', 'CWDC_Cap_vr', ...
    'LITR1C_vr', 'LITR1C_Cap_vr',...
    'LITR2C_vr', 'LITR2C_Cap_vr',...
    'LITR3C_vr', 'LITR3C_Cap_vr',...
    'SOIL1C_vr', 'SOIL1C_Capvr',...
    'SOIL2C_vr', 'SOIL2C_Capvr',...
    'SOIL3C_vr', 'SOIL3C_Capvr',...
    'SOILC_vr'...
    };
plot_unit = 'Pg C';

for ivar = 1:length(var_list)
    var_name = var_list{ivar};
    
    nc_var_info = ncinfo([input_data_path, 'TRENDY2020_S3_CO2ClimateLUC_Matrix.clm2.h0.', var_name, '.170001-201912.nc'], var_name);
    var_dimension = nc_var_info.Dimensions;
    var_unit = struct2table(nc_var_info.Attributes);
    var_unit = var_unit.Value{strcmp(var_unit.Name, 'units')};
    
    var_lat_time_series = nan(lat_size, time_size, depth_size);
    var_globe_time_series = nan(time_size, depth_size);
    
    for ilayer = 1 : depth_size
        disp(['processing var ', var_name, ' layer ', num2str(ilayer)]);
        var_data_origin = ncread([input_data_path, 'TRENDY2020_S3_CO2ClimateLUC_Matrix.clm2.h0.', var_name, '.170001-201912.nc'], var_name, [1, 1, ilayer, 1], [lon_size, lat_size, 1, time_size]);
        var_data_origin = reshape(var_data_origin, [lon_size, lat_size, time_size]).*dz(ilayer); % convert unit from gc/m3 to gc/m2
        
        % calculate stock sum
        var_lat_time_series(:, :, ilayer) = reshape(sum(var_data_origin.*repmat(grid_area_real, [1, 1, time_size]), 1, 'omitnan'), [lat_size, time_size])/10^15; % unit: Pg
        var_globe_time_series(:, ilayer) = reshape(sum(var_data_origin.*repmat(grid_area_real, [1, 1, time_size]), [1, 2], 'omitnan'), [time_size, 1])/10^15;
    end
    
    save([output_data_path, 'lat_time_series_', var_name, '.mat'], 'var_lat_time_series')
    save([output_data_path, 'var_globe_time_series_', var_name, '.mat'], 'var_globe_time_series')
    
    var_lat_time_series_tot = sum(var_lat_time_series, 3, 'omitnan');
    var_globe_time_series_tot = sum(var_globe_time_series, 2, 'omitnan');
    
    var_lat_time_series_diff = var_lat_time_series_tot - repmat(var_lat_time_series_tot(:, 1), [1, time_size]);
    var_globe_time_series_diff = var_globe_time_series_tot - var_globe_time_series_tot(1);
    
    fun_plot(var_name, plot_unit, fig_data_path, time_series, year_start, year_end, lat_grid, var_globe_time_series_diff, var_lat_time_series_diff)
end
