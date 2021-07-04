close all;
clear;
clc;
%% load nc data
input_data_path = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/TRENDY_CLM5/input_data/';
output_data_path = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/TRENDY_CLM5/output_data/';
fig_data_path = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/TRENDY_CLM5/figures/';

% grid info
lon_size = 288;
lat_size = 192;
depth_size = 20;
time_size = 3840;

% load grid info
load([output_data_path, 'sample_grid_info.mat']);
%% processing
var_name_list = {'GPP', 'TOTSOMC', 'SOIL3C_vr'};

for ivar = 1:length(var_name_list)
    var_name = var_name_list{ivar};
    nc_var_info = ncinfo([input_data_path, 'TRENDY2020_S3_CO2ClimateLUC_Matrix.clm2.h0.', var_name, '.170001-201912.nc'], var_name);
    var_dimension = length(nc_var_info.Size);
    
    for igrid = 1:length(select_sample_index(:, 1))
        disp(['processing var ', var_name, ' grid ', num2str(igrid)]);
        
        grid_lon_loc = select_sample_index(igrid, 2);
        grid_lat_loc = select_sample_index(igrid, 3);
        if var_dimension == 3
            var_data_origin = ncread([input_data_path, 'TRENDY2020_S3_CO2ClimateLUC_Matrix.clm2.h0.', var_name, '.170001-201912.nc'],...
                var_name, [grid_lon_loc, grid_lat_loc, 1], [1, 1, time_size], [1, 1, 1]);
            var_data_grid = reshape(var_data_origin, [1, time_size]);
            
            save([output_data_path, 'grid_', num2str(igrid), '_trendy2020_clm5_s3_matrix_', var_name, '_170001_201912.mat'], 'var_data_grid');
        elseif var_dimension == 4
            var_data_origin = ncread([input_data_path, 'TRENDY2020_S3_CO2ClimateLUC_Matrix.clm2.h0.', var_name, '.170001-201912.nc'],...
                var_name, [grid_lon_loc, grid_lat_loc, 1, 1], [1, 1, depth_size, time_size], [1, 1, 1, 1]);
            var_data_grid = reshape(var_data_origin, [depth_size, time_size]);
            
            save([output_data_path, 'grid_', num2str(igrid), '_trendy2020_clm5_s3_matrix_', var_name, '_170001_201912.mat'], 'var_data_grid');
        end
    end
    
end


disp('end of programme')

