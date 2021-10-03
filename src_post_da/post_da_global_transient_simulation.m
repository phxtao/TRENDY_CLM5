clear
clc;
%%
scenario_name = 'S3_CO2ClimateLUC_Matrix';

worker_num = 5;

random_num = 50;

itest = 5;
% scaling_factor_list = [2.4, 1.4, 0.4, 0.1];
scaling_factor_list = [0.5, 0.2, 0.1, 0.05, 0.01];
data_path = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/';
model_src_dir = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/TRENDY_CLM5/src_clm_cen_global';
cd(model_src_dir)

npara = 26;
nsimu2 = 50000;
%% GR statistics
% for igrid = 1:length(scaling_factor_list)
%     for ipara = 1:npara
%         valid_post_para = nan(nsimu2, 6);
%         
%         for iworker = 1:worker_num
%             post_para = load([data_path, 'trendy_clm5/output_data/da_reconstruct/trendy_clm5_global_reconstruct_test', num2str(itest), '_s3_matrix_', num2str(scaling_factor_list(igrid)), '_worker_', num2str(iworker), '_parameters_keep2.mat']);
%             post_para = post_para.parameters_keep2';
%             valid_num = length(find(isnan(post_para(:, 1)) == 0));
%             
%             if valid_num/nsimu2 > 0. && valid_num/nsimu2 < 1
%                 valid_post_para(1:length(round(valid_num/2):valid_num), iworker) = post_para(round(valid_num/2):valid_num, ipara);
%             end
%         end
%         subplot(4, 7, ipara)
%         plot(valid_post_para)
%     end
% end

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

%% forward simulations
% cesm2 resolution
cesm2_resolution_lat = 180/384;
cesm2_resolution_lon = 360/576;
lon_grid = (-180 + cesm2_resolution_lon/2 : cesm2_resolution_lon : 180 - cesm2_resolution_lon/2)';
lat_grid = (90 - cesm2_resolution_lat/2 : -cesm2_resolution_lat : -90 + cesm2_resolution_lat/2)';

month_num = 12;
soil_cpool_num = 7;
soil_decom_num = 20;

nbedrock = 20;

trendy_year_start = 1700;
trendy_year_end = 2019;
trendy_time_series = (trendy_year_start+1/month_num:1/month_num:(trendy_year_end+1))';
days_per_month = [31, 30, 31, 28, 31, 30, 31, 31, 30, 31, 30, 31];
days_per_month_series = repmat(days_per_month, [1, (trendy_year_end - trendy_year_start + 1)]);

time_loop = 319;
time_gap = 321;
da_year_start = 1701:(time_gap-1):2019;
da_year_end = da_year_start+time_loop-1;
period_num = length(da_year_start);

da_start_loc = (da_year_start - trendy_year_start + 1 - 1)*month_num + 1;
da_end_loc = (da_year_end - trendy_year_start + 1)*month_num;
% loc of selected time period in the whole simulation
selected_loc = nan(length(da_year_start)*time_loop*month_num, 1);
for iperiod = 1:period_num
    selected_loc(time_loop*month_num*(iperiod-1)+1:time_loop*month_num*iperiod, :) = da_start_loc(iperiod):da_end_loc(iperiod);
end
% selected time period in the whole simulation
da_time_series = trendy_time_series(selected_loc);

% land area
load([data_path, 'trendy_clm5/input_data/trendy_clm5_global_s3_matrix/global_trendy2020_clm5_land_area.mat']);
global_land_area = sum(sum(land_area, 'omitnan'), 'omitnan'); % unit m2
% load simulation from TRENDY CLM5
input_var_list = {'ALTMAX', 'FPI', 'NPP', 'TOTVEGC', 'O_SCALAR', 'W_SCALAR', 'T_SCALAR', 'COL_FIRE_CLOSS', 'TOTLITC', 'CWDC'};

for ivar = 1:length(input_var_list)
    load([data_path, 'trendy_clm5/input_data/trendy_clm5_global_s3_matrix/global_trendy2020_clm5_', scenario_name, '_', input_var_list{ivar}, '_170001_201912.mat']);
    area_weighted_mean = area_weighted_sum'/global_land_area;
    eval(['frocing_transient.', input_var_list{ivar}, ' = area_weighted_mean(:, selected_loc);']);
    eval(['frocing_steady_state.', input_var_list{ivar}, ' = area_weighted_mean(:, 1:20*month_num);']);
end
frocing_transient.ALTMAX_LAST_YEAR = [frocing_transient.ALTMAX(1), frocing_transient.ALTMAX(1:end-1)];
frocing_transient.nbedrock = nbedrock;

frocing_steady_state.ALTMAX_LAST_YEAR = [frocing_steady_state.ALTMAX(1), frocing_steady_state.ALTMAX(1:end-1)];
frocing_steady_state.nbedrock = nbedrock;

frocing_transient.litter_fall = frocing_transient.NPP - (frocing_transient.TOTVEGC - [frocing_transient.TOTVEGC(1), frocing_transient.TOTVEGC(1:end-1)])./(days_per_month_series(selected_loc)*24*3600); %unit gc/m2/s
frocing_transient.litter_fall(1) = frocing_transient.litter_fall(13);

frocing_steady_state.litter_fall = frocing_steady_state.NPP - (frocing_steady_state.TOTVEGC - [frocing_steady_state.TOTVEGC(1), frocing_steady_state.TOTVEGC(1:end-1)])./(days_per_month_series(1:20*month_num)*24*3600); %unit gc/m2/s
frocing_steady_state.litter_fall(1) = frocing_steady_state.litter_fall(13);

% original simulated soc, unit gc/m2
obs_tot_soc_transient = load([data_path, 'trendy_clm5/input_data/trendy_clm5_global_s3_matrix/global_trendy2020_clm5_', scenario_name, '_TOTSOMC_170001_201912.mat']);
obs_tot_soc_transient = obs_tot_soc_transient.area_weighted_sum(selected_loc)'/global_land_area;
% original simulated litter, unit gc/m2
obs_tot_litter_transient = load([data_path, 'trendy_clm5/input_data/trendy_clm5_global_s3_matrix/global_trendy2020_clm5_', scenario_name, '_TOTLITC_170001_201912.mat']);
obs_tot_litter_transient = obs_tot_litter_transient.area_weighted_sum(selected_loc)'/global_land_area;
% original simulated cwdc, unit gc/m2
obs_tot_cwd_transient = load([data_path, 'trendy_clm5/input_data/trendy_clm5_global_s3_matrix/global_trendy2020_clm5_', scenario_name, '_CWDC_170001_201912.mat']);
obs_tot_cwd_transient = obs_tot_cwd_transient.area_weighted_sum(selected_loc)'/global_land_area;
% original simulated heterotrophic respiration, unit gc/m2/s
obs_tot_hr_transient = load([data_path, 'trendy_clm5/input_data/trendy_clm5_global_s3_matrix/global_trendy2020_clm5_', scenario_name, '_HR_170001_201912.mat']);
obs_tot_hr_transient = obs_tot_hr_transient.area_weighted_sum(selected_loc)'/global_land_area;

% initial cabron pool sizes at the beginning of transient simulation, unit gc/m3
initial_soc1 = load([data_path, 'trendy_clm5/input_data/trendy_clm5_global_s3_matrix/global_trendy2020_clm5_', scenario_name, '_SOIL1C_vr_170001_201912.mat']);
initial_soc1 = initial_soc1.area_weighted_sum'/global_land_area;
initial_soc2 = load([data_path, 'trendy_clm5/input_data/trendy_clm5_global_s3_matrix/global_trendy2020_clm5_', scenario_name, '_SOIL2C_vr_170001_201912.mat']);
initial_soc2 = initial_soc2.area_weighted_sum'/global_land_area;
initial_soc3 = load([data_path, 'trendy_clm5/input_data/trendy_clm5_global_s3_matrix/global_trendy2020_clm5_', scenario_name, '_SOIL3C_vr_170001_201912.mat']);
initial_soc3 = initial_soc3.area_weighted_sum'/global_land_area;
initial_litter1 = load([data_path, 'trendy_clm5/input_data/trendy_clm5_global_s3_matrix/global_trendy2020_clm5_', scenario_name, '_LITR1C_vr_170001_201912.mat']);
initial_litter1 = initial_litter1.area_weighted_sum'/global_land_area;
initial_litter2 = load([data_path, 'trendy_clm5/input_data/trendy_clm5_global_s3_matrix/global_trendy2020_clm5_', scenario_name, '_LITR2C_vr_170001_201912.mat']);
initial_litter2 = initial_litter2.area_weighted_sum'/global_land_area;
initial_litter3 = load([data_path, 'trendy_clm5/input_data/trendy_clm5_global_s3_matrix/global_trendy2020_clm5_', scenario_name, '_LITR3C_vr_170001_201912.mat']);
initial_litter3 = initial_litter3.area_weighted_sum'/global_land_area;
initial_cwd = load([data_path, 'trendy_clm5/input_data/trendy_clm5_global_s3_matrix/global_trendy2020_clm5_', scenario_name, '_CWDC_vr_170001_201912.mat']);
initial_cwd = initial_cwd.area_weighted_sum'/global_land_area;

initial_cpool_total = [initial_cwd; initial_litter1; initial_litter2; initial_litter3; initial_soc1; initial_soc2; initial_soc3];

initial_cpool = nan(soil_cpool_num*soil_decom_num, length(da_start_loc));
% initial_cpool(:, 2:end) = initial_cpool_total(:, da_start_loc(2:end)-1);
initial_cpool(:, 1:end) = initial_cpool_total(:, da_start_loc(1:end)-1);


% carbon storage capacity, unit gc/m3
cap_cwd = load([data_path, 'trendy_clm5/input_data/trendy_clm5_global_s3_matrix/global_trendy2020_clm5_', scenario_name, '_CWDC_Cap_vr_170001_201912.mat']);
cap_cwd = cap_cwd.area_weighted_sum'/global_land_area;
cap_cwd = sum(cap_cwd.*repmat(dz(1:soil_decom_num), [1, length(trendy_time_series)]), 1); % unit gc/m2

cap_litter1 = load([data_path, 'trendy_clm5/input_data/trendy_clm5_global_s3_matrix/global_trendy2020_clm5_', scenario_name, '_LITR1C_Cap_vr_170001_201912.mat']);
cap_litter1 = cap_litter1.area_weighted_sum'/global_land_area;
cap_litter2 = load([data_path, 'trendy_clm5/input_data/trendy_clm5_global_s3_matrix/global_trendy2020_clm5_', scenario_name, '_LITR2C_Cap_vr_170001_201912.mat']);
cap_litter2 = cap_litter2.area_weighted_sum'/global_land_area;
cap_litter3 = load([data_path, 'trendy_clm5/input_data/trendy_clm5_global_s3_matrix/global_trendy2020_clm5_', scenario_name, '_LITR3C_Cap_vr_170001_201912.mat']);
cap_litter3 = cap_litter3.area_weighted_sum'/global_land_area;
cap_litter = cap_litter1 + cap_litter2 + cap_litter3;
cap_litter = sum(cap_litter.*repmat(dz(1:soil_decom_num), [1, length(trendy_time_series)]), 1); % unit gc/m2

cap_soc1 = load([data_path, 'trendy_clm5/input_data/trendy_clm5_global_s3_matrix/global_trendy2020_clm5_', scenario_name, '_SOIL1C_Cap_vr_170001_201912.mat']);
cap_soc1 = cap_soc1.area_weighted_sum'/global_land_area;
cap_soc2 = load([data_path, 'trendy_clm5/input_data/trendy_clm5_global_s3_matrix/global_trendy2020_clm5_', scenario_name, '_SOIL2C_Cap_vr_170001_201912.mat']);
cap_soc2 = cap_soc2.area_weighted_sum'/global_land_area;
cap_soc3 = load([data_path, 'trendy_clm5/input_data/trendy_clm5_global_s3_matrix/global_trendy2020_clm5_', scenario_name, '_SOIL3C_Cap_vr_170001_201912.mat']);
cap_soc3 = cap_soc3.area_weighted_sum'/global_land_area;
cap_soc = cap_soc1 + cap_soc2 + cap_soc3;
cap_soc = sum(cap_soc.*repmat(dz(1:soil_decom_num), [1, length(trendy_time_series)]), 1); % unit gc/m2

for iworker = 1%:worker_num
    post_para = load([data_path, 'trendy_clm5/output_data/da_reconstruct/trendy_clm5_global_reconstruct_test5_s3_matrix_0.05_worker_', num2str(iworker), '_parameters_keep2.mat']);
    post_para = post_para.parameters_keep2';
    
    valid_num = length(find(isnan(post_para(:, 1)) == 0));
    
    post_para_random = post_para(randsample(round(valid_num/2):valid_num, random_num), :);
    
    ensemble_simu_soc = nan(length(obs_tot_soc_transient), random_num);
    ensemble_simu_litter = nan(length(obs_tot_soc_transient), random_num);
    ensemble_simu_cwd = nan(length(obs_tot_soc_transient), random_num);
    ensemble_simu_hr = nan(length(obs_tot_soc_transient), random_num);
    
    ensemble_simu_cap_soc = nan(length(obs_tot_soc_transient), random_num);
    ensemble_simu_cap_litter = nan(length(obs_tot_soc_transient), random_num);
    ensemble_simu_cap_cwd = nan(length(obs_tot_soc_transient), random_num);
    
    for isimu = 1:random_num
        disp(['processing iworker ', num2str(iworker), ' simu ', num2str(isimu)])
        par_new = post_para_random(isimu, :);
        [post_da_mod_cpools_transient, post_da_mod_tot_soc_transient, post_da_mod_tot_litter_transient, post_da_mod_tot_cwd_transient, post_da_mod_tot_hr_transient, ...
            post_da_cap_soc_transient, post_da_cap_litter_transient, post_da_cap_cwd_transient] = ...
            fun_transient_trendy_simu(par_new, da_year_start, da_year_end, initial_cpool, frocing_steady_state, frocing_transient);
        
        ensemble_simu_soc(:, isimu) = post_da_mod_tot_soc_transient;
        ensemble_simu_litter(:, isimu) = post_da_mod_tot_litter_transient;
        ensemble_simu_cwd(:, isimu) = post_da_mod_tot_cwd_transient;
        ensemble_simu_hr(:, isimu) = post_da_mod_tot_hr_transient;
        
        ensemble_simu_cap_soc(:, isimu) = post_da_cap_soc_transient;
        ensemble_simu_cap_litter(:, isimu) = post_da_cap_litter_transient;
        ensemble_simu_cap_cwd(:, isimu) = post_da_cap_cwd_transient;
    end
    
    ensemble_simu.ensemble_simu_soc = ensemble_simu_soc;
    ensemble_simu.ensemble_simu_litter = ensemble_simu_litter;
    ensemble_simu.ensemble_simu_cwd = ensemble_simu_cwd;
    ensemble_simu.ensemble_simu_hr = ensemble_simu_hr;
    
    ensemble_simu.ensemble_simu_cap_soc = ensemble_simu_cap_soc;
    ensemble_simu.ensemble_simu_cap_litter = ensemble_simu_cap_litter;
    ensemble_simu.ensemble_simu_cap_cwd = ensemble_simu_cap_cwd;

    
    %     save([data_path, 'trendy_clm5/output_data/da_reconstruct/trendy_clm5_s3_point_', num2str(igrid), '_reconstruct_worker_', num2str(iworker), '_ensembe_simu.mat'], 'ensemble_simu');
end

%% obs v.s. modelled transient
% plot(trendy_time_series(selected_loc), cap_soc(selected_loc) - obs_tot_soc_transient, '-k', 'LineWidth', 2)
% hold on
% plot(trendy_time_series(selected_loc), mean(ensemble_simu_cap_soc, 2) - mean(ensemble_simu_soc, 2), '-r', 'LineWidth', 2)
% plot(trendy_time_series(selected_loc), quantile(ensemble_simu_cap_soc, 0.975, 2), '-r', 'LineWidth', 0.5)
% plot(trendy_time_series(selected_loc), quantile(ensemble_simu_cap_soc, 0.025, 2), '-r', 'LineWidth', 0.5)
% title('Total SOC stock (gc/m2)', 'FontSize', 24)
% xlabel('Year', 'FontSize', 20)
% ylabel('Simulation', 'FontSize', 20)
% ax=gca;
% ax.FontSize = 20;
% 
% 
% scatter(cap_soc(selected_loc) - obs_tot_soc_transient, mean(ensemble_simu_cap_soc, 2) - mean(ensemble_simu_soc, 2))

subplot(2, 2, 1)
plot(trendy_time_series(selected_loc), obs_tot_soc_transient, '-k', 'LineWidth', 2)
hold on
plot(trendy_time_series(selected_loc), mean(ensemble_simu_soc, 2), '-r', 'LineWidth', 2)
plot(trendy_time_series(selected_loc), quantile(ensemble_simu_soc, 0.975, 2), '-r', 'LineWidth', 0.5)
plot(trendy_time_series(selected_loc), quantile(ensemble_simu_soc, 0.025, 2), '-r', 'LineWidth', 0.5)
title('Total SOC stock (gc/m2)', 'FontSize', 24)
xlabel('Year', 'FontSize', 20)
ylabel('Simulation', 'FontSize', 20)
ax=gca;
ax.FontSize = 20;

subplot(2, 2, 2)
plot(trendy_time_series(selected_loc), obs_tot_litter_transient, '-k', 'LineWidth', 2)
hold on
plot(trendy_time_series(selected_loc), mean(ensemble_simu_litter, 2), '-r', 'LineWidth', 2)
plot(trendy_time_series(selected_loc), quantile(ensemble_simu_litter, 0.975, 2), '-r', 'LineWidth', 0.5)
plot(trendy_time_series(selected_loc), quantile(ensemble_simu_litter, 0.025, 2), '-r', 'LineWidth', 0.5)
title('Total litter stock (gc/m2)', 'FontSize', 24)
xlabel('Year', 'FontSize', 20)
ylabel('Simulation', 'FontSize', 20)
ax=gca;
ax.FontSize = 20;

subplot(2, 2, 3)
plot(trendy_time_series(selected_loc), obs_tot_cwd_transient, '-k', 'LineWidth', 2)
hold on
plot(trendy_time_series(selected_loc), mean(ensemble_simu_cwd, 2), '-r', 'LineWidth', 2)
plot(trendy_time_series(selected_loc), quantile(ensemble_simu_cwd, 0.975, 2), '-r', 'LineWidth', 0.5)
plot(trendy_time_series(selected_loc), quantile(ensemble_simu_cwd, 0.025, 2), '-r', 'LineWidth', 0.5)
title('Total CWD stock (gc/m2)', 'FontSize', 24)
xlabel('Year', 'FontSize', 20)
ylabel('Simulation', 'FontSize', 20)
ax=gca;
ax.FontSize = 20;

subplot(2, 2, 4)
plot(trendy_time_series(selected_loc), obs_tot_hr_transient, '-k', 'LineWidth', 2)
hold on
plot(trendy_time_series(selected_loc), mean(ensemble_simu_hr, 2), '-r', 'LineWidth', 2)
plot(trendy_time_series(selected_loc), quantile(ensemble_simu_hr, 0.975, 2), '-r', 'LineWidth', 0.5)
plot(trendy_time_series(selected_loc), quantile(ensemble_simu_hr, 0.025, 2), '-r', 'LineWidth', 0.5)
title('Heterotrophic Respiration(gc/m2/s)', 'FontSize', 24)
xlabel('Year', 'FontSize', 20)
ylabel('Simulation', 'FontSize', 20)
ax=gca;
ax.FontSize = 20;

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 25 15];

print([data_path, 'TRENDY_CLM5/figures/post_da_obs_vs_simu_transient_full.jpeg'], '-djpeg', '-r300')
close


%% obs v.s. modelled
r2_soc = 1 - sum((obs_tot_soc_transient - mean(ensemble_simu_soc, 2)').^2)/sum((obs_tot_soc_transient - mean(obs_tot_soc_transient)).^2);
r2_litter = 1 - sum((obs_tot_litter_transient - mean(ensemble_simu_litter, 2)').^2)/sum((obs_tot_litter_transient - mean(obs_tot_litter_transient)).^2);
r2_cwd = 1 - sum((obs_tot_cwd_transient - mean(ensemble_simu_cwd, 2)').^2)/sum((obs_tot_cwd_transient - mean(obs_tot_cwd_transient)).^2);
r2_hr = 1 - sum((obs_tot_hr_transient - mean(ensemble_simu_hr, 2)').^2)/sum((obs_tot_hr_transient - mean(obs_tot_hr_transient)).^2);
disp(['r2: soc ', num2str(r2_soc), ' litter ', num2str(r2_litter), ' cwd ', num2str(r2_cwd), ' hr ', num2str(r2_hr)])


subplot(2, 2, 1)
axes_limit = [min([obs_tot_soc_transient, mean(ensemble_simu_soc, 2)']), max([obs_tot_soc_transient, mean(ensemble_simu_soc, 2)'])];
plot(axes_limit, axes_limit, '-k', 'LineWidth', 2)
hold on
scatter(obs_tot_soc_transient, mean(ensemble_simu_soc, 2), 10, 'red', 'filled')
xlim(axes_limit)
ylim(axes_limit)
title('Total SOC stock (gc/m2)', 'FontSize', 25)
xlabel('Original simulation', 'FontSize', 20)
ylabel('Post DA simulation', 'FontSize', 20)
ax=gca;
ax.FontSize = 20;

subplot(2, 2, 2)
axes_limit = [min([obs_tot_litter_transient, mean(ensemble_simu_litter, 2)']), max([obs_tot_litter_transient, mean(ensemble_simu_litter, 2)'])];
plot(axes_limit, axes_limit, '-k', 'LineWidth', 2)
hold on
scatter(obs_tot_litter_transient, mean(ensemble_simu_litter, 2), 10, 'red', 'filled')
xlim(axes_limit)
ylim(axes_limit)
title('Total litter stock (gc/m2)', 'FontSize', 25)
xlabel('Original simulation', 'FontSize', 20)
ylabel('Post DA simulation', 'FontSize', 20)
ax=gca;
ax.FontSize = 20;

subplot(2, 2, 3)
axes_limit = [min([obs_tot_cwd_transient, mean(ensemble_simu_cwd, 2)']), max([obs_tot_cwd_transient, mean(ensemble_simu_cwd, 2)'])];
plot(axes_limit, axes_limit, '-k', 'LineWidth', 2)
hold on
scatter(obs_tot_cwd_transient, mean(ensemble_simu_cwd, 2), 10, 'red', 'filled')
xlim(axes_limit)
ylim(axes_limit)
title('Total CWD stock (gc/m2)', 'FontSize', 25)
xlabel('Original simulation', 'FontSize', 20)
ylabel('Post DA simulation', 'FontSize', 20)
ax=gca;
ax.FontSize = 20;

subplot(2, 2, 4)
axes_limit = [min([obs_tot_hr_transient, mean(ensemble_simu_hr, 2)']), max([obs_tot_hr_transient, mean(ensemble_simu_hr, 2)'])];
plot(axes_limit, axes_limit, '-k', 'LineWidth', 2)
hold on
scatter(obs_tot_hr_transient, mean(ensemble_simu_hr, 2), 10, 'red', 'filled')
xlim(axes_limit)
ylim(axes_limit)
title('Heterotrophic Respiration(gc/m2/s)', 'FontSize', 24)
xlabel('Original simulation', 'FontSize', 20)
ylabel('Post DA simulation', 'FontSize', 20)
ax=gca;
ax.FontSize = 20;

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 15 15];

print([data_path, 'TRENDY_CLM5/figures/post_da_obs_vs_simu.jpeg'], '-djpeg', '-r300')
close



