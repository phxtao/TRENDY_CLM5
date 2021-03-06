clear
clc;
%%
scenario_name = 's3';

grid_num = 10;
worker_num = 6;

random_num = 50;

data_path = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/';
model_src_dir = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/TRENDY_CLM5/src_clm_cen';
cd(model_src_dir)

%% forward simulations
for igrid = 1:grid_num
    valid_post_para = nan(50000, 22, 6);
    for iworker = 1:worker_num
        post_para = load([data_path, 'trendy_clm5/output_data/da_reconstruct/trendy_clm5_s3_point_', num2str(igrid), '_reconstruct_worker_', num2str(iworker), '_parameters_keep2.mat']);
        post_para = post_para.parameters_keep2';
        valid_num = length(find(isnan(post_para(:, 1)) == 0));
        
        valid_post_para(1:length(round(valid_num/2):valid_num), :, iworker) = post_para(round(valid_num/2):valid_num, :);
    end
    save([data_path, 'trendy_clm5/output_data/da_reconstruct/trendy_clm5_s3_point_', num2str(igrid), '_reconstruct_valid_post_para.mat']);
    
end 


%% forward simulations
for igrid = 1:grid_num
    nbedrock = load([data_path, 'trendy_clm5/input_data/nbedrock_grid.mat']);
    nbedrock = nbedrock.nbedrock(igrid);
    
    month_num = 12;
soil_cpool_num = 7;
soil_decom_num = 20;


% Simulation input
trendy_year_start = 1700;
trendy_year_end = 2019;
trendy_time_series = (trendy_year_start+1/month_num:1/month_num:(trendy_year_end+1))';

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

% load simulation from TRENDY CLM5
input_var_list = {'ALTMAX', 'FPI', 'NPP', 'O_SCALAR', 'W_SCALAR', 'T_SCALAR'};

for ivar = 1:length(input_var_list)
    load([data_path, 'trendy_clm5/input_data/grid_s3/grid_', num2str(igrid), '_trendy2020_clm5_', scenario_name, '_', input_var_list{ivar}, '_170001_201912.mat']);
    
    eval(['frocing_transient.', input_var_list{ivar}, ' = var_data_grid(:, selected_loc);']);
    eval(['frocing_steady_state.', input_var_list{ivar}, ' = var_data_grid(:, 1:20*month_num);']);
end
frocing_transient.ALTMAX_LAST_YEAR = [frocing_transient.ALTMAX(1), frocing_transient.ALTMAX(1:end-1)];
frocing_transient.nbedrock = nbedrock;

frocing_steady_state.ALTMAX_LAST_YEAR = [frocing_steady_state.ALTMAX(1), frocing_steady_state.ALTMAX(1:end-1)];
frocing_steady_state.nbedrock = nbedrock;

% original simulated soc, unit gc/m2
obs_tot_soc_transient = load([data_path, 'trendy_clm5/input_data/grid_s3/grid_', num2str(igrid), '_trendy2020_clm5_', scenario_name, '_TOTSOMC_170001_201912.mat']);
obs_tot_soc_transient = obs_tot_soc_transient.var_data_grid(selected_loc);
% original simulated litter, unit gc/m2
obs_tot_litter_transient = load([data_path, 'trendy_clm5/input_data/grid_s3/grid_', num2str(igrid), '_trendy2020_clm5_', scenario_name, '_TOTLITC_170001_201912.mat']);
obs_tot_litter_transient = obs_tot_litter_transient.var_data_grid(selected_loc);
% original simulated cwdc, unit gc/m2
obs_tot_cwd_transient = load([data_path, 'trendy_clm5/input_data/grid_s3/grid_', num2str(igrid), '_trendy2020_clm5_', scenario_name, '_CWDC_170001_201912.mat']);
obs_tot_cwd_transient = obs_tot_cwd_transient.var_data_grid(selected_loc);
% original simulated heterotrophic respiration, unit gc/m2/s
obs_tot_hr_transient = load([data_path, 'trendy_clm5/input_data/grid_s3/grid_', num2str(igrid), '_trendy2020_clm5_', scenario_name, '_HR_170001_201912.mat']);
obs_tot_hr_transient = obs_tot_hr_transient.var_data_grid(selected_loc);

% initial cabron pool sizes at the beginning of transient simulation, unit gc/m3
initial_soc1 = load([data_path, 'trendy_clm5/input_data/grid_s3/grid_', num2str(igrid), '_trendy2020_clm5_', scenario_name, '_SOIL1C_vr_170001_201912.mat']);
initial_soc1 = initial_soc1.var_data_grid;
initial_soc2 = load([data_path, 'trendy_clm5/input_data/grid_s3/grid_', num2str(igrid), '_trendy2020_clm5_', scenario_name, '_SOIL2C_vr_170001_201912.mat']);
initial_soc2 = initial_soc2.var_data_grid;
initial_soc3 = load([data_path, 'trendy_clm5/input_data/grid_s3/grid_', num2str(igrid), '_trendy2020_clm5_', scenario_name, '_SOIL3C_vr_170001_201912.mat']);
initial_soc3 = initial_soc3.var_data_grid;
initial_litter1 = load([data_path, 'trendy_clm5/input_data/grid_s3/grid_', num2str(igrid), '_trendy2020_clm5_', scenario_name, '_LITR1C_vr_170001_201912.mat']);
initial_litter1 = initial_litter1.var_data_grid;
initial_litter2 = load([data_path, 'trendy_clm5/input_data/grid_s3/grid_', num2str(igrid), '_trendy2020_clm5_', scenario_name, '_LITR2C_vr_170001_201912.mat']);
initial_litter2 = initial_litter2.var_data_grid;
initial_litter3 = load([data_path, 'trendy_clm5/input_data/grid_s3/grid_', num2str(igrid), '_trendy2020_clm5_', scenario_name, '_LITR3C_vr_170001_201912.mat']);
initial_litter3 = initial_litter3.var_data_grid;
initial_cwd = load([data_path, 'trendy_clm5/input_data/grid_s3/grid_', num2str(igrid), '_trendy2020_clm5_', scenario_name, '_CWDC_vr_170001_201912.mat']);
initial_cwd = initial_cwd.var_data_grid;

initial_cpool_total = [initial_cwd; initial_litter1; initial_litter2; initial_litter3; initial_soc1; initial_soc2; initial_soc3];

initial_cpool = nan(soil_cpool_num*soil_decom_num, length(da_start_loc));
% initial_cpool(:, 2:end) = initial_cpool_total(:, da_start_loc(2:end)-1);
initial_cpool(:, 1:end) = initial_cpool_total(:, da_start_loc(1:end)-1);



    for iworker = 1:worker_num
        post_para = load([data_path, 'trendy_clm5/output_data/da_reconstruct/trendy_clm5_s3_point_', num2str(igrid), '_reconstruct_worker_', num2str(iworker), '_parameters_keep2.mat']);
        post_para = post_para.parameters_keep2';
        
        valid_num = length(find(isnan(post_para(:, 1)) == 0));
        
        post_para_random = post_para(randsample(round(valid_num/2):valid_num, random_num), :);
        
        ensemble_simu_soc = nan(length(obs_tot_soc_transient), random_num);
        ensemble_simu_litter = nan(length(obs_tot_soc_transient), random_num);
        ensemble_simu_cwd = nan(length(obs_tot_soc_transient), random_num);
        ensemble_simu_hr = nan(length(obs_tot_soc_transient), random_num);
        
        for isimu = 1:random_num
            disp(['processing grid ', num2str(igrid), ' iworker ', num2str(iworker), ' simu ', num2str(isimu)])
            par_new = post_para_random(isimu, :);
            [post_da_mod_cpools_transient, post_da_mod_tot_soc_transient, post_da_mod_tot_litter_transient, post_da_mod_tot_cwd_transient, post_da_mod_tot_hr_transient] = ...
                fun_transient_trendy_simu(par_new, da_year_start, da_year_end, initial_cpool, frocing_steady_state, frocing_transient);
            ensemble_simu_soc(:, isimu) = post_da_mod_tot_soc_transient;
            ensemble_simu_litter(:, isimu) = post_da_mod_tot_litter_transient;
            ensemble_simu_cwd(:, isimu) = post_da_mod_tot_cwd_transient;
            ensemble_simu_hr(:, isimu) = post_da_mod_tot_hr_transient;
        end
        
        ensemble_simu.ensemble_simu_soc = ensemble_simu_soc;
        ensemble_simu.ensemble_simu_litter = ensemble_simu_litter;
        ensemble_simu.ensemble_simu_cwd = ensemble_simu_cwd;
        ensemble_simu.ensemble_simu_hr = ensemble_simu_hr;
        
        save([data_path, 'trendy_clm5/output_data/da_reconstruct/trendy_clm5_s3_point_', num2str(igrid), '_reconstruct_worker_', num2str(iworker), '_ensembe_simu.mat'], 'ensemble_simu');
    end
end

% %% obs v.s. modelled transient
% subplot(2, 2, 1)
% plot(trendy_time_series, post_da_obs_tot_soc_transient, '-k', 'LineWidth', 2)
% hold on
% plot(trendy_time_series, post_da_mod_tot_soc_transient, '-r', 'LineWidth', 2)
% title('Total SOC stock (gc/m2)', 'FontSize', 24)
% xlabel('Year', 'FontSize', 20)
% ylabel('Simulation', 'FontSize', 20)
% ax=gca;
% ax.FontSize = 20;
%
% subplot(2, 2, 2)
% plot(trendy_time_series, post_da_obs_tot_litter_transient, '-k', 'LineWidth', 2)
% hold on
% plot(trendy_time_series, post_da_mod_tot_litter_transient, '-r', 'LineWidth', 2)
% title('Total litter stock (gc/m2)', 'FontSize', 24)
% xlabel('Year', 'FontSize', 20)
% ylabel('Simulation', 'FontSize', 20)
% ax=gca;
% ax.FontSize = 20;
%
% subplot(2, 2, 3)
% plot(trendy_time_series, post_da_obs_tot_cwd_transient, '-k', 'LineWidth', 2)
% hold on
% plot(trendy_time_series, post_da_mod_tot_cwd_transient, '-r', 'LineWidth', 2)
% title('Total CWD stock (gc/m2)', 'FontSize', 24)
% xlabel('Year', 'FontSize', 20)
% ylabel('Simulation', 'FontSize', 20)
% ax=gca;
% ax.FontSize = 20;
%
% subplot(2, 2, 4)
% plot(trendy_time_series, post_da_obs_tot_hr_transient, '-k', 'LineWidth', 2)
% hold on
% plot(trendy_time_series, post_da_mod_tot_hr_transient, '-r', 'LineWidth', 2)
% title('Heterotrophic Respiration(gc/m2/s)', 'FontSize', 24)
% xlabel('Year', 'FontSize', 20)
% ylabel('Simulation', 'FontSize', 20)
% ax=gca;
% ax.FontSize = 20;
%
% fig = gcf;
% fig.PaperUnits = 'inches';
% fig.PaperPosition = [0 0 25 15];
%
% print([data_path, 'TRENDY_CLM5/figures/post_da_obs_vs_simu_transient_full.jpeg'], '-djpeg', '-r300')
% close
%
%
% %% obs v.s. modelled
% subplot(2, 2, 1)
% axes_limit = [min([post_da_obs_tot_soc_transient, post_da_mod_tot_soc_transient']), max([post_da_obs_tot_soc_transient, post_da_mod_tot_soc_transient'])];
% plot(axes_limit, axes_limit, '-k', 'LineWidth', 2)
% hold on
% scatter(post_da_obs_tot_soc_transient, post_da_mod_tot_soc_transient, 10, 'red', 'filled')
% xlim(axes_limit)
% ylim(axes_limit)
% title('Total SOC stock (gc/m2)', 'FontSize', 25)
% xlabel('Original simulation', 'FontSize', 20)
% ylabel('Post DA simulation', 'FontSize', 20)
% ax=gca;
% ax.FontSize = 20;
%
% subplot(2, 2, 2)
% axes_limit = [min([post_da_obs_tot_litter_transient, post_da_mod_tot_litter_transient']), max([post_da_obs_tot_litter_transient, post_da_mod_tot_litter_transient'])];
% plot(axes_limit, axes_limit, '-k', 'LineWidth', 2)
% hold on
% scatter(post_da_obs_tot_litter_transient, post_da_mod_tot_litter_transient, 10, 'red', 'filled')
% xlim(axes_limit)
% ylim(axes_limit)
% title('Total litter stock (gc/m2)', 'FontSize', 25)
% xlabel('Original simulation', 'FontSize', 20)
% ylabel('Post DA simulation', 'FontSize', 20)
% ax=gca;
% ax.FontSize = 20;
%
% subplot(2, 2, 3)
% axes_limit = [min([post_da_obs_tot_cwd_transient, post_da_mod_tot_cwd_transient']), max([post_da_obs_tot_cwd_transient, post_da_mod_tot_cwd_transient'])];
% plot(axes_limit, axes_limit, '-k', 'LineWidth', 2)
% hold on
% scatter(post_da_obs_tot_cwd_transient, post_da_mod_tot_cwd_transient, 10, 'red', 'filled')
% xlim(axes_limit)
% ylim(axes_limit)
% title('Total CWD stock (gc/m2)', 'FontSize', 25)
% xlabel('Original simulation', 'FontSize', 20)
% ylabel('Post DA simulation', 'FontSize', 20)
% ax=gca;
% ax.FontSize = 20;
%
% subplot(2, 2, 4)
% axes_limit = [min([post_da_obs_tot_hr_transient, post_da_mod_tot_hr_transient']), max([post_da_obs_tot_hr_transient, post_da_mod_tot_hr_transient'])];
% plot(axes_limit, axes_limit, '-k', 'LineWidth', 2)
% hold on
% scatter(post_da_obs_tot_hr_transient, post_da_mod_tot_hr_transient, 10, 'red', 'filled')
% xlim(axes_limit)
% ylim(axes_limit)
% title('Heterotrophic Respiration(gc/m2/s)', 'FontSize', 24)
% xlabel('Original simulation', 'FontSize', 20)
% ylabel('Post DA simulation', 'FontSize', 20)
% ax=gca;
% ax.FontSize = 20;
%
% fig = gcf;
% fig.PaperUnits = 'inches';
% fig.PaperPosition = [0 0 15 15];
%
% print([data_path, 'TRENDY_CLM5/figures/post_da_obs_vs_simu.jpeg'], '-djpeg', '-r300')
% close
%
%
%
