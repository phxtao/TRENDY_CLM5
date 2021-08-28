data_path = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/';
model_src_dir = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/TRENDY_CLM5/src_clm_cen';
cd(model_src_dir)

%% da obs vs mod
subplot(2, 2, 1)
plot(1:length(da_time_series), obs_tot_soc_transient, '-k', 'LineWidth', 2)
hold on
plot(1:length(da_time_series), mod_tot_soc_transient, '-r', 'LineWidth', 2)
xticks([12:12:length(da_time_series)])
xticklabels(string(da_time_series(12:12:length(da_time_series))))
title('Total SOC stock (gc/m2)', 'FontSize', 24)
xlabel('Year', 'FontSize', 20)
ylabel('Simulation', 'FontSize', 20)
ax=gca;
ax.FontSize = 20;

subplot(2, 2, 2)
plot(1:length(da_time_series), obs_tot_litter_transient, '-k', 'LineWidth', 2)
hold on
plot(1:length(da_time_series), mod_tot_litter_transient, '-r', 'LineWidth', 2)
xticks([12:12:length(da_time_series)])
xticklabels(string(da_time_series(12:12:length(da_time_series))))
title('Total litter stock (gc/m2)', 'FontSize', 24)
xlabel('Year', 'FontSize', 20)
ylabel('Simulation', 'FontSize', 20)
ax=gca;
ax.FontSize = 20;

subplot(2, 2, 3)
plot(1:length(da_time_series), obs_tot_cwd_transient, '-k', 'LineWidth', 2)
hold on
plot(1:length(da_time_series), mod_tot_cwd_transient, '-r', 'LineWidth', 2)
xticks([12:12:length(da_time_series)])
xticklabels(string(da_time_series(12:12:length(da_time_series))))
title('Total CWD stock (gc/m2)', 'FontSize', 24)
xlabel('Year', 'FontSize', 20)
ylabel('Simulation', 'FontSize', 20)
ax=gca;
ax.FontSize = 20;

subplot(2, 2, 4)
plot(1:length(da_time_series), obs_tot_hr_transient, '-k', 'LineWidth', 2)
hold on
plot(1:length(da_time_series), mod_tot_hr_transient, '-r', 'LineWidth', 2)
xticks([12:12:length(da_time_series)])
xticklabels(string(da_time_series(12:12:length(da_time_series))))
title('Heterotrophic Respiration(gc/m2/s)', 'FontSize', 24)
xlabel('Year', 'FontSize', 20)
ylabel('Simulation', 'FontSize', 20)
ax=gca;
ax.FontSize = 20;


fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 25 15];

print([data_path, 'TRENDY_CLM5/figures/post_da_obs_vs_simu_transient_da.jpeg'], '-djpeg', '-r300')
close

%% full forward simulation
% load simulation from TRENDY CLM5
input_var_list = {'ALTMAX', 'FPI', 'NPP', 'O_SCALAR', 'W_SCALAR', 'T_SCALAR'};

for ivar = 1:length(input_var_list)
    load([data_path, 'TRENDY_CLM5/input_data/grid_s3/grid_', num2str(grid_num), '_trendy2020_clm5_', scenario_name, '_', input_var_list{ivar}, '_170001_201912.mat']);
    
    eval(['post_da_frocing_transient.', input_var_list{ivar}, ' = var_data_grid;']);
end
post_da_frocing_transient.ALTMAX_LAST_YEAR = [post_da_frocing_transient.ALTMAX(1), post_da_frocing_transient.ALTMAX(1:end-1)];
post_da_frocing_transient.nbedrock = nbedrock;

% original simulated soc, unit gc/m2
post_da_obs_tot_soc_transient = load([data_path, 'TRENDY_CLM5/input_data/grid_s3/grid_', num2str(grid_num), '_trendy2020_clm5_', scenario_name, '_TOTSOMC_170001_201912.mat']);
post_da_obs_tot_soc_transient = post_da_obs_tot_soc_transient.var_data_grid;
% original simulated litter, unit gc/m2
post_da_obs_tot_litter_transient = load([data_path, 'TRENDY_CLM5/input_data/grid_s3/grid_', num2str(grid_num), '_trendy2020_clm5_', scenario_name, '_TOTLITC_170001_201912.mat']);
post_da_obs_tot_litter_transient = post_da_obs_tot_litter_transient.var_data_grid;
% original simulated cwdc, unit gc/m2
post_da_obs_tot_cwd_transient = load([data_path, 'TRENDY_CLM5/input_data/grid_s3/grid_', num2str(grid_num), '_trendy2020_clm5_', scenario_name, '_CWDC_170001_201912.mat']);
post_da_obs_tot_cwd_transient = post_da_obs_tot_cwd_transient.var_data_grid;
% original simulated heterotrophic respiration, unit gc/m2/s
post_da_obs_tot_hr_transient = load([data_path, 'TRENDY_CLM5/input_data/grid_s3/grid_', num2str(grid_num), '_trendy2020_clm5_', scenario_name, '_HR_170001_201912.mat']);
post_da_obs_tot_hr_transient = post_da_obs_tot_hr_transient.var_data_grid;
% initial cabron pool sizes at the beginning of transient simulation, unit gc/m3
post_da_initial_soc1 = load([data_path, 'TRENDY_CLM5/input_data/grid_s3/grid_', num2str(grid_num), '_trendy2020_clm5_', scenario_name, '_SOIL1C_vr_170001_201912.mat']);
post_da_initial_soc1 = post_da_initial_soc1.var_data_grid(:, 12);
post_da_initial_soc2 = load([data_path, 'TRENDY_CLM5/input_data/grid_s3/grid_', num2str(grid_num), '_trendy2020_clm5_', scenario_name, '_SOIL2C_vr_170001_201912.mat']);
post_da_initial_soc2 = post_da_initial_soc2.var_data_grid(:, 12);
post_da_initial_soc3 = load([data_path, 'TRENDY_CLM5/input_data/grid_s3/grid_', num2str(grid_num), '_trendy2020_clm5_', scenario_name, '_SOIL3C_vr_170001_201912.mat']);
post_da_initial_soc3 = post_da_initial_soc3.var_data_grid(:, 12);
post_da_initial_litter1 = load([data_path, 'TRENDY_CLM5/input_data/grid_s3/grid_', num2str(grid_num), '_trendy2020_clm5_', scenario_name, '_LITR1C_vr_170001_201912.mat']);
post_da_initial_litter1 = post_da_initial_litter1.var_data_grid(:, 12);
post_da_initial_litter2 = load([data_path, 'TRENDY_CLM5/input_data/grid_s3/grid_', num2str(grid_num), '_trendy2020_clm5_', scenario_name, '_LITR2C_vr_170001_201912.mat']);
post_da_initial_litter2 = post_da_initial_litter2.var_data_grid(:, 12);
post_da_initial_litter3 = load([data_path, 'TRENDY_CLM5/input_data/grid_s3/grid_', num2str(grid_num), '_trendy2020_clm5_', scenario_name, '_LITR3C_vr_170001_201912.mat']);
post_da_initial_litter3 = post_da_initial_litter3.var_data_grid(:, 12);
post_da_initial_cwd = load([data_path, 'TRENDY_CLM5/input_data/grid_s3/grid_', num2str(grid_num), '_trendy2020_clm5_', scenario_name, '_CWDC_vr_170001_201912.mat']);
post_da_initial_cwd = post_da_initial_cwd.var_data_grid(:, 12);

post_da_initial_cpool = [post_da_initial_cwd; ...
    post_da_initial_litter1; post_da_initial_litter2; post_da_initial_litter3;...
    post_da_initial_soc1; post_da_initial_soc2; post_da_initial_soc3];


[post_da_mod_cpools_transient, post_da_mod_tot_soc_transient, post_da_mod_tot_litter_transient, post_da_mod_tot_cwd_transient, post_da_mod_tot_hr_transient] = ...
    fun_transient_trendy_simu(par_new, 1700, 2019, post_da_initial_cpool, frocing_steady_state, post_da_frocing_transient);

%% obs v.s. modelled transient
subplot(2, 2, 1)
plot(trendy_time_series, post_da_obs_tot_soc_transient, '-k', 'LineWidth', 2)
hold on
plot(trendy_time_series, post_da_mod_tot_soc_transient, '-r', 'LineWidth', 2)
title('Total SOC stock (gc/m2)', 'FontSize', 24)
xlabel('Year', 'FontSize', 20)
ylabel('Simulation', 'FontSize', 20)
ax=gca;
ax.FontSize = 20;

subplot(2, 2, 2)
plot(trendy_time_series, post_da_obs_tot_litter_transient, '-k', 'LineWidth', 2)
hold on
plot(trendy_time_series, post_da_mod_tot_litter_transient, '-r', 'LineWidth', 2)
title('Total litter stock (gc/m2)', 'FontSize', 24)
xlabel('Year', 'FontSize', 20)
ylabel('Simulation', 'FontSize', 20)
ax=gca;
ax.FontSize = 20;

subplot(2, 2, 3)
plot(trendy_time_series, post_da_obs_tot_cwd_transient, '-k', 'LineWidth', 2)
hold on
plot(trendy_time_series, post_da_mod_tot_cwd_transient, '-r', 'LineWidth', 2)
title('Total CWD stock (gc/m2)', 'FontSize', 24)
xlabel('Year', 'FontSize', 20)
ylabel('Simulation', 'FontSize', 20)
ax=gca;
ax.FontSize = 20;

subplot(2, 2, 4)
plot(trendy_time_series, post_da_obs_tot_hr_transient, '-k', 'LineWidth', 2)
hold on
plot(trendy_time_series, post_da_mod_tot_hr_transient, '-r', 'LineWidth', 2)
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
subplot(2, 2, 1)
axes_limit = [min([post_da_obs_tot_soc_transient, post_da_mod_tot_soc_transient']), max([post_da_obs_tot_soc_transient, post_da_mod_tot_soc_transient'])];
plot(axes_limit, axes_limit, '-k', 'LineWidth', 2)
hold on
scatter(post_da_obs_tot_soc_transient, post_da_mod_tot_soc_transient, 10, 'red', 'filled')
xlim(axes_limit)
ylim(axes_limit)
title('Total SOC stock (gc/m2)', 'FontSize', 25)
xlabel('Original simulation', 'FontSize', 20)
ylabel('Post DA simulation', 'FontSize', 20)
ax=gca;
ax.FontSize = 20;

subplot(2, 2, 2)
axes_limit = [min([post_da_obs_tot_litter_transient, post_da_mod_tot_litter_transient']), max([post_da_obs_tot_litter_transient, post_da_mod_tot_litter_transient'])];
plot(axes_limit, axes_limit, '-k', 'LineWidth', 2)
hold on
scatter(post_da_obs_tot_litter_transient, post_da_mod_tot_litter_transient, 10, 'red', 'filled')
xlim(axes_limit)
ylim(axes_limit)
title('Total litter stock (gc/m2)', 'FontSize', 25)
xlabel('Original simulation', 'FontSize', 20)
ylabel('Post DA simulation', 'FontSize', 20)
ax=gca;
ax.FontSize = 20;

subplot(2, 2, 3)
axes_limit = [min([post_da_obs_tot_cwd_transient, post_da_mod_tot_cwd_transient']), max([post_da_obs_tot_cwd_transient, post_da_mod_tot_cwd_transient'])];
plot(axes_limit, axes_limit, '-k', 'LineWidth', 2)
hold on
scatter(post_da_obs_tot_cwd_transient, post_da_mod_tot_cwd_transient, 10, 'red', 'filled')
xlim(axes_limit)
ylim(axes_limit)
title('Total CWD stock (gc/m2)', 'FontSize', 25)
xlabel('Original simulation', 'FontSize', 20)
ylabel('Post DA simulation', 'FontSize', 20)
ax=gca;
ax.FontSize = 20;

subplot(2, 2, 4)
axes_limit = [min([post_da_obs_tot_hr_transient, post_da_mod_tot_hr_transient']), max([post_da_obs_tot_hr_transient, post_da_mod_tot_hr_transient'])];
plot(axes_limit, axes_limit, '-k', 'LineWidth', 2)
hold on
scatter(post_da_obs_tot_hr_transient, post_da_mod_tot_hr_transient, 10, 'red', 'filled')
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



