clear;
clc;
grid_num = 1;
worker_num = 6;
%%

% worker_num = 6;
% delete(gcp('nocreate'));
% parpool(worker_num);

%
warning('off');
format long e;

disp([datestr(now,'HH:MM:SS'), ' programme started'])

scenario_name = 's3';

model_src_dir = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/TRENDY_CLM5/src_clm_cen';
% model_src_dir = '/GFPS8p/cess11/taof/trendy_clm5/src_clm_cen';
cd(model_src_dir)

% mac
data_path = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/';

% server
% data_path = '/GFPS8p/cess11/taof/datahub/';

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

nbedrock = load([data_path, 'trendy_clm5/input_data/nbedrock_grid.mat']);
nbedrock = nbedrock.nbedrock(grid_num);
%% Simulation input
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
    load([data_path, 'trendy_clm5/input_data/grid_s3/grid_', num2str(grid_num), '_trendy2020_clm5_', scenario_name, '_', input_var_list{ivar}, '_170001_201912.mat']);
    
    eval(['frocing_transient.', input_var_list{ivar}, ' = var_data_grid(:, selected_loc);']);
    eval(['frocing_steady_state.', input_var_list{ivar}, ' = var_data_grid(:, 1:20*month_num);']);
end
frocing_transient.ALTMAX_LAST_YEAR = [frocing_transient.ALTMAX(1), frocing_transient.ALTMAX(1:end-1)];
frocing_transient.nbedrock = nbedrock;

frocing_steady_state.ALTMAX_LAST_YEAR = [frocing_steady_state.ALTMAX(1), frocing_steady_state.ALTMAX(1:end-1)];
frocing_steady_state.nbedrock = nbedrock;

% original simulated soc, unit gc/m2
obs_tot_soc_transient = load([data_path, 'trendy_clm5/input_data/grid_s3/grid_', num2str(grid_num), '_trendy2020_clm5_', scenario_name, '_TOTSOMC_170001_201912.mat']);
obs_tot_soc_transient = obs_tot_soc_transient.var_data_grid(selected_loc);
% original simulated litter, unit gc/m2
obs_tot_litter_transient = load([data_path, 'trendy_clm5/input_data/grid_s3/grid_', num2str(grid_num), '_trendy2020_clm5_', scenario_name, '_TOTLITC_170001_201912.mat']);
obs_tot_litter_transient = obs_tot_litter_transient.var_data_grid(selected_loc);
% original simulated cwdc, unit gc/m2
obs_tot_cwd_transient = load([data_path, 'trendy_clm5/input_data/grid_s3/grid_', num2str(grid_num), '_trendy2020_clm5_', scenario_name, '_CWDC_170001_201912.mat']);
obs_tot_cwd_transient = obs_tot_cwd_transient.var_data_grid(selected_loc);
% original simulated heterotrophic respiration, unit gc/m2/s
obs_tot_hr_transient = load([data_path, 'trendy_clm5/input_data/grid_s3/grid_', num2str(grid_num), '_trendy2020_clm5_', scenario_name, '_HR_170001_201912.mat']);
obs_tot_hr_transient = obs_tot_hr_transient.var_data_grid(selected_loc);

% initial cabron pool sizes at the beginning of transient simulation, unit gc/m3
initial_soc1 = load([data_path, 'trendy_clm5/input_data/grid_s3/grid_', num2str(grid_num), '_trendy2020_clm5_', scenario_name, '_SOIL1C_vr_170001_201912.mat']);
initial_soc1 = initial_soc1.var_data_grid;
initial_soc2 = load([data_path, 'trendy_clm5/input_data/grid_s3/grid_', num2str(grid_num), '_trendy2020_clm5_', scenario_name, '_SOIL2C_vr_170001_201912.mat']);
initial_soc2 = initial_soc2.var_data_grid;
initial_soc3 = load([data_path, 'trendy_clm5/input_data/grid_s3/grid_', num2str(grid_num), '_trendy2020_clm5_', scenario_name, '_SOIL3C_vr_170001_201912.mat']);
initial_soc3 = initial_soc3.var_data_grid;
initial_litter1 = load([data_path, 'trendy_clm5/input_data/grid_s3/grid_', num2str(grid_num), '_trendy2020_clm5_', scenario_name, '_LITR1C_vr_170001_201912.mat']);
initial_litter1 = initial_litter1.var_data_grid;
initial_litter2 = load([data_path, 'trendy_clm5/input_data/grid_s3/grid_', num2str(grid_num), '_trendy2020_clm5_', scenario_name, '_LITR2C_vr_170001_201912.mat']);
initial_litter2 = initial_litter2.var_data_grid;
initial_litter3 = load([data_path, 'trendy_clm5/input_data/grid_s3/grid_', num2str(grid_num), '_trendy2020_clm5_', scenario_name, '_LITR3C_vr_170001_201912.mat']);
initial_litter3 = initial_litter3.var_data_grid;
initial_cwd = load([data_path, 'trendy_clm5/input_data/grid_s3/grid_', num2str(grid_num), '_trendy2020_clm5_', scenario_name, '_CWDC_vr_170001_201912.mat']);
initial_cwd = initial_cwd.var_data_grid;

initial_cpool_total = [initial_cwd; initial_litter1; initial_litter2; initial_litter3; initial_soc1; initial_soc2; initial_soc3];

initial_cpool = nan(soil_cpool_num*soil_decom_num, length(da_start_loc));
% initial_cpool(:, 2:end) = initial_cpool_total(:, da_start_loc(2:end)-1);
initial_cpool(:, 1:end) = initial_cpool_total(:, da_start_loc(1:end)-1);
%% Data Assimilation - initiation
% Parameter names and their initial values in MCMC
para_name = {'diffus'; 'cryo';...
    'efolding';...
    'taucwd'; 'taul1'; 'taul2';...
    'tau4s1'; 'tau4s2'; 'tau4s3';...
    'fl1s1'; 'fl2s1'; 'fl3s2'; 'fs1s2'; 'fs1s3'; 'fs2s1'; 'fs2s3'; 'fs3s1'; 'fcwdl2';...
    'beta'; ...
    'p4l1'; 'p4l2'; 'p4l3' ...
    };

% number of paras
npara = length(para_name);
para_min = zeros(npara, 1);
para_max = ones(npara, 1);
% set initial values of parameters
% para0 = []';
% clear warning info
valid_para = 0;
while valid_para == 0
    para0 = rand(npara,1);
    
    if sum(para0(end-3:end)) < 1
        lastwarn('');
        [mod_cpools_transient, mod_tot_soc_transient, mod_tot_litter_transient, mod_tot_cwd_transient, mod_tot_hr_transient] = ...
            fun_transient_trendy_simu(para0, da_year_start, da_year_end, initial_cpool, frocing_steady_state, frocing_transient);
        [~, msgid] = lastwarn;
        
        if isnan(mod_tot_soc_transient) == 0
            %calculating prior cost function value
            cost_old1 = cost_fun_trendy(period_num, obs_tot_soc_transient, obs_tot_litter_transient, obs_tot_cwd_transient, obs_tot_hr_transient, mod_tot_soc_transient, mod_tot_litter_transient, mod_tot_cwd_transient, mod_tot_hr_transient);
            if cost_old1 < Inf && isnan(cost_old1) == 0
                valid_para = 1;
            end
        end
    end
end

%% MCMC
disp([datestr(now,'HH:MM:SS'), ' MCMC started'])
%     try
% Predefine the size of some coefficients
nsimu1 = 10000;
nsimu2 = 50000;

for iworker = 1:worker_num
    rng(iworker);
    is_second_time = 0;
    parameters_rec1 = nan(npara, nsimu1);
    parameters_rec2 = nan(npara, nsimu2);
    parameters_keep1 = nan(npara, nsimu1);
    parameters_keep2 = nan(npara, nsimu2);
    
    cost_record1 = nan(nsimu1, 1);
    cost_record2 = nan(nsimu2, 1);
    cost_keep1 = nan(nsimu1, 1);
    cost_keep2 = nan(nsimu2, 1);
    
    soc_mod_trace = nan(nsimu2, length(da_time_series));
    litter_mod_trace = nan(nsimu2, length(da_time_series));
    cwd_mod_trace = nan(nsimu2, length(da_time_series));
    hr_mod_trace = nan(nsimu2, length(da_time_series));

    litter_mod_trace = nan(nsimu2, length(da_time_series));
    cwd_mod_trace = nan(nsimu2, length(da_time_series));
    hr_mod_trace = nan(nsimu2, length(da_time_series));
    
    %----------------------------------------
    % Part1: MCMC run to calculate parameter covariances
    %----------------------------------------
    counter2 = 0;
    
    soc_mod_trace(:, :) = nan;
    litter_mod_trace(:, :) = nan;
    cwd_mod_trace(:, :) = nan;
    hr_mod_trace(:, :) = nan;
    
    cost_old = cost_old1;
    cost_control = 1;
    allow = 20*ones(npara, 1);
    % allow(strcmp(para_name, 'tau4p')) = 90;
    % give original parameter values to par_old
    par_old = para0;
    % interval width
    diff = para_max - para_min;
    upgrade1 = 0;
    simu = 1;
    counter1 = 0;
    ifbreak1 = 0;
    while simu <= nsimu1
        while (true)
            % generate new parameter values in the interval of
            % [-0.5, 0.5]*diff/allow + par_old
            par_new = par_old+(rand(npara,1)-0.5).*diff./allow;
            % all the new valued parameters should be in the preassumed
            % intervals
            if min(par_new) > 0 && max(par_new) < 1 ...
                    && sum(par_new(end-3:end) < 1)
                break;
            end
        end
        % soc_model simulation
        
        % clear warning info
        lastwarn('');
        [mod_cpools_transient, mod_tot_soc_transient, mod_tot_litter_transient, mod_tot_cwd_transient, mod_tot_hr_transient] = ...
            fun_transient_trendy_simu(par_new, da_year_start, da_year_end, initial_cpool, frocing_steady_state, frocing_transient);
        
        
        % original soc_modelled data
        [~, msgid] = lastwarn;
        if strcmp(msgid, 'MATLAB:singularMatrix') || strcmp(msgid, 'MATLAB:nearlySingularMatrix')
            cost_new = cost_old + 100;
        elseif isnan((sum(mod_tot_soc_transient + mod_tot_litter_transient + mod_tot_cwd_transient + mod_tot_hr_transient))) == 1
            cost_new = cost_old + 100;
        else
            % new cost function
            cost_new = cost_fun_trendy(period_num, obs_tot_soc_transient, obs_tot_litter_transient, obs_tot_cwd_transient, obs_tot_hr_transient, mod_tot_soc_transient, mod_tot_litter_transient, mod_tot_cwd_transient, mod_tot_hr_transient);
            
        end
        
        delta_cost = cost_new-cost_old;
        % to decide whether to accept new parameter values
        if min(1, exp(-delta_cost*log(simu + 2))) > rand
            % update the number of simulation
            upgrade1 = upgrade1 + 1;
            
            disp([datestr(now,'HH:MM:SS'), ' upgrade1 iworker ', num2str(iworker), ': ', num2str(upgrade1), ' out of ', num2str(simu), ' cost: ', num2str(cost_new)]);
            
            % record accepted parameter values
            parameters_keep1(:, upgrade1) = par_new;
            % record accepted values of cost function
            cost_keep1(upgrade1, 1)=cost_new;
            % update the value of parameter values
            par_old = par_new;
            % update old cost function value
            cost_old = cost_new;
        end
        % give parameter values to parameters_rec for covarance matrix,
        % rec: record
        parameters_rec1(:, simu)=par_old;
        cost_record1(simu, 1)=cost_old;
        simu = simu + 1;
        
        if simu == nsimu1 && is_second_time == 0
            if upgrade1 < 1
                counter1 = counter1 + 1;
                simu = 1;
                upgrade1 = 0;
                parameters_keep1(:, :) = nan;
                parameters_rec1(:, :) = nan;
                cost_record1(:, 1) = nan;
                cost_keep1(:, 1) = nan;
                allow = allow*1.1;
            end
        end
        if counter1 == 3
            ifbreak1 = 1;
            is_second_time = 1;
            break
        end
    end
    if ifbreak1 == 1
        continue
    end
    
    %----------------------------------------
    % Part 2: MCMC run with updating covariances
    %----------------------------------------
    sd_controling_factor = 0.6; % 2.4; % default to be 2.4^2
    sd = sd_controling_factor/npara;
    epsilon = 0;
    covars=sd*cov(parameters_rec1(:, 1:nsimu1)') + sd*epsilon*diag(ones(npara, 1));
    % covars=sd*cov(reshape(parameters_keep1(iparallel, :, 1:upgrade1(iparallel, 1)), [npara, upgrade1(iparallel, 1)])') + sd*epsilon*diag(ones(npara, 1));
    % update the value of upgrade and simu for new task
    upgrade2 = 0;
    % the very first cost functino value in simulation
    cost_old = min(cost_old1, cost_old*100);
    simu = 1;
    ifbreak2 = 0;
    ifbreak3 = 0;
    while simu <= nsimu2
        monitor_time_2 = tic;
        while (true)
            % mointor the consumed time at the initial stage (if too long i.e. 10min, then back to the first proposal step)
            if upgrade2 == 0
                consume_time_2 = toc(monitor_time_2);
                if consume_time_2 > 1*60
                    ifbreak2 = 1;
                end
            end
            if ifbreak2 == 1
                break
            end
            try
                par_new = mvnrnd(par_old, covars)';
            catch
                ifbreak3 = 1;
                break
            end
            if min(par_new) > 0 && max(par_new) < 1 ...
                    && sum(par_new(end-3:end) < 1)
                break;
            end
        end
        if ifbreak2 == 1 || ifbreak3 == 1
            break
        end
        
        % clear warning info
        lastwarn('');
        [mod_cpools_transient, mod_tot_soc_transient, mod_tot_litter_transient, mod_tot_cwd_transient, mod_tot_hr_transient] = ...
            fun_transient_trendy_simu(par_new, da_year_start, da_year_end, initial_cpool, frocing_steady_state, frocing_transient);
        
        
        % original soc_modelled data
        [~, msgid] = lastwarn;
        if strcmp(msgid, 'MATLAB:singularMatrix') || strcmp(msgid, 'MATLAB:nearlySingularMatrix')
            cost_new = cost_old + 100;
        elseif isnan((sum(mod_tot_soc_transient + mod_tot_litter_transient + mod_tot_cwd_transient + mod_tot_hr_transient))) == 1
            cost_new = cost_old + 100;
        else
            % new cost function
            cost_new = cost_fun_trendy(period_num, obs_tot_soc_transient, obs_tot_litter_transient, obs_tot_cwd_transient, obs_tot_hr_transient, mod_tot_soc_transient, mod_tot_litter_transient, mod_tot_cwd_transient, mod_tot_hr_transient);
            
        end
        
        delta_cost = cost_new-cost_old;
        % to determie whether to accept the newly soc_modelled data
        if min(1, exp(-delta_cost*log(simu + 2))) > rand
            upgrade2 = upgrade2 + 1;
            
            disp([datestr(now,'HH:MM:SS'), ' upgrade2 iworker ', num2str(iworker), ': ', num2str(upgrade2), ' out of ', num2str(simu), ' cost: ', num2str(cost_new)]);
            % record the accepted parameter value
            parameters_keep2(:,upgrade2) = par_new;
            cost_keep2(upgrade2, 1) = cost_new;
            par_old = par_new;
            cost_old = cost_new;
            coef = upgrade2/simu;
            
            soc_mod_trace(upgrade2, :) = mod_tot_soc_transient;
            litter_mod_trace(upgrade2, :) = mod_tot_litter_transient;
            cwd_mod_trace(upgrade2, :) = mod_tot_cwd_transient;
            hr_mod_trace(upgrade2, :) = mod_tot_hr_transient;
            
        end
        parameters_rec2(:, simu) = par_old;
        cost_record2(simu, 1) = cost_old;
        if simu > 4000
            % covars=sd*(cov(reshape(parameters_keep2(iparallel, :, 1:upgrade2(iparallel, 1)), [npara, upgrade2(iparallel, 1)])')) + sd*epsilon*diag(ones(npara, 1));
            covars=sd*(cov(parameters_rec2(:, 1:simu)')) + sd*epsilon*diag(ones(npara, 1));
        end
        simu = simu + 1;
        if mod(simu, 1000) == 0 || simu == nsimu2
            fun_save_da(data_path, grid_num, iworker, soc_mod_trace, litter_mod_trace, cwd_mod_trace, hr_mod_trace, parameters_keep2)
        end
        % to test if the acceptance rate is in a resonable level
        if simu == nsimu2
            if coef < 0 || coef > 1 % coef < 0.1 || coef > 0.5
                ifbreak2 = 1;
                counter2 = counter2 + 1;
                
                if coef < 0.1
                    sd_controling_factor = sd_controling_factor*(1 - 0.3);
                end
                
                if coef > 0.5
                    sd_controling_factor = sd_controling_factor*2;
                end
                
                break
            end
        end
    end
    if counter2 > 5
        error('Error')
    end
end % for iworker = 1:6
disp([datestr(now,'HH:MM:SS'), ' MCMC has been finished']);



