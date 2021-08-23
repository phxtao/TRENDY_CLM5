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
year_end = 1901; %2010;

soc_time_series = (year_start+1/12:1/12:(year_end+1))';

% load simulation from CESM2 (steady state)
load([data_path, 'radiocarbon/demo_harvard_forest/input_data/forcing_steady_state_harvard_forest.mat']);

% load simulation from CESM2 (transient)
load([data_path, 'radiocarbon/demo_harvard_forest/input_data/forcing_transient_harvard_forest.mat']);

obs_tot_soc = reshape(frocing_transient.soc_stock(:, 1:length(year_start:year_end)), [length(soc_time_series), 1]);
%% Data Assimilation - initiation
% Parameter names and their initial values in MCMC
para_name = {'diffus'; 'cryo';...
    'q10';...
    'efolding';...
    'taucwd'; 'taul1'; 'taul2';...
    'tau4s1'; 'tau4s2'; 'tau4s3';...
    'fl1s1'; 'fl2s1'; 'fl3s2'; 'fs1s2'; 'fs1s3'; 'fs2s1'; 'fs2s3'; 'fs3s1'; 'fcwdl2';...
    'w-scaling'; 'beta'};

% number of paras
npara = length(para_name);
para_min = zeros(npara, 1);
para_max = ones(npara, 1);
% set initial values of parameters
para0 = 0.5*ones(npara, 1);
% para0 = []';
% clear warning info
lastwarn('');

% original soc_modelled data
[~, soc_stock_mod, soc_mod] = ...
    fun_forward_simu(para0, year_start, year_end, frocing_steady_state, frocing_transient);

% layer_weighting at different layers of soil profile

%% MCMC
disp([datestr(now,'HH:MM:SS'), ' MCMC started'])
%     try
is_second_time = 0;
% Predefine the size of some coefficients
nsimu1 = 500;
nsimu2 = 1000;
nparallel = 3;

parameters_rec1 = nan(nparallel, npara, nsimu1);
parameters_rec2 = nan(nparallel, npara, nsimu2);
parameters_keep1 = nan(nparallel, npara, nsimu1);
parameters_keep2 = nan(nparallel, npara, nsimu2);

cost_record1 = nan(nsimu1, nparallel);
cost_record2 = nan(nsimu2, nparallel);
cost_keep1 = nan(nsimu1, nparallel);
cost_keep2 = nan(nsimu2, nparallel);

upgrade1 = nan(nparallel, 1);
upgrade2 = nan(nparallel, 1);

soc_mod_trace = nan(nparallel, nsimu2, length(obs_tot_soc));

%calculating prior cost function value
cost_old1 = cost_fun(obs_tot_soc, soc_stock_mod(:, 4));

%----------------------------------------
% Part1: MCMC run to calculate parameter covariances
%----------------------------------------
iparallel = 1;
counter2 = 0;
while iparallel <= nparallel
    parameters_rec1(iparallel, :, :) = nan;
    parameters_rec2(iparallel, :, :) = nan;
    parameters_keep1(iparallel, :, :) = nan;
    parameters_keep2(iparallel, :, :) = nan;
    
    cost_record1(:, iparallel) = nan;
    cost_record2(:, iparallel) = nan;
    cost_keep1(:, iparallel) = nan;
    cost_keep2(:, iparallel) = nan;
    
    soc_mod_trace(iparallel, :, :) = nan;
    
    cost_old = cost_old1;
    cost_control = 1;
    allow = 10*ones(npara, 1);
    % allow(strcmp(para_name, 'tau4p')) = 90;
    % give original parameter values to par_old
    par_old = para0;
    % interval width
    diff = para_max - para_min;
    upgrade1(iparallel, 1) = 0;
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
            if par_new(1)>para_min(1) && par_new(1)<para_max(1)...
                    && par_new(2)>para_min(2) && par_new(2)<para_max(2)...
                    && par_new(3)>para_min(3) && par_new(3)<para_max(3)...
                    && par_new(4)>para_min(4) && par_new(4)<para_max(4)...
                    && par_new(5)>para_min(5) && par_new(5)<para_max(5)...
                    && par_new(6)>para_min(6) && par_new(6)<para_max(6)...
                    && par_new(7)>para_min(7) && par_new(7)<para_max(7)...
                    && par_new(8)>para_min(8) && par_new(8)<para_max(8)...
                    && par_new(9)>para_min(9) && par_new(9)<para_max(9)...
                    && par_new(10)>para_min(10) && par_new(10)<para_max(10)...
                    && par_new(11)>para_min(11) && par_new(11)<para_max(11)...
                    && par_new(12)>para_min(12) && par_new(12)<para_max(12)...
                    && par_new(13)>para_min(13) && par_new(13)<para_max(13)...
                    && par_new(14)>para_min(14) && par_new(14)<para_max(14)...
                    && par_new(15)>para_min(15) && par_new(15)<para_max(15)...
                    && par_new(16)>para_min(16) && par_new(16)<para_max(16)...
                    && par_new(17)>para_min(17) && par_new(17)<para_max(17)...
                    && par_new(18)>para_min(18) && par_new(18)<para_max(18)...
                    && par_new(19)>para_min(19) && par_new(19)<para_max(19)...
                    && par_new(20)>para_min(20) && par_new(20)<para_max(20)...
                    && par_new(21)>para_min(21) && par_new(21)<para_max(21)
                break;
            end
        end
        
        
        % soc_model simulation
        
        % clear warning info
        lastwarn('');
        [~, soc_stock_mod, soc_mod] = ...
            fun_forward_simu(par_new, year_start, year_end, frocing_steady_state, frocing_transient);
        
        % original soc_modelled data
        [~, msgid] = lastwarn;
        if strcmp(msgid, 'MATLAB:singularMatrix') || strcmp(msgid, 'MATLAB:nearlySingularMatrix')
            cost_new = cost_old + 100;
        elseif isnan((sum(sum(soc_stock_mod(:, 4))))) == 1
            cost_new = cost_old + 100;
        else            
            % new cost function
            cost_new = cost_fun(obs_tot_soc, soc_stock_mod(:, 4));
            
        end
        delta_cost = cost_new-cost_old;
        % to decide whether to accept new parameter values
        if min(1, exp(-delta_cost*log(simu + 2))) > rand
            % update the number of simulation
            upgrade1(iparallel, 1) = upgrade1(iparallel, 1) + 1;
            
            disp([datestr(now,'HH:MM:SS'), ' upgrade1 parallel ', num2str(iparallel), ': ', num2str(upgrade1(iparallel, 1)), ' out of ', num2str(simu)]);
            
            % record accepted parameter values
            parameters_keep1(iparallel, :, upgrade1(iparallel, 1)) = par_new;
            % record accepted values of cost function
            cost_keep1(upgrade1(iparallel, 1), iparallel)=cost_new;
            % update the value of parameter values
            par_old = par_new;
            % update old cost function value
            cost_old = cost_new;
        end
        % give parameter values to parameters_rec for covarance matrix,
        % rec: record
        parameters_rec1(iparallel, :, simu)=par_old;
        cost_record1(simu, iparallel)=cost_old;
        simu = simu + 1;
        
        if simu == nsimu1 && is_second_time == 0
            if upgrade1(iparallel, 1) < 10
                counter1 = counter1 + 1;
                simu = 1;
                upgrade1(iparallel, 1) = 0;
                parameters_keep1(iparallel, :, :) = nan;
                parameters_rec1(iparallel, :, :) = nan;
                cost_record1(:, iparallel) = nan;
                cost_keep1(:, iparallel) = nan;
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
    sd_controling_factor = 2.4; % default to be 2.4^2
    sd = sd_controling_factor/npara;
    epsilon = 0;
    covars=sd*cov(reshape(parameters_rec1(iparallel, :, 1:nsimu1), [npara, length(1:nsimu1)])') + sd*epsilon*diag(ones(npara, 1));
    % covars=sd*cov(reshape(parameters_keep1(iparallel, :, 1:upgrade1(iparallel, 1)), [npara, upgrade1(iparallel, 1)])') + sd*epsilon*diag(ones(npara, 1));
    % update the value of upgrade and simu for new task
    upgrade2(iparallel, 1) = 0;
    % the very first cost functino value in simulation
    cost_old = cost_old1;
    simu = 1;
    ifbreak2 = 0;
    ifbreak3 = 0;
    while simu <= nsimu2
        monitor_time_2 = tic;
        while (true)
            % mointor the consumed time at the initial stage (if too long i.e. 10min, then back to the first proposal step)
            if upgrade2(iparallel, 1) == 0
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
            if par_new(1)>para_min(1) && par_new(1)<para_max(1)...
                    && par_new(2)>para_min(2) && par_new(2)<para_max(2)...
                    && par_new(3)>para_min(3) && par_new(3)<para_max(3)...
                    && par_new(4)>para_min(4) && par_new(4)<para_max(4)...
                    && par_new(5)>para_min(5) && par_new(5)<para_max(5)...
                    && par_new(6)>para_min(6) && par_new(6)<para_max(6)...
                    && par_new(7)>para_min(7) && par_new(7)<para_max(7)...
                    && par_new(8)>para_min(8) && par_new(8)<para_max(8)...
                    && par_new(9)>para_min(9) && par_new(9)<para_max(9)...
                    && par_new(10)>para_min(10) && par_new(10)<para_max(10)...
                    && par_new(11)>para_min(11) && par_new(11)<para_max(11)...
                    && par_new(12)>para_min(12) && par_new(12)<para_max(12)...
                    && par_new(13)>para_min(13) && par_new(13)<para_max(13)...
                    && par_new(14)>para_min(14) && par_new(14)<para_max(14)...
                    && par_new(15)>para_min(15) && par_new(15)<para_max(15)...
                    && par_new(16)>para_min(16) && par_new(16)<para_max(16)...
                    && par_new(17)>para_min(17) && par_new(17)<para_max(17)...
                    && par_new(18)>para_min(18) && par_new(18)<para_max(18)...
                    && par_new(19)>para_min(19) && par_new(19)<para_max(19)...
                    && par_new(20)>para_min(20) && par_new(20)<para_max(20)...
                    && par_new(21)>para_min(21) && par_new(21)<para_max(21)
                break;
            end
        end
        if ifbreak2 == 1 || ifbreak3 == 1
            break
        end
        
        % clear warning info
        lastwarn('');
        [~, soc_stock_mod, soc_mod] = ...
            fun_forward_simu(par_new, year_start, year_end, frocing_steady_state, frocing_transient);
        
        % original soc_modelled data
        [~, msgid] = lastwarn;
        if strcmp(msgid, 'MATLAB:singularMatrix') || strcmp(msgid, 'MATLAB:nearlySingularMatrix')
            cost_new = cost_old + 100;
        elseif isnan((sum(sum(soc_stock_mod(:, 4))))) == 1
            cost_new = cost_old + 100;
        else            
            % new cost function
            cost_new = cost_fun(obs_tot_soc, soc_stock_mod(:, 4));
            
        end
        delta_cost = cost_new-cost_old;

        % to determie whether to accept the newly soc_modelled data
        if min(1, exp(-delta_cost*log(simu + 2))) > rand
            upgrade2(iparallel, 1) = upgrade2(iparallel, 1) + 1;
            
            disp([datestr(now,'HH:MM:SS'), ' upgrade2 parallel ', num2str(iparallel), ': ', num2str(upgrade2(iparallel, 1)), ' out of ', num2str(simu)]);
            % record the accepted parameter value
            parameters_keep2(iparallel, :,upgrade2(iparallel, 1)) = par_new;
            cost_keep2(upgrade2(iparallel, 1), iparallel) = cost_new;
            par_old = par_new;
            cost_old = cost_new;
            coef = upgrade2(iparallel, 1)/simu;
            soc_mod_trace(iparallel, upgrade2(iparallel, 1), :) = soc_stock_mod(:, 4);
            
%                         scatter(obs_tot_soc, soc_stock_mod(:, 4))
%                         hold on
%                         plot([min(obs_tot_soc), max(obs_tot_soc)], [min(obs_tot_soc), max(obs_tot_soc)]) 
%                         
%                         plot(obs_tot_soc)
%                         hold on 
%                         plot(soc_stock_mod(:, 4))
            
        end
        parameters_rec2(iparallel, :, simu) = par_old;
        cost_record2(simu, iparallel) = cost_old;
        if simu > 4000
            % covars=sd*(cov(reshape(parameters_keep2(iparallel, :, 1:upgrade2(iparallel, 1)), [npara, upgrade2(iparallel, 1)])')) + sd*epsilon*diag(ones(npara, 1));
            covars=sd*(cov(reshape(parameters_rec2(iparallel, :, 1:simu), [npara, length(1:simu)])')) + sd*epsilon*diag(ones(npara, 1));
        end
        simu = simu + 1;
        if mod(simu, 1000) == 0 || simu == nsimu2
            da_summary.soc_mod_trace = soc_mod_trace;
            da_summary.parameters_keep2 = parameters_keep2;
            
            save([data_path, 'TRENDY_CLM5/output_data/da_reconstruct/clm5_HF_point_reconstruct.mat'], 'da_summary')
        end
        % to test if the acceptance rate is in a resonable level
        if simu == nsimu2
            if coef < 0.1 || coef > 0.5
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
    if ifbreak2 == 0 && ifbreak3 == 0
        iparallel = iparallel + 1;
    end
    if counter2 > 5
        error('Error')
    end
end
disp([datestr(now,'HH:MM:SS'), ' MCMC has been finished']);


%% Gelman-Rubin test
GR = nan(npara, 1);
between_var = nan(npara, 1);
within_var = nan(npara, 1);
for ipara = 1 : npara
    para_mean = nan(nparallel, 1);
    para_sum = nan(nparallel, 1);
    for iparallel = 1 : nparallel
        para_keep_middle = reshape(parameters_keep2(iparallel, 1, :), [nsimu2, 1]);
        valid_num = length(para_keep_middle(~isnan(para_keep_middle), 1));
        para_mean(iparallel, 1) = mean(parameters_keep2(iparallel, ipara, round(valid_num/2)+1:valid_num), 'omitnan');
        para_sum(iparallel, 1) = sum((parameters_keep2(iparallel, ipara, round(valid_num/2)+1:valid_num) - para_mean(iparallel, 1)).^2);
    end
    % calculate between variance and within variance
    between_var = ...
        valid_num/(nparallel-1)*nansum((para_mean - mean(para_mean)).^2);
    within_var = 1/nparallel/(valid_num-1)*sum(para_sum);
    GR(ipara) = sqrt(((within_var*(valid_num-1))/valid_num + between_var/valid_num)/within_var)';
end


%% Coefficient of efficiency
% set the first half number as burn-in
valid_num = nan(nparallel, 1);
soc_mod_parallel = nan(length(obs_tot_soc), nparallel);
soc_std_parallel = nan(length(obs_tot_soc), nparallel);
% col1: correlation, col2: RMSE
stat_parallel = nan(nparallel, 1);
for iparallel = 1 : nparallel
    para_keep_middle = reshape(parameters_keep2(iparallel, 1, :), [nsimu2, 1]);
    % find valid soc_modeled number
    valid_num(iparallel, 1) = length(para_keep_middle(~isnan(para_keep_middle), 1));
    
    soc_mod_parallel(:, iparallel) = mean(reshape(soc_mod_trace(iparallel, round(upgrade2(iparallel, 1)/2) : upgrade2(iparallel, 1), :),...
        [length(round(upgrade2(iparallel, 1)/2) : upgrade2(iparallel, 1)), length(obs_tot_soc)]));
    
    soc_std_parallel(:, iparallel) = std(reshape(soc_mod_trace(iparallel, round(upgrade2(iparallel, 1)/2) : upgrade2(iparallel, 1), :),...
        [length(round(upgrade2(iparallel, 1)/2) : upgrade2(iparallel, 1)), length(obs_tot_soc)]));
    
    SStot = sum((soc_mod_parallel(:, iparallel) - mean(obs_tot_soc)).^2);
    SSres = sum((soc_mod_parallel(:, iparallel) - obs_tot_soc).^2);
    % R squared
    stat_parallel(iparallel, 1) = 1 - SSres/SStot;
end

disp([datestr(now,'HH:MM:SS'), ' Program finished']);


