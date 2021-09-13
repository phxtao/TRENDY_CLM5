function cost = cost_fun_trendy(period_num, obs_tot_soc_transient, obs_tot_litter_transient, obs_tot_cwd_transient, obs_tot_hr_transient, mod_tot_soc_transient, mod_tot_litter_transient, mod_tot_cwd_transient, mod_tot_hr_transient)
weight_soc = ones(length(obs_tot_soc_transient), 1)/length(obs_tot_soc_transient);
weight_litter = ones(length(obs_tot_soc_transient), 1)/length(obs_tot_soc_transient);
weight_cwd = ones(length(obs_tot_soc_transient), 1)/length(obs_tot_soc_transient);
weight_hr = ones(length(obs_tot_soc_transient), 1)/length(obs_tot_soc_transient);

% weight_soc([1, end]) = 0.5;
% weight_litter([1, end]) = 0.5;
% weight_cwd([1, end]) = 0.5;

weight_soc(1) = 1;
weight_litter(1) = 1;
weight_cwd(1) = 1;

weight_total = [weight_soc; weight_litter; weight_cwd; weight_hr];

trans_obs_tot_soc_transient = nan(length(obs_tot_soc_transient), 1);
trans_mod_tot_soc_transient = nan(length(obs_tot_soc_transient), 1);
trans_obs_tot_litter_transient = nan(length(obs_tot_soc_transient), 1);
trans_mod_tot_litter_transient = nan(length(obs_tot_soc_transient), 1);
trans_obs_tot_cwd_transient = nan(length(obs_tot_soc_transient), 1);
trans_mod_tot_cwd_transient = nan(length(obs_tot_soc_transient), 1);
trans_obs_tot_hr_transient = nan(length(obs_tot_soc_transient), 1);
trans_mod_tot_hr_transient = nan(length(obs_tot_soc_transient), 1);

period_time_step = length(obs_tot_soc_transient)/period_num;
for iperiod = 1:period_num
    period_loc = (iperiod-1)*period_time_step+1:iperiod*period_time_step;
    middle_obs_tot_soc_transient = obs_tot_soc_transient(period_loc);
    middle_mod_tot_soc_transient = mod_tot_soc_transient(period_loc);
    middle_obs_tot_litter_transient = obs_tot_litter_transient(period_loc);
    middle_mod_tot_litter_transient = mod_tot_litter_transient(period_loc);
    middle_obs_tot_cwd_transient = obs_tot_cwd_transient(period_loc);
    middle_mod_tot_cwd_transient = mod_tot_cwd_transient(period_loc);
    middle_obs_tot_hr_transient = obs_tot_hr_transient(period_loc);
    middle_mod_tot_hr_transient = mod_tot_hr_transient(period_loc);
    
    trans_obs_tot_soc_transient(period_loc) = (middle_obs_tot_soc_transient - mean(middle_obs_tot_soc_transient))./std(middle_obs_tot_soc_transient);
    trans_mod_tot_soc_transient(period_loc) = (middle_mod_tot_soc_transient - mean(middle_obs_tot_soc_transient))./std(middle_obs_tot_soc_transient);
    
    trans_obs_tot_litter_transient(period_loc) = (middle_obs_tot_litter_transient - mean(middle_obs_tot_litter_transient))./std(middle_obs_tot_litter_transient);
    trans_mod_tot_litter_transient(period_loc) = (middle_mod_tot_litter_transient - mean(middle_obs_tot_litter_transient))./std(middle_obs_tot_litter_transient);
    
    trans_obs_tot_cwd_transient(period_loc) = (middle_obs_tot_cwd_transient - mean(middle_obs_tot_cwd_transient))./std(middle_obs_tot_cwd_transient);
    trans_mod_tot_cwd_transient(period_loc) = (middle_mod_tot_cwd_transient - mean(middle_obs_tot_cwd_transient))./std(middle_obs_tot_cwd_transient);
    
    trans_obs_tot_hr_transient(period_loc) = (middle_obs_tot_hr_transient - mean(middle_obs_tot_hr_transient))./std(middle_obs_tot_hr_transient);
    trans_mod_tot_hr_transient(period_loc) = (middle_mod_tot_hr_transient - mean(middle_obs_tot_hr_transient))./std(middle_obs_tot_hr_transient);
    
end
% trans_obs_total = [trans_obs_tot_soc_transient; trans_obs_tot_litter_transient; trans_obs_tot_cwd_transient; trans_obs_tot_hr_transient];
% trans_mod_total = [trans_mod_tot_soc_transient; trans_mod_tot_litter_transient; trans_mod_tot_cwd_transient; trans_mod_tot_hr_transient];
%
% cost = nansum(weight_total.*(trans_obs_total - trans_mod_total).^2)*100;


unexplained_soc = sum((obs_tot_soc_transient - mod_tot_soc_transient').^2)/sum((obs_tot_soc_transient - mean(obs_tot_soc_transient)).^2);
unexplained_litter = sum((obs_tot_litter_transient - mod_tot_litter_transient').^2)/sum((obs_tot_litter_transient - mean(obs_tot_litter_transient)).^2);
unexplained_cwd = sum((obs_tot_cwd_transient - mod_tot_cwd_transient').^2)/sum((obs_tot_cwd_transient - mean(obs_tot_cwd_transient)).^2);
unexplained_hr = sum((obs_tot_hr_transient - mod_tot_hr_transient').^2)/sum((obs_tot_hr_transient - mean(obs_tot_hr_transient)).^2);
cost = (unexplained_soc + unexplained_litter + unexplained_cwd + unexplained_hr)*100/4;

% plot(trans_obs_tot_soc_transient)
% hold on
% plot(trans_obs_tot_litter_transient)
% hold on
% plot(trans_obs_tot_cwd_transient)
% hold on
% plot(trans_obs_tot_hr_transient)
%

end






