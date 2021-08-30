dev.off()
rm(list = ls())


library(R.matlab)
library(ggplot2)
library(cowplot)
library(jcolors)
library(viridis)

setwd('/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/TRENDY_CLM5/')

####################################################
# default trendy simulation
####################################################
scenario_name = 's3'
grid_num = 10
worker_num = 6

data_path = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/'
year_start = 1700
year_end = 2019
month_num = 12
trendy_time_series = seq(from = (year_start+1/month_num), to = (year_end+1), by = 1/month_num)


igrid = 3
for (igrid in 1:grid_num) {
  # original simulated soc, unit kgc/m2
  obs_tot_soc_transient = readMat(paste(data_path, 'trendy_clm5/input_data/grid_s3/grid_', igrid, '_trendy2020_clm5_', scenario_name, '_TOTSOMC_170001_201912.mat', sep = ''))
  obs_tot_soc_transient = t(obs_tot_soc_transient$var.data.grid)/1000
  # original simulated litter, unit kgc/m2
  obs_tot_litter_transient = readMat(paste(data_path, 'trendy_clm5/input_data/grid_s3/grid_', igrid, '_trendy2020_clm5_', scenario_name, '_TOTLITC_170001_201912.mat', sep = ''))
  obs_tot_litter_transient = t(obs_tot_litter_transient$var.data.grid)/1000
  # original simulated cwdc, unit kgc/m2
  obs_tot_cwd_transient = readMat(paste(data_path, 'trendy_clm5/input_data/grid_s3/grid_', igrid, '_trendy2020_clm5_', scenario_name, '_CWDC_170001_201912.mat', sep = ''))
  obs_tot_cwd_transient = t(obs_tot_cwd_transient$var.data.grid)/1000
  # original simulated heterotrophic respiration, unit gc/m2/year
  obs_tot_hr_transient = readMat(paste(data_path, 'trendy_clm5/input_data/grid_s3/grid_', igrid, '_trendy2020_clm5_', scenario_name, '_HR_170001_201912.mat', sep = ''))
  obs_tot_hr_transient = t(obs_tot_hr_transient$var.data.grid)*365*24*60*60
  
  iworker = 1
  post_da_ensemble_soc = c()
  post_da_ensemble_litter = c()
  post_da_ensemble_cwd = c()
  post_da_ensemble_hr = c()
  for (iworker in 1:worker_num) {
    post_da_ensemble_middle = readMat(paste(data_path, 'trendy_clm5/output_data/da_reconstruct/trendy_clm5_s3_point_', igrid, '_reconstruct_worker_', iworker, '_ensembe_simu.mat', sep = ''))
    post_da_ensemble_middle = post_da_ensemble_middle$ensemble.simu
    post_da_ensemble_soc = cbind(post_da_ensemble_soc, post_da_ensemble_middle[[1]]/1000)
    post_da_ensemble_litter = cbind(post_da_ensemble_litter, post_da_ensemble_middle[[2]]/1000)
    post_da_ensemble_cwd = cbind(post_da_ensemble_cwd, post_da_ensemble_middle[[3]]/1000)
    post_da_ensemble_hr = cbind(post_da_ensemble_hr, post_da_ensemble_middle[[4]]*365*24*60*60)
  }
  
  color_scheme = c('#005AB5', '#DC3220')
  
  ##############plot soc transient
  current_data = data.frame(rbind(cbind(trendy_time_series, 
                                        obs_tot_soc_transient, 
                                        NA, 
                                        NA,
                                        1), 
                                  cbind(trendy_time_series, 
                                        apply(post_da_ensemble_soc, 1, mean, na.rm = TRUE),
                                        apply(post_da_ensemble_soc, 1, quantile, probs = 0.05, na.rm = TRUE),
                                        apply(post_da_ensemble_soc, 1, quantile, probs = 0.95, na.rm = TRUE), 
                                        2)))
  
  colnames(current_data) = c('time', 'simu', 'lower', 'upper', 'source')
  
  p_soc_transient = 
    ggplot() + 
    geom_line(data = current_data, aes(x = time, y = simu, color = as.factor(source)), size = 1.5) +
    geom_ribbon(data = current_data, aes(x = time, y = simu, ymin = lower, ymax = upper, color = NA, fill = as.factor(source)), alpha = 0.3) + 
    scale_color_manual(name = '', values = color_scheme, labels = c('CLM5 default', 'Data assimilation')) +
    scale_fill_manual(name = '', values = color_scheme, labels = c('CLM5 default', 'Data assimilation')) +
    scale_y_continuous(position = 'left') +
    theme_classic() +
    # theme(legend.position = 'None') + 
    theme(legend.justification = c(0, 1), legend.position = c(0, 1), legend.background = element_rect(fill = NA), legend.text.align = 0) +
    theme(legend.text = element_text(size = 30), legend.title = element_text(size = 25))  +
    theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(1, 'inch')) +
    # add title
    labs(title = '', y = expression(paste('SOC stock (kg C m'^'-2', ')', sep = '')), x = 'Time (Year)') +
    # modify the position of title
    # modify the font sizea
    theme(plot.title = element_text(hjust = 0.8, vjust = -40, size = 30)) + 
    theme(axis.title = element_text(size = 32), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1)) +
    # modify the margin
    theme(plot.margin = unit(c(0, 0.4, 0.1, 0.1), 'inch')) +
    theme(axis.text = element_text(size=35)) 
  
  ##############plot soc obs vs mod
  current_data = data.frame(cbind(obs_tot_soc_transient, 
                                  apply(post_da_ensemble_soc, 1, mean, na.rm = TRUE)))
  colnames(current_data) = c('obs', 'mod')
  
  coef_efficiency = 1 - sum((current_data$mod-current_data$obs)**2)/sum((current_data$obs-mean(current_data$obs))**2)
  
  p_soc_regression = 
    ggplot() + 
    geom_abline(slope = 1, intercept = 0, size = 2, color = 'black') + 
    geom_point(data = current_data, aes(x = obs, y = mod), color = 'dark red', shape = 16, size = 1.5) +
    scale_y_continuous(position = 'left') +
    theme_classic() +
    # add title
    labs(title = paste('Explained variation = ', round(coef_efficiency*100, 1), '%', sep = ''), x = expression(paste('SOC default CLM5 (kg C m'^'-2', ')', sep = '')), y = expression(paste('SOC DA CLM5 (kg C m'^'-2', ')', sep = ''))) +
    # modify the position of title
    # modify the font sizea
    theme(plot.title = element_text(hjust = 0.2, vjust = -10, size = 30)) + 
    theme(axis.title = element_text(size = 32), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1)) +
    # modify the margin
    theme(plot.margin = unit(c(0, 0.4, 0.1, 0.1), 'inch')) +
    theme(axis.text = element_text(size=35)) 
  
  
  ##############plot litter transient
  current_data = data.frame(rbind(cbind(trendy_time_series, 
                                        obs_tot_litter_transient, 
                                        NA, 
                                        NA,
                                        1), 
                                  cbind(trendy_time_series, 
                                        apply(post_da_ensemble_litter, 1, mean, na.rm = TRUE),
                                        apply(post_da_ensemble_litter, 1, quantile, probs = 0.05, na.rm = TRUE),
                                        apply(post_da_ensemble_litter, 1, quantile, probs = 0.95, na.rm = TRUE), 
                                        2)))
  
  colnames(current_data) = c('time', 'simu', 'lower', 'upper', 'source')
  
  p_litter_transient =
    ggplot() + 
    geom_line(data = current_data, aes(x = time, y = simu, color = as.factor(source)), size = 1.5) +
    geom_ribbon(data = current_data, aes(x = time, y = simu, ymin = lower, ymax = upper, color = NA, fill = as.factor(source)), alpha = 0.3) + 
    scale_color_manual(name = '', values = color_scheme, labels = c('CLM5 default', 'Data assimilation')) +
    scale_fill_manual(name = '', values = color_scheme, labels = c('CLM5 default', 'Data assimilation')) +
    scale_y_continuous(position = 'left') +
    theme_classic() +
    # theme(legend.position = 'None') + 
    theme(legend.justification = c(0, 1), legend.position = c(0, 1), legend.background = element_rect(fill = NA), legend.text.align = 0) +
    theme(legend.text = element_text(size = 30), legend.title = element_text(size = 25))  +
    theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(1, 'inch')) +
    # add title
    labs(title = '', y = expression(paste('Litter stock (kg C m'^'-2', ')', sep = '')), x = 'Time (Year)') +
    # modify the position of title
    # modify the font sizea
    theme(plot.title = element_text(hjust = 0.8, vjust = -40, size = 30)) + 
    theme(axis.title = element_text(size = 32), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1)) +
    # modify the margin
    theme(plot.margin = unit(c(0, 0.4, 0.1, 0.1), 'inch')) +
    theme(axis.text = element_text(size=35)) 
  
  ##############plot litter obs vs mod
  current_data = data.frame(cbind(obs_tot_litter_transient, 
                                  apply(post_da_ensemble_litter, 1, mean, na.rm = TRUE)))
  colnames(current_data) = c('obs', 'mod')
  
  coef_efficiency = 1 - sum((current_data$mod-current_data$obs)**2)/sum((current_data$obs-mean(current_data$obs))**2)
  
  p_litter_regression =
    ggplot() + 
    geom_abline(slope = 1, intercept = 0, size = 2, color = 'black') + 
    geom_point(data = current_data, aes(x = obs, y = mod), color = 'dark red', shape = 16, size = 1.5) +
    scale_y_continuous(position = 'left') +
    theme_classic() +
    # add title
    labs(title = paste('Explained variation = ', round(coef_efficiency*100, 1), '%', sep = ''), x = expression(paste('Litter default CLM5 (kg C m'^'-2', ')', sep = '')), y = expression(paste('Litter DA CLM5 (kg C m'^'-2', ')', sep = ''))) +
    # modify the position of title
    # modify the font sizea
    theme(plot.title = element_text(hjust = 0.2, vjust = -10, size = 30)) + 
    theme(axis.title = element_text(size = 32), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1)) +
    # modify the margin
    theme(plot.margin = unit(c(0, 0.4, 0.1, 0.1), 'inch')) +
    theme(axis.text = element_text(size=35)) 
  
  
  ##############plot cwd transient
  current_data = data.frame(rbind(cbind(trendy_time_series, 
                                        obs_tot_cwd_transient, 
                                        NA, 
                                        NA,
                                        1), 
                                  cbind(trendy_time_series, 
                                        apply(post_da_ensemble_cwd, 1, mean, na.rm = TRUE),
                                        apply(post_da_ensemble_cwd, 1, quantile, probs = 0.05, na.rm = TRUE),
                                        apply(post_da_ensemble_cwd, 1, quantile, probs = 0.95, na.rm = TRUE), 
                                        2)))
  
  colnames(current_data) = c('time', 'simu', 'lower', 'upper', 'source')
  
  p_cwd_transient =
    ggplot() + 
    geom_line(data = current_data, aes(x = time, y = simu, color = as.factor(source)), size = 1.5) +
    geom_ribbon(data = current_data, aes(x = time, y = simu, ymin = lower, ymax = upper, color = NA, fill = as.factor(source)), alpha = 0.3) + 
    scale_color_manual(name = '', values = color_scheme, labels = c('CLM5 default', 'Data assimilation')) +
    scale_fill_manual(name = '', values = color_scheme, labels = c('CLM5 default', 'Data assimilation')) +
    scale_y_continuous(position = 'left') +
    theme_classic() +
    # theme(legend.position = 'None') + 
    theme(legend.justification = c(0, 1), legend.position = c(0, 1), legend.background = element_rect(fill = NA), legend.text.align = 0) +
    theme(legend.text = element_text(size = 30), legend.title = element_text(size = 25))  +
    theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(1, 'inch')) +
    # add title
    labs(title = '', y = expression(paste('CWD stock (kg C m'^'-2', ')', sep = '')), x = 'Time (Year)') +
    # modify the position of title
    # modify the font sizea
    theme(plot.title = element_text(hjust = 0.8, vjust = -40, size = 30)) + 
    theme(axis.title = element_text(size = 32), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1)) +
    # modify the margin
    theme(plot.margin = unit(c(0, 0.4, 0.1, 0.1), 'inch')) +
    theme(axis.text = element_text(size=35)) 
  
  ##############plot cwd obs vs mod
  current_data = data.frame(cbind(obs_tot_cwd_transient, 
                                  apply(post_da_ensemble_cwd, 1, mean, na.rm = TRUE)))
  colnames(current_data) = c('obs', 'mod')
  
  coef_efficiency = 1 - sum((current_data$mod-current_data$obs)**2)/sum((current_data$obs-mean(current_data$obs))**2)
  
  p_cwd_regression =
    ggplot() + 
    geom_abline(slope = 1, intercept = 0, size = 2, color = 'black') + 
    geom_point(data = current_data, aes(x = obs, y = mod), color = 'dark red', shape = 16, size = 1.5) +
    scale_y_continuous(position = 'left') +
    theme_classic() +
    # add title
    labs(title = paste('Explained variation = ', round(coef_efficiency*100, 1), '%', sep = ''), x = expression(paste('CWD default CLM5 (kg C m'^'-2', ')', sep = '')), y = expression(paste('CWD DA CLM5 (kg C m'^'-2', ')', sep = ''))) +
    # modify the position of title
    # modify the font sizea
    theme(plot.title = element_text(hjust = 0.2, vjust = -10, size = 30)) + 
    theme(axis.title = element_text(size = 32), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1)) +
    # modify the margin
    theme(plot.margin = unit(c(0, 0.4, 0.1, 0.1), 'inch')) +
    theme(axis.text = element_text(size=35)) 
  
  
  ##############plot hr transient
  current_data = data.frame(rbind(cbind(trendy_time_series, 
                                        obs_tot_hr_transient, 
                                        NA, 
                                        NA,
                                        1), 
                                  cbind(trendy_time_series, 
                                        apply(post_da_ensemble_hr, 1, mean, na.rm = TRUE),
                                        apply(post_da_ensemble_hr, 1, quantile, probs = 0.05, na.rm = TRUE),
                                        apply(post_da_ensemble_hr, 1, quantile, probs = 0.95, na.rm = TRUE), 
                                        2)))
  
  colnames(current_data) = c('time', 'simu', 'lower', 'upper', 'source')
  
  p_hr_transient =
    ggplot() + 
    geom_line(data = current_data, aes(x = time, y = simu, color = as.factor(source)), size = 1.5) +
    geom_ribbon(data = current_data, aes(x = time, y = simu, ymin = lower, ymax = upper, color = NA, fill = as.factor(source)), alpha = 0.3) + 
    scale_color_manual(name = '', values = color_scheme, labels = c('CLM5 default', 'Data assimilation')) +
    scale_fill_manual(name = '', values = color_scheme, labels = c('CLM5 default', 'Data assimilation')) +
    scale_y_continuous(position = 'left') +
    theme_classic() +
    # theme(legend.position = 'None') + 
    theme(legend.justification = c(0, 1), legend.position = c(0, 1), legend.background = element_rect(fill = NA), legend.text.align = 0) +
    theme(legend.text = element_text(size = 30), legend.title = element_text(size = 25))  +
    theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(1, 'inch')) +
    # add title
    labs(title = '', y = expression(paste('HR (g C m'^'-2', 'yr'^'-1', ')', sep = '')), x = 'Time (Year)') +
    # modify the position of title
    # modify the font sizea
    theme(plot.title = element_text(hjust = 0.8, vjust = -40, size = 30)) + 
    theme(axis.title = element_text(size = 32), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1)) +
    # modify the margin
    theme(plot.margin = unit(c(0, 0.4, 0.1, 0.1), 'inch')) +
    theme(axis.text = element_text(size=35)) 
  
  ##############plot hr obs vs mod
  current_data = data.frame(cbind(obs_tot_hr_transient, 
                                  apply(post_da_ensemble_hr, 1, mean, na.rm = TRUE)))
  colnames(current_data) = c('obs', 'mod')
  
  coef_efficiency = 1 - sum((current_data$mod-current_data$obs)**2)/sum((current_data$obs-mean(current_data$obs))**2)
  
  p_hr_regression =
    ggplot() + 
    geom_abline(slope = 1, intercept = 0, size = 2, color = 'black') + 
    geom_point(data = current_data, aes(x = obs, y = mod), color = 'dark red', shape = 16, size = 1.5) +
    scale_y_continuous(position = 'left') +
    theme_classic() +
    # add title
    labs(title = paste('Explained variation = ', round(coef_efficiency*100, 1), '%', sep = ''), x =expression(paste('HR default CLM5 (g C m'^'-2', 'yr'^'-1', ')', sep = '')), y = expression(paste('HR DA CLM5 (g C m'^'-2', 'yr'^'-1', ')', sep = ''))) +
    # modify the position of title
    # modify the font sizea
    theme(plot.title = element_text(hjust = 0.2, vjust = -10, size = 30)) + 
    theme(axis.title = element_text(size = 32), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1)) +
    # modify the margin
    theme(plot.margin = unit(c(0, 0.4, 0.1, 0.1), 'inch')) +
    theme(axis.text = element_text(size=35)) 
  
  
  jpeg(paste('./figures/post_da_simu_', igrid, '.jpeg', sep = ''), width = 50, height = 20, units = 'in', res = 300)
  plot_grid(p_soc_transient, p_soc_regression, p_litter_transient, p_litter_regression, 
            p_cwd_transient, p_cwd_regression, p_hr_transient, p_hr_regression, 
            rel_widths = c(3, 2, 3, 2),
            nrow = 2)
  dev.off()
}





