dev.off()
rm(list = ls())


library(R.matlab)
library(ggplot2)
library(cowplot)
library(jcolors)
library(viridis)

library(stableGR)
setwd('/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/TRENDY_CLM5/')

####################################################
# default CLM5 simulation
####################################################
year_start = 1900
year_end = 2010
soc_time_series = seq(from = year_start, to = year_end, by = (year_end-year_start)/((year_end-year_start+1)*12-1))

# simulation from CESM2 (transient)
default_clm5 = readMat('/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/RADIOCARBON/demo_harvard_forest/input_data/forcing_transient_harvard_forest.mat')
default_clm5_soc_stock = as.vector(default_clm5$frocing.transient[[14]])

####################################################
# DA results transient soc
####################################################
simu2_num = 50000
para_names =  c('diffus', 'cryo', 
                'q10', 'efolding', 
                'taucwd', 'taul1', 'taul2', 'tau4s1', 'tau4s2', 'tau4s3', 
                'fs1l1', 'fs1l2', 'fs2l3', 'fs2s1', 'fs3s1', 'fs1s2', 'fs3s2', 'fs1s3', 'fl2cwd', 
                'w-scaling', 'beta')

npara = length(para_names)

da_summary = readMat('/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/TRENDY_CLM5/output_data/da_reconstruct/clm5_HF_point_reconstruct.mat')
soc_mod_trace = da_summary$da.summary[[1]][1, , ]
posterior_para = da_summary$da.summary[[2]][1, , ]

valid_num = length(which(is.na(soc_mod_trace[ , 1]) == 0))
soc_mod_trace = soc_mod_trace[round(valid_num/2):valid_num, ]
posterior_para = posterior_para[ , round(valid_num/2):valid_num]

soc_mod_mean = apply(soc_mod_trace, 2, mean, na.rm = TRUE)
soc_mod_upper = apply(soc_mod_trace, 2, quantile, probs = 0.95, na.rm = TRUE)
soc_mod_lower = apply(soc_mod_trace, 2, quantile, probs = 0.05, na.rm = TRUE)

current_data = data.frame(rbind(cbind(soc_time_series, default_clm5_soc_stock, NA, NA, 1), cbind(soc_time_series, soc_mod_mean, soc_mod_upper, soc_mod_lower, 2)))

colnames(current_data) = c('time', 'soc', 'upper', 'lower', 'source')

color_scheme = c('#005AB5', '#DC3220')

jpeg('./figures/da_transient_soc.jpeg', width = 15, height = 10, units = 'in', res = 300)
ggplot() + 
  geom_line(data = current_data, aes(x = time, y = soc, color = as.factor(source)), size = 1.5) +
  geom_ribbon(data = current_data, aes(x = time, y = soc, ymin = lower, ymax = upper, color = NA, fill = as.factor(source)), alpha = 0.1) + 
  scale_color_manual(name = '', values = color_scheme, labels = c('CLM5 Default', 'Data Assimilation')) +
  scale_fill_manual(name = '', values = color_scheme, labels = c('CLM5 Default', 'Data Assimilation')) +
  scale_y_continuous(position = 'left') +
  theme_classic() +
  # theme(legend.position = 'None') + 
  theme(legend.justification = c(0, 1), legend.position = c(0, 1), legend.background = element_rect(fill = NA), legend.text.align = 0) +
  theme(legend.text = element_text(size = 30), legend.title = element_text(size = 25))  +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(1, 'inch')) +
  # add title
  labs(title = '', y = expression(paste('SOC Stock (g C m'^'-2', ')', sep = '')), x = 'Time (Year)') +
  # modify the position of title
  # modify the font sizea
  theme(plot.title = element_text(hjust = 0.8, vjust = -40, size = 30)) + 
  theme(axis.title = element_text(size = 32), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1)) +
  # modify the margin
  theme(plot.margin = unit(c(0, 0.4, 0.1, 0.1), 'inch')) +
  theme(axis.text = element_text(size=35)) 

dev.off()

####################################################
# DA results posterior para
####################################################
post_para_mean = apply(posterior_para, 1, mean, na.rm = TRUE)
ipara = 1
for (ipara in 1:npara) {
  current_data = data.frame(posterior_para[ipara, ])
  colnames(current_data) = 'para'
  p =
    ggplot() +
    geom_density(data = current_data, aes(x = para), color = 'grey34', size = 2) +
    geom_vline(xintercept = post_para_mean[ipara], size = 2, color = '#005AB5') + 
    
    # scale_x_continuous(limits = c(0, 1)) +
    # scale_color_manual(name = '', values = viridis(worker_num), labels = c()) +
    scale_y_continuous(position = 'left') +
    theme_classic() +
    # theme(legend.position = 'None') +
    theme(legend.justification = c(1, 0.8), legend.position = 'None', legend.background = element_rect(fill = NA), legend.text.align = 0) +
    theme(legend.text = element_text(size = 30), legend.title = element_text(size = 25))  +
    theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(1, 'inch')) +
    # add title
    labs(title = paste(para_names[ipara]), y = 'Density', x = 'Parameter value (scaled)') +
    # modify the position of title
    # modify the font sizea
    theme(plot.title = element_text(hjust = 0.5, vjust = 0, size = 35)) +
    theme(axis.title = element_text(size = 30), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1)) +
    # modify the margin
    theme(plot.margin = unit(c(0, 0.4, 0.1, 0.1), 'inch')) +
    theme(axis.text = element_text(size=30))
  
  eval(parse(text = paste('p', ipara, ' = p', sep = '')))
  
}
jpeg(paste('./figures/para_post_dist.jpeg', sep = ''), width = 42, height = 28, units = 'in', res = 300)
print(plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, nrow = 4))
dev.off()



