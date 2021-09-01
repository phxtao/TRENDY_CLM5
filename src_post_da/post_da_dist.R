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

####################################################
# para convergence
####################################################
para_names =  c('diffus', 'cryo', 
                'efolding', 
                'taucwd', 'taul1', 'taul2', 'tau4s1', 'tau4s2', 'tau4s3', 
                'fs1l1', 'fs1l2', 'fs2l3', 'fs2s1', 'fs3s1', 'fs1s2', 'fs3s2', 'fs1s3', 'fl2cwd', 
                'beta',
                'p4l1', 'p4l2', 'p4l3')

npara = length(para_names)


################ tier individual
igrid = 3
for (igrid in 3) { 
  # Gelman-Rubin test
  GR = array(NA, c(npara, 1))
  between_var = array(NA, c(npara, 1))
  within_var = array(NA, c(npara, 1))
  
  da_results = readMat(paste(data_path, 'trendy_clm5/output_data/da_reconstruct/trendy_clm5_s3_point_', igrid, '_reconstruct_valid_post_para.mat', sep = ''))
  da_results = da_results$valid.post.para
  ipara = 1
  data_dist_tot = c()
  for (ipara in 1:npara) {
    para_mean = array(NA, c(worker_num, 1))
    para_sum = array(NA, c(worker_num, 1))
    
    iworker = 1
    data_dist_middle = c()
    for (iworker in 1:worker_num) {
      opt_para = da_results[ , ipara, iworker]
        num_after_burn_in = length(which(is.na(opt_para) == 0))
        
        para_mean[iworker, 1] = mean(opt_para[1:num_after_burn_in], na.rm = TRUE)
        para_sum[iworker, 1] = sum((opt_para[1:num_after_burn_in] - para_mean[iworker, 1])**2, na.rm = TRUE)
        
        
        data_dist_middle = rbind(data_dist_middle,
                                 cbind(opt_para[1:num_after_burn_in], iworker))
        

    }
    
    data_dist_tot = cbind(data_dist_tot, data_dist_middle[ , 1])
    
    valid_worker_num = length(which(is.na(para_mean) == 0))
    between_var = num_after_burn_in/(valid_worker_num-1)*sum((para_mean - mean(para_mean, na.rm = TRUE))**2, na.rm = TRUE)
    within_var = 1/valid_worker_num/(num_after_burn_in-1)*sum(para_sum, na.rm = TRUE)
    GR[ipara] = sqrt(((within_var*(num_after_burn_in-1))/num_after_burn_in + between_var/num_after_burn_in)/within_var)
    
    print(paste('para ', ipara, ' GR ', GR[ipara]))
  }
  
  data_dist_tot = data.frame(cbind(data_dist_tot, data_dist_middle[ , 2]))
  colnames(data_dist_tot) = c(para_names, 'chain')
  
  
  ipara = 1
  for (ipara in 1:npara) {
    current_data = data_dist_tot[ , c(ipara, npara+1)]
    colnames(current_data) = c('para', 'chain')
    
    p =
      ggplot() +
      geom_density(data = current_data, aes(x = para, group = as.factor(chain), color = as.factor(chain)), color = 'grey34', size = 2) +
      # scale_x_continuous(limits = c(0, 1)) +
      # scale_color_manual(name = '', values = viridis(worker_num), labels = c()) +
      scale_y_continuous(position = 'left') +
      theme_classic() +
      # theme(legend.position = 'None') +
      theme(legend.justification = c(1, 0.8), legend.position = 'None', legend.background = element_rect(fill = NA), legend.text.align = 0) +
      theme(legend.text = element_text(size = 30), legend.title = element_text(size = 25))  +
      theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(1, 'inch')) +
      # add title
      labs(title = paste(para_names[ipara], ' GR = ', round(GR[ipara], 3)), y = 'Density', x = 'Parameter value (scaled)') +
      # modify the position of title
      # modify the font sizea
      theme(plot.title = element_text(hjust = 0.5, vjust = 0, size = 35)) +
      theme(axis.title = element_text(size = 30), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1)) +
      # modify the margin
      theme(plot.margin = unit(c(0, 0.4, 0.1, 0.1), 'inch')) +
      theme(axis.text = element_text(size=30))
    
    eval(parse(text = paste('p', ipara, ' = p', sep = '')))
    
  }
  
  jpeg(paste('./figures/para_post_dist_grid_', igrid, '.jpeg', sep = ''), width = 42, height = 28, units = 'in', res = 300)
  print(plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, nrow = 4))
  dev.off()
}

