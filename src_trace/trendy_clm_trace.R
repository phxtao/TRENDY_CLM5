dev.off()
rm(list = ls())

library(ncdf4)
library(ggplot2)
library(cowplot)

setwd('/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/TRENDY_CLM5')


data_path = '/Users/phoenix/Google_Drive/Tsinghua_Luo/Projects/DATAHUB/TRENDY_CLM5/'
nc_data = nc_open(paste(data_path, 'input_data/TRENDY2020_S3_CO2ClimateLUC_Matrix.clm2.h0.TOTSOMC_1m.170001-201912.nc', sep = ''))
grid_area = ncvar_get(nc_data, 'area')*1000*1000 # convert unit from km^2 to m^2
grid_land_fraction = ncvar_get(nc_data, 'landfrac')
world_land_mask = ncvar_get(nc_data, 'landmask')

var_origin = ncvar_get(nc_data, 'TOTSOMC_1m')

var_lat_time_series = apply(var_origin, c(2, 3), mean, na.rm = TRUE)
var_globe_time_series = apply(var_origin, c(3), mean, na.rm = TRUE)

plot(soc_1m)
