################################# Script for running draft VISIT-e matrix function #################################

# Note: this is a draft approximation of VISIT-e C cycle assuming constant parameters and without considering leaf 
#       processes, C-N coupling, as well as affects of temperature and moisture on turnover rates.

library (raster)
library(ncdf4)
library(sp)
library(gridExtra)
library(ggplot2)
setwd("D:\\TRENDYv9\\DLEM")
setwd("D:\\LC_CMIP6\\8.3\\C4MIP")

######################################## Import CMIP6 outputs ######################################################

# pick up 1 site: 33.3 E, 50.0 N (62.8125 W, 17.5S)
#point<-SpatialPoints(cbind(33.3,50.0), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs"))
point<-SpatialPoints(cbind(-61.8125,-17.5), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs"))


# load all netCDF files in the folder
files<-list.files(pattern="..S2") # list all netCDF files in the folder
(var_names<-substr(files,9,regexpr(".nc",files)-1)) # derive variable names from file names
#(var_names<-substr(files,0,regexpr("_Lmon",files)-1)) 
dat<-data.frame(point)
for (i in 1:length(var_names)){
  r<-stack(files[i])
  r_dat<-extract(r,point) #original:rotate(r)
  dat<-cbind(dat,r_dat[1,])
}

(names(dat)<-c("lon","lat", var_names)) # assign variable names

#calculate wood pool from veg, leaf and root
#dat$cWood=dat$cVeg-dat$cLeaf-dat$cRoot
# correct fluxes from per second to per day
dat$gpp<-dat$gpp*86400
dat$npp<-dat$npp*86400
# dat$nppLeaf<-dat$nppLeaf*86400
# dat$nppWood<-dat$nppWood*86400
# dat$nppRoot<-dat$nppRoot*86400
# dat$fVegLitter<-dat$fVegLitter*86400
# dat$fLitterSoil<-dat$fLitterSoil*86400
# dat$rMaint<-dat$rMaint*86400
# dat$rGrowth<-dat$rGrowth*86400
dat$rh<-dat$rh*86400
dat$CWDC <- dat$CWDC/86400

# explore the data
summary(dat)

#save and load the data
write.csv(dat,"dat_npp_rh.csv")
dat<-read.csv("dat_npp_rh.csv")

########################################## Define matrix simulation function ########################################

# define number of days, years and months
#days = c(31,28,31,30,31,30,31,31,30,31,30,31)
nyears = 320 
tot_len = nyears#*12

## Define matrix simulation function
matrix_simu<-function(pa) {
  # B vector
  beta1=pa[1]; beta2=pa[2] 
  beta3=0.6*(1- beta1- beta2)
  beta4=0.4*(1- beta1- beta2)
  #B = c(beta1, beta2, beta3, 0, 0, 0, 0,0,0)   # allocation
  B = c(beta1, beta2, beta3, beta4, 0, 0, 0,0,0,0,0,0) 
  # A matrix
  # f41 = 1; f52 = 1; f63 = 1; f74 = pa[3]; f75 = pa[4]; f76 = pa[5]; f87 = pa[6]; f98 = pa[7]; f84=pa[8]; f85=pa[9]; f86=pa[10];
  # f78=pa[11]; f79=pa[12]; f97=pa[13]
  # A = c(-1,  0,   0,   0,   0,   0,   0,   0,   0,
  #        0,  -1,   0,   0,   0,   0,   0,   0,   0,
  #        0,   0,  -1,   0,   0,   0,   0,   0,   0,
  #        f41, 0,   0,  -1,   0,   0,   0,   0,   0,
  #        0, f52,   0,   0,  -1,   0,   0,   0,   0,
  #        0,   0, f63,   0,   0,  -1,   0,   0,   0,
  #        0,   0,   0, f74,  f75, f76, -1,   f78, f79,
  #        0,   0,   0, f84,  f85, f86,  f87,  -1, 0,
  #        0,   0,   0,   0,   0,   0,   f97,  f98,  -1 )
  # A = matrix(A, nrow = 9, byrow = TRUE)
  
  f51 = pa[3];   f52 = pa[4];        f53 = pa[5]; f54 = pa[6]; 
  f61 = 1-pa[3]; f62 = 0.0003-pa[4]; f63 = 1-pa[5]; f64 = 0.001-pa[6];
  f75=0.04368; f76=0.0624; f78=0.8; f710=0.368; f711=0.5; f712=0.4;
  f89=0.5;
  f95=0.09632; f96=0.1376;
  f105=0.21; f107=0.6435;f1011=0.3
  f115=0.013; f116=0.64;f1110=0.00015; 
  f127=0.0065; f1210=0.032;
  A = c(-1,  0,   0,   0,   0,   0,   0,   0,   0,  0,  0,  0,  #leaf
        0,  -1,   0,   0,   0,   0,   0,   0,   0,  0,  0,  0,  #sw
        0,   0,  -1,   0,   0,   0,   0,   0,   0,  0,  0,  0,  #fr
        0,   0,   0,  -1,   0,   0,   0,   0,   0,  0,  0,  0,  #cr
        f51, f52,f53, f54, -1,   0,   0,   0,   0,  0,  0,  0,  #AOM1
        f61, f62,f63, f64,  0,  -1,   0,   0,   0,  0,  0,  0,  #AOM2
        0,   0,   0,   0,  f75, f76,  -1, f78,  0,f710,f711,f712,  #SMB1
        0,   0,   0,   0,   0,   0,   0,  -1,  f89, 0,  0,  0,  #SMR
        0,   0,   0,   0,   f95, f96, 0,   0,  -1,  0,  0,  0,  #SMB2
        0,   0,   0,   0,   f105,0,  f107, 0,   0, -1,f1011,0,  #NOM
        0,   0,   0,   0,   f115,f116,0,   0,   0,f1110,-1, 0,  #DOM
        0,   0,   0,   0,   0,   0,  f127, 0,   0,f1210,0, -1) #PSOM
  
  A = matrix(A, nrow = 12, byrow = TRUE)
  # K matrix 
  # temp = c(pa[14],pa[15],pa[16], pa[17],pa[18], pa[19], pa[20], pa[21], pa[22])
  # K = rep(0,81)
  # K = matrix(K, nrow = 9, byrow = TRUE)
  # for (i in 1:9) { K[i,i] = temp[i] }
  clay=0.2028
  temp = c(0.001826,0.000137,0.002740, 0.0000685,0.000608, 0.000148, 0.000458/clay, 
           0.000375, 0.000458,0.010959,0.958904,0.002192)
  K = rep(0,144)
  K = matrix(K, nrow = 12, byrow = TRUE)
  for (i in 1:12) { K[i,i] = temp[i] }
  
  # X vector
  # x=rep(0,tot_len)
  # x_fin=data.frame(x,x,x,x,x,x,x,x,x); 
  # names(x_fin)=c("x1","x2","x3","x4","x5","x6","x7","x8","x9")
  # X vector
  x=rep(0,tot_len)
  x_fin=data.frame(x,x,x,x,x,x,x,x,x,x,x,x); 
  names(x_fin)=c("leaf","sp","froot","croot","AOM1","AOM2","SMB1","SR","SMB2","NOM","DOM","PSOM")
  
  # empty vectors for heterotrophic respiration, litterfall and humus formation
  rh_fin=rep(0,tot_len);
  f_veg_lit_fin=rep(0,tot_len);
  f_lit_soil_fin=rep(0,tot_len)
  rh_fin=rep(0,tot_len);   
  # f_AOM1_in=rep(0,tot_len);   
  # f_AOM2_in=rep(0,tot_len)
  # Initial carbon pool size derived from 1st year outputs where possible
  # x_init = c(dat$cLeaf[1], # leaf
  #            dat$cWood[1], # root
  #            dat$cRoot[1], # wood
  #            pa[23], # leaf litter - unknown
  #            dat$cLitterAbove[1]-pa[23], # wood litter
  #            dat$cLitterBelow[1], # root litter
  #            dat$cSoilFast[1], # soil fast
  #            dat$cSoilMedium[1], # soil medium
  #            dat$cSoilSlow[1])  # soil slow  
  
  ##########all the initial carbon pool are all example value, they are not true value###################
   x_init = c(0.1382,  #leaf
             0.5,     #sapwood
             0.6*(dat2$cVeg[1]/12-0.1382-0.5), #froot
             0.4*(dat2$cVeg[1]-0.1382-0.5), #croot 
             0.15,     #AOM1
             dat2$cLitter[1]/12-0.15, #AOM2
             0.05, #SMB1
             0.2, #SR
             0.05, #SMB2
             0.02, #NOM
             0.5, #DOM
             dat2$cSoil[1]/12-0.05-0.2-0.05-0.02-0.5) # PSOM
  X=x_init   # initialize carbon pools 
 
  # jj=1
  # for (y in 1:nyears){
  #   for (m in 1:12){
  #   npp_in = dat$npp[jj] 
  #   co2_rh = 0; f_veg_lit = 0; f_lit_soil = 0  
  #   # leaf fall starting from critical temperature
  #   # no ts
  #   if (dat$ts[jj]>281) {K[1,1]=0} else {K[1,1]<-pa[14]}
  #   # soil moisture factor ksi: monthly soil moisture data "mrso" multiplied by sensitivity parameters for different pools
  #   ksi=c(1,1,1,pa[24]*dat$mrso[jj],pa[25]*dat$mrso[jj],pa[26]*dat$mrso[jj],pa[27]*dat$mrso[jj],pa[28]*dat$mrso[jj],pa[29]*dat$mrso[jj])
  # 
  #   for (d in 1:days[m]) {
  #     # matrix equation with ksi
  #     X = X + B * npp_in + A %*% K %*% X * ksi
  #     # deriving rh from each litter and soil pool as 1 - sum of transfer coefficients 
  #     co2_rate = c(0,0,0, (1-f74-f84)*K[4,4],(1-f75-f85)*K[5,5],(1-f76-f86)*K[6,6], (1- f87-f97)*K[7,7], (1-f78-f98)*K[8,8], (1-f79)*K[9,9])
  #     co2=sum(co2_rate*X)
  #     co2_rh = co2_rh + co2/days[m]   # monthly average rh
  #     # deriving litterfall
  #     litterfall_rate = c(f41*K[1,1],f52*K[2,2],f63*K[3,3],0,0,0,0,0,0)
  #     litterfall=sum(litterfall_rate*X)
  #     f_veg_lit=f_veg_lit+litterfall/days[m]
  #     # deriving humus formation
  #     humification_rate = c(0,0,0,f74*K[4,4],f75*K[5,5],f76*K[6,6],f84*K[4,4],f85*K[5,5],f86*K[6,6])
  #     humification=sum(litterfall_rate*X)
  #     f_lit_soil=f_lit_soil+humification/days[m]
  #     }
  #   x_fin[jj,]=X
  #   rh_fin[jj]=co2_rh
  #   f_veg_lit_fin[jj]=f_veg_lit
  #   f_lit_soil_fin[jj]=f_lit_soil
  #   jj= jj+1
  #   }
  # } 
  
  jj=1
  for (y in 1:nyears){
    for (m in 1:12){
      npp_in = dat$npp[jj] 
      co2_rh = 0  
      f_veg_lit = 0; 
      #f_lit_soil = 0  
      # leaf fall starting from critical temperature
      # no ts
      #if (dat$ts[jj]>281) {K[1,1]=0} else {K[1,1]<-pa[14]}
      # soil moisture factor ksi: monthly soil moisture data "mrso" multiplied by sensitivity parameters for different pools
      #ksi=c(1,1,1,pa[24]*dat$mrso[jj],pa[25]*dat$mrso[jj],pa[26]*dat$mrso[jj],pa[27]*dat$mrso[jj],pa[28]*dat$mrso[jj],pa[29]*dat$mrso[jj])
      
      for (d in 1:days[m]) {
        # matrix equation with ksi
        X = X + B * npp_in + A %*% K %*% X
        # deriving rh from each litter and soil pool as 1 - sum of transfer coefficients 
        co2_rate = c(0,0,0,0,
                     (1-f75-f95-f105-f115)*K[5,5],
                     (1-f76-f96-f116)*K[6,6],
                     (1-f107-f127)*K[7,7],
                     (1-f78)*K[8,8], 
                     (1-f89)*K[9,9],
                     (1-f710-f1110-f1210)*K[10,10],
                     (1-f711)*K[11,11],
                     (1-f712)*K[12,12])
        co2=sum(co2_rate*X)
        co2_rh = co2_rh + co2/days[m]   # monthly average rh
        # deriving litterfall
        # f_veg_lit=0
        litterfall_rate = c((f51+f61)*K[1,1],(f52+f62)*K[2,2],(f53+f63)*K[3,3],(f54+f64)*K[4,4],0,0,0,0,0,0,0,0)
        litterfall=sum(litterfall_rate*X)
        f_veg_lit=f_veg_lit+litterfall/days[m]
        # # deriving humus formation
        f_lit_soil=0
        humification_rate = c(0,0,0,0,
                              (f75+f95+f105+f115)*K[5,5],
                              (f76+f96+116)*K[6,6],
                              (f107+f127)*K[7,7],
                              f78*K[8,8],
                              f89*K[9,9],
                              (f710+f1110+f1210)*K[10,10],
                              (f711+f1011)*K[11,11],
                              f712*K[10,10])
        humification=sum(litterfall_rate*X)
        f_lit_soil=f_lit_soil+humification/days[m]
      }
      x_fin[jj,]=X
      rh_fin[jj]=co2_rh
      f_veg_lit_fin[jj]=f_veg_lit
      f_lit_soil_fin[jj]=f_lit_soil
      jj= jj+1
    }
  } 
  # outputs: C pools, heterotrophic respiration, litter, soil pool
  result<-list(x_fin, rh_fin,f_veg_lit_fin,f_lit_soil_fin)
  return(result)
}

################################################### Run matrix model ################################################
# paramNames=c("beta1","beta2","f74","f75", "f76", "f87", "f98", "f84", "f85", "f86", "f78", "f79", "f97", 
#              "k1", "k2", "k3", "k4", "k5", "k6", "k7", "k8", "k9", "x4_init", 
#              "sens_moist_x4","sens_moist_x5","sens_moist_x6","sens_moist_x7","sens_moist_x8","sens_moist_x9")
paramNames=c("beta1","beta2","leaf_aom1","sp_aom1", "fr_aom1", "cr_aom1")

# example parameters for running the model
# pa<-c(
#   0.4, # beta 1 - allocation to leaf 
#   0.2, # beta 2 - allocation to wood 
#   0.55, # f74 - fast humus formation from leaf litter
#   0.25, # f75 - fast humus formation from wood litter
#   0.45, # f76 - fast humus formation from root litter
#   0.35, # f87 - transfer from fast SOM to medium SOM
#   0.05, # f98 - transfer from medium SOM to slow SOM
#   0.05, # f84 - medium humus formation from leaf litter
#   0.35, # f85 - medium humus formation from wood litter
#   0.15, # f86 - medium humus formation from root litter
#   0.42, # f78 - transfer from medium SOM to fast SOM  
#   0.45, # f79 - transfer from slow SOM to fast SOM   
#   0.004, # f97 - transfer from fast SOM to slow SOM    
#   1/60, # k1 - turnover rate of leaf pool
#   1/(365*12), # k2 - turnover rate of wood pool
#   1/(365*5), # k3 - turnover rate of root pool
#   1/(365*3), # k4 - turnover rate of leaf litter pool
#   1/(365*7), # k5 - turnover rate of wood litter pool
#   1/(365*4), # k6 - turnover rate of root litter pool
#   0.3/(365*2), # k7 - turnover rate of fast SOM pool
#   0.3/(365*3.5), # k8 - turnover rate of medium SOM pool
#   0.3/(365*80), # k9 - turnover rate of slow SOM pool  
#   0.3, # X4_init - initial size of leaf litter pool
#   0.0003, # sensitivity to moisture of leaf litter
#   0.0003, # sensitivity to moisture of wood litter
#   0.003, # sensitivity to moisture of root litter
#   0.003, # sensitivity to moisture of fast soil
#   0.003, # sensitivity to moisture of medium soil
#   0.003  # sensitivity to moisture of slow soil
# )

###################all are example value, these parameters need to be done in DA#############
pa<-c(
  0.4, # beta 1 - allocation to leaf 
  0.2, # beta 2 - allocation to sw 
  0.25, # f51 - leaf to AOM1
  0.00016, # f52 - sapwood to AOM1
  0.45, # f53 - fine root to AOM1
  0.0007 # f54 - coarse root to AOM1  
)

# test-run the model
test<-matrix_simu(pa)
str(test)
# test$veg=test[[1]]$x1+test[[1]]$x2+test[[1]]$x3+test[[1]]$x4
# test$litter=test[[1]]$x5+test[[1]]$x6
# test$soil=test[[1]]$x7+test[[1]]$x8+test[[1]]$x9+test[[1]]$x10+test[[1]]$x11+test[[1]]$x12

# view results
#View(test[[1]]) # monthly pool sizes
#View(test[[2]]) # monthly rh
#View(test[[3]]) # monthly litterfall
#View(test[[4]]) # monthly humus formation

# Compare modeled results with CMIP6 output for 1st 10 years (blue = modeled, red = CMIP6 output)
# par(mfrow=c(3, 4)) #
# {
# plot(test[[1]]$x1[1:240], type="l", col="red",  xlab="month", 
#       ylim=c(min(min(test[[1]]$x1[1:240]),min(dat$cLeaf[1:240])),max(max(test[[1]]$x1[1:240]),max(dat$cLeaf[1:240]))), 
#       xaxp  = c(0, 120, 10), ylab="cLeaf", main="Leaf Pool modelled vs CMIP6")
# lines(dat$cLeaf[1:240], col="blue")
# 
# plot(test[[1]]$x2[1:240], type="l", col="red", xlab="month", 
#       ylim=c(min(min(test[[1]]$x2[1:240]),min(dat$cWood[1:240])),max(max(test[[1]]$x2[1:240]),max(dat$cWood[1:240]))), 
#       xaxp  = c(0, 120, 10),ylab="cWood", main="Wood Pool modelled vs CMIP6")
# lines(dat$cWood[1:240], col="blue")
# 
# plot(test[[1]]$x3[1:240], type="l", col="red", xlab="month", 
#       ylim=c(min(min(test[[1]]$x3[1:240]),min(dat$cRoot[1:240])),max(max(test[[1]]$x3[1:240]),max(dat$cRoot[1:240]))), 
#       xaxp  = c(0, 120, 10),ylab="cRoot", main="Root Pool modelled vs CMIP6")
# lines(dat$cRoot[1:240], col="blue")
# 
# plot(test[[1]]$x4[1:240]+test[[1]]$x5[1:240], type="l", col="red", 
#       ylim=c(min(min(test[[1]]$x4[1:240]+test[[1]]$x5[1:240]),min(dat$cLitterAbove[1:240])),max(max(test[[1]]$x4[1:240]+test[[1]]$x5[1:240]),max(dat$cLitterAbove[1:240]))), 
#       xlab="month", xaxp  = c(0, 120, 10),ylab="cLitterAbove", main="Above-ground Litter (leaf+wood) modelled vs CMIP6")
# lines(dat$cLitterAbove[1:240], col="blue")
# 
# plot(test[[1]]$x6[1:240], type="l", col="red", 
#       ylim=c(min(min(test[[1]]$x6[1:240]),min(dat$cLitterBelow[1:240])),max(max(test[[1]]$x6[1:240]),max(dat$cLitterBelow[1:240]))), 
#       xlab="month", xaxp  = c(0, 120, 10),ylab="cLitterBelow", main="Below-ground Litter (root) modelled vs CMIP6")
# lines(dat$cLitterBelow[1:240], col="blue")
# 
# plot(test[[1]]$x7[1:240], type="l", col="red", 
#       ylim=c(min(min(test[[1]]$x7[1:240]),min(dat$cSoilFast[1:240])),max(max(test[[1]]$x7[1:240]),max(dat$cSoilFast[1:240]))), 
#       xlab="month", xaxp  = c(0, 120, 10),ylab="cSoilFast", main="Fast SOM pool modelled vs CMIP6")
# lines(dat$cSoilFast[1:240], col="blue")
# 
# plot(test[[1]]$x8[1:240], type="l", col="red",  
#       ylim=c(min(min(test[[1]]$x8[1:240]),min(dat$cSoilMedium[1:240])),max(max(test[[1]]$x8[1:240]),max(dat$cSoilMedium[1:240]))), 
#       xlab="month", xaxp  = c(0, 120, 10),ylab="cSoilMEdium", main="Medium, SOM pool modelled vs CMIP6")
# lines(dat$cSoilMedium[1:240], col="blue")
# 
# plot(test[[1]]$x9[1:240], type="l", col="red",  
#       ylim=c(min(min(test[[1]]$x9[1:240]),min(dat$cSoilSlow[1:240])),max(max(test[[1]]$x9[1:240]),max(dat$cSoilSlow[1:240]))), 
#       xlab="month", xaxp  = c(0, 120, 10),ylab="cSoilSlow", main="Slow SOM pool modelled vs CMIP6")
# lines(dat$cSoilSlow[1:240], col="blue")
# 
# plot(test[[2]][1:240], type="l", col="red",  
#       ylim=c(min(min(test[[2]][1:240]),min(dat$rh[1:240])),max(max(test[[2]][1:240]),max(dat$rh[1:240]))), 
#       xlab="month", xaxp  = c(0, 120, 10),ylab="Rh", main="Heterotrophic respiration modelled vs CMIP6")
# lines(dat$rh[1:240], col="blue")
# 
# plot(test[[3]][1:240], type="l", col="red", 
#       ylim=c(min(min(test[[3]][1:240]),min(dat$fVegLitter[1:240])),max(max(test[[3]][1:240]),max(dat$fVegLitter[1:240]))), 
#       xlab="month", xaxp  = c(0, 120, 10),ylab="fVegLitter", main="Litterfall modelled vs CMIP6")
#   lines(dat$fVegLitter[1:240], col="blue")
# 
# plot(test[[4]][1:240], type="l", col="red", 
#       ylim=c(min(min(test[[4]][1:240]),min(dat$fLitterSoil[1:240])),max(max(test[[4]][1:240]),max(dat$fLitterSoil[1:240]))), 
#       xlab="month", xaxp  = c(0, 120, 10),ylab="fLitterSoil", main="Humus formation modelled vs CMIP6")
#   lines(dat$fLitterSoil[1:240], col="blue")
#   
# plot(2, xlim=c(0,1), ylim=c(0,1), axes = F, main="legend", ylab="")
#   legend(0.1, 0.9, legend=c("CMIP-6 Output", "Modelled"),
#          col=c("blue", "red"), lty=1, cex=1)
# }

#annual vegetation pool
test$veg=test[[1]]$leaf+test[[1]]$sp+test[[1]]$froot+test[[1]]$croot

cveg=NULL
for (i in 0:319) {
  num=sum(test$veg[(12*i+1):(12*i+12)])
  cveg=c(cveg,num)
}

#annual litter pool
clitter=NULL
for (i in 0:319) {
  num=sum(test[[3]][(12*i+1):(12*i+12)])
  clitter=c(clitter,num)
}
#annual soil pool
csoil=NULL
for (i in 0:319) {
  num=sum(test[[4]][(12*i+1):(12*i+12)])
  csoil=c(csoil,num)
}

par(mfrow=c(4, 2),mar) 
{
 plot(cveg, ylim=c(0,max(cveg)),type="l", col="red", xlab="time",ylab="cVeg", main="matrix cVeg vs TRENDYv9")
 plot(dat2$cVeg, col="blue",type='l')

 plot(clitter,ylim=c(0,max(clitter)),type="l", col="red",xlab="time", ylab="cLitter", main="matrix cLitter vs TRENDYv9")
 plot(dat2$cLitter, col="blue",type='l')

 plot(csoil, type="l", col="red",ylim=c(0,max(csoil)),xlab="time", ylab="cSoil", main="matrix cSoil vs TRENDYv9")
 plot(dat2$cSoil, col="blue",type='l')

 plot(test[[2]], type="l", col="red",ylim=c(0,max(test[[2]])),xlab="time", ylab="rh", main="matrix rh vs TRENDYv9")
 plot(dat$rh,col='blue',type='l')
}

#######################################################MCMC###############################################################
# assign min-max values to all parameters

# "beta1","beta2","f74","f75", "f76", "f87", "f98", "f84", "f85", "f86", "f78", "f79", "f97", "k1", "k2", "k3",         "k4",          "k5",     "k6",               "k7",          "k8",         "k9",     "x4_init",  sensitivity to moisture 
# c_min=c(0,  0,   0.1,  0.1,  0.1,  0.1,  0.01,  0.01, 0.1,  0.1,   0.1,  0.1,  0.001, 1/(365*2),  1/(365*60), 1/(365*30), 1/(365*60),   1/(365*10), 1/(365*30),    0.01/(365*2), 0.01/(365*3.5), 0.01/(365*80),  0.2, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001)
# pa=c(0.4, 0.2, 0.55,  0.25, 0.45, 0.35, 0.05,  0.05, 0.35, 0.15, 0.42, 0.45, 0.004,  1/60,  1/(365*12),   1/(365*5),  1/(365*3),    1/(365*7),    1/(365*4),  0.3/(365*2), 0.3/(365*3.5),    0.3/(365*80),     0.3, 0.0003, 0.0003, 0.003, 0.003, 0.003, 0.003)
# c_max=c(1,  1,   0.9,   0.9,  0.9,  0.9,  0.2,  0.2,  0.9,  0.9,  0.9,  0.9,  0.1 ,   1/30,    1/365,   1/(365*0.5),   1/365,     1/(365*0.5),  1/(365*0.5),    1/(365*2),  1/(365*3.5),    1/(365*80),        0.4,  0.003,  0.003,  0.03,  0.03,  0.03,   0.03)

c_min=c(0,  0,   0.1,  0,  0.1,  0)
#pa=c(0.4, 0.2, 0.55,  0.25, 0.45, 0.35, 0.05,  0.05, 0.35, 0.15, 0.42, 0.45, 0.004,  1/60,  1/(365*12),   1/(365*5),  1/(365*3),    1/(365*7),    1/(365*4),  0.3/(365*2), 0.3/(365*3.5),    0.3/(365*80),     0.3, 0.0003, 0.0003, 0.003, 0.003, 0.003, 0.003)
c_max=c(1,  1,   0.9, 0.0003, 1, 0.001)

# set random numbers table to get reproducible results
set.seed(8192021)

paramNum=length(pa)
# number of simulations
nsimu = 2000

upgraded=0
J_last = 400
C_op = pa

C_upgraded = rep(0,paramNum*nsimu)
C_upgraded = matrix(C_upgraded, nrow = paramNum, byrow = TRUE)
J_upgraded = rep(0, nsimu)

isQualified=function(c) {
  flag = T
for (i in 1:paramNum){
  # check for min-max values
  if (c[i] > c_max[i] || c[i] < c_min[i]){
     flag = F
      break}
  # check for total outgoing flux fraction for each pool not exceeding 1
  if(c[1] + c[2] > 1){
   flag = F
   break}
  # if(c[3] + c[8] > 1){
  #   flag = F
  #   break}
  # if(c[4] + c[9] > 1){
  #   flag = F
  #   break}
  # if(c[5] + c[10] > 1){
  #   flag = F
  #   break}
  # if(c[6] + c[13] > 1){
  #   flag = F
  #   break}
  # if(c[7] + c[11] > 1){
  #   flag = F
  #   break}
}
return (flag)
}

GenerateParamValues=function(c_op){
 flag = T
 while (flag){
  c_new = c_op + (runif((paramNum)) - 0.5)*(c_max - c_min)/15.0
  if (isQualified(c_new)){  flag = F }
 }
return (c_new)
}

for (simu in 1:nsimu){
  print (paste0("simulation ",simu))
  c_new = GenerateParamValues(C_op)
 x_simu = matrix_simu(c_new)[[1]]
 rh_simu = matrix_simu(c_new)[[2]]
 # litterfall_simu = matrix_simu(c_new)[[3]]
 # humification_simu = matrix_simu(c_new)[[4]]
 # cost functions
 # J_obj1 = mean (( x_simu$x1 - dat$cLeaf )**2)/(2*var(dat$cLeaf))
 # J_obj2 = mean (( x_simu$x2 - dat$cWood )**2)/(2*var(dat$cWood))
 # J_obj3 = mean (( x_simu$x3 - dat$cRoot )**2)/(2*var(dat$cRoot))
 # J_obj4 = mean (( litterfall_simu - dat$fVegLitter )**2)/(2*var(dat$fVegLitter))
 # J_obj5 = mean (( x_simu$x4+x_simu$x5 - dat$cLitterAbove )**2)/(2*var(dat$cLitterAbove))
 # J_obj6 = mean (( x_simu$x6 - dat$cLitterBelow )**2)/(2*var(dat$cLitterBelow))
 # J_obj7 = mean (( x_simu$x7 - dat$cSoilFast )**2)/(2*var(dat$cSoilFast))
 # J_obj8 = mean (( x_simu$x8 - dat$cSoilMedium )**2)/(2*var(dat$cSoilMedium))
 # J_obj9 = mean (( x_simu$x9 - dat$cSoilSlow )**2)/(2*var(dat$cSoilSlow))
 # J_obj10 = mean (( humification_simu - dat$fLitterSoil )**2)/(2*var(dat$fLitterSoil))
 # J_obj11 = mean (( rh_simu - dat$rh )**2)/(2*var(dat$rh))
 J_obj1 = mean (( x_simu$x1+x_simu$x2+x_simu$x3+x_simu$x4 - dat$cVeg )**2)/(2*var(dat$cVeg))
 J_obj2 = mean (( x_simu$x5+x_simu$x6 - dat$cLitter)**2)/(2*var(dat$cLitter))
 J_obj3 = mean (( x_simu$x7+x_simu$x8+x_simu$x9+x_simu$x10+x_simu$x11+x_simu$x12 - dat$cSoil)**2)/(2*var(dat$cSoil))
 J_obj4 = mean (( rh_simu - dat$rh )**2)/(2*var(dat$rh))
 # total cost function
#J_new= (J_obj1 + J_obj2 + J_obj3 + J_obj4 + J_obj5 + J_obj6 + J_obj7 + J_obj8 + J_obj9 + J_obj10 + J_obj11 )/200
 #J_new= (J_obj1 + J_obj2 + J_obj3 + J_obj4 )/200
 J_new= (J_obj1 + J_obj2 + J_obj3)/200
 delta_J =  J_last - J_new;

 randNum = runif(1)
 if (min(1.0, exp(delta_J)) > randNum) {
   C_op=c_new;
  J_last=J_new;
  C_upgraded[,upgraded]=C_op;
  J_upgraded[upgraded]=J_last; 
  upgraded=upgraded+1 }
}

# results: distributions of parameters and cost functions
df=data.frame(C_upgraded)
df_j=data.frame(J_upgraded)

write.csv(df,'dlem_demo_da_aa.csv')
write.csv(df_j,'dlem_demo_da_j_aa.csv')

# plot distribution of optimized parameters
distr<-t(df) #transpose optimized parameter distributions
distr<-as.data.frame(distr)
names(distr)=paramNames
# { # create density plots of parameter distributions using ggplot2
# p1=ggplot(distr, aes(x=beta1)) + geom_density() + geom_vline(aes(xintercept=median(beta1)),color="blue", linetype="dashed", size=1)
# p2=ggplot(distr, aes(x=beta2)) + geom_density() + geom_vline(aes(xintercept=median(beta2)),color="blue", linetype="dashed", size=1)
# p3=ggplot(distr, aes(x=f74)) + geom_density() + geom_vline(aes(xintercept=median(f74)),color="blue", linetype="dashed", size=1)
# p4=ggplot(distr, aes(x=f75)) + geom_density() + geom_vline(aes(xintercept=median(f75)),color="blue", linetype="dashed", size=1)
# p5=ggplot(distr, aes(x=f76)) + geom_density() + geom_vline(aes(xintercept=median(f76)),color="blue", linetype="dashed", size=1)
# p6=ggplot(distr, aes(x=f87)) + geom_density() + geom_vline(aes(xintercept=median(f87)),color="blue", linetype="dashed", size=1)
# p7=ggplot(distr, aes(x=f98)) + geom_density() + geom_vline(aes(xintercept=median(f98)),color="blue", linetype="dashed", size=1)
# p8=ggplot(distr, aes(x=f84)) + geom_density() + geom_vline(aes(xintercept=median(f84)),color="blue", linetype="dashed", size=1)
# p9=ggplot(distr, aes(x=f85)) + geom_density() + geom_vline(aes(xintercept=median(f85)),color="blue", linetype="dashed", size=1)
# p10=ggplot(distr, aes(x=f86)) + geom_density() + geom_vline(aes(xintercept=median(f86)),color="blue", linetype="dashed", size=1)
# p11=ggplot(distr, aes(x=f78)) + geom_density() + geom_vline(aes(xintercept=median(f79)),color="blue", linetype="dashed", size=1)
# p12=ggplot(distr, aes(x=f79)) + geom_density() + geom_vline(aes(xintercept=median(f76)),color="blue", linetype="dashed", size=1)
# p13=ggplot(distr, aes(x=f97)) + geom_density() + geom_vline(aes(xintercept=median(f97)),color="blue", linetype="dashed", size=1)
# p14=ggplot(distr, aes(x=k1)) + geom_density() + geom_vline(aes(xintercept=median(k1)),color="blue", linetype="dashed", size=1)
# p15=ggplot(distr, aes(x=k2)) + geom_density() + geom_vline(aes(xintercept=median(k2)),color="blue", linetype="dashed", size=1)
# p16=ggplot(distr, aes(x=k3)) + geom_density() + geom_vline(aes(xintercept=median(k3)),color="blue", linetype="dashed", size=1)
# p17=ggplot(distr, aes(x=k4)) + geom_density() + geom_vline(aes(xintercept=median(k4)),color="blue", linetype="dashed", size=1)
# p18=ggplot(distr, aes(x=k5)) + geom_density() + geom_vline(aes(xintercept=median(k5)),color="blue", linetype="dashed", size=1)
# p19=ggplot(distr, aes(x=k6)) + geom_density() + geom_vline(aes(xintercept=median(k6)),color="blue", linetype="dashed", size=1)
# p20=ggplot(distr, aes(x=k7)) + geom_density() + geom_vline(aes(xintercept=median(k7)),color="blue", linetype="dashed", size=1)
# p21=ggplot(distr, aes(x=k8)) + geom_density() + geom_vline(aes(xintercept=median(k8)),color="blue", linetype="dashed", size=1)
# p22=ggplot(distr, aes(x=k9)) + geom_density() + geom_vline(aes(xintercept=median(k9)),color="blue", linetype="dashed", size=1)
# p23=ggplot(distr, aes(x=x4_init)) + geom_density() + geom_vline(aes(xintercept=median(x4_init)),color="blue", linetype="dashed", size=1)
# p24=ggplot(distr, aes(x=sens_moist_x4)) + geom_density() + geom_vline(aes(xintercept=median(sens_moist_x4)),color="blue", linetype="dashed", size=1)
# p25=ggplot(distr, aes(x=sens_moist_x5)) + geom_density() + geom_vline(aes(xintercept=median(sens_moist_x5)),color="blue", linetype="dashed", size=1)
# p26=ggplot(distr, aes(x=sens_moist_x6)) + geom_density() + geom_vline(aes(xintercept=median(sens_moist_x6)),color="blue", linetype="dashed", size=1)
# p27=ggplot(distr, aes(x=sens_moist_x7)) + geom_density() + geom_vline(aes(xintercept=median(sens_moist_x7)),color="blue", linetype="dashed", size=1)
# p28=ggplot(distr, aes(x=sens_moist_x8)) + geom_density() + geom_vline(aes(xintercept=median(sens_moist_x8)),color="blue", linetype="dashed", size=1)
# p29=ggplot(distr, aes(sens_moist_x9)) + geom_density() + geom_vline(aes(xintercept=median(sens_moist_x9)),color="blue", linetype="dashed", size=1)
# 
# # arrange density plots in a 5x4 grid
# grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29, nrow = 6)
# }

{ # create density plots of parameter distributions using ggplot2
  p1=ggplot(distr, aes(x=beta1)) + geom_density() + geom_vline(aes(xintercept=median(beta1)),color="blue", linetype="dashed", size=1)
  p2=ggplot(distr, aes(x=beta2)) + geom_density() + geom_vline(aes(xintercept=median(beta2)),color="blue", linetype="dashed", size=1)
  p3=ggplot(distr, aes(x=leaf_aom1)) + geom_density() + geom_vline(aes(xintercept=median(leaf_aom1)),color="blue", linetype="dashed", size=1)
  p4=ggplot(distr, aes(x=sp_aom1)) + geom_density() + geom_vline(aes(xintercept=median(sp_aom1)),color="blue", linetype="dashed", size=1)
  p5=ggplot(distr, aes(x=fr_aom1)) + geom_density() + geom_vline(aes(xintercept=median(fr_aom1)),color="blue", linetype="dashed", size=1)
  p6=ggplot(distr, aes(x=cr_aom1)) + geom_density() + geom_vline(aes(xintercept=median(cr_aom1)),color="blue", linetype="dashed", size=1)

  # arrange density plots in a 5x4 grid
  grid.arrange(p1,p2,p3,p4,p5,p6, nrow = 3)
}
# set optimized parameters as median values of the distribution
pa_final=rep(0,paramNum)
for (i in 1:paramNum) {pa_final[i]=median(distr[[i]])}
# compare original and optimized parameters
pa
pa_final
# run the simulation with optimized parameters
optimized = matrix_simu(pa_final)

#compare results
# par(mfrow=c(3, 4)) #set plot window to 3x4 grid
# {
#   plot(test[[1]]$x1[1:120], type="l", col="red",  xlab="month", 
#        ylim=c(min(min(test[[1]]$x1[1:120]),min(dat$cLeaf[1:120])),max(max(test[[1]]$x1[1:120]),max(dat$cLeaf[1:120]))), 
#        xaxp  = c(0, 120, 10), ylab="cLeaf", main="Leaf Pool modelled vs CMIP6")
#   lines(dat$cLeaf[1:120], col="blue")
#   lines(optimized[[1]]$x1[1:120], col="green")
#   
#   plot(test[[1]]$x2[1:120], type="l", col="red", xlab="month", 
#        ylim=c(min(min(test[[1]]$x2[1:120]),min(dat$cWood[1:120])),max(max(test[[1]]$x2[1:120]),max(dat$cWood[1:120]))), 
#        xaxp  = c(0, 120, 10),ylab="cWood", main="Wood Pool modelled vs CMIP6")
#   lines(dat$cWood[1:120], col="blue")
#   lines(optimized[[1]]$x2[1:120], col="green")
#   
#   plot(test[[1]]$x3[1:120], type="l", col="red", xlab="month", 
#        ylim=c(min(min(test[[1]]$x3[1:120]),min(dat$cRoot[1:120])),max(max(test[[1]]$x3[1:120]),max(dat$cRoot[1:120]))), 
#        xaxp  = c(0, 120, 10),ylab="cRoot", main="Root Pool modelled vs CMIP6")
#   lines(dat$cRoot[1:120], col="blue")
#   lines(optimized[[1]]$x3[1:120], col="green")
#   
#   plot(test[[1]]$x4[1:120]+test[[1]]$x5[1:120], type="l", col="red", 
#        ylim=c(min(min(test[[1]]$x4[1:120]+test[[1]]$x5[1:120]),min(dat$cLitterAbove[1:120])),max(max(test[[1]]$x4[1:120]+test[[1]]$x5[1:120]),max(dat$cLitterAbove[1:120]))), 
#        xlab="month", xaxp  = c(0, 120, 10),ylab="cLitterAbove", main="Above-ground Litter (leaf+wood) modelled vs CMIP6")
#   lines(dat$cLitterAbove[1:120], col="blue")
#   lines(optimized[[1]]$x4[1:120]+optimized[[1]]$x5[1:120], col="green")
#   
#   plot(test[[1]]$x6[1:120], type="l", col="red", 
#        ylim=c(min(min(test[[1]]$x6[1:120]),min(dat$cLitterBelow[1:120])),max(max(test[[1]]$x6[1:120]),max(dat$cLitterBelow[1:120]))), 
#        xlab="month", xaxp  = c(0, 120, 10),ylab="cLitterBelow", main="Below-ground Litter (root) modelled vs CMIP6")
#   lines(dat$cLitterBelow[1:120], col="blue")
#   lines(optimized[[1]]$x6[1:120], col="green")
#   
#   plot(test[[1]]$x7[1:120], type="l", col="red", 
#        ylim=c(min(min(test[[1]]$x7[1:120]),min(dat$cSoilFast[1:120])),max(max(test[[1]]$x7[1:120]),max(dat$cSoilFast[1:120]))), 
#        xlab="month", xaxp  = c(0, 120, 10),ylab="cSoilFast", main="Fast SOM pool modelled vs CMIP6")
#   lines(dat$cSoilFast[1:120], col="blue")
#   lines(optimized[[1]]$x7[1:120], col="green")
#   
#   plot(test[[1]]$x8[1:120], type="l", col="red",  
#        ylim=c(min(min(test[[1]]$x8[1:120]),min(dat$cSoilMedium[1:120])),max(max(test[[1]]$x8[1:120]),max(dat$cSoilMedium[1:120]))), 
#        xlab="month", xaxp  = c(0, 120, 10),ylab="cSoilMEdium", main="Medium, SOM pool modelled vs CMIP6")
#   lines(dat$cSoilMedium[1:120], col="blue")
#   lines(optimized[[1]]$x8[1:120], col="green")
#   
#   plot(test[[1]]$x9[1:120], type="l", col="red",  
#        ylim=c(min(min(test[[1]]$x9[1:120]),min(dat$cSoilMedium[1:120])),max(max(test[[1]]$x9[1:120]),max(dat$cSoilMedium[1:120]))), 
#        xlab="month", xaxp  = c(0, 120, 10),ylab="cSoilSlow", main="Slow SOM pool modelled vs CMIP6")
#   lines(dat$cSoilSlow[1:120], col="blue")
#   lines(optimized[[1]]$x9[1:120], col="green")
#   
#   plot(test[[2]][1:120], type="l", col="red",  
#        ylim=c(min(min(test[[2]][1:120]),min(dat$rh[1:120])),max(max(test[[2]][1:120]),max(dat$rh[1:120]))), 
#        xlab="month", xaxp  = c(0, 120, 10),ylab="Rh", main="Heterotrophic respiration modelled vs CMIP6")
#   lines(dat$rh[1:120], col="blue")
#   lines(optimized[[2]][1:120], col="green")
#   
#   plot(test[[3]][1:120], type="l", col="red", 
#        ylim=c(min(min(test[[3]][1:120]),min(dat$fVegLitter[1:120])),max(max(test[[3]][1:120]),max(dat$fVegLitter[1:120]))), 
#        xlab="month", xaxp  = c(0, 120, 10),ylab="fVegLitter", main="Litterfall modelled vs CMIP6")
#   lines(dat$fVegLitter[1:120], col="blue")
#   lines(optimized[[3]][1:120], col="green")
#   
#   plot(test[[4]][1:120], type="l", col="red", 
#        ylim=c(min(min(test[[4]][1:120]),min(dat$fLitterSoil[1:120])),max(max(test[[4]][1:120]),max(dat$fLitterSoil[1:120]))), 
#        xlab="month", xaxp  = c(0, 120, 10),ylab="fLitterSoil", main="Humus formation modelled vs CMIP6")
#   lines(dat$fLitterSoil[1:120], col="blue")
#   lines(optimized[[4]][1:120], col="green")
#   
#   plot(2, xlim=c(0,1), ylim=c(0,1), axes = F, main="legend", ylab="")
#   legend(0.1, 0.9, legend=c("CMIP-6 Output", "Modelled Non-optimized", "Modelled Optimized"),
#          col=c("blue", "red", "green"), lty=1, cex=1)
#   
# }
optimized$cveg=optimized[[1]]$x1+optimized[[1]]$x2+optimized[[1]]$x3+optimized[[1]]$x4
optimized$clitter=optimized[[1]]$x5+optimized[[1]]$x6
optimized$csoil=optimized[[1]]$x7+optimized[[1]]$x8+optimized[[1]]$x9+optimized[[1]]$x10+optimized[[1]]$x11+optimized[[1]]$x12

par(mfrow=c(3, 2)) #set plot window to 3x2 grid
{
  plot(test$veg, type="l", col="red",  xlab="month", 
       ylim=c(0,15), ylab="cveg", main=" modelled cveg vs CMIP6")
  lines(dat$cVeg, col="blue")
  lines(optimized$cveg, col="green")
  
  plot(test$litter, type="l", col="red", ylim=c(0,max(test$litter)), 
       xlab="month",ylab="cLitter", main="modelled litter vs CMIP6")
  lines(dat$cLitter, col="blue")
  lines(optimized$clitter, col="green")
  
  plot(test$soil, type="l", col="red", 
       ylim=c(0,20), 
       xlab="month", ylab="cSoil", main="modelled SOC vs CMIP6")
  lines(dat$cSoil, col="blue")
  lines(optimized$csoil, col="green")
  
  plot(test[[2]], type="l", col="red",  
       ylim=c(0,max(test[[2]])), 
       xlab="month", ylab="Rh", main="modelled rh vs CMIP6")
  lines(dat$rh, col="blue")
  lines(optimized[[2]], col="green")
  
  plot(2, xlim=c(0,1), ylim=c(0,1), axes = F, main="legend", ylab="")
  legend(0.1, 0.9, legend=c("CMIP-6 Output", "Modelled Non-optimized", "Modelled Optimized"),
         col=c("blue", "red", "green"), lty=1, cex=1)
  
}
#set plot window back to 1x1
par(mfrow=c(1, 1))
