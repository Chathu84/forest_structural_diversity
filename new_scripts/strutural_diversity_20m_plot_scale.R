
set.seed(123)
library(corrplot)
library(ggplot2)
library(RColorBrewer)

library(GGally)
library(caTools)
#library(randomForest)
library(caret)
library(varSelRF)
library(leaps)
library(MASS)
library(Metrics)
library(mFD)
library(lidR)
library(gstat)
library(terra)
library(raster)


wd = "E:/NEON_Data/Structural_diversity/LiDAR_Data/plot_level_LiDAR/"
setwd(wd)


# create empty data frame
structural_diversity =  data.frame(matrix(ncol = 18, nrow = 0))




# read all file names with .las/.laz extension

file_list = dir(wd, pattern = ".laz", full.names = FALSE, ignore.case = TRUE)
names_to_cloud = substr(file_list, 1, nchar(file_list)-4)

# loop thorugh all plot scale .laz files

for (i in 1:length(names_to_cloud)) #length(names_to_cloud
  
{
  print(i)
  

  split_name = unlist(strsplit(names_to_cloud[i], "[_]"))
  
  options(digits = 11)
  x = as.numeric(split_name[[6]][1])
  y = as.numeric(split_name[[7]][1])
  
  site = split_name[4]
  year_mo = split_name[5]
  
  
  ABBY = lidR::readLAS(file.path(paste0(file_list[i])),filter = "-drop_z_below 0 ")
  

  try (
  data.20m <- 
    clip_rectangle(ABBY, 
                   xleft = (x - 10), ybottom = (y - 10),
                   xright = (x + 10), ytop = (y + 10))
  
  , silent = TRUE)
  
  #when extract lidar data at 20m x 20 x plot level, in some locations there are no points available.
  #I am not sure why. May need to investigate those plots separately. However, in this loop I skipped those plots from further processing.
  
  len_rec = data.20m@header@PHB$`Number of point records`
  print(len_rec)
  
  if (len_rec<= 1 || is.nan(len_rec)== TRUE)
  {
    next
  }
  
  
  # In some cases, there are no ground points availavle to generate the DTM. So, have to skip those plot as well from further processing.
  
  try(
  dtm <- lidR::grid_terrain(las = data.20m,
                            res = 1, # was 0.5 and 0.25 before
                            algorithm = tin())
  , silent = TRUE)
  
  
  # dtm <- grid_terrain(data.20m, 1, kriging(k = 10L)). Use above algorithm to generate dtm.
  
  data.20m <- normalize_height(data.20m, dtm)
  
  # plot(data.200m)
  
  ## chm data includes all heights from 0 to max chm because we are interested in all biomes including shrubs and grass
  
  chm <- grid_canopy(data.20m, res = 1, dsmtin())  
  
  
  pixels_smooth_2 <- 9 #round(((1/chm_res)-1)/2)*1 + 1 # used smoothed parameter as 1. Instead of this I used 9 as the window size
  
  weights <- matrix(1, nrow = pixels_smooth_2, ncol = pixels_smooth_2)
  chm <- terra::focal(x = chm, weights, fun = mean,na.rm=TRUE ) # removed all values of "NA" in the smoothed products
  
  
  # plot(chm)
  ###########################################################################
  ########## All below indices are from https://www.neonscience.org/resources/learning-hub/tutorials/structural-diversity-discrete-return
  ###########################################################################
  
  #MEAN OUTER CANOPY HEIGHT (MOCH)
  #calculate MOCH, the mean CHM height value
  mean.max.canopy.ht <- mean(chm@data@values, na.rm = TRUE) 
  
  #MAX CANOPY HEIGHT
  #calculate HMAX, the maximum CHM height value
  max.canopy.ht <- max(chm@data@values, na.rm=TRUE) 
  
  #RUMPLE
  #calculate rumple, a ratio of outer canopy surface area to 
  #ground surface area (1600 m^2)
  rumple <- rumple_index(chm) 
  
  #TOP RUGOSITY
  #top rugosity, the standard deviation of pixel values in chm and 
  #is a measure of outer canopy roughness
  top.rugosity <- sd(chm@data@values, na.rm = TRUE) 
  
  #DEEP GAPS & DEEP GAP FRACTION
  #number of cells in raster (also area in m2)
  cells <- length(chm@data@values) 
  chm.0 <- chm
  chm.0[is.na(chm.0)] <- 0 #replace NAs with zeros in CHM
  #create variable for the number of deep gaps, 1 m^2 canopy gaps
  zeros <- which(chm.0@data@values == 0) 
  deepgaps <- length(zeros) #number of deep gaps
  #deep gap fraction, the number of deep gaps in the chm relative 
  #to total number of chm pixels
  deepgap.fraction <- deepgaps/cells 
  
  #COVER FRACTION
  #cover fraction, the inverse of deep gap fraction
  cover.fraction <- 1 - deepgap.fraction 
  
  #HEIGHT SD
  #height SD, the standard deviation of height values for all points
  #in the plot point cloud
  vert.sd <- cloud_metrics(data.20m, sd(Z, na.rm = TRUE)) 
  
  #SD of VERTICAL SD of HEIGHT
  #rasterize plot point cloud and calculate the standard deviation 
  #of height values at a resolution of 1 m^2
  sd.1m2 <- grid_metrics(data.20m, sd(Z), 1)
  #standard deviation of the calculated standard deviations 
  #from previous line 
  #This is a measure of internal and external canopy complexity
  sd.sd <- sd(sd.1m2[,3], na.rm = TRUE) 
  
  
  #some of the next few functions won't handle NAs, so we need 
  #to filter these out of a vector of Z points
  Zs <- data.20m@data$Z
  Zs <- Zs[!is.na(Zs)]
  
  #ENTROPY 
  #entropy, quantifies diversity & evenness of point cloud heights 
  #by = 1 partitions point cloud in 1 m tall horizontal slices 
  #ranges from 0-1, with 1 being more evenly distributed points 
  #across the 1 m tall slices 
  entro <- entropy(Zs, by = 1) 
  
  #GAP FRACTION PROFILE 
  #gap fraction profile, assesses the distribution of gaps in the 
  #canopy volume 
  #dz = 1 partitions point cloud in 1 m horizontal slices 
  #z0 is set to a reasonable height based on the age and height of 
  #the study sites 
  gap_frac <- gap_fraction_profile(Zs, dz = 1, z0=3) 
  #defines gap fraction profile as the average gap fraction in each 
  #1 m horizontal slice assessed in the previous line
  GFP.AOP <- mean(gap_frac$gf) 
  
  #VAI
  #leaf area density, assesses leaf area in the canopy volume 
  #k = 0.5 is a standard extinction coefficient for foliage 
  #dz = 1 partitions point cloud in 1 m horizontal slices 
  #z0 is set to the same height as gap fraction profile above
  LADen<-LAD(Zs, dz = 1, k=0.5, z0=3) 
  #vegetation area index, sum of leaf area density values for 
  #all horizontal slices assessed in previous line
  VAI.AOP <- sum(LADen$lad, na.rm=TRUE) 
  
  #VCI
  #vertical complexity index, fixed normalization of entropy 
  #metric calculated above
  #set zmax comofortably above maximum canopy height
  #by = 1 assesses the metric based on 1 m horizontal slices in 
  #the canopy
  VCI.AOP <- VCI(Zs, by = 1, zmax=100) 
  
  
  #OUTPUT CALCULATED METRICS INTO A TABLE
  #creates a dataframe row, out.plot, containing plot descriptors 
  #and calculated metrics
  structural_diversity <- rbind(structural_diversity,data.frame(matrix(c(as.numeric(i), site, year_mo, x, y, mean.max.canopy.ht, max.canopy.ht, 
                                                                         rumple, deepgaps,deepgap.fraction,
                                                                         cover.fraction,top.rugosity, vert.sd, 
                                                                         sd.sd, entro,GFP.AOP, VAI.AOP, VCI.AOP),
                                                                       ncol = 18)) )
    
  
    print(as.character(site))
  

}


colnames(structural_diversity) <- 
  c("raw_id","site", "year_mo", "easting", "northing", "mean.max.canopy.ht.aop",
    "max.canopy.ht.aop", "rumple.aop", "deepgaps.aop",
    "deepgap.fraction.aop","cover.fraction.aop", 
    "top.rugosity.aop", "vert.sd.aop", "sd.sd.aop", 
    "entropy.aop", "GFP.AOP.aop", "VAI.AOP.aop", "VCI.AOP.aop") 

structural_diversity$new_ID = paste(structural_diversity$site,structural_diversity$year_mo,structural_diversity$easting)

# reading Mukund's old data frame to merge other plot level information

old_data = read.csv("E:/NEON_Data/Structural_diversity/structural_metrics_by_plot.csv")

old_data$new_ID = paste(old_data$siteID,old_data$monthyear,old_data$easting)

# find matching rows between two data layers
match_Id = data.frame(match(structural_diversity$new_ID,old_data$new_ID)) # matching fire ID index


pull_Col = data.frame(old_data[match_Id$match.structural_diversity.new_ID..old_data.new_ID.,c(1:13)]) # extract mean unburned biomass

structural_diversity = cbind(pull_Col,structural_diversity)

write.csv(structural_diversity,"E:/NEON_Data/Structural_diversity/struct_diver_20m_plot.csv")
