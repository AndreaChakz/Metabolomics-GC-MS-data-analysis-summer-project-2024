#ANALYSIS OF ALNKANES CDF FILES

#OBTAIN RETENTION INDICES FOR THE ALKANE

#Load all required packages
library(xcms)
library(RAMClustR)
#library(RIAssigner)
library(MSnbase)

# Perform peak detection using CentWave algorithm
# Create an MSnbase OnDiskMSnExp object from the NetCDF file
data_files <- list.files(path = "C://Users//andre//Desktop//SUMMER_COURSE//data_try_2//for_rt", pattern = ".cdf", all.files = TRUE, full.names = TRUE)
raw_data <- readMSData(files = data_files, mode = "onDisk")

# Perform peak picking using CentWave (example parameters)
#snthresh = 15 #or 18 gives 7.9 as first peak
cwp <- CentWaveParam(
  ppm = 25,
  peakwidth = c(5, 20),
  snthresh = 18,
  prefilter = c(3, 100),
  mzCenterFun = "wMean",
  integrate = 1L,
  mzdiff = 12, #previously -0.001,
  fitgauss = FALSE,
  noise = 0,
  verboseColumns = FALSE,
  roiList = list(),
  firstBaselineCheck = TRUE,
  roiScales = numeric(),
  extendLengthMSW = FALSE
)

xset <- findChromPeaks(raw_data, param = cwp)
save(xset, file = "C://Users//andre//Desktop//SUMMER_COURSE//data_try_2//xset_peak_picking_incl_blank.RData")
#load("C://Users//andre//Desktop//SUMMER_COURSE//data_try_2//xset_peak_picking_incl_blank.RData")


#---PREPROCESSING: PEAK GROUPING WITH XCMS--------------------
# Group peaks across samples to take into account slight variations in retention
#times for the same species due to different experimental conditions
# Define sample groups
sample_groups <- factor(c("alkanes", "alkanes2"))

# Assign these groups to the raw_data object
pData(raw_data)$sample_group <- sample_groups
pdp <- PeakDensityParam(sampleGroups = pData(raw_data)$sample_group, bw = 5, minFraction = 0.5)
xset_grouped <- groupChromPeaks(xset, param = pdp)


#---PREPROCESSING: PEAK FILLING WITH XCMS--------------------
#Do peak filling to account for missing peaks
xset_filled <- fillChromPeaks(xset_grouped)


#---DECONVOLUTION WITH RAMCLUSTR--------------------
ramclust_obj <- ramclustR(xset_filled)

alkane_retention_times <- ramclust_obj$clrt
alkane_retention_times <- sort(alkane_retention_times)
alkane_retention_times 


#-----CALCULATE RI FOR ALKANES ----
# Calculate retention indices for alkanes
# Assign known retention indices to each alkane (e.g., C10 = 1000, C11 = 1100, etc.)
num_alkanes <- length(alkane_retention_times)
alkane_indices <- 1100 + 100 * (1:num_alkanes - 1)

# Create a data frame for retention times and retention indices
alkane_data <- data.frame(
  #Alkane = paste0("C", 10:(10 + num_alkanes - 1)),
  rt = alkane_retention_times,
  ri = alkane_indices
)


# Print the alkane data with retention times and indices
print(alkane_data)

#Save alkane retention times and indices
save(alkane_data, file = "C://Users//andre//Desktop//SUMMER_COURSE//data_try_2//alkane_data.RData")
write.csv(alkane_data,"C://Users//andre//Desktop//SUMMER_COURSE//data_try_2//alkane_data.csv", row.names = FALSE)

