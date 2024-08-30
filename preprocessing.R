#----Getting started 
#----Installing packages needed 

#Load all required packages
library(xcms)
library(RAMClustR)
library(MetaboAnalystR)
library(MSnbase)
library(InterpretMSSpectrum) #needed by ramclustr annotate


#---RAW DATA IMPORT--------------------
# Define paths to your CDF files
#This function loads all the CDF files you have in your specified folder, so make sure to specify the right folder and only have those you want analysed 
#Input:path to folder containing CDF files
#Output: list containing paths to all CDF files in that folder
data_files <- list.files(path = "C://Users//andre//Desktop//SUMMER_COURSE//data", pattern = ".cdf", all.files = TRUE, full.names = TRUE)

# Read raw data using readMSData from MSnbase

#OnDisk is a way to be more efficient
raw_data <- readMSData(files = data_files, mode = "onDisk")


#---PREPROCESSING: PEAK PICKING WITH XCMS--------------------
# Perform peak picking using CentWave algorithm
#first set the wavelet parameters 
#(here you can change the parameters to see how the results are affected) 
cwp <- CentWaveParam(ppm = 25,peakwidth = c(5, 20),snthresh = 20,prefilter = c(3, 100),
                     mzCenterFun = "wMean",integrate = 1L,mzdiff = 0.001, fitgauss = FALSE,noise = 10000,
                     verboseColumns = FALSE,roiList = list(),firstBaselineCheck = TRUE,roiScales = numeric(),
                     extendLengthMSW = FALSE)

#then search for peaks in the input data
xset <- findChromPeaks(raw_data, param = cwp)
#Peak picking is slow, we would like to save the output to disk so we don't have
#to redo it every time. Name it according to settings so that you can load your different analysis later
save(xset, file = "C://Users//andre//Desktop//SUMMER_COURSE//data//xset_peak_picking.RData")
load("C://Users//andre//Desktop//SUMMER_COURSE//data//xset_peak_picking.RData")


#---PREPROCESSING: PEAK GROUPING WITH XCMS--------------------
# Group peaks across samples to take into account slight variations in retention
#times for the same species due to different experimental conditions
# Define sample groups
sample_groups <- factor(c("low glucose", "low glucose", "low glucose",
                          "high glucose", "high glucose", "high glucose"))

# Assign these groups to the raw_data object
pData(raw_data)$sample_group <- sample_groups
pdp <- PeakDensityParam(sampleGroups = pData(raw_data)$sample_group, bw = 5, minFraction = 0.5)
xset_grpd <- groupChromPeaks(xset, param = pdp)

#---PREPROCESSING: RETENTION TIME CORRECTION WITH XCMS--------------------
#Do retention time correction
center_sample_index <- 3 #could use any number between 1 and 6 here
obi_params <- ObiwarpParam(binSize = 0.6, response = 1, distFun = "cor", centerSample = center_sample_index)
xset_grpd_rtcorr <- adjustRtime(xset_grpd, param = obi_params)
#It may be of interest to try changing the binSize here 
save(xset_grpd_rtcorr, file = "C://Users//andre//Desktop//SUMMER_COURSE//data//xset_grouped_rt_corrected.RData")
load("C://Users//andre//Desktop//SUMMER_COURSE//data//xset_grouped_rt_corrected.RData")

#---PREPROCESSING: REDO PEAK GROUPING HERE TO FIX ABOVE MENTIONED ISSUE. --------------------
#SEE https://training.galaxyproject.org/training-material/topics/metabolomics/tutorials/gc_ms_with_xcms/tutorial.html#peak-detection-using-xcms
pData(raw_data)$sample_group <- sample_groups
pdp <- PeakDensityParam(sampleGroups = pData(raw_data)$sample_group, bw = 5, minFraction = 0.5)
xset_grpd_rtcorr_grpd <- groupChromPeaks(xset_grpd_rtcorr, param = pdp)


#---PREPROCESSING: PEAK FILLING WITH XCMS--------------------
#Do peak filling to account for missing peaks
#1xset_grpd_rtcorr_grpd_fild <- fillChromPeaks(xset_grpd_rtcorr_grpd)

#1save(xset_grpd_rtcorr_grpd_fild, file = "C://Users//andre//Desktop//SUMMER_COURSE//data//xset_filled_sn5noise10000nofill.RData")
#1load("C://Users//andre//Desktop//SUMMER_COURSE//data//xset_filled_sn20noise10000nofill.RData")


#NOTE! We need to do RI calculation before clustering, because otherwise we are trying
#to add more RI than there are clusters
#NOTE 2: challenge with the above is how to add ri to the xcms object, seems a hassle
#Alternative is to use ramclustr to calibrate ri, export peak table from ramclustr to 
#csv, then import that for matchms and later steps

#---DECONVOLUTION WITH RAMCLUSTR--------------------

#first define experiment (run once, save then load later)
#ExpDes <- defineExperiment(csv = TRUE, force.skip = FALSE)
#save(ExpDes, file = "C://Users//andre//Desktop//SUMMER_COURSE//data//ExpDes.RData")
load("C://Users//andre//Desktop//SUMMER_COURSE//data//ExpDes.RData")


#1rc_obj <- ramclustR(xset_grpd_rtcorr_grpd_fild, ExpDes = ExpDes) #When using peak filling

rc_obj <- ramclustR(xset_grpd_rtcorr_grpd, ExpDes = ExpDes) #When not using peak filling 

#line below is to avoid errors in write.msp. Seems to not be needed now, but if it says something with msint then use it. 
#rc_obj$msmsint <- rc_obj$msint

#---RETENTION INDEX CALCULATION WITH RAMCLUSTR--------------------
rc_obj_1 <- rc.calibrate.ri(ramclustObj = rc_obj, calibrant.data = "C://Users//andre//Desktop//SUMMER_COURSE//data_try_2//alkane_data.csv", poly.order = 3) 

#---PARTIAL ANNOTATION WITH RAMCLUSTR--------------------
#Takes a while to do this step 
rc_obj_2 <- do.findmain(ramclustObj = rc_obj_1)

#line below is to avoid errors in write.msp. Seems to not be needed now, but if it says something with msint then use it. 
rc_obj_2$msmsint <- rc_obj_2$msint

write.msp(ramclustObj = rc_obj_2, one.file = TRUE) #Makes the file called Summer which Python uses (for matchms)
#write.msp(ramclustObj = rc_obj_2, one.file = FALSE) #Makes many folders with compounds for MSFinder 

#export output of above (saves to docs/datasets)
exportDataset(
  ramclustObj = rc_obj_2,
  which.data = "SpecAbund", #you cannot change this name, so instead change names of previous ones 
  label.by = "ann",
  appendFactors = TRUE
)

#You can now go to Spyder and open the code there, it will load the "Summer" file created and 
#identify all the compounds that have been found. It will then save the names to a CSV file called "names"
#which you can take and copy the row of names to insert them in the SpecAbund file created here. 
#In the SpecAbund file you must also add a column that says "Groups" and then LG LG LG HG HG HG 
#This is in order for MetaboAnalystR to know which two groups it should compare. 
#Save the SpecAbund file in the folder that MetaboAnalystR loads it from
