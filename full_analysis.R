#ANALYSIS OF CDF FILES

#Load all required packages
library(xcms)
#library(CAMERA)
library(RAMClustR)
library(MetaboAnalystR)
library(MSnbase)
library(InterpretMSSpectrum) #needed by ramclustr annotate


#---RAW DATA IMPORT--------------------
# Define paths to your CDF files
#Input:path to folder containing CDF files
#Output: list containing paths to all CDF files in that folder
data_files <- list.files(path = "C://Users//andre//Desktop//SUMMER_COURSE//data", pattern = ".cdf", all.files = TRUE, full.names = TRUE)

# Read raw data using readMSData from MSNBase
raw_data <- readMSData(files = data_files, mode = "onDisk")

#---PREPROCESSING: PEAK PICKING WITH XCMS--------------------
# Perform peak picking using CentWave algorithm
#first set the wavelet parameters (try different sntresh and noise) #tryprefilter (3,50)
cwp <- CentWaveParam(ppm = 25,peakwidth = c(5, 20),snthresh = 100,prefilter = c(3, 50),
  mzCenterFun = "wMean",integrate = 1L,mzdiff = 0.001, fitgauss = FALSE,noise = 1000,
  verboseColumns = FALSE,roiList = list(),firstBaselineCheck = TRUE,roiScales = numeric(),
  extendLengthMSW = FALSE)

#then search for peaks in the input data
xset <- findChromPeaks(raw_data, param = cwp)
#Peak picking is slow, we would like to save the output to disk so we don't have
#to redo it every time. Name it according to settings. If you donr need to change it dont redo
save(xset, file = "C://Users//andre//Desktop//SUMMER_COURSE//data//xset_peak_picking_sn100pre50.RData")
load("C://Users//andre//Desktop//SUMMER_COURSE//data//xset_peak_picking_sn100pre50.RData")


#---PREPROCESSING: PEAK GROUPING WITH XCMS--------------------
# Group peaks across samples to take into account slight variations in retention
#times for the same species due to different experimental conditions
# Define sample groups
sample_groups <- factor(c("low glucose", "low glucose", "low glucose",
				 "high glucose", "high glucose", "high glucose"))
                         #"alkanes", "blank"))

# Assign these groups to the raw_data object
pData(raw_data)$sample_group <- sample_groups
pdp <- PeakDensityParam(sampleGroups = pData(raw_data)$sample_group, bw = 5, minFraction = 0.5)
xset_grpd <- groupChromPeaks(xset, param = pdp)

#---PREPROCESSING: RETENTION TIME CORRECTION WITH XCMS--------------------
#Do retention time correction
center_sample_index <- 3 #could use any number between 1 and 6 here
obi_params <- ObiwarpParam(binSize = 0.6, response = 1, distFun = "cor", centerSample = center_sample_index)
xset_grpd_rtcorr <- adjustRtime(xset_grpd, param = obi_params)
save(xset_grpd_rtcorr, file = "C://Users//andre//Desktop//SUMMER_COURSE//data//xset_grouped_rt_corrected_SN100pre50.RData")
load("C://Users//andre//Desktop//SUMMER_COURSE//data//xset_grouped_rt_corrected_sn100pre50.RData")

#---PREPROCESSING: REDO PEAK GROUPING HERE TO FIX ABOVE MENTIONED ISSUE. --------------------
#SEE https://training.galaxyproject.org/training-material/topics/metabolomics/tutorials/gc_ms_with_xcms/tutorial.html#peak-detection-using-xcms
pData(raw_data)$sample_group <- sample_groups
pdp <- PeakDensityParam(sampleGroups = pData(raw_data)$sample_group, bw = 5, minFraction = 0.5)
xset_grpd_rtcorr_grpd <- groupChromPeaks(xset_grpd_rtcorr, param = pdp)



#---PREPROCESSING: PEAK FILLING WITH XCMS--------------------
#Do peak filling to account for missing peaks
xset_grpd_rtcorr_grpd_fild <- fillChromPeaks(xset_grpd_rtcorr_grpd)

save(xset_grpd_rtcorr_grpd_fild, file = "C://Users//andre//Desktop//SUMMER_COURSE//data//xset_filled_SN100pre50.RData")
load("C://Users//andre//Desktop//SUMMER_COURSE//data//xset_filled_SN100pre50.RData")


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


rc_obj <- ramclustR(xset_grpd_rtcorr_grpd_fild, ExpDes = ExpDes) #When using peak filling

#rc_obj <- ramclustR(xset_grpd_rtcorr_grpd, ExpDes = ExpDes) #When not using peak filling 

#line below is to avoid errors in write.msp. Seems to not be needed now, but if it says something with msint then use it. 
#rc_obj$msmsint <- rc_obj$msint

#---RETENTION INDEX CALCULATION WITH RAMCLUSTR--------------------
rc_obj_1 <- rc.calibrate.ri(ramclustObj = rc_obj, calibrant.data = "C://Users//andre//Desktop//SUMMER_COURSE//data_try_2//alkane_data.csv", poly.order = 3) 

#---PARTIAL ANNOTATION WITH RAMCLUSTR--------------------
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


#-------------ONLY needed if MSFinder
save(rc_obj_2, file = "C://Users//andre//Desktop//SUMMER_COURSE//data//ramclust_obj_before_msimport.RData")
load("C://Users//andre//Desktop//SUMMER_COURSE//data//ramclust_obj_before_msimport.RData")

#write.msp(rc_obj_2,"C://Users//andre//Desktop//SUMMER_COURSE//data//sn5noise10.msp") #not needed


#---RUN MSFINDER NOW (IT'S A DIFFERENT PROGRAM) ALTERNATIVELY MATCHMS IN PYTHON--------------------
#:-)

#---IMPORT FORMULAS FROM MSFINDER BACK INTO RAMCLUSTR--------------------
rc_obj_3 <- import.msfinder.formulas(ramclustObj = rc_obj_2) #, msp.dir = "C://Users//andre//Documents//spectra//msfinder_result.msp"
rc_obj_3$msfinder.formula.details
rc_obj_3 <- import.msfinder.mssearch(ramclustObj = rc_obj_2)
rc_obj_3_b <- import.msfinder.structures(ramclustObj = rc_obj_2) #don't know where youd need this


#Save ramclust_obj at this stage
save(rc_obj_3_b, file = "C://Users//andre//Desktop//SUMMER_COURSE//data//ramclust_obj_after_msimport.RData")
load("C://Users//andre//Desktop//SUMMER_COURSE//data//ramclust_obj_after_msimport.RData")

#---USE THE IMPORTED FORMULAS TO ANNOTATE THE RAMCLUSTR OBJECT (PEAKS AFTER DECONVOLUTION/CLUSTERING)--------------------
#Define databases to avoid havoc later
dbs <- c("MoNA", "MINE", "COCONUT", "HMDB", "BLEXP", "UNPD", "MinesNPA") #add more database names here
search.dbs <- dbs #We must do this to prevent the error 'object search.dbs not found'

rc_obj_4 <- annotate(ramclustObj = rc_obj_3_b)

annotation.summary(ramclustObj = rc_obj_4, outfile = "C://Users//andre//Desktop//SUMMER_COURSE//data//annotation_ramclustr.csv")

#---EXPORT ANNOTATED PEAK LIST TO CSV (MEANS COMMA-SEPARATED VALUE, CAN BE READ IN EXCEL)--------------------
exportDataset(
ramclustObj = rc_obj_4,
which.data = "SpecAbund",
label.by = "ann",
appendFactors = TRUE
)
#This spec abund file has names from MSFinder, it may be easier to use specAbund from rc_obj_2 above if you want to manually annotate matchms names to compounds 

#---EXPORT ANNOTATED PEAK LIST TO MSP (NIST FORMAT, CAN BE READ IN WORDPAD AND IS COMPATIBLE WITH METABOANALYSTR)--------------------
write.msp(ramclustObj = rc_obj_4, one.file = FALSE)

#-----------If you used matchms continue here----------

#-------BEGIN METABOANALYSTR ANALYSIS----------------------

#----Statistical--Analysis----
#Find out which metabolites are significantly different between the LG and HG, using T-test
#Make sure your SpecAbund excel file has group names, and names of the compounds. This is all that is used here. 
#following https://www.metaboanalyst.ca/resources/vignettes/Statistical_Analysis_Module.html
# Clean global environment
rm(list = ls())

library(MetaboAnalystR)
#mSet<-InitDataObjects("conc", "stat", FALSE);
mSet<-InitDataObjects("pktable", "pathqea", FALSE); 

# Create mSetObj #try ('conc', 'msetora', FALSE) otherwise cancel it
#mSet<-InitDataObjects("conc", "msetora", FALSE)


mSet<-Read.TextData(mSet, "C://Users//andre//Desktop//SUMMER_COURSE//data//forMetaboanalystR//SpecAbund_sn100pre50.csv", "rowu", "disc"); #Make sure to put the SpecAbund file you want to use in the right folder. 
mSet<-SanityCheckData(mSet);

mSet<-ReplaceMin(mSet);
mSet<-PreparePrenormData(mSet);
mSet<-Normalization(mSet, "NULL", "LogNorm", "MeanCenter", "S10T0", ratio=FALSE, ratioNum=20);


mSet<-PlotNormSummary(mSet, "norm_0_matchms_SN100pre50", format ="png", dpi=72, width=NA); #change names according to settings
mSet<-PlotSampleNormSummary(mSet, "snorm_0_matchms_SN100pre50", format = "png", dpi=72, width=NA);

# Perform T-test (parametric)
mSet<-Ttests.Anal(mSet, nonpar=F, threshp=0.05, paired=FALSE, equal.var=TRUE, "fdr", TRUE)


# Plot of the T-test results
mSet<-PlotTT(mSet, imgName = "matchms_tt_0_SN100pre50", format = "png", dpi = 72, width=NA)

mSet #To see your data 
#----Pathway--Enrichment----

#Load the significantly different metabolites from above and compare to metabolic pathway library
#Following : https://www.metaboanalyst.ca/resources/vignettes/Enrichment_Analysis.html

#List your significant metabolites
#
# Create vector consisting of compounds for enrichment analysis. Change to the compounds you have found with current setting. 

#tmp.vec <- c('glycerol 3-phosphate', 'talose', 'Pyrophosphate', 'L-Pipecolic acid', 'thymol', '5-methoxyindole-3-acetic acid', 'MEPERIDINE', '4-ETHOXYCARBONYL-5-METHYL-6-PHENYL-2-(PARA-TOLYL)PERHYDROPYRROLO(3,4-C)PYRROLE-1,3-DIONE (3A,4-TRANS-6,6A-CIS)', '2-BENZYL-8-METHYL-1,2,3,4-TETRAHYDRO-GAMMA-CARBOLINE')
#tmp.vec <- c('D-glucose', 'Pyrophosphoric acid', '5-methoxyindole-3-acetic acid', '1-Methyladenosine')
#tmp.vec <- c('MEPERIDINE', '4-(2-HYDROXY-5-NITROPHENYL)PYRIMIDINE', 'D-Glucose', 'Pyrophosphoric acid', 'DL-Pipecolic acid', '5-Methoxyindole-3-acetic acid', '3-Phenylpropionic acid', 'L-Serine')
#tmp.vec <- c('myo-Inositol', 'D-glucose', 'Pyrophosphoric acid', '5-methoxyindole-3-acetic acid' ) #S/N 500
#tmp.vec <- c('DIHYDROISOJASMONE','BENZO(B)FLUORENE', 'ISOPRENE', 'Oxidized-adrenal-ferredoxin', '3-Methylbutyl phenylacetate', 'L-Pipecolic acid', 'L-Lysine', '4-(4-HYDROXYPHENYL)BENZOIC ACID','ISOPROPYL CELLOSOLVE', 'rac-Glycerol 3-phosphoate', '(Z)-5-METHOXY-3,5-DIMETHYL-2-HEXENYLTRIMETHYLSILANE', 'talose', 
#'phytosphingosine', '2-Methylpropyl acetate', 'L-Pipecolic acid', 'CIS-9-OCTADECENOIC ACID BUTYL ESTER', 'BUTYLMETHYLSILANE (METHYL-D3)', '2-(3-HYDROXYBUTOXY)-N-(2-(DIMETHYLAMINO)ETHYL)-4-QUINOLINECARBOXAMIDE', 'L-Tyrosine', 
#'METHYL TETRADECANOATE', '4-(2-HYDROXY-5-NITROPHENYL)PYRIMIDINE', 'ETHYLHEXADECYLDIMETHYLAMMONIUM BROMIDE', 'Phensuximide', '3-Methylbutyl phenylacetate', 
#'1,2-DIMETHOXYETHANE', 'ALDRITHIOL-2', '2-hydroxypropanoic acid', '1,6-Anhydro-beta-D-Glucose', '2-OCTENOIC ACID TRIMETHYLSILYL ESTER')
#tmp.vec <- c('Oxidized-adrenal-ferredoxin', 'Glycerol 3-phosphoate', 'talose', 'L-Pipecolic acid', 'Fructose 1,6-bisphosphate', '1,6-Anhydro-beta-D-Glucose', 'beta-Alanine', 'Phensuximide', 'L-Lysine', '(2E)-Decenoyl-ACP')
#tmp.vec <- c('3-Methylbutyl phenylacetate', 'Oxidized-adrenal-ferredoxin', 'Glycerol trihexanoate', '4.4-DDE', 'Pyrophosphate','N-butyl Oleate', 'Glycerol 3-phosphoate', 'L-Pipecolic acid', 'talose', '1,6-Anhydro-beta-D-Glucose', 'CYCLOPENTANE' )
tmp.vec <- c('Oxidized-adrenal-ferredoxin', 'Glycerol 3-phosphoate', 'talose', 'L-Pipecolic acid', 'D-Ribose', '1,6-Anhydro-beta-D-Glucose', '2-hydroxypropanoic acid')


#Set up mSetObj with the list of compounds
mSet<-Setup.MapData(mSet, tmp.vec);

# Cross reference list of compounds against libraries (hmdb, pubchem, chebi, kegg, metlin)
mSet<-CrossReferencing(mSet, "name"); 

mSet$name.map #To view the names found 

# Create the mapping results table
mSet<-CreateMappingResultTable(mSet)

# Input the name of the compound without any matches, only relevant if names not identified
#mSet<-PerformDetailMatch(mSet,'');

# Create list of candidates to replace the compound, only if some were not identified 
#mSet <- GetCandidateList(mSet);

# Set the metabolite filter
mSet<-SetMetabolomeFilter(mSet, F);

# Select metabolite set library, refer to 
mSet<-SetCurrentMsetLib(mSet, "smpdb_pathway", 0);


# Calculate hypergeometric score, results table generated in your working directory

mSet<-CalculateHyperScore(mSet)

# Plot the ORA, bar-graph
mSet<-PlotORA(mSet, "ora_0_sn1000noise10000", "bar", "png", 72, width=NA)



