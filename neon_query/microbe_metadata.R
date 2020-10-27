# Based on Dr. Lee Stanish's code
# https://github.com/lstanish/SummitMicrobes/blob/master/getNEON_soil_microbes_data.R


library(neonUtilities)
library(tidyverse)


##### NOTES #######
# This will need to filter by the sites that we are interested in *BEFORE* 
# we end up pulling the microbiome data 
# qPCR Total Abundances of total archaea, bacteria, and fungi : DP1.10109.001
# Microbe biomass: DP1.10104.001
# Soil Microbe Marker Gene Sequences: DP1.10108.001
# Soil Microbe Metagenome Sequences: DP1.10107.001
# Soil Physical Properties: DP1.10086.001

# Set API key
NEON_TOKEN <- Sys.getenv(x = "NEON_TOKEN")

# Fetch soil microbe marker gene sequence data
marker_genes <- loadByProduct(startdate = "2013-06", enddate = "2019-09",
                            dpID = 'DP1.10108.001', package = 'expanded', 
                            token = NEON_TOKEN, check.size = FALSE, nCores = 15)

#==============================================================================
# Filter Data by DNAsampleID's that match the 16S primers of interest
#==============================================================================

#Make a list to catch output of the multi-step 
    # filtering (cut last 3 tables from original list)
processed_marker_genes <- vector(mode = "list", length = 6)

#filter and store 16S Pcr Amplification data frame
processed_marker_genes[[4]]<- marker_genes$mmg_soilPcrAmplification_16S[marker_genes$mmg_soilPcrAmplification_16S$forwardPrimer == "CCTACGGGNBGCASCAG",]

amp_meta <- c(1, 2, 3, 5, 6)
for(i in amp_meta){
  processed_marker_genes[[i]]<- marker_genes[[i]][marker_genes[[i]]$dnaSampleID
      %in% unique(processed_marker_genes[[4]]$dnaSampleID),]
}
names(processed_marker_genes) <- names(marker_genes)[1:6]

# get metadata for sequencing only done at Battelle 
neon_marker_genes <- vector(mode = "list", length = 6)
for(i in 1:length(neon_marker_genes)){
neon_marker_genes[[i]] <- processed_marker_genes[[i]][processed_marker_genes[[i]]$laboratoryName != "Argonne National Laboratory",]
}
names(neon_marker_genes) <- names(processed_marker_genes)

# filter down to DNA samples with ITS sequences
its_filtered_markers <- vector(mode = "list", length = 6)
for(i in 1:length(its_filtered_markers)){
  its_filtered_markers[[i]] <- neon_marker_genes[[i]][neon_marker_genes[[i]]$dnaSampleID %in% neon_marker_genes$mmg_soilMarkerGeneSequencing_ITS$dnaSampleID,]
}
names(its_filtered_markers) <- names(neon_marker_genes)

#make output directories
if(dir.exists("~/fastq")){
  print("Warning: directory ~/fastq/ and subdirectories already exist!")
}else{
system("mkdir ~/fastq") #fastq folder in home of rstudio user
system("mkdir ~/fastq/its") #its fastq directory
system("mkdir ~/fastq/16s") #16s fastq directory
system("mkdir ~/fastq/18s") #18s fastq directory
}
#write out tables from ITS filtered data
for(i in 1:length(its_filtered_markers)){
  write.csv(its_filtered_markers[[i]],
            file = paste0('~/fastq/', names(its_filtered_markers)[[i]],'.csv'),
            row.names = FALSE)
}
#write out variables file
write.csv(marker_genes$variables_10108, file = "~/fastq/variables.csv", row.names=FALSE)

#download fastq files
zipsByURI(filepath = "/home/rstudio/fastq", savepath = "/home/rstudio/fastq",
          unzip = FALSE, check.size = FALSE, saveZippedFiles = TRUE)

#move files into sub folders
#its mv #1
system('mv *ITS*R[0-9].fastq.tar.gz its')
#16s mv #1
system('mv *16S*R[0-9].fastq.tar.gz 16s')
#its mv #2
system('mv *ITS_R[0-9]_fastq.tar.gz its')
#16s mv #2
system('mv *16S_R[0-9]_fastq.tar.gz 16s')
#its mv #3
system('mv *ITS_R[0-9].fastq.tar.gz its')
#16s mv #3
system('mv *16S_R[0-9].fastq.tar.gz 16s')

#check to see if all the files have matching R1 & R2 files

#16s sanity check
length(list.files(path = "~/fastq/16s/", pattern = "BMI*", full.names = FALSE)) #670
length(list.files(path = "~/fastq/16s/", pattern = "*_R1.*.gz", full.names = FALSE)) #335
length(list.files(path = "~/fastq/16s/", pattern = "*_R2.*.gz", full.names = FALSE)) #335

#its sanity check
length(list.files(path = "~/fastq/its/", pattern = "BMI*", full.names = FALSE)) #1373
length(list.files(path = "~/fastq/its/", pattern = "*_R1.*.gz", full.names = FALSE)) #687
length(list.files(path = "~/fastq/its/", pattern = "*_R2.*.gz", full.names = FALSE)) #686
its_r1<-list.files(path = "~/fastq/its/", pattern = "*_R1.*.gz", full.names = FALSE)
its_r2<-list.files(path = "~/fastq/its/", pattern = "*_R2.*.gz", full.names = FALSE)

its_r1_test <- gsub("(_R).*", replacement = "\\1", x = its_r1)
its_r2_test <- gsub("(_R).*", replacement = "\\1", x = its_r2)
setdiff(its_r1_test,its_r2_test)
#[1] "BMI_B69RN_ITS_R1" does not have a mattching R2 file
#executed: rm BMI_B69RN_ITS_R1_fastq.tar.gz in ~/fastq/its/

#==================================================================================
# Get other microbial and soil metadata
#==================================================================================

#get microbial biomass metadata: DP1.10104.001; lipid analysis
microbial_biomass <- loadByProduct(startdate = "2013-06", enddate = "2019-09",
                             dpID = 'DP1.10104.001', package = 'expanded', 
                         token = NEON_TOKEN, check.size = FALSE, nCores = 15)

if(dir.exists("~/neon_amplicon/biomass/")){
  print("Warning: microbial biomass directory already exists!")
}else{
system('mkdir ~/neon_amplicon/biomass')
}
for(i in 1:length(microbial_biomass)){
  write.csv(microbial_biomass[[i]],
            file = paste0('~/neon_amplicon/biomass/', names(microbial_biomass)[[i]],'.csv'),
            row.names = FALSE)
}


#get relative abundances of Archaea, Fungi, and Bacteria: DP1.10109.001; qPCR
qpcr_abundances <- loadByProduct(startdate = "2013-06", enddate = "2019-09",
                      dpID = 'DP1.10109.001', package = 'expanded', 
                      token = NEON_TOKEN, check.size = FALSE, nCores = 15)


#get soil physical properties DP1.10086.001; biogeochemical measurements
soil_properties <- loadByProduct(startdate = "2013-06", enddate = "2019-09",
                                 dpID = 'DP1.10086.001', package = 'expanded', 
                                 token = NEON_TOKEN, check.size = FALSE, nCores = 15)





################################################################################
#                          _ _               
#                         | | |              
#      ___  __ _ _ __   __| | |__   _____  __
#    / __|/ _` | '_ \ / _` | '_ \ / _ \ \/ /
#    \__ \ (_| | | | | (_| | |_) | (_) >  < 
#    |___/\__,_|_| |_|\__,_|_.__/ \___/_/\_\
#                                        
################################################################################
#==================================================================================
# 16S rRNA gene PCR primer diagnostics
#==================================================================================

#get unique forward primers
unique(unique(marker_genes$mmg_soilPcrAmplification_16S$forwardPrimer))
# "GTGYCAGCMGCCGCGGTAA" is the modified earth microbiome project forward primer Parada, et al. 2016

#make counts by 16S rRNA gene primer set, label them, barplot to see distribution of data 
length(marker_genes$mmg_soilPcrAmplification_16S$forwardPrimer)
#count each unique primer in the dataset
forward_primer_dist <- as.data.frame(table(as.vector(marker_genes$mmg_soilPcrAmplification_16S$forwardPrimer)))
primer_set <- c("Pro341Fi","Mislabeled_ITS", "New_EMP")

#bind forward primer counts into a new dataframe to plot
f_prime <- cbind(forward_primer_dist, primer_set)

#ggplot call for primer data
p<-ggplot(data=f_prime, aes(x=primer_set, y=Freq)) +
  geom_bar(stat = "identity", fill="steelblue")+
  theme_minimal()+theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab(NULL)+labs(title = "Forward Primer Distribution")
p

#filter 16S data for CCTACGGGNBGCASCAG Pro341Fi, which targets the V3-V4 regions
#of the 16S rRNA gene
# ITS data only has one primer set so the rest of the data in the list should be
# filtered by the DNA sample ID's left in the 16S dataframe

#==============================================================================
#   Code Debugging Section
#==============================================================================

# debugging function to count data elements
n2tab_count <-function(n){
  df <- as.data.frame(table(as.vector(n)))
  return(df)}

# count each DNA sample debug
dna_sample_debug <- vector(mode = "list", length = length(its_filtered_markers))

for(i in 1:length(its_filtered_markers)){
  dna_sample_debug[[i]] <- n2tab_count(its_filtered_markers[[i]]$dnaSampleID)
}
names(dna_sample_debug) <- names(its_filtered_markers)

#sanity check for duplicates

# mmg_soilDnaExtraction duplicates
length(which(dna_sample_debug[[1]]$Freq > 1)) # 3 
length(unique(its_filtered_markers[[1]]$collectDate))
# mmg_soilMarkerGeneSequencing_16S duplicates
length(which(dna_sample_debug[[2]]$Freq > 1)) # 54 

write.table(x = its_filtered_markers[[1]], 
            file = "~/neon_amplicon/amplicon_sites_ver1.txt", sep = "\t")

# mmg_soilMarkerGeneSequencing_ITS
length(which(dna_sample_debug[[3]]$Freq > 1)) # 2

# mmg_soilPcrAmplification_16S
length(which(dna_sample_debug[[4]]$Freq > 1)) # 385

# mmg_soilPcrAmplification_ITS
length(which(dna_sample_debug[[5]]$Freq > 1)) # 794

# mmg_soilRawDataFiles - less meaningful due to paired end reads
length(which(dna_sample_debug[[6]]$Freq > 2)) # 6147


#===============================================================================
# Lee Stanish's Code Base Starts Here
#===============================================================================
## Grab soil microbe data ##
L1mic <- loadByProduct(startdate = "2016-09", enddate = "2017-01", dpID = 'DP1.10108.001', 
                       package = 'expanded', check.size = FALSE)
L1mic.dna <- L1mic$mmg_soilDnaExtraction   # read in soilDnaExtraction L1 data

length(grep("marker gene|marker gene and metagenomics",
            processed_marker_genes$mmg_soilDnaExtraction$sequenceAnalysisType))
unique(processed_marker_genes$mmg_soilDnaExtraction$sequenceAnalysisType)

length(unique(processed_marker_genes$mmg_soilDnaExtraction$collectDate))
length(unique(processed_marker_genes$mmg_soilMarkerGeneSequencing_16S$collectDate))
length(unique(processed_marker_genes$mmg_soilMarkerGeneSequencing_ITS$collectDate))
# need to remove Argonne sequencing
# need to filter again by ITS collection date and plotID's



marker_genes_dna <- marker_genes$mmg_soilDnaExtraction   # read in soilDnaExtraction L1 data



#grep filter sequence analysis by all types
L1mic.dna <- marker_genes_dna[grep("marker gene|marker gene and metagenomics", marker_genes_dna$sequenceAnalysisType),]
#make dnaSampleID upper case
L1mic.dna$dnaSampleID <- toupper(L1mic.dna$dnaSampleID)

# 16S sequencing metadata
L1mmg16S <- L1mic$mmg_soilMarkerGeneSequencing_16S   # read in marker gene sequencing 16S L1 data
L1mmg16S$dnaSampleID <- toupper(L1mmg16S$dnaSampleID)

# ITS sequencing metadata
L1mmgITS <- L1mic$mmg_soilMarkerGeneSequencing_ITS   # read in marker gene sequencing ITS L1 data
L1mmgITS$dnaSampleID <- toupper(L1mmgITS$dnaSampleID)

# 16S rawDataFiles metadata - this contains URL links to sequence data
L1mmgRaw <- L1mic$mmg_soilRawDataFiles   # read in soilDnaExtraction L1 data
L1mmgRaw16S <- L1mmgRaw[grep('16S', L1mmgRaw$rawDataFileName), ]
L1mmgRaw16S <- L1mmgRaw16S[!duplicated(L1mmgRaw16S$dnaSampleID), ]

# ITS rawDataFiles metadata - this contains URL links to sequence data
L1mmgRawITS <- L1mmgRaw[grep('ITS', L1mmgRaw$rawDataFileName), ]
L1mmgRawITS <- L1mmgRawITS[!duplicated(L1mmgRawITS$dnaSampleID), ]

# variables file - needed to use zipsByURI
varFile <- L1mic$variables

# export data
targetGene <- '16S'  # change to ITS if you want the ITS data instead

if(!dir.exists(paste0(outDir, 'mmg/')) ) {
  dir.create(paste0(outDir, 'mmg/'))
}

write.csv(L1sls.scc, paste0(outDir, 'mmg/', "soilFieldData.csv"), row.names=FALSE)
write.csv(L1.sls.bgc, paste0(outDir, 'mmg/', "soilBGCData.csv"), row.names=FALSE)
write.csv(L1.sls.sm, paste0(outDir, 'mmg/', "soilMoistureData.csv"), row.names=FALSE)
write.csv(L1.sls.ph, paste0(outDir, 'mmg/', "soilpHData.csv"), row.names=FALSE)
write.csv(L1.sls.mg, paste0(outDir,'mmg/', "soilmetagenomicsPoolingData.csv"), row.names=FALSE)
write.csv(L1mic.dna, paste0(outDir,'mmg/', "soilDNAextractionData.csv"), row.names=FALSE)
write.csv(L1mmgITS, paste0(outDir, 'mmg/', "soilITSmetadata.csv"), row.names=FALSE)
if(targetGene=="16S") {
  write.csv(L1mmgRaw16S[1,], paste0(outDir, 'mmg/', "mmg_soilrawDataFiles.csv"), row.names=FALSE)
} else {
  write.csv(L1mmgRawITS, paste0(outDir, 'mmg/', "mmg_soilrawDataFiles.csv"), row.names=FALSE)
}
write.csv(varFile, paste0(outDir,'mmg/', "variables.csv"), row.names=FALSE)


# Download sequence data (lots of storage space needed!)
rawFile <- paste0(outDir, 'mmg/')
zipsByURI(filepath = rawFile, savepath = outDir, unzip = FALSE, saveZippedFiles = TRUE)