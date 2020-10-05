# Based on Dr. Lee Stanish's code
# https://github.com/lstanish/SummitMicrobes/blob/master/getNEON_soil_microbes_data.R


library(neonUtilities)
# Set output directory, relative path to Rproj
outDir <- "~/NEON_metagenomics/test_out/"

##### NOTES #######
# This will need to filter by the sites that we are interested in *BEFORE* 
# we end up pulling the microbiome data 
# qPCR Total Abundances of total archaea, bacteria, and fungi : DP1.10109.001
# Microbe biomass: DP1.10104.001
# Soil Microbe Marker Gene Sequences: DP1.10108.001
# Soil Microbe Metagenome Sequences: DP1.10107.001
# Soil Phyical Properties: DP1.10086.001

# Set API key
NEON_TOKEN <- Sys.getenv(x = "NEON_TOKEN")
#set timezone 
#options(readr.default_locale=readr::locale(tz="USA/Mountain"))

# Fetch metadata for shotgun metagenomes
#sg_met <- loadByProduct(startdate = "2013-06", enddate = "2019-09",
#                       dpID = 'DP1.10107.001', package = 'expanded',
#                       token = NEON_TOKEN, check.size = FALSE, nCores = 3)

# Fetch soil microbe marker gene sequence data
marker_genes <- loadByProduct(startdate = "2013-06", enddate = "2019-09",
                            dpID = 'DP1.10108.001', package = 'expanded', 
                            token = NEON_TOKEN, check.size = FALSE, nCores = 15)

unique(unique(marker_genes$mmg_soilPcrAmplification_16S$forwardPrimer))
# "GTGYCAGCMGCCGCGGTAA" is the modified earth microbiome project forward primer Parada, et al. 2016
length(marker_genes$mmg_soilPcrAmplification_16S$forwardPrimer)
forward_primer_dist <- as.data.frame(table(as.vector(marker_genes$mmg_soilPcrAmplification_16S$forwardPrimer)))
primer_set <- c("Pro341Fi","Mislabeled_ITS", "New_EMP")
f_prime <- cbind(forward_primer_dist, primer_set)

p<-ggplot(data=f_prime, aes(x=primer_set, y=Freq)) +
  geom_bar(stat = "identity", fill="steelblue")+
  theme_minimal()+theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab(NULL)+labs(title = "Forward Primer Distribution")
p

length(marker_genes$mmg_soilPcrAmplification_16S$reversePrimer)
reverse_primer_dist <- as.data.frame(table(as.vector(marker_genes$mmg_soilPcrAmplification_16S$reversePrimer)))


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