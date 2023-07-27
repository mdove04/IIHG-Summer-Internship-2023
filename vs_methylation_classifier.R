## Methylation classifier implementation
library(minfi)
library(data.table)
library(caret)
library(kernlab)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(tidyverse)

## Preprocessing of data #################
##########################################
# process data - data set to be analyzed #
baseDir <- "/Users/makayladove/Downloads/IDAT FILES SCHWANNOMAS"

# load array data
targets <- read.metharray.sheet(baseDir)
RGSet <- read.metharray.exp(targets=targets)
RGSet@annotation = c(array= "IlluminaHumanMethylationEPIC", annotation = "ilm10b2.hg19")
phenodata <- pData(RGSet)

# normalize array data
manifest <- getManifest(RGSet)
MSet.funnorm <- preprocessFunnorm(RGSet)

# filter data based on common SNPs, X/Y chromosome, cross reactive probes, detection p-value
# remove SNPs
mSetSqFlt = dropLociWithSnps(MSet.funnorm)

#remove X/Y chromsome
annEpic = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
keep = !(featureNames(mSetSqFlt) %in% annEpic$Name[annEpic$chr %in% c("chrX","chrY")])
mSetSqFlt = mSetSqFlt[keep,]

# remove probes that do not meet mean detection p-value < 0.05 (SKIPPED)
detP = detectionP(RGSet)
keep = colMeans(detP) < 0.05

detP = detP[match(featureNames(mSetSqFlt),rownames(detP)),]
keep = rowSums(detP < 0.05) == ncol(mSetSqFlt)
mSetSqFlt.2 = mSetSqFlt[keep,]

#get beta values (WITHOUT REMOVED PROBES)
bVals = getBeta(mSetSqFlt)
mVals = getM(mSetSqFlt)

##########################################
# Apply the classifier 
# 0 = Neural Crest Schwannoma
# 1 = Immune Enriched Schwannoma
load("/Users/makayladove/Downloads/svm_linear_classifier.Rds")

testLines <- t(bVals)
rownames(testLines) <- as.character(phenodata[rownames(testLines),"sample"])
p <- predict(svm_Linear, newdata=testLines)

# Output sample ID's and methylation classification 
rownames(testLines)
p


#Create table of classifications
s_classifications <- data.frame(rownames(testLines), p)
s_classifications <- s_classifications %>%
  rename("sample"= "rownames.testLines.") %>%
  rename("classification" = "p")
s_classifications$classification <- as.character(s_classifications$classification)
s_classifications['classification'][s_classifications['classification'] == '0'] <- 'Neural Crest'
s_classifications['classification'][s_classifications['classification'] == '1'] <- 'Immune Enriched'
View(s_classifications)



write.csv(s_classifications, "schwannoma_classifications.csv")





