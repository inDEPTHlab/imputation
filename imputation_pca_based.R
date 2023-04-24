# Load packages
library(psych) #For descriptives
library(foreign) #To read spss files
library(missMDA) #For PCA based imputation

# Load ADHD symptoms
adhd.data <- read.spss("data/CHILDCBCL_6_incl_Tscores_20201111.sav", to.data.frame=T)
adhd.data <- adhd.data[c("IDC","sum_att_5")]

# Load general data
general.data <- read.spss("data/CHILD-ALLGENERALDATA_24102022.sav", to.data.frame=T)

# Merge ADHD and general data
adhd_general.data <- merge(adhd.data, general.data, by = "IDC", all = T)

# Load maternal smoking
smoking.data <- read.spss("data/130815_ GEDRAGSGROEP Maternal smoking pregnancy.sav", to.data.frame=T)

# Merge maternal smoking
adhd_general_smoking.data <- merge(adhd_general.data, smoking.data, by = "IDM", all = T)

# Load combined cell count data
load("/home/599032/methylationR2/Datasets/celltypeCounts_birth_Combined.RData")

# Merge cell count data
adhd_general_smoking_counts.data <- merge(adhd_general_smoking.data, wbc_birth_Combined, by = "IDC", all = T)

# Merge sample plate and ID
load("~/alspac/GenR_METH_JointedQC_with_ALSPAC/Data/jointQC_samplesfile_180529.RData")
# Filter for birth
samples_birth.data <- samples[samples$Period == "Birth", ]

# Reduce data to only children with DNA methylation data
perinatal.data <- merge(adhd_general_smoking_counts.data, samples_birth.data, by = "IDC")

# Reduce data to only relevant variables and complete ADHD 
perinatal.data <- perinatal.data[!is.na(perinatal.data$sum_att_5),c("IDC", "Sample_ID", "sum_att_5", "GENDER", "GESTBIR", "msmoke", "EDUCM", "CD8T", "NK", "CD4T", "Bcell", "Gran", "Mono", "nRBC", "WEIGHT", "BMI_0", "BMI_1")]

# Convert education to continuous
perinatal.data$EDUCM <- as.numeric(perinatal.data$EDUCM)

# Check data
str(perinatal.data)
describe(perinatal.data)

# Variables to be included in imputation
variable_names <- c("sum_att_5", "GENDER", "GESTBIR", "msmoke", "EDUCM", "CD8T", "NK", "CD4T", "Bcell", "Gran", "Mono", "nRBC", "WEIGHT", "BMI_0", "BMI_1")

# Maximum number of PCs, which can be fitted (number of variables - 1)
ncp.max <- length(variable_names)-1

# Use Kfold cross-validation to determine optimum number of componencts
set.seed(20230423)
nb <- estim_ncpFAMD(perinatal.data[,variable_names], ncp.max = ncp.max)
nb$ncp

# Save elbow plot
png(file="figures/elbow_plot.png")
plot(0:ncp.max, nb$criterion, xlab = "nb dim", ylab = "MSEP")
dev.off()

# Perform PCA based imputation using the optimum number of components determined by CV
set.seed(20230423)
perinatal.imputePCA <- imputeFAMD(perinatal.data[,variable_names], ncp = nb$ncp)

# Recover the IDs, which were not included in imputation
perinatal_imputed.data <- as.data.frame(perinatal.imputePCA$completeObs)
perinatal_imputed.data <- cbind(perinatal.data[c("IDC", "Sample_ID")],perinatal_imputed.data)

# Compare associations with and without imputation
summary(lm(scale(sum_att_5) ~ scale(BMI_0) + GENDER + CD8T + NK + CD4T + Bcell + Gran + Mono + nRBC + EDUCM, data = perinatal.data))
summary(lm(scale(sum_att_5) ~ scale(BMI_0) + GENDER + CD8T + NK + CD4T + Bcell + Gran + Mono + nRBC + EDUCM, data = perinatal_imputed.data))

