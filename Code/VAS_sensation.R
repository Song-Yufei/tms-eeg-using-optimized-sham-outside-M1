# this script perform Wilcoxon signed-rank test on the VAS data

#### install package ####
# List of packages to install
packages_to_install <- c(
  "foreign", "ggplot2", "tidyverse", "psych", "PairedData", 
  "reshape", "dplyr", "ggpubr", "pastecs", "car", "plyr", 
  "rstatix", "coin", "rcompanion", "DescTools", "ggthemes", "reshape2"
)

# Install and load the packages
install.packages(packages_to_install, dependencies = TRUE)

library(foreign)
library(ggplot2)
library(tidyverse)
library(psych)
library(PairedData)
library(reshape)
library(dplyr)
library(ggpubr)
library(pastecs)
library(car)
library(plyr)
library(rstatix)
library(coin)
library(rcompanion)
library(DescTools)
library(ggthemes) 
library(reshape2)

#### Read VAS ####
setwd("")
VASData = read.csv("VAS_sensation.csv", header=TRUE, sep = ",", row.names=1)
VASData = VASData[,-49]

# create participant IDs
VASData$participant <- paste0("p", 1:24)

# define prefixes and columns
prefixes <- c("R_sma_", "S_sma_", "R_mpfc_", "S_mpfc_", "R_ag_", "S_ag_")
columns <- c("AudiIntense", "ScapIntense", "ScapArea", "Pain")

# take the mean of two replications
# Compute means for each combination of prefix and column 
for (prefix in prefixes) {
  for (column in columns) {
    col1 <- paste0(prefix, column, "_1")
    col2 <- paste0(prefix, column, "_2")
    mean_col_name <- paste0(prefix, column)
    VASData[[mean_col_name]] <- (VASData[[col1]] + VASData[[col2]]) / 2
  }
}

# Remove the original columns with numeric suffixes
VASData <- VASData[, -grep("_\\d$", names(VASData))]

# Convert the 'participant' column to a factor
VASData$participant <- factor(VASData$participant)

# View the structure of the VASData
str(VASData)

# wide to long format
vasdata <- melt(VASData, id.vars = "participant", 
                measure.vars = c("R_sma_AudiIntense", "R_sma_ScapIntense", "R_sma_ScapArea", "R_sma_Pain",
                                 "S_sma_AudiIntense", "S_sma_ScapIntense", "S_sma_ScapArea", "S_sma_Pain",
                                 "R_mpfc_AudiIntense", "R_mpfc_ScapIntense", "R_mpfc_ScapArea", "R_mpfc_Pain",
                                 "S_mpfc_AudiIntense", "S_mpfc_ScapIntense", "S_mpfc_ScapArea", "S_mpfc_Pain",
                                 "R_ag_AudiIntense", "R_ag_ScapIntense", "R_ag_ScapArea", "R_ag_Pain",
                                 "S_ag_AudiIntense", "S_ag_ScapIntense", "S_ag_ScapArea", "S_ag_Pain"),
                value.name = "score", variable.name = "groups")
# Extract and insert categorical variables
vasdata$Targets <- factor(gl(3, 192, labels = c("sma", "mpfc", "ag")))
vasdata$Conditions <- factor(gl(2, 96, labels = c("Active", "Sham")))
vasdata$Items <- factor(gl(4, 24, labels = c("AudiIntense", "ScapIntense", "ScapArea", "Pain")))

levels(vasdata$Targets)
levels(vasdata$Conditions)
levels(vasdata$Items)

str(vasdata)

#### SMA####
df.sma = subset(vasdata,vasdata$Targets=="sma")

#Boxplot
ggboxplot(df.sma, x = "Items", y = "score",
          color = "Conditions", palette = "jco",
          add = "jitter",shape = "Items",
          bxp.errorbar = TRUE,
          bxp.errorbar.width = 0.4) +
         labs(title="VAS_SMA", y = "VAS Score")+ 
         coord_cartesian(ylim = c(0, 10)) +
         scale_y_continuous(breaks = seq(0, 10, 2)) +
         theme_tufte(base_size = 12, ticks = TRUE)

#Statistics (WILCOXON SIGNED RANK SUM TEST)
res <- wilcox.test(VASData$R_sma_AudiIntense, VASData$S_sma_AudiIntense,
                   paired = TRUE, exact =  FALSE)
print(res) #p-value = 0.1109

res <- wilcox.test(VASData$R_sma_ScapIntense, VASData$S_sma_ScapIntense,
                   paired = TRUE, exact =  FALSE)
print(res) #p-value = 0.08608

res <- wilcox.test(VASData$R_sma_ScapArea, VASData$S_sma_ScapArea,
                   paired = TRUE, exact =  FALSE)
print(res) #p-value = 0.7777

res <- wilcox.test(VASData$R_sma_Pain, VASData$S_sma_Pain,
                   paired = TRUE, exact =  FALSE)
print(res) #p-value = 0.9041

#### mPFC ####
df.mpfc = subset(vasdata,vasdata$Targets=="mpfc")

#Boxplot
ggboxplot(df.mpfc, x = "Items", y = "score",
          color = "Conditions", palette = "jco",
          add = "jitter",shape = "Items",
          bxp.errorbar = TRUE,
          bxp.errorbar.width = 0.4) +
  labs(title="VAS_mPFC", y = "VAS Score")+ 
  coord_cartesian(ylim = c(0, 10)) +
  scale_y_continuous(breaks = seq(0, 10, 2)) +
  theme_tufte(base_size = 12, ticks = TRUE)

#Statistics
res <- wilcox.test(VASData$R_mpfc_AudiIntense, VASData$S_mpfc_AudiIntense,
                   paired = TRUE, exact =  FALSE)
print(res) #p-value = 0.07548

res <- wilcox.test(VASData$R_mpfc_ScapIntense, VASData$S_mpfc_ScapIntense,
                   paired = TRUE, exact =  FALSE)
print(res) #p-value = 0.001119

res <- wilcox.test(VASData$R_mpfc_ScapArea, VASData$S_mpfc_ScapArea,
                   paired = TRUE, exact =  FALSE)
print(res) #p-value = 0.02615

res <- wilcox.test(VASData$R_mpfc_Pain, VASData$S_mpfc_Pain,
                   paired = TRUE, exact =  FALSE)
print(res) #p-value = 0.008139

#### AG ####
df.ag = subset(vasdata,vasdata$Targets=="ag")

#Boxplot
ggboxplot(df.ag, x = "Items", y = "score",
          color = "Conditions", palette = "jco",
          add = "jitter",shape = "Items",
          bxp.errorbar = TRUE,
          bxp.errorbar.width = 0.4) +
  labs(title="VAS_AG", y = "VAS Score")+ 
  coord_cartesian(ylim = c(0, 10)) +
  scale_y_continuous(breaks = seq(0, 10, 2)) +
  theme_tufte(base_size = 12, ticks = TRUE)

#Statistics
res <- wilcox.test(VASData$R_ag_AudiIntense, VASData$S_ag_AudiIntense,
                   paired = TRUE, exact =  FALSE)
print(res) #p-value = 0.3932

res <- wilcox.test(VASData$R_ag_ScapIntense, VASData$S_ag_ScapIntense,
                   paired = TRUE, exact =  FALSE)
print(res) #p-value = 0.005254

res <- wilcox.test(VASData$R_ag_ScapArea, VASData$S_ag_ScapArea,
                   paired = TRUE, exact =  FALSE)
print(res) #p-value = 0.6587

res <- wilcox.test(VASData$R_ag_Pain, VASData$S_ag_Pain,
                   paired = TRUE, exact =  FALSE)
print(res) #p-value = 0.09595
