##%######################################################%##
#                                                          #
####               R package requirements               ####
#                                                          #
##%######################################################%##

# This script details all of the datasets and R packages required to 
# run the analyses associated with the paper.

outDir <- "0_data_package_requirements/"

if(!dir.exists(outDir)) dir.create(outDir)

sink(paste0(outDir,"log.txt"))

t.start <- Sys.time()

print(t.start)

############ 1. R packages required from github ############

# A number of R packages have been developed by Dr Tim Newbold for use when
# analysing the PREDICTS database. Code to download and install these 
# packages from Github is detailed here.

# the devtools library is required to download packages from Github
library(devtools)

# The library predictsFunctions includes various functions for organising, analysing
# and plotting outputs from analyses of the PREDICTS database. 
install_github(repo = "timnewbold/predicts-demo",subdir = "predictsFunctions")

# The StatisticalModels package includes various functions for analysing the
# PREDICTS database
install_github(repo = "timnewbold/StatisticalModels")

############ 2. Data to be downloaded ############

# Datasets used in this analysis can be downloaded from the following webpages.

# PREDICTs data - Biodiversity samples
# https://data.nhm.ac.uk/dataset/the-2016-release-of-the-predicts-database


# CRU data - TMP and TMX data for mean and max monthly temperatures
# https://crudata.uea.ac.uk/cru/data/hrg/

# Download and unzip into a 0_data folder the "cru_ts4.03.1901.2018.tmp.dat.nc.gz" from  "tmp" folder.


# Datasets downloaded from the above links should be saved into a directory named 0_data
datadir <- "0_data/"

if(!dir.exists(datadir)) dir.create(datadir)


# To see a list of all R packages used in this work, please run the following code from a 
# working directory containing all scripts:
library(renv) 
pckgs <- unique(dependencies(getwd()))
unique(pckgs$Package)

sessionInfo()

t.end <- Sys.time()

print(round(t.end - t.start,0))

sink()
