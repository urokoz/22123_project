library('tidyverse')
library('broom')
library('GSVA')
library(utils)
library('estimate')
source("R/99_func_file.R")


source("R/01_load.R")
#pick pipeline
#CIT
source("R/02_CIT_signature.R")
source("R/03_CIT_signature1.R")
source("R/04_CIT_signature3.R")
source("R/05_purity.R")


#GBM
source("R/GBM_signatures1.R")
source("R/GBM_signatures2.R")
source("R/GBM_signatures3.R")
source("R/purity_GBM.R")
