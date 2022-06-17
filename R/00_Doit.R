library('tidyverse')
library('broom')
library('GSVA')
source("R/99_func_file.R")

#run project pipeline
source("R/01_load.R") #define dataset to run
source("R/02_CIT_signature1.R")
source("R/03_CIT_signature2.R")
source("R/04_CIT_signature3.R")
source("R/05_purity.R")