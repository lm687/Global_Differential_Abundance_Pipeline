rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(ggplot2)
library(gridExtra)
library(dplyr)
library(readxl)

source("../3_analysis/helper/pcawg.colour.palette.R")

## moved to summary_TMB_PCAWG

