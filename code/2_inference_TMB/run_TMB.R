## Creating dataset to assess the performance of models for parameter recovery

rm(list = ls())

library(TMB)
library(optparse)

source("1_create_ROO/roo_functions.R")
source("2_inference/helper/helper_DA_stan.R") ## for normalise_rw
source("2_inference_TMB/mm_multinomial/helper_functions.R")
source("2_inference_TMB/helper_TMB.R")

option_list = list(
  make_option(c("--model"), type="character", default=NA,
              help="Which model to use for inference", metavar="character"),
  make_option(c("--input"), type="character", default=NA,
              help="Input file with dataset (RDS)", metavar="character"),
  make_option(c("--output"), type="character", default=NA,
              help="Output file in which to write the results of the inference (RDS file)", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if(opt$model == "fullREM"){
  TMB::compile("2_inference_TMB/mm_multinomial/fullRE_ME_multinomial.cpp",  "-std=gnu++17")
  dyn.load(dynlib("2_inference_TMB/mm_multinomial/fullRE_ME_multinomial"))
  mod_model_name = "fullRE_M"
}else if(opt$model == "fullREDM"){
  TMB::compile("2_inference_TMB/mm_multinomial/fullRE_ME_dirichletmultinomial.cpp",  "-std=gnu++17")
  dyn.load(dynlib("2_inference_TMB/mm_multinomial/fullRE_ME_dirichletmultinomial"))
  mod_model_name = "fullRE_DM"
}else if(opt$model == "diagREDM"){
  TMB::compile("2_inference_TMB/mm_multinomial/diagRE_ME_dirichletmultinomial.cpp",  "-std=gnu++17")
  dyn.load(dynlib("2_inference_TMB/mm_multinomial/diagRE_ME_dirichletmultinomial"))
  mod_model_name = "diagRE_DM"
}else{
  stop('Specifiy a valid <model>')
}

cat(opt$input)
dataset = load_PCAWG(ct = opt$input, typedata = "signatures", simulation = T, path_to_data = NA)
# results_inference = try(wrapper_run_TMB(opt$input, model = mod_model_name, typedata = "simulation", simulation = TRUE))
results_inference = try(wrapper_run_TMB(object = dataset, sort_columns = T, model = mod_model_name))

saveRDS(object = results_inference, file = opt$output)
