## replot the output of analyse_inference_simulations_integrate.R, if more plots are needed
## this way the p-values for competing methods do not need to be re-computed

rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(ggplot2)
library(ggrepel)
library(dplyr)

source("helper_model_assessment.R")
source("../../../2_inference_TMB/helper_TMB.R")
source("../../../1_create_ROO/roo_functions.R")

multiple_runs <- T
if(multiple_runs){
  flder_out <- "../../../../results/results_TMB/simulated_datasets/mixed_effects_models_multiple/"
  flder_in <- "../../../../data/assessing_models_simulation/inference_results/TMB/nlminb/summaries_multiple/"
}else{
  flder_out <- "../../../../results/results_TMB/simulated_datasets/mixed_effects_models/"
  flder_in <- "../../../../data/assessing_models_simulation/inference_results/TMB/nlminb/summaries/"
}

# datasetgeneration=['GenerationMixturePCAWG', #'GenerationMixturefewersignaturesPCAWG', \
#                    'GenerationMixturefewersignaturespairedPCAWG', \
#                    'GenerationMixturefewersignaturespairedstomachPCAWG',\
#                    'GenerationMixturefewersignaturespairedKidneyRCCPCAWG',\

generation <- c('GenerationMixturefewersignaturespairedKidneyRCCPCAWG')
generation <- c('GenerationMixturefewersignaturespairedPCAWG')
generation <- c('GenerationMixturefewersignaturespairedstomachPCAWG')
generation <- c('GenerationMixturefewersignaturesPCAWG')

a <- readRDS(paste0("../../../../data/assessing_models_simulation/summaries_synthetic_DA/", generation, ".RDS"))


varying_betashape <-give_accuracies_with_varying_var(var = 'beta_gamma_shape',
                                                     datasets_arg = a$datasets,
                                                     pvals_data_frame_arg = a$pvals_data_frame)

if((generation %in% c("GenerationMixturePCAWG", "GenerationMixturefewersignaturesPCAWG", "GenerationMixturefewersignaturespairedPCAWG")) | grepl('GenerationMixturefewersignaturespaired', generation) ){
  varying_betashape$beta_gamma_shape <- signif(varying_betashape$beta_gamma_shape, 2)
}

varying_betashape$model <- gsub("pvals_", "", varying_betashape$model)

.xx <- varying_betashape[which(varying_betashape$beta_gamma_shape == max(varying_betashape$beta_gamma_shape)),]
## sometimes there are problems with numeric precision. if that is the case, select as the last value only one of the selected rows
.xx <- .xx[!duplicated(.xx$model),]

if(generation == 'GenerationMixturefewersignaturespairedKidneyRCCPCAWG'){
  title = 'Kidney-RCC (paired)'
}else if( generation == 'GenerationMixturefewersignaturespairedPCAWG'){
  title = 'Liver-HCC (paired)'
}else if ( generation  == 'GenerationMixturefewersignaturespairedstomachPCAWG'){
  title = 'Stomach-AdenoCa (paired)'
}else if( generation == 'GenerationMixturefewersignaturesPCAWG'){
  title = 'Liver-HCC (not paired)'
}

ggplot(varying_betashape, aes(x=factor(beta_gamma_shape), y = Accuracy, col=model, group=model, label=model,
                              lty=model%in% c('fullREM', 'fullREDMSL', 'diagREDMSL', 'diagREDM')))+
  geom_point()+geom_line()+theme_bw()+
  scale_color_manual(values=colours_models2)+labs(col='', x='Percentage of mixture')+guides(col='none', lty='none')+
  ggtitle(title)+
  geom_label_repel(data = .xx,
                   max.overlaps = Inf, aes(x=factor(max(varying_betashape$beta_gamma_shape))), direction = "y",
                   nudge_x=max(varying_betashape$beta_gamma_shape)+2,
                   size=2.5)+
  coord_cartesian(xlim = c(0, length(unique(varying_betashape$beta_gamma_shape))*1.5))+
  theme_bw()+theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust=1))
ggsave(paste0(flder_out, generation, "/summaries/accuracy_with_betagammashape_palette2_factorv2.pdf"),
       height = 3.0, width = 4.0)

