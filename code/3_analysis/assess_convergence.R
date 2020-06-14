#setwd(dirname(rstudioapi::getSourceEditorContext()$path)); setwd("../")

source("2_inference/helper/monitornew.R")
source("2_inference/helper/monitorplot.R")

library(rstan)
library(bayesplot)
library(ggplot2)
library(dplyr)
library(reshape2)
library(patchwork)
library(latex2exp)
library(bayesplot)
#library(rstanarm)
library(plyr)
library(cowplot)
library(optparse)

option_list = list(
  make_option(c("--file_posterior"), type="character", default=NA, 
              help="File with the posterior, with directory included", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# opt=list(); opt$file_posterior = "../data/inference/M_6000_CNS-PiloAstro_signatures.Rdata"

opts = opt
load(opt$file_posterior)
opt = c(opt, opts)
name = gsub("[.]", "_", gsub("-", "_", gsub(".Rdata", "", basename(opt$file_posterior))))

folder_pdfs = "../results/convergence/pdfs/"
pdf_name = paste0(folder_pdfs, "convergence_report_", name, ".tex")
path_figures = "../results/convergence/figures/"

day = paste0(strsplit(as.character(file.info(opt$file_posterior)$ctime), ' ')[[1]][1])
fit = extract(fit_stan)
posterior2 <- as.matrix(fit_stan)
posterior2_no_theta = posterior2[,!grepl('theta', colnames(posterior2))]

color_scheme_set("mix-blue-pink")
filename_betachains = paste0(path_figures, name, "beta_chains.png")
# p <- bayesplot::mcmc_trace(posterior2,  pars = colnames(posterior2)[grep("beta", colnames(posterior2))],
#                            facet_args = list(nrow = 1, labeller = label_parsed))+ theme(text = element_text(size=10))
# p + facet_text(size = 6)
# ggsave(filename_betachains, width = 7, height = 1.5, device = png())


# rhat =  rstan::Rhat(posterior2)
# rhat2 =  rstan::Rhat(posterior2_no_theta)
# rhat_per_par = apply(posterior2_no_theta, 2, rstan::Rhat)
essb = ess_bulk(posterior2)
esst = rstan::ess_tail(posterior2)

mon2 <- monitor_extra(fit_stan)

which_min_ess <- which.min(mon2[!grepl('theta', rownames(mon2)), 'tailseff'])

# print(summary(fit_stan)$summary)
# mu_tau_summary <- summary(fit_stan, probs = c(0.1, 0.9))$summary

# color_scheme_set("green")
# mcmc_acf(fit_stan, pars = c("beta[1,1]"))
mcmc_rhat(apply(posterior2_no_theta, 2, rhat))
filename_mcmc_rhat = paste0(path_figures, name, "mcmc_rhat.png")
ggsave(filename_mcmc_rhat, height = 2, width = 8)
## Plotting

filename_local_ess_min = paste0(path_figures, name, "local_ess_min.png")
plot_local_ess(fit = fit_stan, par = which_min_ess, nalpha = 20)
ggsave(filename_local_ess_min, height = 4, width = 4)

filename_local_ess_beta = paste0(path_figures, name, "local_ess_beta.png")
a = lapply(which(grepl('beta',  rownames(mon2))), function(i) plot_local_ess(fit = fit_stan, par = i, nalpha = 20))
plot_grid(plotlist = a, nrow = 2)
ggsave(filename_local_ess_beta, height = 2, width = 10)


filename_quantile_ess = paste0(path_figures, name, "quantile_ess.png")
plot_quantile_ess(fit = fit_stan, par = which_min_ess, nalpha = 40)
ggsave(filename_quantile_ess, height = 4, width = 4)

filename_change_ess = paste0(path_figures, name, "quantile_ess_min.png")
plot_change_ess(fit = fit_stan, par = which_min_ess)
ggsave(filename_change_ess, height = 4, width = 5)

samp <- as.array(fit_stan)
xmin <- "beta[1,1]" #paste0("x[", which_min_ess, "]")

#mcmc_hist_r_scale(samp[, , xmin])

#plot_ess(mon2) 
#plot_rhat(mon2)

#pairs(fit_stan)

maxrhat = max(apply(posterior2_no_theta, 2, rhat))
## Write to TeX

write_append = paste0("\\documentclass{article}\\usepackage{graphicx}\\date{\\today}\\setlength{\\textwidth}{7.1in}\\setlength{\\oddsidemargin}{-0.4in}\\setlength{\\textheight}{9.6in}\\setlength{\\topmargin}{-0.8in}\\begin{document}\\today\\section*{ Convergence report for ", gsub("_", " ", name), "} ")
write_append = paste0(write_append,  paste0("Date of inference: ", day, "\\\\"))
write_append = paste0(write_append, paste0("\\textbf{Explanation of parameters} Rhat indicates the ratio of variance within a change over variance pooling all chains. \\emph{We recommend running at least four chains by default and only using the sample if R-hat is less than 1.05.} ESS bulk/tail indicate the bulk/tail effective sample size estimate, i.e. show the sampling efficiency of mean and median estimates. \\emph{We recommend running at least four chains by default and only using the sample if R-hat is less than 1.05. Both bulk-ESS and tail-ESS should be at least 100 (approximately) per Markov Chain in order to be reliable and indicate that estimates of respective posterior quantiles are reliable.}"))
write_append = paste0(write_append, paste0("\\begin{itemize}\\item Max Rhat: ", maxrhat, "\n"))
write_append = paste0(write_append, paste0("\\item ess bulk: ", essb, "\n"))
write_append = paste0(write_append, paste0("\\item ess tail: ", esst, "\\end{itemize}\n"))
# write_append = paste0(write_append, paste0("\\subsection*{Chains for $\\beta$}\\includegraphics[width=\\textwidth]{../figures/", basename(filename_betachains), "}\n"))
write_append = paste0(write_append, paste0("\\begin{minipage}{\\textwidth}\\includegraphics[width=\\textwidth]{../figures/", basename(filename_mcmc_rhat), "}\\end{minipage}\\\\"))
write_append = paste0(write_append, paste0("In the plots below ESS should have high values\\\\"))
write_append = paste0(write_append, paste0("\\begin{minipage}{.3\\textwidth}\\includegraphics[width=\\textwidth]{../figures/", basename(filename_local_ess_min), "}\\end{minipage}"))
write_append = paste0(write_append, paste0("\\begin{minipage}{.3\\textwidth}\\includegraphics[width=\\textwidth]{../figures/", basename(filename_quantile_ess), "}\\end{minipage}"))
write_append = paste0(write_append, paste0("\\begin{minipage}{.36\\textwidth}\\includegraphics[width=\\textwidth]{../figures/", basename(filename_change_ess), "}\\end{minipage}\n"))
write_append = paste0(write_append, paste0("\\subsection*{Local ESS for $\\beta$}\\includegraphics[width=\\textwidth]{../figures/", basename(filename_local_ess_beta), "}\n"))
# write_append = paste0(write_append, "\n\n", xtable::xtable(mon), "\n\n")
write_append = paste0(write_append, paste0( "\\end{document}"))
write(file = pdf_name, x = write_append)
print(getwd())
print(pdf_name)
print( basename(pdf_name) )

system(paste0("cd ", folder_pdfs, " pdflatex ", basename(pdf_name)))
