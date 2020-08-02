setwd("/Users/morril01/Documents/PhD/GlobalDA/code/")
library(dplyr)
library(gridExtra)
source("3_analysis/helper/helper_plots.R")
folder_objects="../data/roo/"

## Read all the ROO files which contain exposures in two groups
fles_in = list.files(folder_objects, full.names=TRUE)
roo_files = sapply(fles_in, readRDS)

## Determine the cancer types and the type of mutation fetaure under consideration
types = do.call('rbind', sapply(gsub("_ROO.RDS", "", basename(fles_in)), strsplit, split = '_'))
types = data.frame(types, stringsAsFactors=FALSE)
colnames(types) = c('Cancer_type', 'Type_substitution')

## Functions
give_filename = function(strct, strfeat){
  paste0("../data/roo//", strct, '_', strfeat, "_ROO.RDS")
}

save_plots = function(names_roo_files, slot_name){
  sapply(names_roo_files, function(roo_file_name){
    i = roo_files[[roo_file_name]]
    png_name = paste0("../results/reports/figure_roo_summary/", paste0( paste0(gsub("[.]", "_", gsub(".RDS", "", basename(roo_file_name))), slot_name, collapse = "")), ".png")
    if(typeof(i) == "logical"){ if(is.na(i)){
      png(png_name, width=6, height=3, unit="in", res=300)
      par(mfrow=c(1,2))
      plot(0,0)
      plot(0,0)
      dev.off()
    }}else{
      x = slot(i, slot_name)
      if(length(x[[1]]) == 0){
        png(png_name, width=6, height=3, unit="in", res=300)
        par(mfrow=c(1,2))
        plot(0,0,)
        plot(0,0)
        dev.off()
      }else{
        png(png_name, width=6, height=3, unit="in", res=300)
        par(mfrow=c(1,2))
        image(t(as(x[[1]], 'matrix')))
        image(t(as(x[[2]], 'matrix')))
        dev.off()
      }
    }
    png_name
  })
}
## Read in and save as table
table_characteristics = lapply(unique(types$Cancer_type), function(i){  data.frame(t(sapply(unique(types$Type_substitution), function(j) {
    .x = roo_files[[(give_filename(i,j))]]
    if(typeof(.x) != 'S4'){
      if(is.na(.x)){
      a = rep('-', 4)
      }
    }else{
      a = c(nrow(.x@count_matrices_all[[1]]), ncol(.x@count_matrices_all[[1]]), paste0(apply(.x@count_matrices_all[[1]], 2, sum), collapse = ", "),
            paste0(apply(.x@count_matrices_all[[2]], 2, sum), collapse = ", "))
    }
    a = sapply(a, function(i) if(nchar(i)>50){paste0(substr(i, start =1, stop=77), '...')}else{i})
    names(a) = c('Number of samples', 'Number of features', 'Sum of features for group 1', 'Sum of features for group 2')
    a
    })), stringsAsFactors = FALSE)
})

table_characteristics = do.call('rbind', lapply(1:length(table_characteristics), function(i) cbind(unique(types$Cancer_type)[i], rownames(table_characteristics[[i]]), table_characteristics[[i]])))
rownames(table_characteristics) = NULL
colnames(table_characteristics)[1:2] = c('Cancer_type', 'Feature_type')
table_characteristics$Number.of.samples = as.numeric(table_characteristics$Number.of.samples)
table_characteristics$Number.of.features = as.numeric(table_characteristics$Number.of.features)

### Save table in TeX
table_characteristics = table_characteristics %>% group_by(Feature_type) %>% dplyr::mutate(Number_samples_total=sum(Number.of.samples, na.rm = TRUE))
filename = "../results/reports/report_PCAWG.tex"

## Save heatmaps with exposures in groups, and append it to the report
pngnames1 = save_plots(names(roo_files), "count_matrices_all")
pngnames2 = save_plots(names(roo_files)[grepl("signatures", names(roo_files))], "count_matrices_active")
pngs = as.character(sort(c(pngnames1, pngnames2)))
write_pdfs = sapply(as.character(sort(unique(table_characteristics$Cancer_type))), function(ct){
  
  paste0(c("\\begin{figure}", sapply(paste0("figure_roo_summary/", basename(pngs[grepl(pattern = ct, x = pngs)])), function(i) paste0("\\begin{minipage}{.24\\textwidth}\\includegraphics[width=\\textwidth]{", i, "}\\end{minipage}")),
         "\\caption{", ct, "}\\end{figure}"), collapse = "")
})


## Write TeX
write_append = paste0("\\documentclass{article}\\usepackage{longtable}\\usepackage{graphicx}\\usepackage{array}\\date{\\today}\\newcolumntype{L}[1]{>{\\raggedright\\let\\newline\\\\\\arraybackslash\\hspace{0pt}}m{#1}}\\newcolumntype{C}[1]{>{\\centering\\let\\newline\\\\\\arraybackslash\\hspace{0pt}}m{#1}}\\newcolumntype{R}[1]{>{\\raggedleft\\let\\newline\\\\\\arraybackslash\\hspace{0pt}}m{#1}}\\setlength{\\textwidth}{7.1in}\\setlength{\\oddsidemargin}{-0.4in}\\setlength{\\textheight}{9.6in}\\setlength{\\topmargin}{-0.8in}\\begin{document}\\today\\section{Report of PCAWG objects}\\tiny ")
write_append = paste0(write_append, "\n\n", print(xtable::xtable(table_characteristics)), "\n\n")
# write_append = gsub("\\begin\\{table\\}[ht]\n\\centering\n\\begin\\{tabular\\}\\{rllrrllr\\}",
#      "\\begin\\{longtable\\}{R\\{.2in\\}L\\{1in\\}L\\{1in\\}R\\{.6in\\}R\\{.6in\\}L\\{1in\\}L\\{1in\\}R\\{1in\\}\\}",
#      write_append)

write_append = gsub("rllrrllr", "R{.2in}L{1in}L{.5in}R{.2in}R{.2in}L{1in}L{1in}R{1in}", write_append)

write_append = gsub("begin\\{table\\}", "begin\\{longtable\\}", write_append)
write_append = gsub("\\\\centering", "", write_append)
# write_append = gsub("\\n\\n\\\\begin\\{", "", write_append)
write_append = gsub("end\\{table\\}", "end\\{longtable\\}", write_append)
write_append = gsub("\\[ht\\]\n\n", "", write_append)
write_append = gsub("\\\\begin\\{tabular\\}", "", write_append)
write_append = gsub("\\\\end\\{tabular\\}", "", write_append)

substr(write_append, 1, 1000)
write_append = paste0(write_append, "\\clearpage", paste0(write_pdfs, collapse = ""))
write_append = paste0(write_append, "\\end{document}")

write(file = filename, x = write_append)


## Save a pdf with all barplots for all cancer types
## for each cancer type, plot three side-by-side barplots
roo_files
createBarplot()

pdf("../results/reports/barplots_PCAWG.pdf", width = 20)
for(cancer_type in sort(unique(types$Cancer_type))){
  plts = c()
  for(type_substitution in types %>% filter(Cancer_type == cancer_type) %>% select(Type_substitution) %>% unlist){
    x = roo_files[[paste0(folder_objects, '/', cancer_type, '_', type_substitution, '_ROO.RDS')]]
    if(is.na(x)){
      plts1 = list(ggplot(data.frame()) + geom_point() + xlim(0, 10) + ylim(0, 100)+ggtitle(paste0(cancer_type, ' ', type_substitution))) ## empty plot
    }else{
      plts1 = lapply(wrapper_barplot_ROO(x, remove_labels = TRUE), function(i) i+ggtitle(paste0(cancer_type, ' ', type_substitution))+theme(legend.position = "none"))
      plts1 = c(plts1, lapply(wrapper_barplot_ROO(x, bool_normalised = TRUE, remove_labels = TRUE), function(i) i+ theme(legend.position = "none")))
    }
    plts = c(plts, plts1)
  }  
  do.call("grid.arrange", c(plts, nrow=2, as.table = FALSE)) ## as.table = FALSE [saves it plots it by row
}
dev.off()

