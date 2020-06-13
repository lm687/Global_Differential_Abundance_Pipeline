
setwd("/Users/morril01/Documents/PhD/GlobalDA/code/")
library(optparse)
library(dplyr)

option_list = list(
  make_option(c("--file_ROO"), type="character", default=NA, 
              help="File with the posterior, with directory included", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

opt = list(); opt$file_ROO = "../data/roo/Biliary-AdenoCA_nucleotidesubstitution1_ROO.RDS"

roo_file = readRDS(opt$file_ROO)
fles_in = list.files("../data/roo/", full.names=TRUE)
roo_files = sapply(fles_in, readRDS)

nrow(roo_file@count_matrices_all[[1]])
roo_file@count_matrices_active

types = do.call('rbind', sapply(gsub("_ROO.RDS", "", basename(fles_in)), strsplit, split = '_'))
types = data.frame(types, stringsAsFactors=FALSE)
colnames(types) = c('Cancer_type', 'Type_substitution')


give_filename = function(strct, strfeat){
  paste0("../data/roo//", strct, '_', strfeat, "_ROO.RDS")
}

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

### why are there so few samples?
table_characteristics = table_characteristics %>% group_by(Feature_type) %>% dplyr::mutate(sum_all_ct=sum(Number.of.samples, na.rm = TRUE))
filename = "../results/reports/report_PCAWG.tex"

write_append = paste0("\\documentclass{article}\\usepackage{longtable}\\usepackage{graphicx}\\usepackage{array}\\date{\\today}\\newcolumntype{L}[1]{>{\\raggedright\\let\\newline\\\\\\arraybackslash\\hspace{0pt}}m{#1}}\\newcolumntype{C}[1]{>{\\centering\\let\\newline\\\\\\arraybackslash\\hspace{0pt}}m{#1}}\\newcolumntype{R}[1]{>{\\raggedleft\\let\\newline\\\\\\arraybackslash\\hspace{0pt}}m{#1}}\\setlength{\\textwidth}{7.1in}\\setlength{\\oddsidemargin}{-0.4in}\\setlength{\\textheight}{9.6in}\\setlength{\\topmargin}{-0.8in}\\begin{document}\\today\\section{Report of PCAWG objects}\\tiny ")
write_append = paste0(write_append, "\n\n", print(xtable::xtable(table_characteristics)), "\n\n")
# write_append = gsub("\\begin\\{table\\}[ht]\n\\centering\n\\begin\\{tabular\\}\\{rllrrllr\\}",
#      "\\begin\\{longtable\\}{R\\{.2in\\}L\\{1in\\}L\\{1in\\}R\\{.6in\\}R\\{.6in\\}L\\{1in\\}L\\{1in\\}R\\{1in\\}\\}",
#      write_append)

write_append = gsub("rllrrllr", "R{.2in}L{1in}L{1in}R{.6in}R{.6in}L{1in}L{1in}R{1in}", write_append)

write_append = gsub("begin\\{table\\}", "begin\\{longtable\\}", write_append)
write_append = gsub("\\\\centering", "", write_append)
write_append = gsub("\\n\\n\\\\begin\\{", "", write_append)
write_append = gsub("end\\{table\\}", "end\\{longtable\\}", write_append)
write_append = paste0(write_append, "\\end{document}")
write(file = filename, x = write_append)


sapply(names(roo_files), function(roo_file_name){
  i = roo_files[[roo_file_name]]
  png_name = paste0("../results/reports/figure_roo_summary/", gsub(".RDS", ".png", basename(roo_file_name)))
  if(typeof(i) == "logical"){ if(is.na(i)){
    png(png_name, width=6, height=3, unit="in", res=300)
    par(mfrow=c(1,2))
    plot(0,0, main=roo_file_name)
    plot(0,0, main=roo_file_name)
    dev.off()
  }}else{
    x = i@count_matrices_all
    png(png_name, width=6, height=3, unit="in", res=300)
    par(mfrow=c(1,2))
    image(as(x[[1]], 'matrix'), main=roo_file_name)
    image(as(x[[2]], 'matrix'), main=roo_file_name)
    dev.off()
  }
})




