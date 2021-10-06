library(readxl)
library(gridExtra)
gerstung_changing_sigs <- readxl::read_excel("/Users/morril01/Documents/PhD/GlobalDA/data/restricted/pcawg/41586_2019_1907_MOESM9_ESM.xlsx", sheet = "signatures-changing")
gerstung_constant_sigs <- readxl::read_excel("/Users/morril01/Documents/PhD/GlobalDA/data/restricted/pcawg/41586_2019_1907_MOESM9_ESM.xlsx", sheet = "signatures-constant")
gerstung_changing_sigs
gerstung_constant_sigs

gerstung_changing_sigs_earlylate <- gerstung_changing_sigs[gerstung_changing_sigs$time == "early-late",]
gerstung_changing_sigs_clonalsubclonal <- gerstung_changing_sigs[gerstung_changing_sigs$time == "clonal-subclonal",]

df_changes_el <- gerstung_changing_sigs_earlylate %>% group_by(signature) %>% dplyr::summarise(median(mean_change))
df_changes_cs <- gerstung_changing_sigs_clonalsubclonal %>% group_by(signature) %>% dplyr::summarise(median(mean_change))


grid.arrange(ggplot(gerstung_changing_sigs_earlylate, aes(x=factor(signature, levels=df_changes_el$signature[order(df_changes_el$`median(mean_change)`)]),
                                   y=mean_change, group=signature,col=histology_abbreviation))+geom_boxplot()+geom_jitter()+
               geom_hline(yintercept = 0, lty='dashed')+
  theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust=1))+guides(col=FALSE)+theme_bw(),
ggplot(gerstung_changing_sigs_clonalsubclonal, aes(x=factor(signature, levels=df_changes_cs$signature[order(df_changes_cs$`median(mean_change)`)]),
                                             y=mean_change, group=signature,col=histology_abbreviation))+geom_hline(yintercept = 0, lty='dashed')+
  geom_boxplot()+  geom_jitter()+
  theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust=1))+guides(col=FALSE)+theme_bw(),
nrow=2)
