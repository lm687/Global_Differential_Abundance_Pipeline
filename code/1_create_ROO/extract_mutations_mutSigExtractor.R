if(debug_bool){
  setwd("~/Documents/PhD/GlobalDA/code/")
  opt=list();
  opt$input_files='../data/restricted/pcawg/pcawg_restricted_snv_counts/00b9d0e6-69dc-4345-bffd-ce32880c8eef ../data/restricted/pcawg/pcawg_restricted_snv_counts/02917220-6a7a-46a1-8656-907e96bef88e ../data/restricted/pcawg/pcawg_restricted_snv_counts/03ad38a6-0902-4aaa-84a3-91ea88fa9883 ../data/restricted/pcawg/pcawg_restricted_snv_counts/05616329-e7ba-4efd-87b1-d79cd0f7af3d ../data/restricted/pcawg/pcawg_restricted_snv_counts/068f4f69-d2fe-4f25-912e-ca7d4623efb6 ../data/restricted/pcawg/pcawg_restricted_snv_counts/0e7f46ca-6f5c-4538-b6d6-00af65d57fcf ebe0ed67-2d3f-45cd-8f9b-4912595b16a0 ../data/restricted/pcawg/pcawg_restricted_snv_counts/fa676301-902f-473f-8313-5bff34ae549a ../data/restricted/pcawg/pcawg_restricted_snv_counts/fce8d8c6-f2a0-43a8-9a7a-b9c519a6686c'
  opt$cancer_type=" Lymph-BNHL";
}else{
  option_list = list(
    make_option(c("--cancer_type"), type="character", default=NA, 
                help="", metavar="character"),
    make_option(c("--output"), type="character", default=NA, 
                help="", metavar="character"),
    make_option(c("--input_files"), type="character", default=NA, 
                help="", metavar="character")
  );
  
  opt_parser = OptionParser(option_list=option_list);
  opt = parse_args(opt_parser);
}

flder_clonal = "../data/restricted/pcawg/consensus_subclonal_reconstruction_20170325/"
flder_vcf = "../data/restricted/pcawg/pcawg_restricted_snv_counts/"
fles_vcf = list.files(flder_vcf)
fles_cancer_type = paste0(basename(strsplit(opt$input_files, " ")[[1]]), "_cluster_assignments.txt.gz")


# install.packages("devtools")
# devtools::install_github('https://github.com/UMCUGenetics/mutSigExtractor/')
library(mutSigExtractor)
# library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)

#' https://github.com/UMCUGenetics/mutSigExtractor

##' Problem: the vcf from PCAWG is from their own aasembly hs37d5.
##' Looking at genome.fa, it looks as though it's GRCh37 (i.e. hg19)
##' mutSigExtractor uses as a default genome hg19, or we can also use hg38, but that throws an error


## read ct and information

## read vcf
# list.files("../../../data/restricted/pcawg/pcawg_restricted_snv/")
# vcf_snv <- "../../../data/restricted/pcawg/pcawg_restricted_snv/0a6be23a-d5a0-4e95-ada2-a61b2b5d9485.consensus.20160830.somatic.snv_mnv.vcf"

## read ccf

contexts_snv_hg38 <- extractSigsSnv(vcf.file=vcf_snv, output='contexts',
                                    ref.genome=BSgenome.Hsapiens.UCSC.hg38)
contexts_snv <- extractSigsSnv(vcf.file=vcf_snv, output='contexts',
                               ref.genome=BSgenome.Hsapiens.UCSC.hg19)
head(contexts_snv)

sigs_snv <- fitToSignatures(
  mut.context.counts=contexts_snv[,1], 
  signature.profiles=SBS_SIGNATURE_PROFILES_V3
)
head(sigs_snv)

# 0a6be23a-d5a0-4e95-ada2-a61b2b5d9485	Prost-AdenoCA	TRUE	TRUE

my_ct <- load_PCAWG(ct = "Prost-AdenoCA", typedata = "signatures", path_to_data = "../../../data/")
my_ct <- colSums(my_ct$Y[(grepl("0a6be23a-d5a0-4e95-ada2-a61b2b5d9485", rownames(my_ct$Y))),])

comparison <- cbind(my_ct, sigs_snv[match(names(my_ct), names(sigs_snv))])
plot(comparison)


