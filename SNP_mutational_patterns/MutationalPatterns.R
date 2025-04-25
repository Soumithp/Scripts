#module load R/3.3.2-gccmkl

library(MutationalPatterns)
library(BSgenome)
library(NMF)
library(lsa)
library(gridExtra)


args = commandArgs(trailingOnly=TRUE)
if (length(args) != 4 ) {
  stop("Usage: Rscript --vanilla MutationalPatterns.R <vcf_location> <sample_name> <ref_code> <cosmic_sig>", call.=FALSE)
} else if (length(args) == 4) {
  vcf_location <- args[1]
  sample_name <- args[2]
  ref_code <- args[3]
  cosmic_sig <- args[4]
}



mutationSignatureContribution <- function(vcf_location, sample_name, ref_code = "hg19", sp_url, sig_version='v3.4' )
{
    ref_genome <- paste0("BSgenome.Hsapiens.UCSC.",ref_code)
    #ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
    library(ref_genome, character.only = TRUE)
    
    # vcfs <- read_vcfs_as_granges( as.character(meta$loc), as.character(meta$sample), genome = "hg19")
    vcfs <- read_vcfs_as_granges( c(as.character(vcf_location)), c(as.character(sample_name)), genome = ref_code)
    auto <- extractSeqlevelsByGroup(species="Homo_sapiens",style="UCSC",group="auto")
    vcfs <- lapply(vcfs, function(x) keepSeqlevels(x, auto, pruning.mode = "coarse"))
    
    # count mutation type occurrences
    type_occurrences <- mut_type_occurrences(vcfs, ref_genome)
    # plot mutation spectrum
    # plot_spectrum(type_occurrences, CT = TRUE)
    # 96 mutation profile
    mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
    # plot_96_profile(mut_mat)
    # avoid 0 in matrix
    mut_mat <- mut_mat + 1e-4
    #mut_mat2 = apply(mut_mat, 1, function(x) exp( sum(log(x))/length(x) ) )
    #mut_mat2 = apply(mut_mat, 1, median )
    #mut_mat2 = as.matrix(mut_mat2)
    #plot_96_profile(mut_mat2)
    
    SUBSTITUTIONS <- c( "[C>A]", "[C>G]", "[C>T]", "[T>A]", "[T>C]", "[T>G]" )
    substitution = rep(SUBSTITUTIONS, each = 16)
    mut_type = paste0( substr(rownames(mut_mat),1,1), substitution, substr(rownames(mut_mat),3,3) )
    
    #sig_dir <- "/project/shared/DSSR/s229294/hoshida_lab_rnaseq/rnaseq_SNP_calling/workflow/output/mut_sig/"
    #sig_version <- "v3.4"
    #if(grepl('hg19',ref_genome)) sp_url <- paste0(sig_dir,"/COSMIC_",sig_version,"_SBS_GRCh37.txt")
    cancer_signatures <- read.table(sp_url, sep = "\t", header = TRUE)
    rownames(cancer_signatures) = as.character(cancer_signatures[,1])
    cancer_signatures = as.matrix(cancer_signatures[,-1])
    cancer_signatures = cancer_signatures[ match(mut_type,rownames(cancer_signatures)), ]
    mut_type == rownames(cancer_signatures)
    
    fit_res <- fit_to_signatures(mut_mat, cancer_signatures)
    select <- which(rowSums(fit_res$contribution) > 0)
    #plot_contribution(fit_res$contribution[select,], cancer_signatures[,select], coord_flip = FALSE,  mode = "relative")
    
    #contribution = data.frame( signature=rownames(fit_res$contribution[select,]), fit_res$contribution[select,] )
    contribution <- data.frame( signature=rownames(fit_res$contribution), fit_res$contribution )
    signature_contribution_name <- paste0(sample_name,"_COSMIC_Mutational_Signatures_",sig_version,"_sampleLevelContribution.txt")
    write.table(contribution,signature_contribution_name,sep="\t",col.names=T,row.names=F,quote=F)
    
}

mutationSignatureContribution(
  vcf_location = vcf_location,
  sample_name = sample_name,
  ref_code = ref_code,
  sp_url = cosmic_sig
)
