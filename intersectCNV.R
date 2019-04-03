#!/usr/bin/env Rscript
library(data.table)
library(dplyr)
library(stringr)

## functions
splitV <- function(dt){
  dt[, Alter := sapply(V68, function(x){
    tmp <- unlist(strsplit(x, ":"))
    tmp <- unlist(strsplit(tmp[3], ","))
    return(as.numeric(tmp[2]))
  })]
  dt[, Depth := sapply(V68, function(x){
    tmp <- unlist(strsplit(x, ":"))
    return(as.numeric(tmp[2]))
  })]
  dt[, MAF := Alter/Depth]
  
  ## filters
  ## remove low allelic fractions (germline has this step later as well at .25)
  # dt <- dt[MAF > 0.1 | Gene.refGene == "KRAS"]
  dt <- dt[MAF >= 0.05]
  ## scale by log2depth and remove too those too low
  dt[, log2Depth := log2(Depth)]
  dt[, Depth.z := scale(log2Depth), by="Chr"]
  # dt <- dt[Depth.z > -1.5 | Gene.refGene == "KRAS"]
  # dt <- dt[Depth.z > -2]
  # ## optional (remove singletons)
  # dt <- dt[Alter > 1 | Gene.refGene == "KRAS"]
  # dt <- dt[Depth.z > -1.75]
  return(dt)
}

dbFilter <- function(dt){
  ## Exac
  dt[, c("ExAC_ALL", "ExAC_AFR", "ExAC_AMR", "ExAC_EAS", "ExAC_FIN", "ExAC_NFE", "ExAC_OTH", "ExAC_SAS") := .(as.numeric(ExAC_ALL), as.numeric(ExAC_AFR), as.numeric(ExAC_AMR), as.numeric(ExAC_EAS), as.numeric(ExAC_FIN), as.numeric(ExAC_NFE), as.numeric(ExAC_OTH), as.numeric(ExAC_SAS))]
  dt[is.na(ExAC_ALL), ExAC_ALL := 0]
  dt[, ExacMax := pmax(ExAC_ALL, ExAC_AFR, ExAC_AMR, ExAC_EAS, ExAC_FIN, ExAC_NFE, ExAC_OTH, ExAC_SAS, na.rm=TRUE)]
  dt <- dt[ExAC_ALL < 0.01 & ExacMax < 0.01]
  dt[, c("Polyphen2_HDIV_score", "Polyphen2_HVAR_score") := .(as.numeric(Polyphen2_HDIV_score), as.numeric(Polyphen2_HVAR_score))]
  
  dt[, c("sift", "poly", "lrt", "taste", "provean", "vest", "cadd", "fathmm", "gerp", "siphy") := .(ifelse(SIFT_score == ".", NA, ifelse(as.numeric(SIFT_score) < 0.03, 1, 0)),
                                                                                                    ifelse(Polyphen2_HDIV_score == ".", NA, ifelse(as.numeric(Polyphen2_HDIV_score) > 0.95 , 1, 0)),
                                                                                                    ifelse(LRT_score == ".", NA, ifelse(as.numeric(LRT_score) < 0.005, 1, 0)),
                                                                                                    ifelse(MutationTaster_score == ".", NA, ifelse(as.numeric(MutationTaster_score) > 0.995, 1, 0)),
                                                                                                    ifelse(PROVEAN_score == ".", NA, ifelse(as.numeric(PROVEAN_score) < -2.4, 1, 0)),
                                                                                                    ifelse(VEST3_score == ".", NA, ifelse(as.numeric(VEST3_score) > 0.75, 1, 0)),
                                                                                                    ifelse(CADD_phred == ".", NA, ifelse(as.numeric(CADD_phred) > 22, 1, 0)),
                                                                                                    ifelse(`fathmm-MKL_coding_score` == ".", NA, ifelse(as.numeric(`fathmm-MKL_coding_score`) > 0.92, 1, 0)),
                                                                                                    ifelse(`GERP++_RS` == ".", NA, ifelse(as.numeric(`GERP++_RS`) > 4.4, 1, 0)),
                                                                                                    ifelse(SiPhy_29way_logOdds == ".", NA, ifelse(as.numeric(SiPhy_29way_logOdds) > 14, 1, 0))
  )]
  dt[, dbRatio := (sum(sift, poly, lrt, taste, provean, vest, cadd, fathmm, gerp, siphy, na.rm = TRUE)/sum(!is.na(sift), !is.na(poly), !is.na(lrt), !is.na(taste), !is.na(provean), !is.na(vest), !is.na(cadd), !is.na(fathmm), !is.na(gerp), !is.na(siphy))), by=c("Chr", "Start")]
  dt <- dt[dbRatio > 0.5]
  
  return(dt)
}

cosmicCounts <- function(x){
  ## removes ID=COSMXXXX
  prelist <- unlist(strsplit(x, ";"))
  preoccur <- unlist(strsplit(prelist[2], "="))
  occur <- unlist(strsplit(preoccur[2], ","))
  occur <- gsub(".*[(]haematopoietic_and_lymphoid_tissue[)]", "0", occur)
  occur <- lapply(occur, function(x){gsub("[(].*[)]", "", x)})
  occur <- as.numeric(occur)
  return(sum(occur))
}

## arguments 
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = basename(args[1]) %>% str_replace(replacement = "", pattern = ".txt") %>% paste0(".filteredAnnotation.txt")
}
dn <- dirname(args[1])

## read in files
dt <- fread(args[1], sep="\t") %>% splitV %>% dbFilter
dt <- dt[Func.refGene == "exonic" & !is.element(ExonicFunc.refGene, c("synonymous SNV"))]
dt[, cosm := sapply(cosmic86, cosmicCounts)]

write.table(dt,paste0(dn,"/",args[2]), sep="\t", quote = FALSE, row.names = FALSE)

