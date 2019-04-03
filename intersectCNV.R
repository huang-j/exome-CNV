#!/usr/bin/env Rscript
library(data.table)
library(dplyr)
library(stringr)
# library(IRanges)

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0 | length(args)==1) {
  stop("At least two argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==2) {
  # default output file
  args[3] = "~/rsrch/PDAC_genelist_hg19.bed"
}
# args=c("~/rsrch/CNV/BW13_T_reseq_S4.umi/BW13_T_reseq_S4.umi.called.seg",
#        "~/Downloads/crev2-dp10-family1-paired-events.txt",
#        "~/rsrch/PDAC_genelist_hg19.bed")
dn <- dirname(args[1])

samplename <- unlist(str_split(args[1], "/"))
samplename <- unlist(str_split(samplename[length(samplename)], "[.]"))
samplename <- samplename[1]
if(samplename == "GV79-PBMC_S5"){
  samplename <- "GV79-T_S6"
}
## genelist
print("cnvGL")
cnvGL <- fread(args[3], sep="\t", col.names = c("Chr", "Start", "End", "Gene stable ID", "Gene", "Transcript stable ID"))
setkey(cnvGL, "Chr", "Start", "End")
## arg 1 is Segmentation
print("gseg")
gseg <- fread(args[1], sep="\t", col.names=c("Chr", "Start", "End", "Markers", "Copy.Ratio", "g.Event"), skip = "CONTIG")
setkey(gseg, "Chr", "Start", "End")
## arg 2 is haploseqstuff
print("hap")
hap <- fread(args[2], sep="\t", col.names=c("Chr", "Start", "End", "EVENT_STATE", "NUM_INFORMATIVE_MARKERS", "MEAN_POSTPROB", "PHASE_CONCORDANCE",	"sample_id",	"baf",	"pval",	"mean_markers",	"sample_baf",	"baf_dev", "event"))
setkey(hap, "Chr", "Start", "End")
hap <- hap[grepl(samplename, sample_id)]

test <- foverlaps(gseg, cnvGL, type="any")
colnames(test) <- c("Chr","Start","End","Gene stable ID","Gene","Transcript stable ID","g.Start","g.End","Markers","Copy.Ratio","g.Event")
test2 <- foverlaps(test, hap, by.x=c("Chr", "g.Start", "g.End"), type="any")
# test2 <- test2[!is.na(Gene)]
test2 <- test2[, .(sample_id, Chr, Gene, i.Start, i.End, g.Start, g.End, Markers, Copy.Ratio, g.Event, Start, End,  MEAN_POSTPROB,  pval, NUM_INFORMATIVE_MARKERS, event, EVENT_STATE, PHASE_CONCORDANCE, baf,mean_markers, sample_baf, baf_dev, `Gene stable ID`, `Transcript stable ID`)]
## if Haplo doesn't call it, then use g.Event (compared to copy ratio) if it does, then we need to run a comparison between the calls
##            The comparison should be as follows:
##                check if the calls are the same: if yes, then go with call
##                    check pvalue
##                if not we check to see if copy ratio. If ratio doesn't match call then we use ratio to determine the call
test2[, sample_id := eval(samplename)]
callTest <- function(MEAN_POSTPROB, g.Event, pval, event, Copy.Ratio, NUM_INFORMATIVE_MARKERS){
  checkCopyRatio <- function(post_prob, Copy.Ratio){
    if (!post_prob) {
      cr <- 0.6
    } else if(post_prob >= 0.9) {
      cr <- 0.4
    } else if(post_prob >= 0.85) {
      cr <- 0.45
    } else if(post_prob >= 0.8) {
      cr <- 0.5
    } else { return("0") }
    if (Copy.Ratio > cr) {
      return("+")
    } else if ( Copy.Ratio < (-1 * (cr - .05 )) ){
      return("-")
    } else if (abs(Copy.Ratio) < 0.15){
      return("0")
    } else {
      return("undetermined")
    }
  }
  ## check markers first
  ## is haploseq called?
  if(is.na(MEAN_POSTPROB)){
    checkCopyRatio(post_prob = FALSE, Copy.Ratio)
  } else if (NUM_INFORMATIVE_MARKERS > 40) {
    ## 
    # if ( ( g.Event == "+" & event == "gain") | ( g.Event == "-" & event == "loss" )  ) {
    #   if (pval < 0.05){
    #     
    #   } else { return("0")}
    # } else if ( g.Event == "0" ) {
    #   if(event == "cnloh"){ return("cnloh") }
    #   else if (MEAN_POSTPROB > 0.8 & pval < 0.05) { return("cnloh")}
    #   else { checkCopyRatio(g.Event, Copy.Ratio) }
    # } else { checkCopyRatio(g.Event, Copy.Ratio) }
    call <- checkCopyRatio(post_prob = MEAN_POSTPROB, Copy.Ratio)
    if (call == "0") { 
      if(event == "cnloh"){ return("cnloh") }
      else if (MEAN_POSTPROB > 0.8) { return("cnloh") }
      else { return(call) }
    } else { return(call) }
    
  } else {return("0")}
  
}
test2[, Call := mapply(callTest, as.numeric(MEAN_POSTPROB), g.Event, pval, event, as.numeric(Copy.Ratio), as.numeric(NUM_INFORMATIVE_MARKERS)), by=.(Chr, i.Start)]

write.table(test2[Markers > 20], paste0(dn,"/",samplename,".hap.annotated.segs.tsv"), sep="\t", quote = FALSE, row.names = FALSE)
