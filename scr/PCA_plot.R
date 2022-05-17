library(tidyverse)
library(bracer)
library(ggh4x)

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
 stop('specify five arguments: <eigenVecFilesPath> <method> <outputFile> ', call.=FALSE)
}

all = data.frame()
fileList=glob(args[1], engine = "r")
for (f in fileList){
myd<-read.table(f, sep = "\t", header = T, comment = '')
all <-rbind(all,myd)
}

all$chrom<-factor(all$chrom, levels=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"))
all$length<- factor(all$length, levels=c("whole_chromosome", "p_arm", "q_arm"))
mycol<-c("#FF6B6B", "#FFD93D", "#6BCB77", "#4D96FF")
all %>% mutate(acro = ifelse(chrom == "chr13" | chrom == "chr14" | chrom == "chr15" | chrom == "chr21" | chrom == "chr22", "yes", "no")) %>% filter(method == paste(args[2])) %>% ggplot(aes(PC1, PC2, color = Superpopulation)) + geom_point(size = 1, alpha = 0.6) + facet_nested(chrom~length + assembly) + theme_bw() + theme(legend.text=element_text(size=13), axis.text.x=element_text(size=9), axis.text.y = element_text(size = 9) ,axis.title=element_text(size = 12), strip.text.x = element_text(size = 13), strip.text.y = element_text(size = 12) ,legend.title=element_blank()) + scale_color_manual(values = mycol)
ggsave(args[3], width = 8, heigh = 25)