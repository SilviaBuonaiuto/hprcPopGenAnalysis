library(tidyverse)
library(bracer)
library(ggh4x)
library(gridExtra)
library (grid)

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
all$length<- factor(all$length, levels=c("whole_chromosome", "q_arm", "p_arm"))
mycol<-c("#FBCB0A", "#F32424", "#5FD068", "#612897")
p1<-all %>% mutate(acro = ifelse(chrom == "chr13" | chrom == "chr14" | chrom == "chr15" | chrom == "chr21" | chrom == "chr22", "yes", "no")) %>% filter(method == paste(args[2])) %>% filter(chrom == "chr1" | chrom == "chr2"| chrom =="chr3"| chrom=="chr4"| chrom=="chr5" |chrom=="chr6"| chrom=="chr7"| chrom=="chr8"| chrom=="chr9"| chrom=="chr10"| chrom=="chr11") %>% ggplot(aes(PC1, PC2, color = Superpopulation)) + geom_point(size = 1, alpha = 0.6) + facet_nested(chrom~length + assembly) + theme_bw() + theme(legend.text=element_text(size=13), axis.text.x=element_text(size=9), axis.text.y = element_text(size = 9) ,axis.title=element_text(size = 12), strip.text.x = element_text(size = 10), strip.text.y = element_text(size = 11) ,legend.title=element_blank()) + scale_color_manual(values = mycol) + xlim(-0.6, 0.2)
p2<-all %>% mutate(acro = ifelse(chrom == "chr13" | chrom == "chr14" | chrom == "chr15" | chrom == "chr21" | chrom == "chr22", "yes", "no")) %>% filter(method == "PGGB") %>% filter(chrom == "chr12"| chrom =="chr13"| chrom=="chr14"| chrom=="chr15" |chrom=="chr16"| chrom=="chr17"| chrom=="chr18"| chrom=="chr19"| chrom=="chr20"| chrom=="chr21"| chrom == "chr22") %>% ggplot(aes(PC1, PC2, color = Superpopulation)) + geom_point(size = 1, alpha = 0.6) + facet_nested(chrom~length + assembly) + theme_bw() + theme(legend.text=element_text(size=13), axis.text.x=element_text(size=9), axis.text.y = element_text(size = 9) ,axis.title=element_text(size = 12), strip.text.x = element_text(size = 10), strip.text.y = element_text(size = 11) ,legend.title=element_blank()) + scale_color_manual(values = mycol) + guides(color = guide_legend(direction = "horizontal")) +  xlim(-0.6, 0.2)

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
mylegend<-g_legend(p2)
p3 <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"), p2 + theme(legend.position="none"), nrow=1), mylegend,nrow = 2, heights=c(10, 1))

ggsave(args[3], plot = p3, width = 10, heigh = 10)
