library(tidyverse)
library(ggh4x)
library(bracer)

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
 stop('specify five arguments: <filePath> <pathToPlot> ', call.=FALSE)
}

allClusters = data.frame()
fileList=glob(args[1], engine = "r")
for (f in fileList){
myd<-read.table(f, sep = "\t", header = T)
allClusters <-rbind(allClusters, myd)
}

allClusters$chr<-factor(allClusters$chr, levels=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"))
allClusters$length<- factor(allClusters$length, levels=c("whole_chromosome", "p_arm", "q_arm"))
allClusters$method<- factor(allClusters$method, levels=c("PGGB", "Minigraph_Cactus"))
mycolor=c("#11468F","#DA1212", "#FDCA40")
box.colors <- c( "#DDDDDD",  "#EEEEEE")

allClusters %>% mutate(num = ifelse(chr == "chr1", 1, ifelse(chr == "chr2", 2, ifelse(chr == "chr3", 3, ifelse(chr == "chr4", 4, ifelse(chr == "chr5", 5, ifelse(chr ==  "chr6", 6, ifelse(chr == "chr7", 7, ifelse(chr == "chr8", 8, ifelse(chr == "chr9", 9, ifelse(chr == "chr10", 10, ifelse(chr == "chr11", 11, ifelse(chr == "chr12", 12, ifelse(chr == "chr13",13, ifelse(chr == "chr14", 14, ifelse(chr == "chr15", 15, ifelse(chr == "chr16",16, ifelse(chr == "chr17",17, ifelse(chr == "chr18", 18, ifelse(chr == "chr19", 19, ifelse(chr == "chr20", 20, ifelse(chr == "chr21", 21, 22)))))))))))))))))))))) %>% mutate(acro=ifelse(chr == "chr13" | chr == "chr14" | chr == "chr15" | chr == "chr21" | chr == "chr22", "y", "n")) %>% ggplot(aes(num, clusters, color = length))+ geom_rect(aes(xmin=num-0.5, xmax=num+0.5 , ymin=0, ymax=max(clusters) , fill = acro), alpha=0.2, color=NA ) + scale_fill_manual(values = box.colors) + geom_point(size = 3, alpha =0.7) +scale_x_continuous(breaks = seq(from = 1, to = 22, by = 1))  + facet_nested_wrap( ~assembly + method, nrow = 4, strip.position = "left") + theme_bw() + scale_y_continuous(position = "right") + xlab("") + theme(legend.text=element_text(size=13), axis.text=element_text(size=13), axis.title=element_text(size = 12), strip.text.y = element_text(size = 13), legend.title=element_blank()) + scale_color_manual(values = mycolor)

ggsave(args[2], width = 16, heigh = 8)