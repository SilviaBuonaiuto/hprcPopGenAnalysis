library(tidyverse)
library(factoextra)
library(bracer)
library(magrittr)

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
 stop('specify five arguments: <filePath> <length> <assembly> <method> <outputFile> ', call.=FALSE)
}

##### Calculate optimal cluster
k_opt<- function(df) {
	gap_stat<-fviz_nbclust(df, kmeans, method = 'gap_stat')
	firstlocalmax_position <- which(diff(gap_stat$data$gap)<0)[1]
	firstlocalmax_gapvalue <-gap_stat$data[firstlocalmax_position, ]$gap
	firstlocalmax_SE.sim <-gap_stat$data[firstlocalmax_position, ]$SE.sim
	gapmin <- firstlocalmax_gapvalue - firstlocalmax_SE.sim
	gapmax <- firstlocalmax_gapvalue + firstlocalmax_SE.sim
	optimal_k<-gap_stat$data[1:firstlocalmax_position,]  %>% filter( between (gap, gapmin, gapmax )) %>% select(clusters) %>% filter(clusters == min(as.numeric(clusters)))
	return(optimal_k)
}


#### Create df with optimal cluster number for each chromosome and concatenate
new = data.frame()
fileList=glob(args[1], engine = "r")
for (f in fileList){
myd<-read.table(f, sep = "\t", header = T, comment = '')
mydB<- myd %>% select(PC1, PC2)
mydO<- k_opt(mydB)
chrom = strsplit(basename(f), ".", fixed = T) %>% sapply(extract2, 3)
all<-mydO %>% mutate(chr = paste(chrom), length = paste(args[2]), assembly = paste(args[3]), method = paste(args[4]))
new <-rbind(new,all)
}

new %>% write.table(args[5], sep = "\t", quote = F, col.names = T, row.names = F)