library(splitstackshape)
library(GenomicRanges)
library(ggplot2)
library(plyr)
library(RColorBrewer)
library(ggsignif)

#### Set Working Directory, Input Files & Variables #####################################################################
argv = commandArgs(trailingOnly=TRUE)
C1 = argv[1L]
C2 = argv[2L]
domains.file.C1 = argv[3L] #domains of the condition 1 (HiC-Bench) 
domains.file.C2 = argv[4L] #domains of the condition 2 (HiC-Bench)
ratio.scores.file = argv[5L] #ratio scores file (HiC-Bench)

padj.cut=0.05 #padjusted cutoff to determine boundaries with significant insulation changes
lfc.cut=0.5 #absolute log2FC cutoff of the Mean Ratio Boundary Score (condition1/condition2) to determine boundaries with significant insulation changes 

#### Set Conditions and Replicates ###########################################################################

## Get the Exact Sample Labels ##
ratio.scores.tab <- read.delim(ratio.scores.file,header = T,stringsAsFactors = F)
labels=colnames(ratio.scores.tab)[2:dim(ratio.scores.tab)[2]]

conditions=c(C1,C2) #Pick the 2 conditions. The order is: condition 1, Condition 2. This determines the way the log2FC will be computed [log2(condition1/condition2)]

#### START ###################################################################################################
samples=labels[grep(paste0(conditions[1],"|",conditions[2]),labels)]
samples.C1 = samples[grep(conditions[1],samples)]
samples.C2 = samples[grep(conditions[2],samples)]

num.reps=length(samples)/2
n=num.reps-1

options("scipen"=100, "digits"=4)

#### Get Boundaries in each condition and merge them all in one table #######################################
boundaries.merged = data.frame()

group=1
for (domains.i.path in c(domains.file.C1,domains.file.C2)){
  print(paste0("Proccessing: ", conditions[group]))
  
  ## Get Domains (HiC-Bench) ##
  print("Getting domains... ")
  
  domains <- read.delim(domains.i.path,header = F,stringsAsFactors = F)
  domains$V2 <- domains$V2+1
  names(domains)  <- c("chr", "start", "end")
  domains <- domains[!is.na(domains$chr),]
  domains <- domains[order(domains$chr,domains$start),]
  dim(domains)
  
  ## Get Boundaries ##
  print("Getting boundaries... ")
  boundaries.i <- data.frame()
  
  for (i in 1:(nrow(domains)-1)){
    
    if (domains$chr[i+1] == domains$chr[i]){
      temp <- data.frame(chr=character(1),start=numeric(1),end=numeric(1),boundary_ID=character(1),stringsAsFactors = F)
      temp$chr <- domains$chr[i]
      temp$start <- domains$end[i]
      temp$end <- domains$start[i+1]
      boundaries.i <- rbind(boundaries.i,temp)
    }
  }
  boundaries.i$boundary_ID <- paste(boundaries.i$chr,as.numeric(boundaries.i$start),as.numeric(boundaries.i$end),sep = ":")
  boundaries.i$condition <- conditions[group]
  group=group+1
  
  print(paste0("Number of boundaries: ",dim(boundaries.i)[1]))
  boundaries.merged=rbind(boundaries.merged,boundaries.i)
  
}
print(paste0("Number of boundaries after merging: ",dim(boundaries.merged)[1]))
boundaries.merged$condition[duplicated(boundaries.merged$boundary_ID)] <- "BOTH"

#Remove Duplicated Boundaries (present in both conditions)#
boundaries.merged = boundaries.merged[!duplicated(boundaries.merged$boundary_ID),]
print(paste0("Number of boundaries after removing duplicated boundaries: ",dim(boundaries.merged)[1]))

#### Create Master Boundaries Dataframe w/ 1 ins.score column by replicate ###################################

master.df <- as.data.frame(matrix(nrow=length(boundaries.merged$chr),ncol=length(conditions)+8))
colnames(master.df) <- c(colnames(boundaries.merged),conditions,"log2FC","pvalue","padjusted")
master.df[,1:5] <- boundaries.merged
master.df.GR <-  makeGRangesFromDataFrame(master.df,keep.extra.columns = T)

#### Prepare Ratio Insulation Scores Table ###################################################################
colnames(ratio.scores.tab)[1] <- "boundary.ID"
ratio.scores.tab = subset(ratio.scores.tab, select = c("boundary.ID",samples)) #Remove the columns of the unused conditions 
ratio.scores.tab$boundary.ID <- gsub(ratio.scores.tab$boundary.ID,pattern = "-",replacement = ":")
ratio.scores.tab = cSplit(ratio.scores.tab,splitCols = "boundary.ID",sep = ":",drop = F)
colnames(ratio.scores.tab)[(2+2*num.reps):(2+2*num.reps+2)] <- c("chr","start","end")
ratio.scores.tab$C1.mean <- rowMeans(subset(ratio.scores.tab, select = samples.C1)) #Get mean ratio score of each bin-score in condition1
ratio.scores.tab$C2.mean <- rowMeans(subset(ratio.scores.tab, select = samples.C2)) #Get mean ratio score of each bin-score in condition2
ratio.scores.tab.GR <-  makeGRangesFromDataFrame(ratio.scores.tab,keep.extra.columns = T)

#### Calculate Mean (Ratio) Insulation Score for each boundary for each replicate and perform a TTEST #########
print("Compute Mean (Ratio) Insulation Scores and perform a differential analysis (ttest)...")
min.ovl <- 30000

for (i in 1:dim(master.df)[1]){ ## MASTER LOOP ##
if(i %% 500==0) {print(paste0(i," of ",dim(master.df)[1], " boundaries processed..."))}

ovl.Bound.scores <- subsetByOverlaps(ratio.scores.tab.GR,master.df.GR[i],minoverlap = min.ovl) #Get all the ins.scores bins that fall in the boundary

if (length(ovl.Bound.scores) > 0){  #Check if at least 1 ratio score falls in the boundary coordinates

ovl.Bound.scores.df <- as.data.frame(ovl.Bound.scores)
master.df[i,6] <- mean(colMeans(ovl.Bound.scores.df[,7:(7+n)])) #Compute the Mean Boundary Ratio Score in the condition 1
master.df[i,7] <- mean(colMeans(ovl.Bound.scores.df[,(7+n+1):(7+2*n+1)]))  #Compute the Mean Boundary Ratio Score in the condition 2
master.df$log2FC[i] <- log2(master.df[i,6]/master.df[i,7]) #Compute the log2 fold change (condition1/condition2)

check=try(t.test(x=ovl.Bound.scores.df$C1.mean,y=ovl.Bound.scores.df$C2.mean,paired = T),silent = T) #Check if we have enough data to perform a ttest

if(length(check) != 1){
ttest = t.test(x=ovl.Bound.scores.df$C1.mean,y=ovl.Bound.scores.df$C2.mean,paired = F) #Perform a pooled paired-ttest
master.df$pvalue[i] <- ttest$p.value

} else {master.df$pvalue[i] <- NA}
} else {master.df[i,5:dim(master.df)[2]] <- NA}

master.df$padjusted <- p.adjust(master.df$pvalue,method="BH" ,n=dim(master.df[!is.na(master.df$pvalue),])[1]) #Compute FDR values
master.df[,6:dim(master.df)[2]] <- apply(master.df[,6:dim(master.df)[2]],2,as.numeric )
}
print("Done")

#### Clean Master Dataframe ###############################################################################
master.df.clean <- master.df[!is.na(master.df$padjusted),]
master.df.clean <- master.df.clean[!is.na(master.df.clean$log2FC),]
master.df.clean <- master.df.clean[!is.na(master.df.clean$chr),]

#### Volcano Plot (FDR/MBS) ###############################################################################

## Set cutoffs ## (Use your own criteria)
padj.cut=0.05
lfc.cut=0.5

master.df.clean$insulation.status <- ""
master.df.clean$insulation.status[master.df.clean$padjusted < padj.cut & master.df.clean$log2FC > lfc.cut] <- "increased"
master.df.clean$insulation.status[master.df.clean$padjusted < padj.cut & master.df.clean$log2FC < -lfc.cut] <- "decreased"
master.df.clean$insulation.status[master.df.clean$padjusted >= padj.cut | abs(master.df.clean$log2FC) <= lfc.cut] <- "stable"

write.table(master.df.clean, file = paste0("diff.boundaries_master_",conditions[1],"_vs_",conditions[2],".tsv"), col.names = T, quote = F, sep = "\t")

## Get Boundaries with significant changes ##
boundaries_up <- master.df.clean[master.df.clean$padjusted < padj.cut & master.df.clean$log2FC > lfc.cut,]
boundaries_down <- master.df.clean[master.df.clean$padjusted < padj.cut & master.df.clean$log2FC < -lfc.cut,]

print(paste0("Number of boundaries with significant increased insulation: ", dim(boundaries_up)[1]))
print(paste0("Number of boundaries with  significant decreased insulation: ", dim(boundaries_down)[1]))

write.table(boundaries_up[,1:3], file = paste0("diff.boundaries_increased_",conditions[1],"_vs_",conditions[2],".bed"),row.names = F, col.names = F, quote = F, sep = "\t")
write.table(boundaries_down[,1:3], file = paste0("diff.boundary_decreased_",conditions[1],"_vs_",conditions[2],".bed"),row.names = F, col.names = F, quote = F, sep = "\t")

## Get Volcano Plot ##
min.lfc=min(boundaries_down$log2FC,na.rm = T)*1.1
max.lfc=max(boundaries_up$log2FC,na.rm = T)*1.1
max.padj=-log10(min(c(boundaries_down$padjusted,boundaries_up$padjusted),na.rm = T))*1.05

png(file=paste0("FDR_vs_MBS_volcano_",conditions[1],"_vs_",conditions[2],".png"),width=2048, height=2048, pointsize=100)
plot(master.df.clean$log2FC, -log10(as.numeric(master.df.clean$padjusted)),main=paste0("Cutoffs: log2FC=",lfc.cut," ; padjusted=",padj.cut),xlab = "log2FC(MBS.Ratio)",cex.main=0.5,cex.lab=0.7, ylab="-log10(FDR)",xlim=c(min.lfc,max.lfc),ylim=c(0,max.padj),col="grey",pch=16)
box(lwd=4)
axis(side = 1, lwd = 4)
axis(side = 2, lwd = 4)
points(boundaries_down$log2FC[boundaries_down$log2FC <= -lfc.cut], -log10(as.numeric(boundaries_down$padjusted[boundaries_down$log2FC <= -lfc.cut])) ,col="blue",pch=20,cex=1.5)
points(boundaries_up$log2FC[boundaries_up$log2FC >= lfc.cut], -log10(as.numeric(boundaries_up$padjusted[boundaries_up$log2FC >= lfc.cut])),col="darkred",pch=20,cex=1.5)
abline(h=c(-log10(padj.cut)), col="red", lwd=4, lty=2)
abline(v=lfc.cut,col="red",lwd=4,lty=2) 
abline(v=-lfc.cut,col="red",lwd=4,lty=2) 
dev.off()

## Boxplot [log2FC (MBS.ratio) by insulation status] ##

master.df.clean$condition <- as.factor(master.df.clean$condition)

p_meds <- ddply(master.df.clean, .(insulation.status), summarise, med = median(log2FC))
p_meds$med <- round(p_meds$med,digits = 3)

pdf(file=paste0("boxplot_MBS.ratio_",conditions[1],"_vs_",conditions[2],".pdf"),width=14, height=7)

ggplot(master.df.clean, aes(x=insulation.status, y=log2FC, group=insulation.status, fill=insulation.status))+
xlab(label = "Insulation Status") + 
geom_boxplot(width=0.8,outlier.size = 0.1,alpha=0.6)+
ylab(label = paste0("log2FC (MBS.Ratio) [",conditions[1]," / ",conditions[2],"]"))+
ylim(c(min(master.df.clean$log2FC)*1.1,max(master.df.clean$log2FC)*1.1))+
geom_text(data = p_meds, aes(x = insulation.status, y = med, label = med), size = 2.8, vjust = -1)+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(hjust = 0.5))+
geom_signif(comparisons = list(c("increased","decreased"),c("increased","stable"),c("decreased","stable")), map_signif_level=TRUE,test = 'wilcox.test',na.rm = T,y_position = max(master.df.clean$log2FC)*1.1)

dev.off()

