#statistical exploration of the data geenrated from PPI networks

# Load required R packages
library(tidyverse)
library(rstatix)
library(ggpubr)
library(effectsize)


#pairwise comparison
pwc <- df_adj_transformed %>% pairwise_t_test(entropy_500 ~ superkingdom, p.adjust.method = "bonferroni")
pwc

cohens_d(df_adj_transformed$entropy_500[which(df_adj_transformed$superkingdom=="Bacteria")], df_adj_transformed$entropy_500[which(df_adj_transformed$superkingdom=="Archaea")])
cohens_d(df_adj_transformed$entropy_500[which(df_adj_transformed$superkingdom=="Bacteria")], df_adj_transformed$entropy_500[which(df_adj_transformed$superkingdom=="Eukaryota")])
cohens_d(df_adj_transformed$entropy_500[which(df_adj_transformed$superkingdom=="Archaea")], df_adj_transformed$entropy_500[which(df_adj_transformed$superkingdom=="Eukaryota")])

res.aov <- df_adj_transformed %>% anova_test(entropy_500 ~ superkingdom)
df_adj_transformed$logentropy500=log(df_adj_transformed$entropy_500)
colorss=c("indianred3","gray70","dodgerblue3")
res.aov
pwc <- pwc %>% add_xy_position(x = "Superkingdom")
ggboxplot(df_adj_transformed, x = "superkingdom", y = "logentropy500", fill=colorss) +
  stat_pvalue_manual(pwc, label = "p.adj", tip.length = 0) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )

colorss=c("gray70","indianred3","dodgerblue3")
boxplot(df_adj_transformed$entropy_500~df_adj_transformed$superkingdom, col=colorss)


#phylogenetic placement analysis
library(ape)
library(phytools)
library(dplyr)
library(ggplot2)
library(ggtree)

#Load Tree and Data
#setwd("~/hugtree")
ht <- read.tree("~/Documents/Research/Manuscripts/PNAS eukaryotes/scripts/HugTree.newick")
#ce <- read.csv("ce_evo.csv",stringsAsFactors = F)

df_new3=df_new2[which(df_new2$entropy_500>0.7),]
df_new3=df_new3[which(df_new3$Domain=="Bacteria"),]

###########################################
# Check Coverage of Hug Tree W/ CE Data ###
###########################################

# How much overlap between datasets?
ce_names <- gsub(" ","_",df_new3$Name)
tips <- ht$tip.label
x <- lapply(ce_names,grepl,x=tips) %>% lapply(sum) %>% unlist()
table(x>0)

# xt <- lapply(ce_names,grep,x=tips) %>% unlist()

#Side by side Hug tree and Hug tree pruned to only CE spp
pdf("HugTreeVsPrunedTree.pdf",width=8,height=5)
par(mfrow=c(1,2))
plot(ht,type="fan",show.tip.label = F)
plot(keep.tip(ht,
              unlist(lapply(ce_names,grep,x=ht$tip.label))),
     type="fan",show.tip.label = F)
par(mfrow=c(1,1))
dev.off()
library(ggtree)
#plot which species on hug tree have CE data with ggtree
y <- ht$tip.label %in% ht$tip.label[unlist(lapply(ce_names,grep,x=ht$tip.label))]
d <- data.frame(id=ht$tip.label,HasCE=y)
p <- ggtree(ht,layout="circular")

# Format data for visualizing on tips
ht_trimmed <- keep.tip(ht,unlist(lapply(ce_names,grep,x=ht$tip.label)))
ind <- lapply(ce_names,grep,x=ht_trimmed$tip.label)
names(ind) <- as.character(1:length(ind))
ind <- unlist(ind)
match_ce <- data.frame(Tip=ht_trimmed$tip.label[ind],
                       Name=ce_names[as.numeric(names(ind))],
                       stringsAsFactors = F)
rownames(df_new2) <- gsub(" ","_",df_new2$Name)
d <- data.frame(id=match_ce$Tip,
                SpectralEntropy=df_new3[match_ce$Name,"entropy_500"],
                EIMicro=df_new3[match_ce$Name,"effectiveness"],
                Domain=df_new3[match_ce$Name,"Domain"])

#plot EI Micro 
#p <- ggtree(ht_trimmed,layout="daylight", color="grey39", size=0.3)
p <- ggtree(ht_trimmed,layout="circular", color="grey39", size=0.3)
p <- p %<+% d + geom_tippoint(aes(color=SpectralEntropy,shape=Domain)) + scale_color_distiller(palette = "Spectral")
pdf("Entropy500bOnHugTree.pdf",width=8,height=8)
plot(p)
dev.off()


#caculating GC vals
#biomartr, magrittr
library(R.utils)
organism = namepb$ID

for (i in 1200:length(namepb$ID)){
  getGenome(db = "refseq", namepb$ID[i], reference = FALSE, release = NULL, gunzip = FALSE, path = file.path("_ncbi_downloads", "genomes"))
}

dataFiles <- lapply(Sys.glob("doc_*_db_refseq_summary_statistics.tsv"), read.table)

library(plyr)
library(readr)


genomedata= matrix(0, nrow=length(dataFiles), ncol= 9)
colnames(genomedata)= c("ID", "genome_length_mbp", "50_mbp", "n_seqs", "n_nnn", "rel_nnn", "genome_entropy", "n_gc", "rel_gc")
genomedata=as.data.frame(genomedata)

for (i in 1:length(dataFiles)) {
  genomedata$ID[i]= as.character(dataFiles[[i]]$V1[2])
  genomedata$genome_length_mbp[i]= as.character(dataFiles[[i]]$V2[2])
  genomedata$`50_mbp`[i]= as.character(dataFiles[[i]]$V3[2])
  genomedata$n_seqs[i]= as.character(dataFiles[[i]]$V4[2])
  genomedata$n_nnn[i]= as.character(dataFiles[[i]]$V5[2])
  genomedata$rel_nnn[i]= as.character(dataFiles[[i]]$V6[2])
  genomedata$genome_entropy[i]= as.character(dataFiles[[i]]$V7[2])
  genomedata$n_gc[i]= as.character(dataFiles[[i]]$V8[2])
  genomedata$rel_gc[i]= as.character(dataFiles[[i]]$V9[2])
}

namepb = namep[which(namep$Domain=="Bacteria"),]

genomedata$ID= as.integer(genomedata$ID)
genomedata$`50_mbp`= as.integer(genomedata$`50_mbp`)
genomedata$genome_length_mbp= as.integer(genomedata$genome_length_mbp)


library(dplyr)
genomeica <- full_join(genomedata, namepb, by = 'ID')


library(R.utils)
filegen = list.files(path="genomes", pattern="*_genomic_refseq.fna.gz", full.names=T)

filegen2 = list.files(path="genomes", pattern="*_genomic_refseq.fna.gz", full.names=F)

GCval=NULL

for (i in 1:length(filegen)){
  nam=paste(filegen[i])
  nam0=paste(filegen2[i])
  NN=gunzip(nam, remove=F)
  NNseq=read.fasta(NN)[[1]]
  nam2=strsplit(nam0, "_")[[1]][1]
  
  each=c(nam2,GC(NNseq))
  GCval=rbind(GCval,each)
  print(i)
}


