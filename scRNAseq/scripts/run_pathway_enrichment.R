#' ============================================================================
#' KEGG Pathway Enrichment Visualization
#' ============================================================================
#'
#' Manuscript: Multi-omics Reveals How Diet and Exercise Improve Peripheral Neuropathy
#' Authors:    Eid SA, Guo K, Noureldein MH, Jang D-G, Allouch AM, Miller CM, Kiriluk CP,
#   Lentz W, Boutary S, Mule JJ, Pennathur S, Bhatt DK, Feldman EL, Hur J
#'
#' Contact:
#   Dr. Kai Guo        <guokai8@gmail.com>
#   Dr. Junguk Hur     <junguk.hur@med.UND.edu>
#'
#' Description:
#'   Visualizes KEGG pathway enrichment results as bubble plots for Schwann cell
#   subtypes (mySC, nmSC, ImmSC) and other cell types. Bubble size reflects rich
#   factor and color indicates -log10(p-value).
#'
#' Input:
#'   KEGG enrichment CSV files (*KEGG_enrichment_pval.csv) per cell type
#'
#' Output:
#'   PDF enrichment bubble plots per cell type and comparison
#'
#' Related figures:
#'   Fig. 3C-F, Fig. S5-S7, Tables S5-S7
#'
#' ============================================================================
setwd("Desktop/Project/Step/final/enrichmentp/")
library(tidyverse)
filenames<-list.files(pattern=".*HFD.*KEGG_enrichment_pval.csv")
d<-lapply(filenames, function(x)read.csv(x,row.names=1))
names(d)<-sub('_KEGG_enrichment_pval.csv','',filenames)
####mySC
mySC<-do.call(rbind,d[grep('mySC',names(d))])
mySC$Group<-sub('\\_.*','',rownames(mySC))
ggplot(mySC,aes(Group,Term,color=-log10(Pvalue),size=Significant/Annotated))+geom_point()+scale_color_gradient(low="pink",high="red")+theme_minimal()+theme(axis.text.x=element_text(angle=90,vjust = 0.5,hjust = 1,size=4),axis.text.y=element_text(size=4))+xlab("")+labs(size="RichFactor")
dev.print(pdf,file="mySC_enrichment_p.pdf")
######
nmSC<-do.call(rbind,d[grep('nmSC',names(d))])
nmSC$Group<-sub('\\_.*','',rownames(nmSC))
ggplot(nmSC,aes(Group,Term,color=-log10(Pvalue),size=Significant/Annotated))+geom_point()+scale_color_gradient(low="pink",high="red")+
  theme_minimal()+theme(axis.text.x=element_text(angle=90,vjust = 0.5,hjust = 1),axis.text.y=element_text(size=10))+xlab("")+labs(size="RichFactor")+ylab("")
dev.print(pdf,file="nmSC_enrichment_p.pdf")
#########
Mac<-do.call(rbind,d[grep('Mac',names(d))])
Mac$Group<-sub('\\_.*','',rownames(Mac))
ggplot(Mac,aes(Group,Term,color=-log10(Pvalue),size=Significant/Annotated))+geom_point()+scale_color_gradient(low="pink",high="red")+
  theme_minimal()+theme(axis.text.x=element_text(angle=90,vjust = 0.5,hjust = 1),axis.text.y=element_text(size=9))+xlab("")+labs(size="RichFactor")+ylab("")
dev.print(pdf,file="Mac_enrichment_p.pdf")
########
ImmSC<-do.call(rbind,d[grep('ImmSC',names(d))])
ImmSC$Group<-sub('\\_.*','',rownames(ImmSC))
ggplot(ImmSC,aes(Group,Term,color=-log10(Pvalue),size=Significant/Annotated))+geom_point()+scale_color_gradient(low="pink",high="red")+
  theme_minimal()+theme(axis.text.x=element_text(angle=90,vjust = 0.5,hjust = 1),axis.text.y=element_text(size=9))+xlab("")+labs(size="RichFactor")+ylab("")
dev.print(pdf,file="ImmSC_enrichment_p.pdf")
#######
Perineurial<-do.call(rbind,d[grep('Perineurial',names(d))])
Perineurial$Group<-sub('\\_.*','',rownames(Perineurial))
ggplot(Perineurial,aes(Group,Term,color=-log10(Pvalue),size=Significant/Annotated))+geom_point()+scale_color_gradient(low="pink",high="red")+
  theme_minimal()+theme(axis.text.x=element_text(angle=90,vjust = 0.5,hjust = 1),axis.text.y=element_text(size=9))+xlab("")+labs(size="RichFactor")+ylab("")
dev.print(pdf,file="Perineurial_enrichment_p.pdf")
########
#######
EndoFib1<-do.call(rbind,d[grep('EndoFib1',names(d))])
EndoFib1$Group<-sub('\\_.*','',rownames(EndoFib1))
ggplot(EndoFib1,aes(Group,Term,color=-log10(Pvalue),size=Significant/Annotated))+geom_point()+scale_color_gradient(low="pink",high="red")+
  theme_minimal()+theme(axis.text.x=element_text(angle=90,vjust = 0.5,hjust = 1),axis.text.y=element_text(size=9))+xlab("")+labs(size="RichFactor")+ylab("")
dev.print(pdf,file="EndoFib1_enrichment_p.pdf")
########
EndoFib2<-do.call(rbind,d[grep('EndoFib2',names(d))])
EndoFib2$Group<-sub('\\_.*','',rownames(EndoFib2))
ggplot(EndoFib2,aes(Group,Term,color=-log10(Pvalue),size=Significant/Annotated))+geom_point()+scale_color_gradient(low="pink",high="red")+
  theme_minimal()+theme(axis.text.x=element_text(angle=90,vjust = 0.5,hjust = 1),axis.text.y=element_text(size=9))+xlab("")+labs(size="RichFactor")+ylab("")
dev.print(pdf,file="EndoFib2_enrichment_p.pdf")
########
SMC<-do.call(rbind,d[grep('SMC',names(d))])
SMC$Group<-sub('\\_.*','',rownames(SMC))
ggplot(SMC,aes(Group,Term,color=-log10(Pvalue),size=Significant/Annotated))+geom_point()+scale_color_gradient(low="pink",high="red")+
  theme_minimal()+theme(axis.text.x=element_text(angle=90,vjust = 0.5,hjust = 1),axis.text.y=element_text(size=3))+xlab("")+labs(size="RichFactor")+ylab("")
dev.print(pdf,file="SMC_enrichment_p.pdf")
######
########
VSMC<-do.call(rbind,d[grep('VSMC',names(d))])
VSMC$Group<-sub('\\_.*','',rownames(VSMC))
ggplot(VSMC,aes(Group,Term,color=-log10(Pvalue),size=Significant/Annotated))+geom_point()+scale_color_gradient(low="pink",high="red")+
  theme_minimal()+theme(axis.text.x=element_text(angle=90,vjust = 0.5,hjust = 1),axis.text.y=element_text(size=5))+xlab("")+labs(size="RichFactor")+ylab("")
dev.print(pdf,file="VSMC_enrichment_p.pdf")
######EpiFib
EpiFib<-do.call(rbind,d[grep('EpiFib',names(d))])
EpiFib$Group<-sub('\\_.*','',rownames(EpiFib))
ggplot(EpiFib,aes(Group,Term,color=-log10(Pvalue),size=Significant/Annotated))+geom_point()+scale_color_gradient(low="pink",high="red")+
  theme_minimal()+theme(axis.text.x=element_text(angle=90,vjust = 0.5,hjust = 1),axis.text.y=element_text(size=5))+xlab("")+labs(size="RichFactor")+ylab("")
dev.print(pdf,file="EpiFib_enrichment_p.pdf")
####Pericytes
Pericytes<-do.call(rbind,d[grep('Pericytes',names(d))])
Pericytes$Group<-sub('\\_.*','',rownames(Pericytes))
ggplot(Pericytes,aes(Group,Term,color=-log10(Pvalue),size=Significant/Annotated))+geom_point()+scale_color_gradient(low="pink",high="red")+
  theme_minimal()+theme(axis.text.x=element_text(angle=90,vjust = 0.5,hjust = 1),axis.text.y=element_text(size=8))+xlab("")+labs(size="RichFactor")+ylab("")
dev.print(pdf,file="Pericytes_enrichment_p.pdf")
#######
SC3<-do.call(rbind,d[grep('SC3',names(d))])
SC3$Group<-sub('\\_.*','',rownames(SC3))
ggplot(SC3,aes(Group,Term,color=-log10(Pvalue),size=Significant/Annotated))+geom_point()+scale_color_gradient(low="pink",high="red")+
  theme_minimal()+theme(axis.text.x=element_text(angle=90,vjust = 0.5,hjust = 1),axis.text.y=element_text(size=6))+xlab("")+labs(size="RichFactor")+ylab("")
dev.print(pdf,file="SC3_enrichment_p.pdf")
#####Endo
Endo<-do.call(rbind,d[grep('Endo',names(d))])
Endo$Group<-sub('\\_.*','',rownames(Endo))
ggplot(Endo,aes(Group,Term,color=-log10(Pvalue),size=Significant/Annotated))+geom_point()+scale_color_gradient(low="pink",high="red")+
  theme_minimal()+theme(axis.text.x=element_text(angle=90,vjust = 0.5,hjust = 1),axis.text.y=element_text(size=6))+xlab("")+labs(size="RichFactor")+ylab("")
dev.print(pdf,file="Endo_enrichment_p.pdf")


