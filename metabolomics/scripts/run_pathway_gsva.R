#' ============================================================================
#' GSVA Pathway Scoring for Metabolomics
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
#'   GSVA scoring of metabolic pathways followed by Wilcoxon testing and PLS-DA.
#   Generates VIP-scored significant pathways across all pairwise comparisons
#   using Gaussian kernel CDF with minimum pathway size of 2.
#'
#' Input:
#'   data_prep.rdata (preprocessed metabolomics data), KEGG pathway annotations
#'
#' Output:
#'   Pathway statistical results, Venn diagrams, PLS-DA models, VIP scores
#   with fold-changes
#'
#' Related figures:
#'   Fig. 4A-B, Fig. S10
#'
#' ============================================================================
library(scGSVA)
load("data_prep.rdata")
annot<-path
annot[,3]<-annot[,2]
ddp<-lapply(dd, function(x)scgsva(as.matrix(t(x)),annot = annot,kcdf="Gaussian",min.sz = 2))
ddpd<-lapply(ddp, function(x)x@gsva)
####
library(rstatix)
library(tidyverse)
ddpd$Acyl<-ddpd$Acyl[.1:2]
scnwp<-ddpd$SCN%>%rownames_to_column(var="id")%>%gather(met,val,-id)%>%mutate(Group=metas[id,"Group"])%>%
  group_by(met)%>%wilcox_test(val~Group)%>%adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")
acywp<-ddpd$Acyl%>%rownames_to_column(var="id")%>%gather(met,val,-id)%>%mutate(Group=metaa[id,"Group"])%>%
  group_by(met)%>%wilcox_test(val~Group)%>%adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")
untarwp<-ddpd$Untargeted%>%rownames_to_column(var="id")%>%gather(met,val,-id)%>%mutate(Group=metau[id,"Group"])%>%
  group_by(met)%>%wilcox_test(val~Group)%>%adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")
polarwp<-ddpd$Polar%>%rownames_to_column(var="id")%>%gather(met,val,-id)%>%mutate(Group=metap[id,"Group"])%>%
  group_by(met)%>%wilcox_test(val~Group)%>%adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")
#####################################
write.csv(scnwp,file="SCN_wilox_all_pathway.csv")
write.csv(acywp,file="Acyl_wilox_pathway.csv")
write.csv(untarwp,file="Untargeted_wilox_all_pathway.csv")
write.csv(polarwp,file="Polar_wilox_pathway.csv")
######
make.comparison<-function(group,ref=NULL){
  if(is.null(ref)){
    ref=group
  }
  tmp<-expand.grid(group,ref,stringsAsFactors = F)
  tmp<-tmp[tmp$Var1!=tmp$Var2,]
  tmp <- tmp[!duplicated(t(apply(tmp, 1, sort))), ]
  return(tmp)
}
comp<-make.comparison(unique(metas$Group),ref=c("WT_18WK","HFD_18WK","WT_26WK","HFD_26WK"))
#####
scnmp<-ddpd$SCN%>%rownames_to_column(var="id")%>%gather(met,val,-id)%>%mutate(Group=metas[id,"Group"])%>%group_by(met,Group)%>%summarise(mu=mean(val))%>%spread(Group,mu)
acymp<-ddpd$Acyl%>%rownames_to_column(var="id")%>%gather(met,val,-id)%>%mutate(Group=metaa[id,"Group"])%>%group_by(met,Group)%>%summarise(mu=mean(val))%>%spread(Group,mu)
polarmp<-ddpd$Polar%>%rownames_to_column(var="id")%>%gather(met,val,-id)%>%mutate(Group=metap[id,"Group"])%>%group_by(met,Group)%>%summarise(mu=mean(val))%>%spread(Group,mu)
untmp<-ddpd$Untargeted%>%rownames_to_column(var="id")%>%gather(met,val,-id)%>%mutate(Group=metau[id,"Group"])%>%group_by(met,Group)%>%summarise(mu=mean(val))%>%spread(Group,mu)
#######
scnfp<-cbind(scnmp,do.call(cbind,apply(comp, 1, function(x){
  nam<-paste0(x[1],"vs",x[2])
  res<-scnmp[,x[1]]-scnmp[,x[2]]
  names(res)<-nam
  return(res)
})))
####
acyfp<-cbind(acymp,do.call(cbind,apply(comp, 1, function(x){
  nam<-paste0(x[1],"vs",x[2])
  res<-log2(acymp[,x[1]]/acymp[,x[2]])
  names(res)<-nam
  return(res)
})))
###
polarfp<-cbind(polarmp,do.call(cbind,apply(comp, 1, function(x){
  nam<-paste0(x[1],"vs",x[2])
  res<-log2(polarmp[,x[1]]/polarmp[,x[2]])
  names(res)<-nam
  return(res)
})))
####
untfp<-cbind(untmp,do.call(cbind,apply(comp, 1, function(x){
  nam<-paste0(x[1],"vs",x[2])
  res<-untmp[,x[1]]-untmp[,x[2]]
  names(res)<-nam
  return(res)
})))
########
scnprp<-scnwp%>%mutate(comparision=paste0(group1,"vs",group2))%>%filter(group2%in%c("WT_18WK","HFD_18WK","WT_26WK","HFD_26WK"))%>%left_join(scnfp,by=c("met"="met"))
acyrp<-acywp%>%mutate(comparision=paste0(group1,"vs",group2))%>%filter(group2%in%c("WT_18WK","HFD_18WK","WT_26WK","HFD_26WK"))%>%left_join(acyfp,by=c("met"="met"))
polarrp<-polarwp%>%mutate(comparision=paste0(group1,"vs",group2))%>%filter(group2%in%c("WT_18WK","HFD_18WK","WT_26WK","HFD_26WK"))%>%left_join(polarfp,by=c("met"="met"))
untrp<-untarwp%>%mutate(comparision=paste0(group1,"vs",group2))%>%filter(group2%in%c("WT_18WK","HFD_18WK","WT_26WK","HFD_26WK"))%>%left_join(untfp,by=c("met"="met"))
#####
write.csv(scnprp,file="SCN_wilox_all_pathway.csv")
write.csv(acyrp,file="Acyl_wilox_pathway.csv")
write.csv(polarrp,file="Untargeted_wilox_all_pathway.csv")
write.csv(untrp,file="Polar_wilox_pathway.csv")
#########
library(pheatmap)
pheatmap(t(ddpd$SCN[,unique(subset(scnwp,p<0.05)$met)]),scale="row",annotation_col = metas[rownames(ddpd$SCN),"Group",drop=F],cluster_cols = F,show_colnames = F,
         color = colorRampPalette(c("darkblue","white","red"))(128),border_color = "white",cellwidth = 6,cellheight = 7,fontsize_row = 6)
dev.print(pdf,file="SCN_wilcox_p_all_pathway.pdf")
pheatmap(t(ddpd$Acyl[,unique(subset(acywp,p<0.05)$met),drop=F]),scale="row",annotation_col = metaa[rownames(ddpd$Acyl),"Group",drop=F],cluster_cols = F,cluster_rows = F,show_colnames = F,
         color = colorRampPalette(c("darkblue","white","red"))(128),border_color = "white",cellwidth = 8,cellheight = 10)
dev.print(pdf,file="Acyl_wilcox_p.pdf")
###
pheatmap(t(ddpd$Polar[,unique(subset(polarwp,p<0.05)$met)]),scale="row",annotation_col = metap[rownames(ddpd$Polar),"Group",drop=F],cluster_cols = F,show_colnames = F,
         color = colorRampPalette(c("darkblue","white","red"))(128),border_color = "white",cellwidth = 8,cellheight = 10)
dev.print(pdf,file="Polar_wilcox_p_pathway.pdf")
###
pheatmap(t(ddpd$Untargeted[,unique(subset(untarwp,p<0.05)$met)]),scale="row",annotation_col = metau[rownames(ddpd$Untargeted),"Group",drop=F],cluster_cols = F,show_colnames = F,
         color = colorRampPalette(c("darkblue","white","red"))(128),border_color = "white",cellwidth = 6,cellheight = 7,fontsize_row = 8)
dev.print(pdf,file="Untargeted_wilcox_p_pathway.pdf")
###########

scnrs<-split(scnprp,scnprp$comparision)
acyrs<-split(acyrp,acyrp$comparision)
polarrs<-split(polarrp,polarrp$comparision)
untrs<-split(untrp,untrp$comparision)
#####
acyrd<-lapply(acyrs, function(x)as.data.frame(x))
polarrd<-lapply(polarrs, function(x)as.data.frame(x))
untrd<-lapply(untrs, function(x)as.data.frame(x))
###
scnrg<-lapply(scnrs, function(x)subset(x,p<0.05)$met)
acyrg<-lapply(acyrs, function(x)subset(x,p<0.05)$met)
polarrg<-lapply(polarrs, function(x)subset(x,p<0.05)$met)
untrg<-lapply(untrs, function(x)subset(x,p<0.05)$met)
#####
library(VennDetail)
cont<-c("HFD_18WKvsWT_18WK","HFD_26WKvsWT_26WK","DR_26wkvsHFD_26WK","EX_26WKvsHFD_26WK","DREX_26WKvsHFD_26WK")
plot(venndetail(untrg[cont]))
write.csv(getFeature(venndetail(untrg[cont]),rlist = untrd[cont],userowname = F,gind=rep(1,5)),file="untargeted_venndetail.csv")
dev.print(pdf,file="Untargeted_venn.pdf")
plot(venndetail(scnrg[cont]))
dev.print(pdf,file="SCN_venn.pdf")
###
plot(venndetail(acyrg[cont]))
dev.print(pdf,file="Acyr_venn.pdf")
write.csv(getFeature(venndetail(acyrg[cont]),rlist = acyrd[cont],userowname = F,gind=rep(1,5)),file="Acyr_venndetail.csv")

####No
source("../../omics.r")
plun<-make.PLSDA.model(metau$Group,ddpd$Untargeted)
plotIndiv(plun$plsda,ellipse = F, style = "ggplot2",title = "Untargted",group=metau$Group,legend=T,ind.names = T)
dev.print(pdf,file="untarget_plsda_name.pdf")
plotIndiv(plun$plsda,ellipse = F, style = "ggplot2",title = "Untargted",group=metau$Group,legend=T,ind.names = F)
dev.print(pdf,file="untarget_plsda_wo_name.pdf")
plotIndiv(plun$plsda,ellipse = T, style = "ggplot2",title = "Untargted",group=metau$Group,legend=T,ind.names = F)
dev.print(pdf,file="untarget_plsda.pdf")
###
vipun<-plun$vip
plot.vip2(vipun)
dev.print(pdf,file="plsda_vip_untargeted.pdf")
#### only two pathways
plac<-make.PLSDA.model(metaa$Group,ddpd$Acyl)
plotIndiv(plac$plsda,ellipse = F, style = "ggplot2",title = "Acyl",group=metaa$Group,legend=T,ind.names = T)
dev.print(pdf,file="Acyl_plsda_name.pdf")
plotIndiv(plac$plsda,ellipse = F, style = "ggplot2",title = "Acyl",group=metaa$Group,legend=T,ind.names = F)
dev.print(pdf,file="Acyl_plsda_wo_name.pdf")
plotIndiv(plac$plsda,ellipse = T, style = "ggplot2",title = "Acyl",group=metaa$Group,legend=T,ind.names = F)
dev.print(pdf,file="Acyl_plsda.pdf")
###
vipac<-plac$vip
plot.vip2(vipac)
dev.print(pdf,file="plsda_vip_Acyl.pdf")
#####
plsc<-make.PLSDA.model(metas$Group,ddpd$SCN)
plotIndiv(plsc$plsda,ellipse = F, style = "ggplot2",title = "SCN",group=metas$Group,legend=T,ind.names = T)
dev.print(pdf,file="SCN_plsda_name.pdf")
plotIndiv(plsc$plsda,ellipse = F, style = "ggplot2",title = "SCN",group=metas$Group,legend=T,ind.names = F)
dev.print(pdf,file="SCN_plsda_wo_name.pdf")
plotIndiv(plsc$plsda,ellipse = T, style = "ggplot2",title = "SCN",group=metas$Group,legend=T,ind.names = F)
dev.print(pdf,file="SCN_plsda.pdf")
###
vipsc<-plsc$vip
plot.vip2(vipsc)
dev.print(pdf,file="plsda_vip_SCN.pdf")
#####
plsp<-make.PLSDA.model(metap$Group,ddpd$Polar,cv = F)
plotIndiv(plsp$plsda,ellipse = F, style = "ggplot2",title = "Polar",group=metap$Group,legend=T,ind.names = T)
dev.print(pdf,file="Polar_plsda_name.pdf")
plotIndiv(plsp$plsda,ellipse = F, style = "ggplot2",title = "Polar",group=metap$Group,legend=T,ind.names = F)
dev.print(pdf,file="Polar_plsda_wo_name.pdf")
plotIndiv(plsp$plsda,ellipse = T, style = "ggplot2",title = "Polar",group=metap$Group,legend=T,ind.names = F)
dev.print(pdf,file="Polar_plsda.pdf")
###
vipp<-plsp$vip
plot.vip2(vipp)
dev.print(pdf,file="plsda_vip_Polar.pdf")
#####
plsdacom<-function(meta,group,data,cv = F,ncomp = 4,log = T){
  mm<-subset(meta,Group%in%group)
  id<-rownames(mm)
  dd<-as.matrix(data[id,])
  plsp<-make.PLSDA.model(mm$Group,dd,cv=cv,ncomp=ncomp,log=log)
  return(list(plsp=plsp,group=mm$Group))
}
####
compp<-comp[c(1,12,16:18),]
###
scnv<-list()
for(i in 1:nrow(compp)){
  group<-compp[i,]
  tmp<-plsdacom(metas,group=group,ddpd$SCN,cv=F,ncomp=4,log=F)
  nam<-paste0("SCN_",group[1],"vs",group[2])
  scnv[[nam]]<-tmp
}
#####
acyv<-list()
for(i in 1:nrow(compp)){
  group<-compp[i,]
  tmp<-plsdacom(metaa,group=group,ddpd$Acyl,cv=F,ncomp=4,log=F)
  nam<-paste0("Acyl_",group[1],"vs",group[2])
  acyv[[nam]]<-tmp
}
#
ploarv<-list()
for(i in 1:nrow(compp)){
  group<-compp[i,]
  tmp<-plsdacom(metap,group=group,ddpd$Polar,cv=F,ncomp=4,log=F)
  nam<-paste0("Polar_",group[1],"vs",group[2])
  ploarv[[nam]]<-tmp
}
#untarget log =T
untv<-list()
for(i in 1:nrow(compp)){
  group<-compp[i,]
  tmp<-plsdacom(metau,group=group,ddpd$Untargeted,cv=F,ncomp=4,log=F)
  nam<-paste0("Untarget_",group[1],"vs",group[2])
  untv[[nam]]<-tmp
}
####
######
write.vip2<-function(data,name,cutoff=1){
  data<-as.data.frame(data)
  data$species=rownames(data)
  dd<-data%>%filter(comp1>=cutoff)%>%arrange(desc(comp1))
  write.csv(dd,file=paste(name,"_VIP_cut_pathway.csv",sep=""),row.names=F)
}
for(i in names(scnv)[c(1,4)]){
  tmp<-scnv[[i]]
  plotIndiv(tmp$plsp$plsda,ellipse = F, style = "ggplot2",title = "SCN",group=tmp$group,legend=T,ind.names = T)
  ggsave(file=paste0(i,"_plsda_name.pdf"),width=5.5,height=4.2)
  plotIndiv(tmp$plsp$plsda,ellipse = T, style = "ggplot2",title = "SCN",group=tmp$group,legend=T,ind.names = F)
  ggsave(file=paste0(i,"_plsda_wo_name_pathway.pdf"),width=5.5,height=4.2)
  vipt<-tmp$plsp$vip
  plot.vip2(vipt)
  ggsave(file=paste0(i,"_plsda_vip_fc_pathway.pdf"),width=12,height=5)
  write.vip2(vipt,name=i,cutoff=1)
}
###
 ###
for(i in names(ploarv)){
  tmp<-ploarv[[i]]
  plotIndiv(tmp$plsp$plsda,ellipse = F, style = "ggplot2",title = "Polar",group=tmp$group,legend=T,ind.names = T)
  ggsave(file=paste0(i,"_plsda_name_pathway.pdf"),width=5.5,height=4.2)
  plotIndiv(tmp$plsp$plsda,ellipse = T, style = "ggplot2",title = "Polar",group=tmp$group,legend=T,ind.names = F)
  ggsave(file=paste0(i,"_plsda_wo_name_pathway.pdf"),width=5.5,height=4.2)
  vipt<-tmp$plsp$vip
  plot.vip2(vipt,fc = tmp$plsp$fc)
  ggsave(file=paste0(i,"_plsda_vip_fc_pathway.pdf"),width=8,height=6)
  write.vip2(vipt,name=i,cutoff=1)
}
#######
for(i in names(untv)){
  tmp<-untv[[i]]
  plotIndiv(tmp$plsp$plsda,ellipse = F, style = "ggplot2",title = "Untargeted",group=tmp$group,legend=T,ind.names = T)
  ggsave(file=paste0(i,"_plsda_name_pathway.pdf"),width=5.5,height=4.2)
  plotIndiv(tmp$plsp$plsda,ellipse = T, style = "ggplot2",title = "Untargeted",group=tmp$group,legend=T,ind.names = F)
  ggsave(file=paste0(i,"_plsda_wo_name_pathway.pdf"),width=5.5,height=4.2)
  vipt<-tmp$plsp$vip
  plot.vip2(vipt,fc = tmp$plsp$fc)
  ggsave(file=paste0(i,"_plsda_vip_fc_pathway.pdf"),width=8,height=6)
  write.vip2(vipt,name=i,cutoff=1)
}



