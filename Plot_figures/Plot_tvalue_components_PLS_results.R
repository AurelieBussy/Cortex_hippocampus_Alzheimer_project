#### library modules ####
library("RMINC")
library("lmerTest")
library("corrplot")
library("rlang")
library("ggplot2")
library("Hmisk") 
library("reshape2")
library("plyr")
library("gplots")
library("Matrix")
library("table1")
library("magrittr")
library("ggpubr)"
library("viridis")
library("effects")
library("splines")
library("lme4")
library("radarchart")
library("fmsb")
library("ggradar")
library("devtools")
library("ComplexHeatmap")
library("grid")
library("gridExtra")
library("ComplexHeatmap")
library("grid")
library("ggpubr")
library("circlize")
library("matrixStats")
library("lsa")
library("ggridges")
library("hrbrthemes")
library("forcats")
library("scatterplot3d")
library("dplyr")
library("survival")
library("survminer")
library("psych")   
library("ggcorrplot")  
library("vcd")
library("ggthemes")
library("tidyr")
library("scales")     
library("missForest")
library("SpatioTemporal")
library("extrafont")
library("BSDA")
library("viridisLite")
library("tidyverse")
library("rstatix")
library("ggpubr")
library("R.matlab")
library("RColorBrewer")


        # From NMF cortex
demographics_NMF_CIVET<-read.csv('NMF_outputs/Cortex/k10/demographics_and_nmfweights_k10.csv')
demographics_NMF_CIVET$Group<-factor(demographics_NMF_CIVET$Group, levels=c("HC","FAMHX","MCI","AD"))

        # From NMF hippocampus
demographics_NMF_hipp<-read.csv('NMF_outputs/Hippocampus/demographics_and_nmfweights_k4.csv')
demographics_NMF_hipp$Group<-factor(demographics_NMF_hipp$Group, levels=c("HC","FAMHX","MCI","AD"))

        #define number of each component/metric
metrics_CIVET=c("CT","SA","T1","T2")
nmetric_CIVET=4
ncomp_CIVET=10
ncol_CIVET=ncol(demographics_NMF_CIVET)
metrics_hipp=c("DBM","T1","T2")
nmetric_hipp=3
ncomp_hipp=4
ncol_hipp=ncol(demographics_NMF_hipp)

        # select only NMF weights hipp
components_hipp<-demographics_NMF_hipp[,(ncol_hipp-ncomp_hipp*nmetric_hipp+1):ncol_hipp]

        # make sure we are using the same participants from hipp et cortex
demographics_NMF_hipp_ID<-as.data.frame(demographics_NMF_hipp$ID)
colnames(demographics_NMF_hipp_ID)<-"ID"
demographics_CIVET_bis<-merge(demographics_NMF_hipp_ID,demographics_NMF_CIVET,by="ID")

        # select only NMF weights cortex
components_CIVET<-demographics_CIVET_bis[,(ncol_CIVET-ncomp_CIVET*nmetric_CIVET+1):ncol_CIVET]
demographics_all<-cbind(demographics_NMF_hipp[,1:5],components_CIVET,components_hipp)

        # prepare proper rowname order for plots      
rowname_CIVET<-paste(rep(metrics_CIVET,ncomp_CIVET),paste(rep(1:ncomp_CIVET, each=nmetric_CIVET)),sep="_CIVET_")
ideal_rowname_CIVET<-paste(rep(metrics_CIVET, each=ncomp_CIVET),paste(rep(1:ncomp_CIVET, nmetric_CIVET)),sep="_CIVET_")
rowname_hipp<-paste(rep(paste(rep("HIPP_C", ncomp_hipp),c(seq( from = 1, to = ncomp_hipp )),sep=""), each=nmetric_hipp),rep(metrics_hipp,ncomp_hipp),sep="_")
ideal_rowname_hipp<-paste(rep(rep(paste(rep("HIPP_C", ncomp_hipp),c(seq( from = 1, to = ncomp_hipp )),sep="")),nmetric_hipp),rep(metrics_hipp,each=ncomp_hipp),sep="_")
metrics_all_CIVET<-rep(metrics_CIVET,each=ncomp_CIVET)
metrics_all_hipp<-rep(metrics_hipp,each=ncomp_hipp)
ideal_rowname<-c(ideal_rowname_CIVET,ideal_rowname_hipp)
rowname<-c(rowname_hipp,rowname_CIVET)
metrics_all<-c(metrics_all_CIVET,metrics_all_hipp)

        # load PLS outputs csv values
dec = "." 
tvalue<-read.csv("PLS_outputs/tvalues_PLS_no_info.csv",header=F, sep=" ")

ideal_rowname<-c(ideal_rowname_CIVET,ideal_rowname_hipp)
rowname<-c(rowname_CIVET,rowname_hipp)
metrics_all<-c(metrics_all_CIVET,metrics_all_hipp)
metrics_all_diff<-c(paste(metrics_all_CIVET, "CTX",sep="_"),paste(metrics_all_hipp, "hipp",sep="_"))

        # add 4th column with component's names
tvalue[,4]<-as.data.frame(rowname)
colnames(tvalue)<-c("tvalue_LV1","tvalue_LV2","tvalue_LV3","rowname")

        # reorder to match proper order
tvalue_sorted<-tvalue[match(ideal_rowname, tvalue$rowname),]
tvalue_sorted$rowname <- factor(tvalue_sorted$rowname, levels = tvalue_sorted$rowname)

        # binarize tvalue 0=not significant 1= significant
for (i in 1:nrow(tvalue_sorted)){
  if (abs(tvalue_sorted[i,1])<1.96){tvalue_sorted[i,5]<-0}else{tvalue_sorted[i,5]<-1}
  if (abs(tvalue_sorted[i,2])<1.96){tvalue_sorted[i,6]<-0}else{tvalue_sorted[i,6]<-1}
  if (abs(tvalue_sorted[i,3])<1.96){tvalue_sorted[i,7]<-0}else{tvalue_sorted[i,7]<-1}
}

tvalue_sorted$V5<-factor(tvalue_sorted$V5)
tvalue_sorted$V6<-factor(tvalue_sorted$V6)
tvalue_sorted$V7<-factor(tvalue_sorted$V7)
tvalue_sorted$V8<-metrics_all
tvalue_sorted$V12<-metrics_all_diff

        # create new column LV1=#9, LV2=#10 and LV3=#11 where if signficant = metric name, not significant=0 (for plotting later)
for (i in 1:nrow(tvalue_sorted)){
  if(tvalue_sorted[i,5]==0) {tvalue_sorted$V9[i]<-0} else {tvalue_sorted$V9[i]<-paste(metrics_all[i],tvalue_sorted$V5[i],sep="")}
  if(tvalue_sorted[i,6]==0) {tvalue_sorted$V10[i]<-0} else {tvalue_sorted$V10[i]<-paste(metrics_all[i],tvalue_sorted$V6[i],sep="")}
  if(tvalue_sorted[i,7]==0) {tvalue_sorted$V11[i]<-0} else {tvalue_sorted$V11[i]<-paste(metrics_all[i],tvalue_sorted$V7[i],sep="")}
}

        # Plot tvalues for LV1 in the right order and with each bar colored depending of metric and significance
tvalue_plot_LV1<- ggplot(tvalue_sorted,aes(x=as.numeric(as.character(-tvalue_sorted$tvalue_LV1)),y=fct_rev(rowname),fill=V9))+scale_y_discrete()+
  annotate("rect",xmin = -Inf, xmax = Inf, ymin = 0.5, ymax = 4.5, fill = '#2F435A', alpha = 0.2) +
  annotate("rect",xmin = -Inf, xmax = Inf, ymin = 4.5, ymax = 8.5, fill = '#39918C', alpha = 0.2) +
  annotate("rect",xmin = -Inf, xmax = Inf, ymin = 8.5, ymax = 12.5, fill = '#E0CFA7', alpha = 0.2) +
  annotate("rect",xmin = -Inf, xmax = Inf, ymin = 12.5, ymax = 22.5, fill = '#2F435A', alpha = 0.2) +
  annotate("rect",xmin = -Inf, xmax = Inf, ymin = 22.5, ymax = 32.5, fill = '#39918C', alpha = 0.2) +
  annotate("rect",xmin = -Inf, xmax = Inf, ymin = 32.5, ymax = 42.5, fill = '#D0B49F', alpha = 0.2) +
  annotate("rect",xmin = -Inf, xmax = Inf, ymin = 42.5, ymax = 52.5, fill = '#AB6B51', alpha = 0.2) +
  annotate("rect",ymin = c(0.5,4.5,8.5,12.5,22.5,32.5,42.5), ymax = c(4.5,8.5,12.5,22.5,32.5,42.5,52.5), xmin = -Inf, xmax = -8, fill = c('#2F435A','#39918C','#E0CFA7','#2F435A','#39918C','#D0B49F','#AB6B51'), alpha = 1) +
  geom_bar(stat='identity',colour="black",size=0.1)+
  theme_classic()+theme(legend.position="none")+xlim(-8, 11) + 
  geom_vline(xintercept = c(-1.96,1.96),color="black")+
  geom_vline(xintercept = c(-2.58,2.58),color="grey50")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "bold"))+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "bold"))+
  scale_fill_manual(values=c("white","#AB6B51","#E0CFA7","#D0B49F","#39918C" ,"#2F435A"))

ggsave(tvalue_plot_LV1, file="Figures_final/tvalue_plot_HIPP_CIVET_PLS_LV1.png", width=5, height=9)

        # Plot tvalues for LV2 in the right order and with each bar colored depending of metric and significance

        tvalue_plot_LV2<-ggplot(tvalue_sorted,aes(x=tvalue_LV2,y=fct_rev(rowname),fill=V10))+
  scale_y_discrete()+
  annotate("rect",xmin = -Inf, xmax = Inf, ymin = 0.5, ymax = 4.5, fill = '#2F435A', alpha = 0.2) +
  annotate("rect",xmin = -Inf, xmax = Inf, ymin = 4.5, ymax = 8.5, fill = '#39918C', alpha = 0.2) +
  annotate("rect",xmin = -Inf, xmax = Inf, ymin = 8.5, ymax = 12.5, fill = '#E0CFA7', alpha = 0.2) +
  annotate("rect",xmin = -Inf, xmax = Inf, ymin = 12.5, ymax = 22.5, fill = '#2F435A', alpha = 0.2) +
  annotate("rect",xmin = -Inf, xmax = Inf, ymin = 22.5, ymax = 32.5, fill = '#39918C', alpha = 0.2) +
  annotate("rect",xmin = -Inf, xmax = Inf, ymin = 32.5, ymax = 42.5, fill = '#D0B49F', alpha = 0.2) +
  annotate("rect",xmin = -Inf, xmax = Inf, ymin = 42.5, ymax = 52.5, fill = '#AB6B51', alpha = 0.2) +
  annotate("rect",ymin = c(0.5,4.5,8.5,12.5,22.5,32.5,42.5), ymax = c(4.5,8.5,12.5,22.5,32.5,42.5,52.5), xmin = -Inf, xmax = -8, fill = c('#2F435A','#39918C','#E0CFA7','#2F435A','#39918C','#D0B49F','#AB6B51'), alpha = 1) +
  geom_bar(stat='identity', color="black", size=0.1)+
  theme_classic()+theme(legend.position="none")+  xlim(-8, 11) + 
  geom_vline(xintercept = c(-1.96,1.96),color="black")+
  geom_vline(xintercept = c(-2.58,2.58),color="grey50")+
  theme(axis.title.x = element_blank(), axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "bold"))+
  theme(axis.title.y = element_blank(), axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "bold"))+
  scale_fill_manual(values=c("white","#AB6B51","#E0CFA7","#D0B49F","#39918C" ,"#2F435A"))

ggsave(tvalue_plot_LV2, file="Figures_final/tvalue_plot_HIPP_CIVET_PLS_LV2.png", width=5, height=9)
        
        # Plot tvalues for LV3 in the right order and with each bar colored depending of metric and significance

tvalue_plot_LV3<-ggplot(tvalue_sorted,aes(x=tvalue_LV3,y=fct_rev(rowname),fill=V11))+
  scale_y_discrete()+
  annotate("rect",xmin = -Inf, xmax = Inf, ymin = 0.5, ymax = 4.5, fill = '#2F435A', alpha = 0.2) +
  annotate("rect",xmin = -Inf, xmax = Inf, ymin = 4.5, ymax = 8.5, fill = '#39918C', alpha = 0.2) +
  annotate("rect",xmin = -Inf, xmax = Inf, ymin = 8.5, ymax = 12.5, fill = '#E0CFA7', alpha = 0.2) +
  annotate("rect",xmin = -Inf, xmax = Inf, ymin = 12.5, ymax = 22.5, fill = '#2F435A', alpha = 0.2) +
  annotate("rect",xmin = -Inf, xmax = Inf, ymin = 22.5, ymax = 32.5, fill = '#39918C', alpha = 0.2) +
  annotate("rect",xmin = -Inf, xmax = Inf, ymin = 32.5, ymax = 42.5, fill = '#D0B49F', alpha = 0.2) +
  annotate("rect",xmin = -Inf, xmax = Inf, ymin = 42.5, ymax = 52.5, fill = '#AB6B51', alpha = 0.2) +
  annotate("rect",ymin = c(0.5,4.5,8.5,12.5,22.5,32.5,42.5), ymax = c(4.5,8.5,12.5,22.5,32.5,42.5,52.5), xmin = -Inf, xmax = -8, fill = c('#2F435A','#39918C','#E0CFA7','#2F435A','#39918C','#D0B49F','#AB6B51'), alpha = 1) +
  geom_bar(stat='identity', color="black", size=0.1)+
  theme_classic()+theme(legend.position="none")+
  xlim(-8, 11) + 
  geom_vline(xintercept = c(-1.96,1.96),color="black")+
  geom_vline(xintercept = c(-2.58,2.58),color="grey50")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "bold"))+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "bold"))+
  scale_fill_manual(values=c("white","#E0CFA7","#D0B49F","#39918C" ,"#2F435A"))

ggsave(tvalue_plot_LV3, file="Figures_final/tvalue_plot_HIPP_CIVET_PLS_LV3.png", width=5, height=9)
