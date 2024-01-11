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

        # load PLS correlation oututs
demog<-read.csv("PLS_outputs/PLS_behaviour_LV123.csv",header=T,sep=",")
demog$Demographics<-factor(demog$Demographics, levels = demog$Demographics)

        # reorganise demographics to have order you want
Demographics<-(demog$Demographics)
Demographics_order<-as.data.frame(c( "alcohol" ,"current_smoking","past_smoking", "drugs","high_BP","high_chol" , "diabetes",
                                     "MOCA","RBANS_total","imm_mem" , "del_mem","attention" , "language","visuospatial" ,
                                     "depression" ,"anxiety",
                                     "hearing_problem","concussion","TIA","brain_injury","headaches_migraines","seizures",
                                     "heart_disease" ,"liver_disease" , "kidney_disease", "thyroid_disease","cancer",
                                     "arthritis","neck_back_prob","allergies"))
colnames(Demographics_order)<-"Demographics"

demog_sorted<-merge(Demographics_order,demog, by="Demographics", sort = F)
demog_sorted$Demographics <- factor(demog_sorted$Demographics, levels = demog_sorted$Demographics)

        # plot correlation (for LV1 I added a '-' to corr_LV1 and llcorr_LV1; same as without any "-"; it's just a preference of plot orientation)
demog_plot_LV1<-ggplot(demog_sorted,aes(x=-corr_LV1,y=fct_rev(Demographics),fill=Demographics))+
  xlim(-0.7, 0.7) +  scale_y_discrete()+
  annotate("rect",xmin = -Inf, xmax = Inf, ymin = 0.5, ymax = 14.5, fill = '#533440', alpha = 0.2) +
  annotate("rect",xmin = -Inf, xmax = Inf, ymin = 14.5, ymax = 16.5, fill = '#A47786', alpha = 0.2) +
  annotate("rect",xmin = -Inf, xmax = Inf, ymin = 16.5, ymax = 23.5, fill = '#E4D4C8', alpha = 0.2) +
  annotate("rect",xmin = -Inf, xmax = Inf, ymin = 23.5, ymax = 30.5, fill = '#9DB6CC', alpha = 0.2) +
  annotate("rect",ymin = c(0.5,14.5,16.5,23.5), ymax = c(14.5,16.5,23.5,30.5), xmin = -Inf, xmax = -0.7, fill = c('#533440','#A47786','#E4D4C8','#9DB6CC'), alpha = 1) +
  geom_bar(stat='identity',colour="black",size=0.1)+
  geom_errorbar(aes(xmin=-llcorr_LV1, xmax=-ulcorr_LV1), width=.2,
                position=position_dodge(.9)) +
  theme_classic()+theme(legend.position="none")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "bold"))+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "bold"))+
  scale_fill_manual(values=c("white","white","#9DB6CC","white","#9DB6CC","#9DB6CC","white",
                             "#E4D4C8","#E4D4C8","#E4D4C8","#E4D4C8","#E4D4C8","#E4D4C8","#E4D4C8",
                             "white","#A47786",
                             "white","white","white","white","white","white",
                             "white","white","white","white","white","white","white","white","white"))

ggsave(demog_plot_LV1, file="demog_PLS_LV1.png", width=6, height=9)

demog_plot_LV2<-ggplot(demog_sorted,aes(x=corr_LV2,y=fct_rev(Demographics),fill=Demographics))+
  scale_y_discrete()+
  annotate("rect",xmin = -Inf, xmax = Inf, ymin = 0.5, ymax = 14.5, fill = '#533440', alpha = 0.2) +
  annotate("rect",xmin = -Inf, xmax = Inf, ymin = 14.5, ymax = 16.5, fill = '#A47786', alpha = 0.2) +
  annotate("rect",xmin = -Inf, xmax = Inf, ymin = 16.5, ymax = 23.5, fill = '#E4D4C8', alpha = 0.2) +
  annotate("rect",xmin = -Inf, xmax = Inf, ymin = 23.5, ymax = 30.5, fill = '#9DB6CC', alpha = 0.2) +
  annotate("rect",ymin = c(0.5,14.5,16.5,23.5), ymax = c(14.5,16.5,23.5,30.5), xmin = -Inf, xmax = -0.7, fill = c('#533440','#A47786','#E4D4C8','#9DB6CC'), alpha = 1) +
  geom_bar(stat='identity',colour="black",size=0.1)+
  geom_errorbar(aes(xmin=llcorr_LV2, xmax=ulcorr_LV2), width=.2,
                position=position_dodge(.9)) +xlim(-0.7, 0.7) + 
  theme_classic()+theme(legend.position="none")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "bold"))+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "bold"))+
  scale_fill_manual(values=c("white","white","white","white","white","white","#9DB6CC","white","white",
                             "white","white","white","white","white","white","white",
                             "white","white",
                             "white","#533440","white","white","white","white",
                             "#533440","white","white","#533440","#533440","white"))

ggsave(demog_plot_LV2, file="demog_PLS_LV2.png", width=6, height=9)

demog_plot_LV3<-ggplot(demog_sorted,aes(x=corr_LV3,y=fct_rev(Demographics),fill=Demographics))+
  scale_y_discrete()+
  annotate("rect",xmin = -Inf, xmax = Inf, ymin = 0.5, ymax = 14.5, fill = '#533440', alpha = 0.2) +
  annotate("rect",xmin = -Inf, xmax = Inf, ymin = 14.5, ymax = 16.5, fill = '#A47786', alpha = 0.2) +
  annotate("rect",xmin = -Inf, xmax = Inf, ymin = 16.5, ymax = 23.5, fill = '#E4D4C8', alpha = 0.2) +
  annotate("rect",xmin = -Inf, xmax = Inf, ymin = 23.5, ymax = 30.5, fill = '#9DB6CC', alpha = 0.2) +
  annotate("rect",ymin = c(0.5,14.5,16.5,23.5), ymax = c(14.5,16.5,23.5,30.5), xmin = -Inf, xmax = -0.7, fill = c('#533440','#A47786','#E4D4C8','#9DB6CC'), alpha = 1) +
  geom_bar(stat='identity',colour="black",size=0.1)+
  geom_errorbar(aes(xmin=llcorr_LV3, xmax=ulcorr_LV3), width=.2,
                position=position_dodge(.9)) +
  xlim(-0.7, 0.7) + 
  theme_classic()+theme(legend.position="none")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "bold"))+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "bold"))+
  scale_fill_manual(values=c("white","white","white","white","white","white","white","white","white",
                             "white","white","white","white","white",
                             "#A47786","white",                            
                             "white","white", "white","white","white","white","#533440","#533440",
                             "white","white","white","white","white","white"))

ggsave(demog_plot_LV3, file="demog_PLS_LV3.png", width=6, height=9)
