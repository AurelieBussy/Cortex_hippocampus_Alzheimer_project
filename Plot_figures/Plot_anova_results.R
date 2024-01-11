#####
### run the other R scripts from repository to have necessary data loaded (Plot_correlation_PLS_results.R and Plot_BSR_components_PLS_results.R)

        # create empty matrices of size (52=number of components; 6=number of variables of interest from model)
matrix_padj<-data.frame(matrix(ncol=52,nrow=6))
matrix_diff<-data.frame(matrix(ncol=52,nrow=6))
matrix_lwr<-data.frame(matrix(ncol=52,nrow=6))
matrix_upr<-data.frame(matrix(ncol=52,nrow=6))

        # save Group pvalue, lower, upper, and diff values from ancova)
for (i in 1:52){
  model=lm(all[,i]~Group+Age,demographics_NMF)
  matrix_padj[,i]<-TukeyHSD(x=aov(model), 'Group', conf.level=0.95)$Group[,4]
  matrix_diff[,i]<-TukeyHSD(x=aov(model), 'Group', conf.level=0.95)$Group[,1]
  matrix_lwr[,i]<-TukeyHSD(x=aov(model), 'Group', conf.level=0.95)$Group[,2]
  matrix_upr[,i]<-TukeyHSD(x=aov(model), 'Group', conf.level=0.95)$Group[,3]
}

        # name groups of all matrics (check order; depends on which group is the reference in your model)
Group<-c("HC-FAMHX","MCI-FAMHX","AD-FAMHX","MCI-HC","AD-HC","AD-MCI")
matrix_padj[,53]<-Group
matrix_diff[,53]<-Group
matrix_lwr[,53]<-Group
matrix_upr[,53]<-Group

rownames(matrix_padj)<-c("HC-FAMHX","MCI-FAMHX","AD-FAMHX","MCI-HC","AD-HC","AD-MCI")
colnames(matrix_padj)<-c(rowname_hipp,rowname_CIVET,"Group")
rownames(matrix_diff)<-c("HC-FAMHX","MCI-FAMHX","AD-FAMHX","MCI-HC","AD-HC","AD-MCI")
colnames(matrix_diff)<-c(rowname_hipp,rowname_CIVET,"Group")
rownames(matrix_lwr)<-c("HC-FAMHX","MCI-FAMHX","AD-FAMHX","MCI-HC","AD-HC","AD-MCI")
colnames(matrix_lwr)<-c(rowname_hipp,rowname_CIVET,"Group")
rownames(matrix_upr)<-c("HC-FAMHX","MCI-FAMHX","AD-FAMHX","MCI-HC","AD-HC","AD-MCI")
colnames(matrix_upr)<-c(rowname_hipp,rowname_CIVET,"Group")

        # stack all pvalues to perform FDR correction on the entirety of pvalues
newmydata<-stack(matrix_padj[,1:52]) # stack pvalues
datafdr<-newmydata[1]
datafdrr<-as.numeric(unlist(datafdr))
datacorrected<-p.adjust(datafdrr,method = "fdr") # perform FDR correction
datacorrected<-data.frame(datacorrected)
Unstack<-cbind(datacorrected,newmydata[,2])
pvaluecorrected<-unstack(Unstack) # unstack to original format
pvalueTRUE_FALSE<-pvaluecorrected<=0.05 # matrix FALSE=not-sign, TRUE=sign
rownames(pvaluecorrected)<-rownames(matrix_padj) # renames matrices
rownames(pvalueTRUE_FALSE)<-rownames(matrix_padj)

matrix_diff$Group<-as.factor(matrix_diff$Group)
matrix_upr$Group<-as.factor(matrix_upr$Group)
matrix_lwr$Group<-as.factor(matrix_lwr$Group)
matrix_padj$Group<-as.factor(matrix_padj$Group)

        # sort matrices in order you want in fig (bigger effect 1st, smaller effects last)
Group_ideal<-c("AD-FAMHX","AD-HC","AD-MCI","MCI-FAMHX","MCI-HC","HC-FAMHX")
matrix_diff_sorted<-matrix_diff[c(3 ,5 ,6 ,2, 4 ,1),]
matrix_upr_sorted<-matrix_upr[c(3 ,5 ,6 ,2, 4 ,1),]
matrix_lwr_sorted<-matrix_lwr[c(3 ,5 ,6 ,2, 4 ,1),]
matrix_padj_sorted<-matrix_padj[c(3 ,5 ,6 ,2, 4 ,1),]
pvalueTRUE_FALSE_sorted<-pvalueTRUE_FALSE[c(3 ,5 ,6 ,2, 4 ,1),]

        # create matrices for colors
color_values<-data.frame(matrix(ncol=6,nrow=52))
fill_values<-data.frame(matrix(ncol=6,nrow=52))
matrix_diff_sorted$Group<-factor(matrix_diff_sorted$Group, levels = matrix_diff_sorted$Group)

        # loop to assign color depending on significance (grey=not-sign, color= sign)
for (i in 1:52){
  p_truth<-pvalueTRUE_FALSE_sorted[,i]
  if (p_truth[1]=="FALSE"){color_values[i,1]<-"grey50"}else{color_values[i,1]<-"#39918C"}
  if (p_truth[2]=="FALSE"){color_values[i,2]<-"grey50"}else{color_values[i,2]<-"#2F435A"}
  if (p_truth[3]=="FALSE"){color_values[i,3]<-"grey50"}else{color_values[i,3]<-"#9DB6CC"}
  if (p_truth[4]=="FALSE"){color_values[i,4]<-"grey50"}else{color_values[i,4]<-"#E4D4C8"}
  if (p_truth[5]=="FALSE"){color_values[i,5]<-"grey50"}else{color_values[i,5]<-"#A47786"}
  if (p_truth[6]=="FALSE"){color_values[i,6]<-"grey50"}else{color_values[i,6]<-"#533440"}
  if (p_truth[1]=="FALSE"){fill_values[i,1]<-"grey50"}else{fill_values[i,1]<-"#39918C"}
  if (p_truth[2]=="FALSE"){fill_values[i,2]<-"grey50"}else{fill_values[i,2]<-"#2F435A"}
  if (p_truth[3]=="FALSE"){fill_values[i,3]<-"grey50"}else{fill_values[i,3]<-"#9DB6CC"}
  if (p_truth[4]=="FALSE"){fill_values[i,4]<-"grey50"}else{fill_values[i,4]<-"#E4D4C8"}
  if (p_truth[5]=="FALSE"){fill_values[i,5]<-"grey50"}else{fill_values[i,5]<-"#A47786"}
  if (p_truth[6]=="FALSE"){fill_values[i,6]<-"grey50"}else{fill_values[i,6]<-"#533440"}
}
colnames(fill_values)<-c("AD-FAMHX","AD-HC","AD-MCI","MCI-FAMHX","MCI-HC","FAMHX-HC")
rownames(fill_values)<-c(rowname_hipp,rowname_CIVET)
colnames(color_values)<-c("AD-FAMHX","AD-HC","AD-MCI","MCI-FAMHX","MCI-HC","FAMHX-HC")
rownames(color_values)<-c(rowname_hipp,rowname_CIVET)

        # create plots (one plot per metric and per component --> not ideal visualization; see other versions below)
myplots <- vector('list', ncol(matrix_diff_sorted)-1)

i=1
for (i in 1:52){
  myplots[[i]] <- local({
    i<-i
    p1<-ggplot(matrix_diff_sorted,aes(x=as.matrix(matrix_diff_sorted[,i]),y=fct_rev(Group),color=Group, fill=Group))+
      scale_y_discrete()+
      annotate("rect",xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 1.5, fill = '#533440', alpha = 0.3) +
      annotate("rect",xmin = -Inf, xmax = Inf, ymin = 1.5, ymax = 2.5, fill = '#A47786', alpha = 0.3) +
      annotate("rect",xmin = -Inf, xmax = Inf, ymin = 2.5, ymax = 3.5, fill = '#E4D4C8', alpha = 0.3) +
      annotate("rect",xmin = -Inf, xmax = Inf, ymin = 3.5, ymax = 4.5, fill = '#9DB6CC', alpha = 0.3) +
      annotate("rect",xmin = -Inf, xmax = Inf, ymin = 4.5, ymax = 5.5, fill = '#2F435A', alpha = 0.3) +
      annotate("rect",xmin = -Inf, xmax = Inf, ymin = 5.5, ymax = Inf, fill = '#39918C', alpha = 0.3) +
      geom_vline(xintercept = c(0),color="grey30")+
      geom_errorbar(aes(xmin=matrix_lwr_sorted[,i], xmax=matrix_upr_sorted[,i]), width=.2,
                    position=position_dodge(.9)) + #xlim(-35, 55) + 
      geom_point(aes(fill=Group),pch=21, colour="black", size=4)+
      theme_classic()+theme(legend.position="none")+
      theme(axis.title.x = element_blank(),
            axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "bold"))+
      theme(axis.title.y = element_blank(),
            axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "bold"))+
      scale_color_manual(values=as.character(color_values[i,]))+
      scale_fill_manual(values=as.character(fill_values[i,]))
    print(i)
    print(p1)
  })}


      # other configuration of plots (better for plot)

matrix_padj_sorted_t<-as.data.frame(t(matrix_padj_sorted))[1:52,]
matrix_diff_sorted_t<-as.data.frame(t(matrix_diff_sorted))[1:52,]
matrix_upr_sorted_t<-as.data.frame(t(matrix_upr_sorted))[1:52,]
matrix_lwr_sorted_t<-as.data.frame(t(matrix_lwr_sorted))[1:52,]
pvalueTRUE_FALSE_sorted_t<-as.data.frame(t(pvalueTRUE_FALSE_sorted))[1:52,]
        # reorganise matrices 
names_hipp<-paste(paste("HIPP_",paste(rep(1:ncomp_hipp, each=nmetric_hipp)),sep="C"),rep(metrics_hipp,ncomp_hipp),sep="_")
names_CIVET<-paste(rep(metrics_CIVET,ncomp_CIVET),paste(rep(1:ncomp_CIVET, each=nmetric_CIVET)),sep="_CIVET_")
names<-c(names_hipp,names_CIVET)

ideal_rowname_hipp<-paste(paste("HIPP_",paste(rep(1:ncomp_hipp, nmetric_hipp)),sep="C"),rep(metrics_hipp, each=ncomp_hipp),sep="_")
ideal_rowname_CIVET<-paste(rep(metrics_CIVET, each=ncomp_CIVET),paste(rep(1:ncomp_CIVET, nmetric_CIVET)),sep="_CIVET_")
ideal_names<-c(ideal_rowname_CIVET,ideal_rowname_hipp)

ideal_metric<-c(rep(metrics_CIVET, each=ncomp_CIVET),rep(metrics_hipp, each=ncomp_hipp))

matrix_padj_sorted_t_sorted<-matrix_padj_sorted_t[match(ideal_names, names),]
matrix_diff_sorted_t_sorted<-matrix_diff_sorted_t[match(ideal_names, names),]
matrix_upr_sorted_t_sorted<-matrix_upr_sorted_t[match(ideal_names, names),]
matrix_lwr_sorted_t_sorted<-matrix_lwr_sorted_t[match(ideal_names, names),]
pvalueTRUE_FALSE_sorted_t_sorted<-pvalueTRUE_FALSE_sorted_t[match(ideal_names, names),]

matrix_diff_sorted_t_sorted[,7]<-ideal_names
matrix_lwr_sorted_t_sorted[,7]<-ideal_names
matrix_upr_sorted_t_sorted[,7]<-ideal_names
matrix_padj_sorted_t_sorted[,7]<-ideal_names
pvalueTRUE_FALSE_sorted_t_sorted[,7]<-ideal_names
pvalueTRUE_FALSE_sorted_t_sorted[,8]<-ideal_metric

colnames(matrix_diff_sorted_t_sorted)<-c("AD-FAMHX","AD-HC","AD-MCI","MCI-FAMHX","MCI-HC","FAMHX-HC","names")
colnames(matrix_lwr_sorted_t_sorted)<-c("AD-FAMHX","AD-HC","AD-MCI","MCI-FAMHX","MCI-HC","FAMHX-HC","names")
colnames(matrix_upr_sorted_t_sorted)<-c("AD-FAMHX","AD-HC","AD-MCI","MCI-FAMHX","MCI-HC","FAMHX-HC","names")
colnames(matrix_padj_sorted_t_sorted)<-c("AD-FAMHX","AD-HC","AD-MCI","MCI-FAMHX","MCI-HC","FAMHX-HC","names")
colnames(pvalueTRUE_FALSE_sorted_t_sorted)<-c("AD-FAMHX","AD-HC","AD-MCI","MCI-FAMHX","MCI-HC","FAMHX-HC","names","metric")

matrix_diff_sorted_t_sorted$names<-factor(matrix_diff_sorted_t_sorted$names, levels = matrix_diff_sorted_t_sorted$names)
matrix_lwr_sorted_t_sorted$names<-factor(matrix_lwr_sorted_t_sorted$names, levels = matrix_lwr_sorted_t_sorted$names)
matrix_upr_sorted_t_sorted$names<-factor(matrix_upr_sorted_t_sorted$names, levels = matrix_upr_sorted_t_sorted$names)
matrix_padj_sorted_t_sorted$names<-factor(matrix_padj_sorted_t_sorted$names, levels = matrix_padj_sorted_t_sorted$names)
pvalueTRUE_FALSE_sorted_t_sorted$names<-factor(pvalueTRUE_FALSE_sorted_t_sorted$names, levels = pvalueTRUE_FALSE_sorted_t_sorted$names)

    # loop to automatically assign color based on pvalues

color_values<-data.frame(matrix(ncol=6,nrow=52))
fill_values<-data.frame(matrix(ncol=6,nrow=52))

i=1 
j=1
for (i in 1:52){
  for (j in 1:6){
  p_truth<-pvalueTRUE_FALSE_sorted_t_sorted[i,]
  if (p_truth[j]=="FALSE"){color_values[i,j]<-"white"}
  if (p_truth[j]=="TRUE" & p_truth[8]=="CT"){color_values[i,j]<-"#AB6B51"}
  if (p_truth[j]=="TRUE" & p_truth[8]=="SA"){color_values[i,j]<-"#D0B49F"}
  if (p_truth[j]=="TRUE" & p_truth[8]=="T1"){color_values[i,j]<-"#39918C"}
  if (p_truth[j]=="TRUE" & p_truth[8]=="T2"){color_values[i,j]<-"#2F435A"}
  if (p_truth[j]=="TRUE" & p_truth[8]=="DBM"){color_values[i,j]<-"#E0CFA7"}
  
  }
}

      # rename again 
colnames(fill_values)<-c("AD-FAMHX","AD-HC","AD-MCI","MCI-FAMHX","MCI-HC","FAMHX-HC")
rownames(fill_values)<-ideal_names
colnames(color_values)<-c("AD-FAMHX","AD-HC","AD-MCI","MCI-FAMHX","MCI-HC","FAMHX-HC")
rownames(color_values)<-ideal_names
fill_values[,7]<-ideal_names
color_values[,7]<-ideal_names
colnames(fill_values)<-c("AD-FAMHX","AD-HC","AD-MCI","MCI-FAMHX","MCI-HC","FAMHX-HC","names")
colnames(color_values)<-c("AD-FAMHX","AD-HC","AD-MCI","MCI-FAMHX","MCI-HC","FAMHX-HC","names")

fill_values$names<-factor(fill_values$names, levels = fill_values$names)
color_values$names<-factor(color_values$names, levels = color_values$names)

      # loop to create one plot per group comparison --> much better)
i=1
for (i in 1:6){
  myplots[[i]] <- local({
    i<-i
    p1<-ggplot(matrix_diff_sorted_t_sorted,aes(y=as.numeric(as.character(matrix_diff_sorted_t_sorted[,i])),x=factor(names,levels=unique(names))),color=factor(names,levels=unique(names)), fill=factor(names,levels=unique(names)))+
      scale_x_discrete()+
      annotate("rect",ymin = -Inf, ymax = Inf, xmin = -Inf, xmax = 10.5, fill = '#AB6B51', alpha = 0.2) +
      annotate("rect",ymin = -Inf, ymax = Inf, xmin = 10.5, xmax = 20.5, fill = '#D0B49F', alpha = 0.2) +
      annotate("rect",ymin = -Inf, ymax = Inf, xmin = 20.5, xmax = 30.5, fill = '#39918C', alpha = 0.2) +
      annotate("rect",ymin = -Inf, ymax = Inf, xmin = 30.5, xmax = 40.5, fill = '#2F435A', alpha = 0.2) +
      annotate("rect",ymin = -Inf, ymax = Inf, xmin = 40.5, xmax = 44.5, fill = '#E0CFA7', alpha = 0.2) +
      annotate("rect",ymin = -Inf, ymax = Inf, xmin = 44.5, xmax = 48.5, fill = '#39918C', alpha = 0.2) +
      annotate("rect",ymin = -Inf, ymax = Inf, xmin = 48.5, xmax = Inf, fill = '#2F435A', alpha = 0.2) +
      geom_hline(yintercept = c(0),color="grey30")+
      geom_errorbar(aes(ymin=as.numeric(as.character(matrix_lwr_sorted_t_sorted[,i])), ymax=as.numeric(as.character(matrix_upr_sorted_t_sorted[,i]))), width=.3,
                    position=position_dodge(.9)) + ylim(-90, 55) + 
      geom_point(aes(fill=names),pch=21, colour="black", size=6)+
      theme_classic()+theme(legend.position="none")+
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank())+
      theme(axis.title.y = element_blank(),
            axis.text.y = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = .5, face = "bold"))+
      scale_color_manual(values=as.character(color_values[,i]))+
      scale_fill_manual(values=as.character(color_values[,i]))
    print(i)
    print(p1)
  })}

myplots[[1]]
myplots[[2]]
myplots[[3]]
myplots[[4]]
myplots[[5]]
myplots[[6]]

p1 <- ggplot() + theme_void()+ theme(panel.background = element_rect(colour="white"))
plots<-ggarrange(myplots[[1]],p1,myplots[[2]],p1,myplots[[3]],p1,myplots[[5]],p1,myplots[[4]],p1,myplots[[6]],nrow=11,ncol=1,heights = c(1,0.2,1,0.2,1,0.2,1,0.2,1,0.2,1))
ggsave(plots, file="/data/chamal/projects/aurelie/T2star_MP2RAGE/NEW/Project3/Figure_anova.png", width=15, height=14)

