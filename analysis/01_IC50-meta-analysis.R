# -------------------------------------------------------------------------
#' 
#' IC50 meta analysis using data from the Stanford database
#' 
# -------------------------------------------------------------------------


if (run_stanford_db_meta_regression==FALSE){
  load("raw-data/MetaAnalysisCalibration_Stanford_20230221_h09_m28_s22.RData")
  Means_table$Antibodies <- tolower(Means_table$Antibodies)
} else if (run_stanford_db_meta_regression==TRUE) {
  
  
  y = (CombinedData$IC50_scale) #LogN Transformed data
  ylog=-y
  Threshold=(10000)
  
  
  
  cens = 1*(CombinedData$IC50.Cens==">") #0 if real, 1 if it's censored (ie, 0 values)
  cluster = match(CombinedData$ReferenceAssay,sort(unique(CombinedData$ReferenceAssay))) #ID for each patient
  
  
  ListOfGroups=sort(unique(CombinedData$Group))
  TableOfGroups<-unique(CombinedData[,c("Antibodies","VariantCategory","Group")])[match(ListOfGroups,unique(CombinedData[,c("Antibodies","VariantCategory","Group")])$Group),]
  fixed_effect_Group=match(CombinedData$Group,ListOfGroups)
  
  # fixed_effect_mAb = match(CombinedData$Antibodies,unique(CombinedData$Antibodies)) #ID for each patient
  # fixed_effect_Variant = match(CombinedData$VariantCategory,unique(CombinedData$VariantCategory))-1 #ID for each patient
  id = unique(cluster)
  
  
  #Initialize design matrix for the random effects
  #A is for the random intercept
  # A2 <- matrix(0,length(cluster),ncol=length(unique(MetaBA1_2_reshaped_noNA_limitAb$variants_cat)))
  # A2<-matrix(0,length(cluster),ncol=3)
  # A2[,1]<-1
  
  A <- matrix(1, nrow=NROW(cluster), ncol=1)
  
  # for (i in 1:length(cluster)){
  #   if (fixed_effect_Variant[i]!=0) {
  #     A2[i,1+fixed_effect_Variant[i]]=cluster[i]
  #   }
  # }
  
  #B is for the random slope
  # M<-A*day #to tell lmec that we need random slopes
  # B <- matrix(NA, nrow=NROW(cluster), ncol=1)
  B <- A
  # B2 <- A2
  # B[, 2] <- M
  
  # fixed_effect_matrix<-matrix(0,length(cluster),ncol=length(unique(MetaBA1_2_reshaped_noNA_limitAb$ab_name))+length(unique(MetaBA1_2_reshaped_noNA_limitAb$variants_cat)))
  
  fixed_effect_matrix0<-matrix(0,length(cluster),ncol=length(ListOfGroups))
  # fixed_effect_matrix0[,1]<-1
  
  
  
  
  for (i in 1:length(cluster)){
    
    fixed_effect_matrix0[i,fixed_effect_Group[i]]=1
    
  }
  

X0 = fixed_effect_matrix0

fit3_nocens = lme(-IC50_scale~Group+0,random=(~1|ReferenceAssay),data=CombinedData)

fit3 = lmec(yL=ylog,cens=cens, X=X0, Z=B, cluster=cluster, method='ML',maxstep =10000,init=list("beta"=fixed.effects(fit3_nocens),"bi"=t(random.effects(fit3_nocens)[,1])))
# save(fit3,file=paste("metafit_",format(Sys.time(), "%Y%m%d_h%H_m%M_s%S"),".RData",sep=""))

parameters=TableOfGroups
parameters$log_IC50=-fit3$beta
parameters$IC50=exp(parameters$log_IC50)

sizeR=nrow(fit3$varFix)
sizeC=ncol(fit3$varFix)
contrast=matrix(0,nrow=sizeR,ncol=sizeC)
diag(contrast)=1

# parameters_SE=diag(sqrt(t(contrast) %*% (fit3$varFix %*% contrast)))
parameters_SE=sqrt(diag(fit3$varFix))

parameters$log_IC50_SE=parameters_SE

Means_table=parameters
Means_table$IC50_Lower=exp((Means_table$log_IC50)-1.96*Means_table$log_IC50_SE)
Means_table$IC50_Upper=exp((Means_table$log_IC50)+1.96*Means_table$log_IC50_SE)
Means_table$IC50[Means_table$Antibodies=="Convalescent"]=1/Means_table$IC50[Means_table$Antibodies=="Convalescent"]
Means_table$IC50_Lower[Means_table$Antibodies=="Convalescent"]=1/Means_table$IC50_Lower[Means_table$Antibodies=="Convalescent"]
Means_table$IC50_Upper[Means_table$Antibodies=="Convalescent"]=1/Means_table$IC50_Upper[Means_table$Antibodies=="Convalescent"]




Means_table$text=paste(ifelse(Means_table$IC50>10000,">10000",signif(Means_table$IC50,3)),
                       "\n(95% CI: ",
                       ifelse(Means_table$IC50_Lower>10000,">10000",signif(Means_table$IC50_Lower,3)),
                       " - ",
                       ifelse(Means_table$IC50_Upper>10000,">10000",signif(Means_table$IC50_Upper,3)),")",sep='')


outputtable<-acast(Means_table,Antibodies~VariantCategory,function(x){paste(x,collapse="")},value.var="text")

write.csv(outputtable,"output/IC50_Table_Stanford_v2.csv")

save(Means_table,file=paste("output/MetaAnalysisCalibration_Stanford_",format(Sys.time(), "%Y%m%d_h%H_m%M_s%S"),".RData",sep=""))

Means_table$Antibodies <- tolower(Means_table$Antibodies)


plot_ic50_regression<-ggplot(CombinedData,aes(x=tolower(Antibodies),y=IC50,color=(Antibodies=="Convalescent")))+
  # geom_boxplot(outlier.shape = NA)+
  geom_point(position=position_dodge2(width=0.3))+
  geom_point(data=Means_table,aes(y=pmin(10000,IC50),shape=(IC50<=10000),color=(Antibodies=="convalescent")),fill="black",size=3,stroke=1)+
  geom_errorbar(data=Means_table[Means_table$IC50<10000,],aes(ymin=pmin(10000,IC50_Lower),ymax=pmin(10000,IC50_Upper),color=(Antibodies=="convalescent")))+
  scale_y_log10(breaks=c(0.1,1,10,100,1000,10000),labels=c("0.1","1","10","100","1,000","10,000"))+
  scale_shape_manual(values=c(16,1),guide="none")+
  scale_color_manual(values=c("black","dodgerblue"),guide="none")+
  coord_cartesian(ylim=c(0.1,10000))+
  theme_linedraw() +
  labs(y=("<span style='font-size:14pt;'><span style='color:#000000;'>IC50 (ng/ml)</span>
       <span style='color:#1E90FF;'><br>(or Neutralization titer in the case of the convalescent plasma)</span></span>"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(colour = "gray95"),
        panel.grid.minor = element_line(colour = "gray95"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        # panel.spacing.y = unit(-2.5,"lines"),
        strip.background = element_rect(fill="white",color=NA),
        strip.text = element_text(colour="black"),
        legend.spacing.y = unit(0, 'cm'),
        legend.key = element_rect(size = 1),
        legend.key.size = unit(0.45, "cm"),
        legend.margin = margin(t=0.8,b=-0.4,unit="cm"),
        legend.background = element_blank(),
        axis.text.x=element_text(angle=45,hjust=1),
        axis.title.x=element_blank(),
        axis.title.y=element_markdown())+
  facet_wrap(~VariantCategory,ncol=3)
plot_ic50_regression

ggsave(plot_ic50_regression,file="output/FigS1_StanfordMetaAnalysis.pdf",width=14,height=18)

  # cleanup:
  rm(y,logy,Threshold,cens,cluster,ListOfGroups,TableOfGroups,fixed_effect_Group,id,A,B,fixed_effect_matrix0,i,
     X0,fit3_nocens,fit3,parameters,sizeR,sizeC,contrast,parameters_SE,outputtable)

}