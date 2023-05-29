# -------------------------------------------------------------------------
#' 
#' Process data from Stanford database on IC50 and convalescent plasma studies
#' 
# -------------------------------------------------------------------------


# load data ---------------------------------------------------------------

# load concentration data 
Data_mAb<-read.csv("raw-data/Stanford-approved mAbs and adintrevimab -2023-01-17_Updated.csv")
Data_CP<-read.csv("raw-data/Stanford-CP-2023-01-17.csv")


#organise and join mAb and convalescent plasma data from stanford database.
mAbSTUDYLIST<-unique(Data_mAb[,c(1,2)])
Data_mAb_WTonly_full<-Data_mAb[,c("Reference","Assay","Antibodies","Control","Control..IC50.Cmp","Control..IC50.GeoMean")]
Data_mAb_WTonly_full<-Data_mAb_WTonly_full[!is.na(Data_mAb_WTonly_full$Control..IC50.GeoMean),]
Data_mAb_WTonly<-unique(Data_mAb_WTonly_full)
colnames(Data_mAb_WTonly)[colnames(Data_mAb_WTonly)=="Control"]="Variant"
colnames(Data_mAb_WTonly)[colnames(Data_mAb_WTonly)=="Control..IC50.Cmp"]="IC50.Cens"
colnames(Data_mAb_WTonly)[colnames(Data_mAb_WTonly)=="Control..IC50.GeoMean"]="IC50"

Data_mAb_Varaints_only_full<-Data_mAb[,c("Reference","Assay","Antibodies","Variant","Potency..IC50.Cmp","Potency..IC50.GeoMean")]
Data_mAb_Varaints_only_full<-Data_mAb_Varaints_only_full[!is.na(Data_mAb_Varaints_only_full$Potency..IC50.GeoMean),]
Data_mAb_Varaints_only<-unique(Data_mAb_Varaints_only_full)
colnames(Data_mAb_Varaints_only)[colnames(Data_mAb_Varaints_only)=="Potency..IC50.Cmp"]="IC50.Cens"
colnames(Data_mAb_Varaints_only)[colnames(Data_mAb_Varaints_only)=="Potency..IC50.GeoMean"]="IC50"


CPSTUDYLIST<-unique(Data_CP[,c(1,2)])

StudyInmAb<-CPSTUDYLIST$Reference[CPSTUDYLIST$Reference %in% mAbSTUDYLIST$Reference]

## Take only WT 1 month after infection (trhee categories, 1m, 2-6 m, >6m)
Data_CP_WTonly_full<-Data_CP[(Data_CP$Reference %in% StudyInmAb) & Data_CP$Infection..CP.=="Wild Type" & Data_CP$Months=="1m",]
Data_CP_WTonly_full<-Data_CP_WTonly_full[,c("Reference","Assay","Control","Control..NT50.Cmp","Control..NT50.GeoMean")]
Data_CP_WTonly_full$Antibodies<-"Convalescent"
Data_CP_WTonly_full<-Data_CP_WTonly_full[!is.na(Data_CP_WTonly_full$Control..NT50.GeoMean),]
Data_CP_WTonly<-unique(Data_CP_WTonly_full)
colnames(Data_CP_WTonly)[colnames(Data_CP_WTonly)=="Control"]="Variant"
colnames(Data_CP_WTonly)[colnames(Data_CP_WTonly)=="Control..NT50.Cmp"]="IC50.Cens"
colnames(Data_CP_WTonly)[colnames(Data_CP_WTonly)=="Control..NT50.GeoMean"]="IC50"


CombinedDataFull=rbind(Data_mAb_WTonly,Data_mAb_Varaints_only,Data_CP_WTonly)
SplitUpVariant=str_split_fixed(CombinedDataFull$Variant,"/",3)

CombinedDataFull$VariantCategory_High=SplitUpVariant[,1]
CombinedDataFull$VariantCategory_High[grepl("Wild Type",CombinedDataFull$VariantCategory_High)]="Wild Type"
CombinedDataFull$VariantCategory=CombinedDataFull$VariantCategory_High
CombinedDataFull$VariantCategory[CombinedDataFull$VariantCategory=="Omicron"]=paste(CombinedDataFull$VariantCategory_High[CombinedDataFull$VariantCategory_High=="Omicron"],SplitUpVariant[CombinedDataFull$VariantCategory_High=="Omicron",2],sep="/")


##### Remove all censored data with censoring LOD below 10,000 ng/ml
CombinedData<-CombinedDataFull[!(CombinedDataFull$Antibodies!="Convalescent" &
                                   CombinedDataFull$IC50.Cens==">" &
                                   CombinedDataFull$IC50<10000),]

CombinedData<-unique(CombinedData)

###Remove Variants except in list
VariantList=c("Wild Type","Alpha","Beta","Gamma","Delta","Omicron")
CombinedData<-CombinedData[CombinedData$VariantCategory_High %in% VariantList,]

####Remove variants with random additional mutations
CombinedData<-CombinedData[(CombinedData$Variant %in% unique(CombinedData$VariantCategory)) | (CombinedData$VariantCategory_High=="Wild Type"),]


###Remove mAbs not in studies (Prophylaxis and therapies manuscripts)
mAbIncludeList=c("Convalescent",
                 "Bamlanivimab + Etesevimab","Bamlanivimab","Etesevimab",
                 "Casirivimab + Imdevimab","Casirivimab","Imdevimab",
                 "Cilgavimab + Tixagevimab","Cilgavimab","Tixagevimab",
                 "Amubarvimab + Romlusevimab","Amubarvimab","Romlusevimab",
                 "Adintrevimab",
                 "Regdanvimab",
                 "Sotrovimab",
                 "Bebtelovimab",
                 "Bamlanivimab + Bebtelovimab + Etesevimab")

CombinedData<-CombinedData[CombinedData$Antibodies %in% mAbIncludeList,]


###Exclude variant/mAb combinations where there was less than 3 studies.
CombinedData$Group<-paste(CombinedData$Antibodies,CombinedData$VariantCategory)
CombinedDataFull$Group<-paste(CombinedDataFull$Antibodies,CombinedDataFull$VariantCategory)
CombinedDataFulltemp<-CombinedDataFull[CombinedDataFull$Antibodies %in% mAbIncludeList,]
CombinedDataFulltemp<-CombinedDataFulltemp[CombinedDataFulltemp$VariantCategory_High %in% VariantList,]
CombinedDataFulltemp<-CombinedDataFulltemp[(CombinedDataFulltemp$Variant %in% unique(CombinedDataFulltemp$VariantCategory)) | (CombinedDataFulltemp$VariantCategory_High=="Wild Type"),]
VariantsPerStudy<-unique(CombinedDataFulltemp[,c("Reference","Group")])
VariantListwithMoreThan1Study<-names(table(VariantsPerStudy$Group)[table(VariantsPerStudy$Group)>2])

CombinedData<-CombinedData[CombinedData$Group %in% VariantListwithMoreThan1Study,]


#### All values above 10,000 ng/ml censor at 10,000
CombinedData$IC50.Cens[CombinedData$Antibodies!="Convalescent" & CombinedData$IC50>=10000]=">"
CombinedData$IC50[CombinedData$Antibodies!="Convalescent" & CombinedData$IC50>=10000]=10000

CombinedData$IC50_scale=CombinedData$IC50
CombinedData$IC50_scale[CombinedData$Antibodies=="Convalescent"]=1/(CombinedData$IC50[CombinedData$Antibodies=="Convalescent"])
CombinedData$IC50_scale=log(CombinedData$IC50_scale)

###Remove all zeros
CombinedData<-CombinedData[!(CombinedData$IC50==0),]

CombinedData$ReferenceAssay=paste(CombinedData$Reference,CombinedData$Assay)

