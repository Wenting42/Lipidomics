#####
#LiPa Viewer_ManuelePiccolis_HSPH_2016_ver.0.6
#######
#Libraries
library(reshape2)
library(ggplot2)
library(grid)
library(shiny)
library(plotrix)
library(png)
library(scales)
#####
#Functions
WNormalize<-function(DFnorm,SumCdata,Sta,nxsat,nxsat2){
  #SumCdata=SumCD
  #DFnorm=Ori
  #print("Please open the file with Standards")
  #Sta=read.delim(file.choose(), header=TRUE)
  a=merge.data.frame(Sta,DFnorm,by.x = "LipidMolecule",by.y = "LipidGroup")
  CorretionValues=a[c(2,5:length(a[1,]))]
  Correctiondf=merge.data.frame(CorretionValues, SumCdata, by.x = "Class.x",by.y="Group.1")
  Scorr=as.matrix(Correctiondf[nxsat])
  Sori=as.matrix(Correctiondf[nxsat2])
  CorValues=matrix(nrow=length(Scorr[,1]))
  CorValues=CorValues[,-1]
  i=1
  while (i<=(length(Scorr[1,]))){
    x=Scorr[,1]/Scorr[,i]
    CorValues=cbind(CorValues,x)
    i=i+1
  }
  NormValues=cbind(Correctiondf[1],as.data.frame(Sori*CorValues))
  colnames(NormValues)<-colnames(SumCdata)
  remove=c(match(NormValues$Group.1,SumCdata$Group.1))
  SumCdata2=SumCdata[-remove,]
  
  NORM=rbind(SumCdata2,NormValues)
  return(NORM)
}
WNormalizeIndividual<-function(Ori,Sdata,Sta,nxsat,nxsat2){
  #print("Please open the file with Standards")
  #Sta=read.delim(file.choose(), header=TRUE)
  a=merge.data.frame(Sta,Ori,by.x = "LipidMolecule",by.y = "LipidGroup")
  CorretionValues=a[c(2,5:length(a[1,]))]
  Correctiondf=merge.data.frame(CorretionValues, Sdata, by.x = "Class.x",by.y="Class")
  Scorr=as.matrix(Correctiondf[nxsat])
  Sori=as.matrix(Correctiondf[nxsat2+1])
  CorValues=matrix(nrow=length(Scorr[,1]))
  CorValues=CorValues[,-1]
  i=1
  while (i<=(length(Scorr[1,]))){
    x=Scorr[,1]/Scorr[,i]
    CorValues=cbind(CorValues,x)
    i=i+1
  }
  NormValues=cbind(Correctiondf[c((2+length(nxsat)),1)],as.data.frame(Sori*CorValues))
  colnames(NormValues)<-colnames(Sdata)
  
  remove=which(Sdata$LipidGroup %in% NormValues$LipidGroup)
  SumCdata2=Sdata[-remove,]
  NORM=rbind(SumCdata2,NormValues)
  return(NORM)
}
WFindpmol<-function(DFnorm,Sta){
  ##DFnorm is the original lisp of lipids from which the function estract and calculate intensity corresponding to 1 pmol
  #DFnorm=Ori
  #print("Please open the file with Standards")
  #Sta=read.delim(file.choose(), header=TRUE)
  a=merge.data.frame(Sta,DFnorm,by.x = "LipidMolecule",by.y = "LipidGroup")
  CorretionValues=a[c(5:length(a[1,]))]
  #Correctiondf=merge.data.frame(CorretionValues, SumCdata, by.x = "Class.x",by.y="Group.1")
  #Scorr=as.matrix(Correctiondf[2:16])
  #Sori=as.matrix(Correctiondf[17:31])
  CorValues=matrix(nrow=length(CorretionValues[,1]))
  CorValues=CorValues[,-1]
  i=1
  while (i<=length(CorretionValues[1,])){
    x=CorretionValues[,i]/a$pmol
    CorValues=cbind(CorValues,x)
    i=i+1
  }
  NormValues=cbind(a[2],as.data.frame(CorValues))
  return(NormValues)
}
WGetSaturated<-function(oriclean){
includeS1="\\(..:[0]\\)|\\(..:[0][a-z]\\)|\\([a-z]..:[0]\\)"
Satdata=oriclean[grepl(includeS1,cleandata$FA1),]
includeS2="\\(..:[123456789]\\)"
Satdata2=Satdata[!grepl(includeS2,Satdata$FA2),]
Satdata3=Satdata2[!grepl(includeS2,Satdata2$FA3),]
return(Satdata3)
}
WGetSingleSaturated<-function(oriclean){
  includeS1="\\(..:[0]\\)|\\(..:[0][a-z]\\)|\\([a-z]..:[0]\\)"
  Satdata=oriclean[grepl(includeS1,cleandata$FA1),]
  includeS3="."
  SatdataM=Satdata[grepl(includeS3,Satdata$FA2),]
  includeS2="\\(..:[0]\\)"
  Satdata2=SatdataM[!grepl(includeS2,SatdataM$FA2),]
  SatdataA=oriclean[!grepl(includeS1,cleandata$FA1),]
  SatdataB=SatdataA[grepl(includeS2,SatdataA$FA2),]
  final=rbind(Satdata2,SatdataB)
  return(final)
}
WCleandata<-function(data){
  nonrej<-data[data$Rej.<1,]
  includeL="\\([123456][02468]:.\\)|\\([123456][02468]:.[a-z]\\)|\\([a-z][123456][02468]:.\\)"
  cleandata=nonrej[grepl(includeL,nonrej$FA1),]
  Oddchain=nonrej[!grepl(includeL,nonrej$FA1),]
  excludeL2="\\(.[13579]:.\\)"
  cleandata2=cleandata[!grepl(excludeL2,cleandata$FA2),]
  cleandata3=cleandata2[!grepl(excludeL2,cleandata2$FA3),]
  return(cleandata3)
}
WTransp_Mean<-function(SumCdata,nrep,nsamp){
  #####Get SD and mean and prepare for plot
  #SumCdata=s[-2]
  tdata=t(SumCdata[-1])
  colnames(tdata)=(SumCdata[,1])
  tdata=as.data.frame(tdata)
  labels=c(LETTERS[seq(from = 1, to = 15 )])
  labels2=c(rep(labels[1:nsamp],each=nrep))
  
  factors=as.data.frame(labels2)
  colnames(factors)="Samples"
  d1<-cbind(factors,tdata)
  
  meanall=aggregate(x = d1[2:length(d1)],by = list(d1$Samples),mean)
  sdall=aggregate(x = d1[2:length(d1)],by = list(d1$Samples),sd)
  alld=cbind(meanall,sdall)
  
  mdata <- melt(meanall, id=c("Group.1"))
  msd <- melt(sdall, id=c("Group.1"))
  rec=cbind(mdata,msd[3])
  colnames(rec)=c("Sample","Lipid","Mean","Sd")
  
  a=which(rec$Sample=="A")
  b=rep(x = a,times = 1,each = length(meanall$Group.1))
  mean_C=rec[b,]
  nrec=cbind(rec,mean_C$Mean)
  return(rec)
}
#PLOT FUNCTIONS
mybarplot_Individual<-function(recIND3){
#bar plot for individual lipids
dodge = position_dodge(width=0.9)
ggplot(data=recIND3,aes(y=Mean,x=Lipid,fill=Sample))+
  geom_bar(position = position_dodge(), 
           stat="identity") +
  geom_errorbar(aes(ymin=Mean-Sd,ymax=Mean+Sd),
                position=dodge,width=0.1,size=0.3)+
  theme_bw() +
  theme(
    #legend.position=c(0.1,0.78),
    axis.text.x = element_text(hjust = 1,color="Black", size=12, angle=45),
    axis.text.y = element_text(color="Black", size=12)
  )+
  xlab("Identified Lipids") +
  ylab("Area Peak") 
  #scale_fill_manual(values=c("Gray", "Red", "Blue","Orange","orangered1"),name="Legend",labels=c("ControlKD","Gpat3KD+G3/4i","GPAT3KD","GPAT3KD+G3/4i","Palmitate0.2mM +Oleate 0.05mM"))
}
mybarplot_all<-function(rec){
ggplot(rec, aes(Sample, Mean)) +
  geom_bar(aes(fill = Sample), position = "dodge", stat="identity")+
  facet_wrap(~Lipid,ncol=8,scales = "free")+
  theme_bw() +
  theme(axis.text.x=element_blank(),strip.background = element_blank()) +
  geom_errorbar(aes(ymin=Mean-Sd, ymax=Mean+Sd), width=.3, color="black")
  #scale_fill_manual(values=c("Gray", "Red", "Blue","Orange","orangered1"),name="Legend",labels=c("ControlKD","GPAT3KD+G3/4i","ControlKD+Palmitate","GPAT3KD+G3/4i+Palmitate","Palmitate0.2mM +Oleate 0.05mM"))
}
circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}
#plot for percentage pmol
pmolplot<-function(pmol_rec){
  ggplot(pmol_rec, aes(Lipid, Mean)) +
  geom_bar(aes(fill = Lipid), position = "dodge", stat="identity")+
  theme_bw() +
  theme(axis.text.x=element_blank(),strip.background = element_blank(), strip.text = element_blank()) +
  facet_wrap(~Sample,nrow=1,)+
  geom_errorbar(aes(ymin=Mean-Sd, ymax=Mean+Sd), width=.15, color="black")+
  ggtitle("% of the Total") + 
    theme(plot.title = element_text(size = 40,lineheight=0.8, face="bold"))+
    ylab("% of Total") +
    scale_y_continuous(labels=percent)
  #scale_fill_manual(values=c("Gray", "Red", "Blue","Orange","orangered1"))
}
cleanandmerge<-function(){
  print("select the negative mode file")
  NEG=read.delim2(file.choose(),sep="\t",row.names = NULL)
  start=(which(NEG[1]=="Rej."))
  del=c(1:(start-1))
  NEG2=NEG[-del,]
  colnames(NEG2)<-c(as.character(NEG[start,]))
  
  
  print("select the positive mode file")
  which(colnames(NEG) %in% "Area.c.1.")
  
  
  POS=read.delim(file.choose(), header=TRUE)
 
  
  
}
#####
#Inputs
#Number of replicates and samples and counters
nrep<-as.numeric(readline("How many replicates?  --> "))
#nrep=2
nsamp<-as.numeric(readline("How many conditions? NB max numer is 15  -->  "))
#nsamp=7
print("Open cleaned and combined file from Lipid Search")
data<-read.csv(file.choose(), header=TRUE)
print("Open Standard file")
Sta=read.delim(file.choose(), header=TRUE)
#Sta=read.delim("Standards_to_use_pmol.txt", header=TRUE)


#Calculates column numbers
start=(which(colnames(data) == "ARatio.s1.c."))
end=(which(colnames(data) %in% "ARatio.s1.c.")+(nsamp*nrep-1))
col_inc=(start:end)
totsamp=end-start
nxsat=(2:(totsamp+2))
nxsat2=(2:(totsamp+2))+totsamp+1


#Analysis_Normalization Areas of all lipids
cleandata<-WCleandata(data)
Sdata=cleandata[c(3,4,col_inc)]
SumCD<-aggregate(Sdata[-c(1,2)],by = list(Sdata$Class),FUN = sum)
nonrej<-data[data$Rej.<1,]
Ori=nonrej[c(3:4,col_inc)]
#WNormalize is only for aggregated lipid species and not for individual ones
Normdf=WNormalize(Ori,SumCD,Sta,nxsat,nxsat2)
rec=WTransp_Mean(Normdf,nrep,nsamp)
#only if I want fewer samples
#rec3=subset.data.frame(x = rec, rec$Sample == ("C") |rec$Sample == ("S1"))
#mybarplot_all(rec)


#Analysis_Normalization Areas of all lipids as ratio
Ratios=as.matrix(Normdf[-1])
Ratios2=Ratios/Ratios[,1]
Ratios3=cbind(Normdf[1],as.data.frame(Ratios2),row.names=NULL)
Ratios4=WTransp_Mean(Ratios3,nrep,nsamp)
#mybarplot_all(Ratios4)


#Analysis_normalization Individaul species (eg.PA=17)
NormdfIND=WNormalizeIndividual(Ori,Sdata,Sta,nxsat,nxsat2)
Normdfsplit=split.data.frame(NormdfIND,NormdfIND$Class)
s=as.data.frame(Normdfsplit[1])
recIND=WTransp_Mean(s[,-2],nrep,nsamp)
#recIND3=subset.data.frame(x = recIND, recIND$Sample == ("C") |recIND$Sample == ("S1") |recIND$Sample == ("S2"))
#mybarplot_Individual(recIND)


#Analysis Double Saturation
Saturation<-WGetSaturated(cleandata)
Sdata2=Saturation[c(3,4,col_inc)]
SumCD2<-aggregate(Sdata2[-c(1,2)],by =list(Sdata2$Class),FUN = sum)
SatNormdf=WNormalize(Ori,SumCD2,Sta,nxsat,nxsat2)
#Satrec=WTransp_Mean(SumCD2)
PerSat=merge.data.frame(Normdf,SatNormdf,by.x = "Group.1",by.y = "Group.1")
PerSatTable=cbind(PerSat[1],(PerSat[nxsat2]/PerSat[nxsat]))
Satrec=WTransp_Mean(PerSatTable,nrep,nsamp)
Satrec3=subset.data.frame(x = Satrec, Satrec$Sample == ("C") |Satrec$Sample == ("S1") |Satrec$Sample == ("S2"))
#mybarplot_all(Satrec)


#Analysis Single Saturation
S_Saturation<-WGetSingleSaturated(cleandata)
Sdata3=S_Saturation[c(3,4,col_inc)]
SumCD3<-aggregate(Sdata3[-c(1,2)],by =list(Sdata3$Class),FUN = sum)
Sat_SNormdf=WNormalize(Ori,SumCD3,Sta,nxsat,nxsat2)
#Satrec=WTransp_Mean(SumCD2)
PerSat_S=merge.data.frame(Normdf,Sat_SNormdf,by.x = "Group.1",by.y = "Group.1")
PerSat_STable=cbind(PerSat_S[1],(PerSat_S[nxsat2]/PerSat_S[nxsat]))
Sat_Srec=WTransp_Mean(PerSat_STable,nrep,nsamp)
#mybarplot_all(Sat_Srec)


#Analysis in absolute amount (pmol)
corrPmol<-WFindpmol(Ori,Sta)
colnames(corrPmol)<-colnames(Normdf)
Pmoldf=merge.data.frame(corrPmol, Normdf, by.x = "Group.1",by.y="Group.1")
unitpmol=as.matrix(Pmoldf[nxsat])
tocorrectOri=as.matrix(Pmoldf[nxsat2])
LinpMol=cbind(Pmoldf[1],as.data.frame(tocorrectOri/unitpmol))
ABSpmol_rec=WTransp_Mean(LinpMol,nrep,nsamp)
#mybarplot_all(ABSpmol_rec)


#Analysis pmol in percentage %
Mpmol=c()
i=2
while (i<=(length(LinpMol[1,]))){
  x=sum(LinpMol[,i])
  Mpmol=cbind(Mpmol,x)
  i=i+1
}
LinpMol2=LinpMol
LinpMol2=LinpMol2[-1]
Mfinal=matrix(ncol=length(LinpMol2[1,]))
Mfinal=Mfinal[-1,]
i=1
while (i<=(length(LinpMol2[,1]))){
  x=LinpMol2[i,]/Mpmol[1,]
  Mfinal=rbind(Mfinal,x)
  i=i+1
}
Perpmol=cbind(LinpMol[1],as.data.frame(Mfinal))
#in percentage pmol
pmol_rec=WTransp_Mean(Perpmol,nrep,nsamp)
#pmolplot(pmol_rec)



#Summary file of all Analyses combined 
Graphic=merge(rec,Satrec, by=c("Sample","Lipid"),all.x=TRUE)
Graphic2=merge(Graphic,ABSpmol_rec, by=c("Sample","Lipid"),all.x=TRUE)
Graphic3=merge(Graphic2,Ratios4, by=c("Sample","Lipid"),all.x=TRUE)
colnames(Graphic3)<-c("Sample","Lipid",
                      "norm_AVG","norm_SD",
                      "sat_AVG","sat_SD",
                      "pmol_AVG","pmol_SD",
                      "Ratio_AVG","Ratio_SD"
                      )
Gmap=read.delim(file="Def_Graphic_Map.txt", header=TRUE)
Graphic4=merge(Gmap,Graphic3,by.x="Lipid",by.y="Lipid",all.x=TRUE)
#Save files to local folder
#write.table(rec, file="1.Areas.txt",sep="\t",quote = FALSE,row.names = FALSE,col.names = TRUE)
#write.table(Satrec, file="2.Saturation.txt",sep="\t",quote = FALSE,row.names = FALSE,col.names = TRUE)
write.table(Graphic4, file="Summary Analysis.txt",sep="\t",quote = FALSE,row.names = FALSE,col.names = TRUE)


##PLOTs
mybarplot_all(rec)
mybarplot_all(Satrec)
mybarplot_all(Sat_Srec)
mybarplot_all(Ratios4)
mybarplot_all(ABSpmol_rec)
pmolplot(pmol_rec)
mybarplot_Individual(recIND)


#Lipa Plot Viewer
ima <- readPNG("lipid_map.png")
ggplot(Graphic4, aes(X, Y),na.rm = TRUE)+
  annotation_custom(rasterGrob(ima,hjust=0.482, vjust=0.545))+
  geom_point(aes(colour = sat_AVG,size=Ratio_AVG), na.rm = TRUE)+
  #geom_point(aes(colour = Mean.y,size=Mean.y), na.rm = TRUE)+
  #geom_point(shape=21, colour = "black", size = 6,fill = "gray", stroke = 1.5)+
  #geom_text(hjust = -1,vjust = -1)+
  scale_size_continuous(range = c(1,10))+
  facet_wrap(~ Sample) +
  xlim(0, 300)+
  ylim(0, 240)+
  coord_fixed(ratio = 1)+
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  scale_colour_gradient(low="Green", high="red")


#Annexes
#Asolutepmol=WTransp_Mean(LinpMol,nrep,nsamp)