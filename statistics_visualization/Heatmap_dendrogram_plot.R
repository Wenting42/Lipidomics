library(tidyverse)
library(broom)
library(readxl)
library (ggdendro)
library (plotly)
library(heatmap.plus)
library(RColorBrewer)
library(gplots)


#input
file_name1="Interaction_statistics.xlsx"
file_name2="filtered_lipidomics_cosolidate_KWcurated.xlsx"


#read the file and change the first row to ID
df1<- 
  read_excel(file_name1) %>% 
  as.data.frame() 
colnames(df1)[1] <- 'ID'

df2<- 
  read_excel(file_name2) %>% 
  as.data.frame()
colnames(df2)[1] <- 'ID'


#Combine files

df <- merge(df2, df1, by='ID')

# only cut off p >0.05
df_sig <- subset(df,adj.P.Val <= 0.05, select=ID:MainArea.s22.)

#remove column1 ID, and copy it to the colnames.
rownames(df_sig) <- df_sig[,1]
df_sig <- df_sig[,-1]


#transpose and scaling
df_t <- t(df_sig) %>% as.data.frame()

df_scaled <- as.matrix(scale(df_t))


#subset teh data

df_neutral <- df_scaled[,grep(('[TD]G'), colnames(df_scaled))]
df_phospho <- df_scaled[,grep(('(P[ACEGSI]|LP[CEGSI]|CL)'), colnames(df_scaled))]
df_spingo <- df_scaled[,grep(('(Cer|CerG1|SM|ChE)'), colnames(df_scaled))]


Condition_colors <-  unlist(lapply(rownames(df_sig), function(x){
  if (grepl('(TG|DG)', x)) '#FFFF00' #yellow
  else if (grepl('(P[ACEGSI]|LP[CEGSI]|CL)', x)) '#006400'
  else '#9933FF'
  
  
}))

Condition_colors
#input <- as.matrix(t(df_sig))
#input <- as.matrix(df_sig)
dev.off()
graphics.off()



pdf(file='E:/Dropbox/Farese Walther Lab/Collaboration/20180826 Nina-DGAT2/Output/Heatmap/interactionlipids.pdf')

#heatmap.2(input, trace = 'none', density='none', col=redblue(100), cexRow = 1, cexCol = 0.2, margins = c(20, 13),
          #ColSideColors = Condition_colors, scale = "column",
          #hclust=function(x) hclust(x, method='average'))


heatmap.2(df_scaled, trace = 'none', density='none', col=bluered(100), cexRow = 1, cexCol = 0.2, margins = c(20, 13),
          ColSideColors = Condition_colors, 
          hclust=function(x) hclust(x, method='average'))

legend(0.75, 0.95, legend = c('Neutral lipids', 'Glycophospholipids', 'Sphingolipids $ sterols'), 
       fill=c('#FFFF00', '#006400','#9933FF'), cex=0.5)

dev.off()

#heatmap(df_scaled, Colv=F, scale='none')





