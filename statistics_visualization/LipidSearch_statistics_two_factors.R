

#libraries

library(readxl)
library(limma)
library(tidyverse)

#file
file_name="filtered_lipidomics_cosolidate_KWcurated.xlsx"

#read the file
df  <- 
  read_excel(file_name) %>% 
  as.data.frame()

#cosolidate rows with same name by sum


#fix row name
rownames(df)<-df[,1]
df <- df[,-1] %>% as.data.frame()

#remove the rows containing negative values(backgroud)
df <- df[!rowSums(df <= 0), ] 


#log2 transformation
df <- log(df,2)


#design matrix
sample <- c(rep('wt_chow',4), rep('ko_chow',5), 
            rep('wt_fpc',4), rep('ko_fpc',5))
sample <- factor(sample)

design.mat <- model.matrix(~0+sample)

colnames(design.mat) <- levels(sample)

#contrast matrix

contrast.mat <- makeContrasts(
  Chow = "ko_chow - wt_chow",
  Fpc = "ko_fpc - wt_fpc",
  Dietwt ="wt_fpc - wt_chow",
  Diff = "ko_fpc - wt_fpc - ko_chow + wt_chow",
  levels = design.mat
)

#anotherway of doing this, using classical 
#TWO factor ANOVA using model matrix

#genotype <- factor(rep(c('wt','ko', 'wt','ko'), c(4,5,4,5)))
#diet <- factor(rep(c('chow', 'fpc'), c(9,9)))
#design.mat <- model.matrix(~0+genotype*diet)



#Limma package
fit <- lmFit(df, design=design.mat)
fit2 <- contrasts.fit(fit,contrast.mat)
fit3 <- eBayes(fit2)

toptable_OUTPUT1 <- topTable(fit3, coef="Chow", adjust.method = 'fdr',
                             lfc=0, number=nrow(df)) %>% as.data.frame()

toptable_OUTPUT2 <- topTable(fit3, coef="Fpc", adjust.method = 'fdr',
                             lfc=0, number=nrow(df)) %>% as.data.frame()

toptable_OUTPUT3 <- topTable(fit3, coef="Diff", adjust.method = 'fdr',
                             lfc=0, number=nrow(df)) %>% as.data.frame()

toptable_OUTPUT0 <- topTable(fit3, adjust.method = 'fdr',
                             lfc=0, number=nrow(df)) %>% as.data.frame()


toptable_OUTPUT4 <- topTable(fit3, coef="Dietwt", adjust.method = 'fdr',
                             lfc=0, number=nrow(df)) %>% as.data.frame()

#write the output
write.csv (x=toptable_OUTPUT0, file=sprintf("ANOVA_statistics.csv"))

write.csv (x=toptable_OUTPUT1, file=sprintf("Chow_statistics.csv"))

write.csv (x=toptable_OUTPUT2, file=sprintf("Fpc_statistics.csv"))

write.csv (x=toptable_OUTPUT3, file=sprintf("Interaction_statistics.csv"))

write.csv (x=toptable_OUTPUT4, file=sprintf("Diet_wt_statistics.csv"))



#test

deg.ANOVA<- toptable_OUTPUT11 <- topTable(fit3, adjust.method = 'fdr', p.value=0.05, lfc=log2(1.5),
                                         number=nrow(df)) %>% as.data.frame()

deg.chow<- toptable_OUTPUT11 <- topTable(fit3, coef="Chow", adjust.method = 'fdr', p.value=0.05, lfc=log2(1.5),
                             number=nrow(df)) %>% as.data.frame()

deg.fpc<- toptable_OUTPUT12 <- topTable(fit3, coef="Fpc", adjust.method = 'fdr', p.value=0.05, lfc=log2(1.5),
                            number=nrow(df)) %>% as.data.frame()

deg.diff<- toptable_OUTPUT13 <- topTable(fit3, coef="Diff", adjust.method = 'fdr', p.value=0.05, lfc=log2(1.5),
                             number=nrow(df)) %>% as.data.frame()

deg.dietwt<- toptable_OUTPUT14 <- topTable(fit3, coef="Dietwt", adjust.method = 'fdr', p.value=0.05, lfc=log2(1.5),
                                         number=nrow(df)) %>% as.data.frame()


dim(deg.ANOVA)

dim(deg.chow)

dim(deg.fpc)

dim(deg.diff)
