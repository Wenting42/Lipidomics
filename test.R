#rm(list=ls())  ##### empty the environment first. If you need to keep former variables in workspace, 
#just comment this part
#### Please install the lacking packages ------> install.packages("package")
library("FactoMineR")
library("tidyverse")
library("ggplot2")
library("magrittr")
library("ggrepel")
library("reshape2")
library("stargazer")
library("Hmisc")
library("limma")
library("factoextra")
library("scales")
library("RColorBrewer")
library("stringr")
library("readxl")

# Read in file, and store the data into tible structure for the tidyverse to apply
lipidomics <- read_csv(file.choose())

# reformat the colnames for better dealing

#colnames(lipidomics) <- make.names(colnames(lipidomics))  

message("If your samples contains more control groups, you need to remove them for further analysis.")

############################################### 
# interactive way

extracted.Grade <- as.character(readline("extracted.Grade -----> "))
unextracted.Grade <- as.character(readline("unextracted.Grade -----> "))

message("Please input the extrated sample names, e.g. APValue[s17]")
extracted.P.No <- as.character(readline("extracted.P.No -----> "))




message("Please input the unextrated sample names, e.g. APValue[s18]")
unextracted.P.No <- as.character(readline("unextracted.P.No -----> "))

###### if the file has extracted and unextracted samples
### either -contains(colnumber of these samples) or -contains(name of the samples) // contains(pattern)
# eg. 
# extracted.P.No <- "APValue[s17]"   #### if you want interactive way, please comment this line with "#"
# 
# unextracted.P.No <- "APValue[s18]"  #### if you want interactive way, please comment this line with "#"



# select and count the sample grades you want, filter the p value which less than 0.001
lipidCount <- lipidomics %>% 
  select(contains("Grade"), -contains("Grade[c]" ), -contains(extracted.Grade), -contains(unextracted.Grade),
         contains("APValue"),-contains(extracted.P.No), -contains(unextracted.P.No))  %>% 
  transmute(A=rowSums(.=="A"), B=rowSums(.=="B"), 
            C=rowSums(.=="C"), D=rowSums(.=="D"), 
            No.grade=rowSums(.=="-"), APvalue.001=rowSums(. <= "0.001"))

# add this count into original strucutes
lipidSelect <- lipidomics %>% bind_cols(lipidCount)

# Filter the dataset based on your criteria  
# flexible parameter 4 which depends on the experiment for total number of A and B 
filtered.lipidomics <- lipidSelect %>% 
  rowwise() %>% 
  filter( Rej == 0 &
            sum(A, B) >= 4 &
            APvalue.001 >= 4)

# Check how many lipids per class that passed the filtering
#summary(as.data.frame(filtered_lipidomics$Class), maxsum=nrow(filtered_lipidomics))
# by(filtered_lipidomics, filtered_lipidomics$Class, summary)

describe(filtered.lipidomics$Class)
# stargazer(summary(as.data.frame(filtered_lipidomics$Class), maxsum=nrow(filtered_lipidomics)), type="latex",align=TRUE)
#aa <- filtered_lipidomics %>% rowwise %>% filter(Class=="PS")


# QC PLOT 1 - Abundance vs. retention time [requires ggplot2]
##############################################################################################
# interactive way
### e.g. abundanceSample -----> "log10(MainArea.s1._DMSO)"  ## don't forget the quote and log10
##############################################################
message("Please input the MainArea of the sample name to check its abundance vs. retention time, e.g. `MainArea[s1]`. 
        Please note that ` ` is the backtick which is same position as tilde")

abundance.sample <- readline("abundance.sample -----> ")


abundance.sample <- "`MainArea[s4]`"   # for running convinence 

ggplot(data=filtered.lipidomics, aes_string(x = sprintf("log10(%s)", abundance.sample), y = "BaseRt")) +
  geom_point() +
  theme_bw() +
  facet_grid(.~Class) +
  labs("Abundance VS. Retention time", x="lipid class (log10(MainArea))", y="Retention time")
message("Please input the name you want to store for the graph. e.g. retention.pdf")
plot.name <- readline("QC plot name: ")
ggsave(filename = plot.name, path = 'plot/', device = "pdf")

############################################################################







##### if you don't want interactive way, just uncomment below part, fill the x value and 
##### comment the interactive part using -----> command/Ctrl + shift + c (mac/windows, linux)
#
# ggplot(data=filtered_lipidomics, aes(x = log10(MainArea.s1._DMSO), y = BaseRt)) +
#   geom_point() +
#   theme_bw() +
#   facet_grid(.~Class)
#
#############################################################



# making group
message("Please input the info of the experiment groups below and end with 'Enter'
        Please Only input the sample groups for analysis.
        e.g. Group number: 2   
        Name of Group 1: wild type
        Sample names of Group 1: s1 s2 s3 s4
        Name of Group 2: knockout
        Sample names of Group 2: s4 s5 s6 s7")

ngroups <- readline("Group number: ")

inputGroups <- function(n){
  group.info <- c()
  sample.info <- c()
  info <- list()
  for(i in 1:n){
    ms1 <- paste("Name of Group ", i, ": ")
    group.info[i] <- readline(prompt= ms1)
    ms2 <- paste("Sample names of Group ", i, ": ")
    sample.info[i] <- readline(prompt=ms2)
    info[[group.info[i]]] <- sample.info[i]
    
  }
  return(info)
}


total.info <- inputGroups(ngroups)

# Grep index and sample names 

grepIndex <- function(sample, data){
  sample.names <- list()
  group.names <- c()
  col.index <- list()
  for( i in 1:length(names(sample))){
    group.names[i] <- names(sample)[i]
    group.sample.names <- sample[[i]] %>% strsplit(., "\\s+") %>% unlist()
    
    sample.names[[group.names[i]]] <- paste("MainArea[", group.sample.names, "]", sep="")
    
    col.index[[group.names[i]]] <- which(colnames(data) %in% sample.names[[i]])
    
  }
  info <- rbind(sample.names, col.index)
  return(info)
}

# store the sample and index in a list
sample.info <- grepIndex(total.info, filtered.lipidomics) 

group.names <- dimnames(sample.info)[[2]]



# Preparation for pair-wise correlations
inf2NA <- function(x) { x[is.infinite(x)] <- NA; x }

panel.cor <- function(x, y, digits = 2, cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  # correlation coefficient
  r <- cor(inf2NA(x), inf2NA(y), use = "pairwise.complete.obs", method = "pearson")
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste("r= ", txt, sep = "")
  text(0.5, 0.6, txt)
  
  # p-value calculation
  p <- cor.test(x, y)$p.value
  txt2 <- format(c(p, 0.123456789), digits = digits)[1]
  txt2 <- paste("p= ", txt2, sep = "")
  if(is.na(p<0.01)) txt2 <- paste("p= ", "<0.01", sep = "")
  text(0.5, 0.4, txt2)
}

panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "#4393C3", ...)
}




# QC PLOT 2 - Pair-wise correlation between replicates
index <- 1:length(sample.info)
sample.index <- index[index %% 2 ==0]


for(i in 1:length(sample.index)){
  range <- sample.index[i]
  plot.name <- paste("pairs.plot.", i, ".pdf",sep="")
  path <- file.path("plot/", plot.name)
  pdf(file=path)
  pairs(log10(filtered.lipidomics[, sample.info[[range]]]), 
        lower.panel=panel.smooth, diag.panel=panel.hist, 
        upper.panel = panel.cor, xlim=c(4,12), ylim=c(4,12))
  dev.off()
}


index.list <- index[!index %% 2 ==0]
sample.list <- c()
for(i in length(index.list):1){
  j <- index.list[i]
  sample.list <- c(sample.info[[j]], sample.list)
}

group.repeats <- c()
for(i in length(group.names):1){
  k <- index.list[i]
  len <- length(sample.info[[k]])
  repeats <- rep(group.names[i], len)
  group.repeats <- c(repeats, group.repeats)
}


# Formatting the table for PCA
filtered.lipidomics.PCA <-  filtered.lipidomics %>% 
  select(sample.list) %>% 
  t %>% 
  as.data.frame()

colnames(filtered.lipidomics.PCA) <- filtered.lipidomics$LipidMolec                           

log2.filtered.lipidomics.PCA <- log((filtered.lipidomics.PCA+1), 2)



############# naming the group 
# group.repeats <- c()
# times <- c()
# repeats <- list()
# for (i in 1:ngroups){
#   group.names[i] <- dimnames(sample.info)[[2]][i]
#   range <- index.list[i]
#   times <- length(sample.info[[range]])
#   print(times)
#   repeats[i] <- c(rep(group.names[i], times)) 
#   a <- repeats[i] %>% unlist()
#   group.repeats <- c(group.repeats, a) 
# }
#group.repeats <- c(rep("H7.G1" , 7), rep("H15.G2" , 7), rep("L19.G3", 3), rep("L21.G4", 3))

filtered.lipidomics.PCA$Group <- group.repeats
log2.filtered.lipidomics.PCA$Group <- group.repeats











# Perform PCA [requires FactoMineR]
res.pca <-  PCA(filtered.lipidomics.PCA, scale.unit=TRUE, ncp=5, quali.sup=ncol(filtered.lipidomics.PCA), graph=T)

concat <-  cbind.data.frame(filtered.lipidomics.PCA[, ncol(filtered.lipidomics.PCA)], res.pca$ind$coord)

ellipse.coord = coord.ellipse(concat, bary=TRUE)
pdf(file="plot/pca_wo_s16_s3.pdf")
plot.PCA(res.pca, habillage=ncol(filtered.lipidomics.PCA), ellipse=ellipse.coord, cex=0.8, label="all")
dev.off()

pdf(file="plot/pca.pdf")
plot.PCA(res.pca, habillage=ncol(filtered.lipidomics.PCA), ellipse=ellipse.coord, cex=0.8, label="quali")
dev.off()


# # proportion of variance explained by components
# pdf(file="plot/scree.plot.pdf")
# fviz_screeplot(res.pca, ncp=10, ylim=c(0, 45), xlab="Principle Components")
# dev.off()


# Output the filtered lipidome
write.csv(filtered.lipidomics, "data/filtered.lipidomics.csv")




### making data for barplot graph


############ check if need to remove the standards for these part 

# standards <- c("TG(17:1/17:1/17:1)","PS(17:0/20:4)","PI(17:0/20:4)","PG(17:0/14:1)","PE(17:0/14:1)","LPS(17:1)","LPE(17:1)","LPE(17:1)","LPC(17:1)","DG(19:0/19:0)","CerG1(d18:1/12:0)","Cer(d18:1/17:0)","CL(14:1/14:1/15:1/14:1)")


filtered.groups <- filtered.lipidomics %>% select(Class, sample.list) 
filtered.groups <- aggregate(.~Class, data=filtered.groups, FUN=sum)

data <- filtered.groups[,-1] %>% t %>% data.frame(group.repeats, .)
names(data) <- c("experiment.group",filtered.groups$Class) 



data.sd <- aggregate(.~experiment.group, data=data, function(x) sd(x)) 
data.sd <- data.sd %>% gather( class, sd, -experiment.group)
  
  
# b <- data %>% group_by(experiment.group) %>% mutate_if(is.numeric, sd)
# 
# sd <- data %>% group_by(experiment.group) %>% mutate(.funs=)
# 
# kk <- data %>% group_by(experiment.group) %>% summarise_all(.funs=sum)
# data.sd <- lapply(kk[,-1], function(x)sd(unlist(x)))
# #sd <- data.frame(data.sd, experiment.group="sd")
# #dd <-  data.sd %>% unlist() %>% as.numeric()
# dd <- rbind(as.numeric(unlist(data.sd)), names(data.sd), deparse.level=0) %>% t() %>% as.data.frame()
# colnames(dd) <- c("sd", "class")
# dd$sd <- as.numeric(as.character(dd$sd))
# dd$class <- as.character(dd$class)

wide.data <- data  %>% group_by(experiment.group) %>% summarise_all(.funs=sum) 
# 
data.long <- wide.data %>% gather(class, value, -experiment.group)

data.long <- bind_cols(data.long, sd=data.sd$sd)

ggplot(data.long, aes(x=class, fill=experiment.group, y=value)) + 
  geom_bar(position="dodge", stat="identity") +
  #ylim(10^5, Inf) +
  geom_errorbar(aes(ymin=value, ymax=value+sd), width=.3,
                position=position_dodge(.9)) +
  theme_bw() +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
               labels = trans_format("log10", math_format(10^.x))) +
  # scale_y_continuous(labels = scales::scientific) +
  labs(ylab="Main Area", fill="Groups") +
  ggtitle("Total lipid classes") +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_text(angle=45, hjust=1))+
  coord_cartesian(ylim=c(1e5,1e13))



message("Please input the name you want to store for the graph. e.g. retention.pdf")
name.class  <- readline("Total lipid classes plot name: ")
ggsave(filename = name.class, path = 'plot/', device = "pdf")





###### for each class
message("How many classes you want to visualize for barplots?")
ntimes <- readline("Please input the numbers of classes you want to visualize: ")



message("which class you want to visualize for barplot?")
pick.class <- as.character(readline("Please input the class for barplot visualization: "))


filtered.class <- filtered.lipidomics %>% 
  select(Class, LipidMolec, sample.list) %>% 
  filter(grepl(pick.class, Class))

class.data <- aggregate(.~Class+LipidMolec, data=filtered.class, FUN=sum) %>% select(-Class)
data.wide <- data.frame(t(class.data[,-1]), group.repeats)
#molec.names <- class.data[,1] %>% unlist() %>% str_extract_all(., "\\(\\d*.*\\)") %>% str_remove_all(., "[\\(\\)]") 

names(data.wide) <- c(class.data[,1], "experiment.group")


class.sd <- aggregate(.~experiment.group, data=data.wide, function(x) sd(x)) 
class.sd <- class.sd %>% gather( Aceyl, sd, -experiment.group)


class.long <- data.wide %>% group_by(experiment.group) %>% summarise_all(funs(sum)) 
# long.sd <- apply(class.long[,-1], 2, sd)
# 
# df <- data.frame(as.numeric(unlist(long.sd)), names(long.sd))
# colnames(df) <- c("sd", "Aceyl")
# df$Aceyl <- as.character(df$Aceyl)

class.data.long <- class.long %>% gather(Aceyl, AreaValue, -experiment.group)

class.data.long <- bind_cols(class.data.long, sd=class.sd$sd)


# class.data.long <- full_join(class.data.long, df, by="Aceyl")


ggplot(class.data.long, aes(x=Aceyl, fill=experiment.group, y=AreaValue)) + 
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin=AreaValue, ymax=AreaValue+sd), width=.8,
                position=position_dodge(.9)) +
  theme_bw() +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),      ###### base of 10
               labels = trans_format("log10", math_format(10^.x))) +
  #scale_y_continuous(labels = scales::scientific)+                      ###### scientific notation
  labs(ylab="Main Area", fill="Groups", xlab="AT") +
  ggtitle(pick.class) +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle=45, hjust=1)) +
  xlab("Aceyl composition") +
  coord_cartesian(ylim=c(1e5,1e11))


message("Please input the name you want to store for the graph. e.g. TG.pdf")
class.type <- readline("The lipid class plot name: ")
ggsave(filename = class.type, path = 'plot/', device = "pdf")







# original version without error bars
# filtered.groups <- filtered.lipidomics %>% select(Class, sample.list) 
# filtered.groups <- aggregate(.~Class, data=filtered.groups, FUN=sum)
# 
# data <- filtered.groups[,-1] %>% t %>% data.frame(., group.repeats)
# names(data) <- c(filtered.groups$Class, "experiment.group") 
# 
# data <- data  %>% group_by(experiment.group) %>% summarise_all(.funs=sum) 
# 
# data.long <- data %>% gather(class, value, -experiment.group)
# 
# ggplot(data.long, aes(x=class, fill=experiment.group, y=value)) + geom_bar(position="dodge", stat="identity") +
#   theme_bw() +
#   scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
#                 labels = trans_format("log10", math_format(10^.x))) +
#   #scale_y_continuous(labels = scales::scientific) +
#   labs(ylab="Main Area", fill="Groups") +
#   ggtitle("Total lipid classes") +
#   theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle=45, hjust=1))
# 
# message("Please input the name you want to store for the graph. e.g. retention.pdf")
# name.class  <- readline("Total lipid classes plot name: ")
# ggsave(filename = name.class, path = 'plot/', device = "pdf")
# 
# 



###### for each class
# message("How many classes you want to visualize for barplots?")
# ntimes <- readline("Please input the numbers of classes you want to visualize: ")
# 
# 
# 
# message("which class you want to visualize for barplot?")
# pick.class <- as.character(readline("Please input the class for barplot visualization: "))
# 
# filtered.class <- filtered.lipidomics %>% 
#   select(Class, LipidMolec, sample.list) %>% 
#   filter(grepl(pick.class, Class))
# 
# class.data <- aggregate(.~Class+LipidMolec, data=filtered.class, FUN=sum) %>% select(-Class)
# data.wide <- data.frame(t(class.data[,-1]), group.repeats)
# molec.names <- class.data[,1] %>% unlist() %>% str_extract_all(., "\\(\\d*.*\\)") %>% str_remove_all(., "[\\(\\)]") 
# 
# names(data.wide) <- c(molec.names, "experiment.group")
# class.long <- data.wide %>% group_by(experiment.group) %>% summarise_all(funs(sum)) 
# class.data.long <- class.long %>% gather(Aceyl, AreaValue, -experiment.group)
# 
# ggplot(class.data.long, aes(x=Aceyl, fill=experiment.group, y=AreaValue)) + 
#   geom_bar(position="dodge", stat="identity") +
#   theme_bw() +
#   # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),      ###### base of 10
#   #              labels = trans_format("log10", math_format(10^.x))) +
#   scale_y_continuous(labels = scales::scientific)+                      ###### scientific notation
#   labs(ylab="Main Area", fill="Groups", xlab="AT") +
#   ggtitle(pick.class) +
#   theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle=45, hjust=1)) +
#   xlab("Aceyl composition")
# 
# message("Please input the name you want to store for the graph. e.g. TG.pdf")
# class.type <- readline("The lipid class plot name: ")
# ggsave(filename = class.type, path = 'plot/', device = "pdf")
# 
# 
# 



##### STATISTICS WITH LIMMA FOR INDIVIDUAL LIPIDS #####

# Create a design matrix 

samples <- factor(group.repeats)

design <- model.matrix(~0+samples)

colnames(design) <- levels(samples)

# Reformat the filtered_lipidomics_PCA, consolidate duplicated lipids and log2+1 transform [requires magrittr, dplyr]
filtered.lipids <- filtered.lipidomics %>% select(Class, LipidMolec, sample.list)

filtered.lipids <- aggregate(. ~LipidMolec+Class, data=filtered.lipids, FUN=sum)

rownames(filtered.lipids) <- filtered.lipids$LipidMolec

filtered.lipids <- filtered.lipids %>% select(-LipidMolec, -Class)

log2.filtered.lipids <- log((filtered.lipids+1), 2)

filtered.lipids %>% write_csv('data/filtered.lipids.csv')

# Fit model and extract contrasts [requires limma]
fit <- lmFit(log2.filtered.lipids, design)

###################################################################################################
message("Please note that you NEED make contrast groups manually to Compare the difference between/among groups.
        e.g. compare B+C against A: A-(B+C)/2,
        A against B:  A-B,
        A+B against C+D: (A+B)/2-(C+D)/2, etc.")

cont.matrix <- makeContrasts(
  G2_vs_G1_H = H15.G2-H7.G1,
  G3_vs_G4_L = L19.G3- L21.G4,
  L_H   = (L19.G3- L21.G4)-(H15.G2-H7.G1),
  levels=design)
##################################################################################################
cont.matrix <- makeContrasts(
  G2_G1 = H15.G2-H7.G1,
  G2_G3 = H15.G2-L19.G3,
  levels=design
)

# fit3 <- eBayes(fit)
# anova.groups <- toptable_OUTPUT0 <- topTable(fit3, adjust.method = 'fdr', 
#                                              lfc=0, number=nrow(filtered.lipids)) %>% as.data.frame()

fit2 <- contrasts.fit(fit, cont.matrix)

fit2 <- eBayes(fit2)

#G2.3 <- topTable(fit2, coef="G2_G3",adjust.method = 'fdr', 
#                 lfc=0, number=nrow(filtered.lipids)) %>% as.data.frame()



# Use topTable to obtain significantly altered lipids
# Diff <- topTable(fit2, coef="Diff",adjust.method = 'fdr', p.value=0.05, lfc=0,
#                  number=nrow(df)) %>% as.data.frame()
# 
# refed_vs_fasted_in_WT <- topTable(fit2, coef="GPAT3_4i.vs.DMSO.WT", adjust.method = "fdr",
#                                   lfc=0, number=nrow(filtered.lipids))
# 
# refed_vs_fasted_in_KO <- topTable(fit2, coef="DGAT2_GPAT3_4i.vs.DGAT2i.KO", adjust.method = "fdr",
#                                   lfc=0, number=nrow(filtered.lipids))
# 
# ko_vs_wt_in_fasted <- topTable(fit2, coef="DGAT2i.ko.vs.DMSO.fasted", adjust.method = "fdr",
#                                lfc=0, number=nrow(filtered.lipids))
# 
# ko_vs_wt_in_refed <- topTable(fit2, coef="ko_vs_wt_in_refed", adjust.method = "fdr",
#                               lfc=0, number=nrow(filtered.lipids))

anova.groups <- toptable_OUTPUT0 <- topTable(fit2, adjust.method = 'fdr',
                                             lfc=0, number=nrow(filtered.lipids)) %>% as.data.frame()
# G2.G1 <- toptable_OUTPUT1 <- topTable(fit2, coef="G2_G1", adjust.method = 'fdr',
#                                         lfc=0, number=nrow(filtered.lipids)) %>% as.data.frame()
# 
# G2.G3 <- toptable_OUTPUT2 <- topTable(fit2, coef="G2_G3", adjust.method = 'fdr',
#                                         lfc=0, number=nrow(filtered.lipids)) %>% as.data.frame()
# 


###################################
G2.G1.H <- toptable_OUTPUT1 <- topTable(fit2, coef="G2_vs_G1_H", adjust.method = 'fdr',
                                        lfc=0, number=nrow(filtered.lipids)) %>% as.data.frame()

G3.G4.L <- toptable_OUTPUT2 <- topTable(fit2, coef="G3_vs_G4_L", adjust.method = 'fdr',
                                        lfc=0, number=nrow(filtered.lipids)) %>% as.data.frame()

HL.dif <- toptable_OUTPUT3 <- topTable(fit2, coef="L_H", adjust.method = 'fdr',
                                       lfc=0, number=nrow(filtered.lipids)) %>% as.data.frame()


# diet.wt <- toptable_OUTPUT4 <- topTable(fit2, coef="Dietwt", adjust.method = 'fdr',
#                              lfc=0, number=nrow(filtered.lipids)) %>% as.data.frame()
# 
# 
# ko.wt.diff <- toptable_OUTPUT5 <- topTable(fit2, coef="Diff", adjust.method = 'fdr',
#                              lfc=0, number=nrow(filtered.lipids)) %>% as.data.frame()
# 




write.csv (x=toptable_OUTPUT0, file=sprintf("data/ANOVA_statistics.csv"))

write.csv (x=toptable_OUTPUT1, file=sprintf("data/G2.G1_statistics.csv"))

write.csv (x=toptable_OUTPUT2, file=sprintf("data/G3.G4_statistics.csv"))

write.csv (x=toptable_OUTPUT3, file=sprintf("data/H.L_statistics.csv"))

# write.csv (x=toptable_OUTPUT4, file=sprintf("data/Diet_wt_statistics.csv"))
# 
# write.csv (x=toptable_OUTPUT5, file=sprintf("data/Interaction_statistics.csv"))

#fold.change <- as.numeric(readline("Please input the threshold for fold change: "))

test.ANOVA<- toptable_OUTPUT0.1 <- topTable(fit2, adjust.method = 'fdr', p.value=0.05, lfc=log2(1.5),
                                            number=nrow(filtered.lipids)) %>% as.data.frame()
test.G2.G1.dif <- toptable_OUTPUT1.1 <- topTable(fit2, coef="G2_vs_G1_H", adjust.method = 'fdr',
                                                 p.value=0.05, lfc=log2(1.5), number=nrow(filtered.lipids)) %>% as.data.frame()

test.G3.G4.dif <- toptable_OUTPUT1.2 <- topTable(fit2, coef="G3_vs_G4_L", adjust.method = 'fdr',
                                                 p.value=0.05, lfc=log2(1.5), number=nrow(filtered.lipids)) %>% as.data.frame()

test.H.L.diff <- toptable_OUTPUT1.3 <- topTable(fit2, coef="L_H", adjust.method = 'fdr',
                                                p.value=0.05, lfc=log2(1.5), number=nrow(filtered.lipids)) %>% as.data.frame()


# test.diet.wt <- toptable_OUTPUT1.4 <- topTable(fit2, coef="Dietwt", adjust.method = 'fdr',
#                                         p.value=0.05, lfc=log2(fold.change), number=nrow(filtered.lipids)) %>% as.data.frame()
# 
# 
# test.ko.wt.diff <- toptable_OUTPUT1.5 <- topTable(fit2, coef="Diff", adjust.method = 'fdr',
#                                                   p.value=0.05, lfc=log2(fold.change), number=nrow(filtered.lipids)) %>% as.data.frame()


write.csv (x=toptable_OUTPUT0.1, file=sprintf("data/test.ANOVA_statistics.csv"))

write.csv (x=toptable_OUTPUT1.1, file=sprintf("data/test.G2.G1_statistics.csv"))

write.csv (x=toptable_OUTPUT1.2, file=sprintf("data/test.G3.G4_statistics.csv"))

write.csv (x=toptable_OUTPUT1.3, file=sprintf("data/test.H.L_statistics.csv"))

# write.csv (x=toptable_OUTPUT1.4, file=sprintf("data/test.Diet_wt_statistics.csv"))
# 
# write.csv (x=toptable_OUTPUT1.5, file=sprintf("data/test.Interaction_statistics.csv"))








#fold.change <- as.numeric(readline("Please input the fold change: "))

input <- G3.G4.L

input$sig <- factor(input$adj.P.Val< 0.05 & abs(input$logFC) > 1.5)
input <- input %>% group_by(sig) 
sig.lipids <- sum(input$adj.P.Val< 0.05 & abs(input$logFC) > 1.5)
message("When the tests' q value treshold is 0.05, and the fold change threshold is 1.5. The number of lipids which are statistically significant are: ", sig.lipids)



points <- input %>% 
  rownames_to_column('lipid') %>% 
  filter(abs(logFC) > 1.5 & adj.P.Val< 0.05 ) %>% 
  column_to_rownames('lipid')

volc = ggplot(input, aes(logFC, -log10(adj.P.Val))) + #volcanoplot with log2Foldchange versus pvalue
  geom_rect(aes(xmin = -1.55, xmax = 1.5, ymin = -Inf, ymax = Inf),
            fill = "grey90")+
  geom_vline(xintercept=0, linetype="solid", colour="grey60", size=0.4) +
  #Add a horizontal line for P-value cut-off
  geom_hline(yintercept=-log10(0.05), linetype="solid", colour="grey60", size=0.4)+
  geom_hline(yintercept=-log10(0.01), linetype="solid", colour="grey60", size=0.4)+
  #scale_y_continuous(breaks=seq(,by=0.5))+  # Ticks from 0-10, every .25
  geom_point(aes(col=sig, size=AveExpr)) + #add points colored by significance, point size by average value
  scale_size_continuous(range = c(0.25,2.5))+ #add size limit
  scale_color_manual(values=c("dodgerblue3", "darkorange3", 'chartreuse4', 'grey45')) + 
  #ggtitle ('Total Lipid species') +
  #xlim(c(-2.7, 2.7)) + ylim(c(0, 2)) + #set the x,y axis limits.
  xlab("Log2 fold change ") + ylab("-Log10(q value)") + 
  theme_bw()+ # remove background
  theme(line=element_blank())+
  geom_text_repel(data=points, aes(label=rownames(points)), size=2.5)#adding text for the FDR points.

volc
#Add a vertical line for fold change cut-offs
plot.name <- readline("Volcano plot: ")
ggsave(filename = plot.name, path = 'plot/', device = "pdf")



