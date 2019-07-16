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
message("Please input the extrated sample names, e.g. APValue[s17]")
extracted.P.No <- as.character(readline("extracted.P.No -----> "))

message("Please input the unextrated sample names, e.g. APValue[s18]")
unextracted.P.No <- as.character(readline("unextracted.P.No -----> "))

###### if the file has extracted and unextracted samples
### either -contains(colnumber of these samples) or -contains(name of the samples) // contains(pattern)
# eg. 
extracted.P.No <- "APValue[s17]"   #### if you want interactive way, please comment this line with "#"

unextracted.P.No <- "APValue[s18]"  #### if you want interactive way, please comment this line with "#"



# select and count the sample grades you want, filter the p value which less than 0.001
lipidCount <- lipidomics %>% 
  select(contains("Grade"), -contains("Grade.c" ), -contains("extracted"), 
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


abundance.sample <- "`MainArea[s1]`"   # for running convinence 

ggplot(data=filtered.lipidomics, aes_string(x = sprintf("log10(%s)", abundance.sample), y = "BaseRt")) +
  geom_point() +
  theme_bw() +
  facet_grid(.~Class) +
  labs("Abundance VS. Retention time", x="lipid class (log10(MainArea))", y="Retention time")
message("Please input the name you want to store for the graph. e.g. retention.pdf")
name <- readline(as.character())
ggsave(name)

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
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}




# QC PLOT 2 - Pair-wise correlation between replicates
index <- 1:length(sample.info)
sample.index <- index[index %% 2 ==0]


for(i in 1:length(sample.index)){
  range <- sample.index[i]
  plot.name <- paste("pairs.plot.", i, ".png",sep="")
  png(plot.name)
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





# Formatting the table for PCA
filtered.lipidomics.PCA <-  filtered.lipidomics %>% 
                            select(sample.list) %>% 
                            t %>% 
                            as.data.frame()

colnames(filtered.lipidomics.PCA) <- filtered.lipidomics$LipidMolec                           

log2.filtered.lipidomics.PCA <- log((filtered.lipidomics.PCA+1), 2)



############# naming the group 
group.repeats <- c()
times <- c()

for (i in ngroups:1){
  group.names[i] <- dimnames(sample.info)[[2]][i]
  times[i] <- length(sample.info[[i]])
  repeats <- c(rep(group.names[i], times[i]))
 group.repeats <- c(repeats, group.repeats)
}

filtered.lipidomics.PCA$Group <- group.repeats
log2.filtered.lipidomics.PCA$Group <- group.repeats













# Perform PCA [requires FactoMineR]
res.pca <-  PCA(filtered.lipidomics.PCA, scale.unit=TRUE, ncp=5, quali.sup=ncol(filtered.lipidomics.PCA), graph=T)

concat <-  cbind.data.frame(filtered.lipidomics.PCA[, ncol(filtered.lipidomics.PCA)], res.pca$ind$coord)

ellipse.coord = coord.ellipse(concat, bary=TRUE)

plot.PCA(res.pca, habillage=ncol(filtered.lipidomics.PCA), ellipse=ellipse.coord, cex=0.8, label="all")

plot.PCA(res.pca, habillage=ncol(filtered.lipidomics.PCA), ellipse=ellipse.coord, cex=0.8, label="quali", legend = "FALSE")



# proportion of variance explained by components
fviz_screeplot(res.pca, ncp=10, ylim=c(0, 45), xlab="Principle Components")



# Output the filtered lipidome
write.csv(filtered.lipidomics, "filtered.lipidomics.csv")


# comparison expression level total lipid classes of all experimental condition using AUC and volcano plots
# main.area <- filtered.lipidomics %>% 
#   group_by(Class, Group) %>% 
#   select(contains("MainArea[s")) %>%  
#   summarise_all(.funs=sum) %>% 
#   mutate(area.value=Reduce("+", .[2:ncol(.)])) %>% ungroup() %>% 
#   select(Class, ncol(.)) %>% view()
# 
# 
# 
#    
      
### making data for barplot graph


############ check if need to remove the standards for these part 

# standards <- c("TG(17:1/17:1/17:1)","PS(17:0/20:4)","PI(17:0/20:4)","PG(17:0/14:1)","PE(17:0/14:1)","LPS(17:1)","LPE(17:1)","LPE(17:1)","LPC(17:1)","DG(19:0/19:0)","CerG1(d18:1/12:0)","Cer(d18:1/17:0)","CL(14:1/14:1/15:1/14:1)")



# filtered.groups <- filtered.lipidomics %>% select(Class, sample.list) %>% as.data.frame()
#       
# filtered.groups <- filtered.groups %>% group_by(Class) %>% summarise_all(funs(sum))  %>% ungroup() %>% t
# 
# data <- filtered.groups[-1, ] %>% as.data.frame()
# 
# colnames(data) <- filtered.groups[1,] 
# 
# data <- data.frame(data, experiment.group=group.repeats, stringsAsFactors=FALSE) %>% 
#         rownames_to_column("rowname") %>% 
#         mutate_if(., is.factor, ~as.numeric(levels(.x))[.x]) %>% 
#         column_to_rownames("rowname")

filtered.groups <- filtered.lipidomics %>% select(Class, sample.list) 
filtered.groups <- aggregate(.~Class, data=filtered.groups, FUN=sum)

data <- filtered.groups[,-1] %>% t %>% data.frame(., group.repeats)
names(data) <- c(filtered.groups$Class, "experiment.group") 

data <- data  %>% group_by(experiment.group) %>% summarise_all(funs(sum)) 


data.long <- data %>% gather(class, value, -experiment.group)

ggplot(data.long, aes(x=class, fill=experiment.group, y=value)) + geom_bar(position="dodge", stat="identity") +
  theme_bw() +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  labs(ylab="Main Area", fill="Groups") +
  ggtitle("Total lipid classes") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle=45, hjust=1))
  
 







 
 
  

###### for each class
message("How many classes you want to visualize for barplots?")
ntimes <- readline("Please input the numbers of classes you want to visualize: ")




message("which class you want to visualize for barplot?")
pick.class <- as.character(readline("Please input the class for barplot visualization: "))
each.data <- filtered.lipidomics %>% 
              select(Class, LipidMolec, sample.list) %>% 
              filter(grepl(pick.class, Class)) %>% 
              select(-Class) %>% 
              t 

filtered.class <- filtered.lipidomics %>% 
  select(Class, LipidMolec, sample.list) %>% 
  filter(grepl(pick.class, Class))

class_data <- aggregate(.~Class+LipidMolec, data=filtered.class, FUN=sum) %>% select(-Class)
data_wide <- data.frame(t(class_data[,-1]), group.repeats)
colnames(data_wide) <- c(class_data$LipidMolec, "experiment.group")
class_long <- data_wide %>% group_by(experiment.group) %>% summarise_all(funs(sum)) 
class.data.long <- class_long %>% gather(Aceyl, AreaValue, -experiment.group)

# data.wide <- each.data[-1,] %>% as.data.frame()

# indx <- sapply(data.wide, is.factor)
# data.wide[indx] <- lapply(data.wide[indx], function(x) as.numeric(as.character(x)))
# 
# class.data <- data.frame(data.wide, experiment.group=group.repeats, stringsAsFactors=FALSE) %>%
#   rownames_to_column("rowname") %>%
#   mutate_if(., is.factor, ~as.numeric(levels(.x))[.x]) %>%
#   column_to_rownames("rowname")
# 
# 
# 
# class.data <- data.frame(data.wide, group.repeats)
#col.names <- each.data[1,] %>% str_extract_all(., "\\(\\d*.*\\)") %>% str_remove_all(., "[\\(\\)]")
#col.names <- c(col.names, "experiment.group")



# class.data <- class.data %>% group_by(group.repeats) %>% summarise_all(funs(sum))
# names(class.data) <- c("experiment.group", each.data[1,])
# duplicate.samples.indx <- which(duplicated(names(class.data))) 
# duplicate.samples <- names(class.data[duplicate.samples.indx])
# duplicate.indx.all <- lapply(duplicate.samples, function(x) which(names(class.data)==x)) %>% unlist()
# 
# 
# for(i in 1:length(duplicate.samples)){
#   indx <- which(names(class.data)==duplicate.samples[i]) 
#   names[i] <- duplicate.samples[i]
#   sample <- rowSums(class.data[,indx])
#   class.data <- cbind(class.data, sample) 
# }
# class.data <- class.data[,-duplicate.indx.all]
# 
# names(class.data)[(ncol(class.data)-length(duplicate.samples)+1):ncol(class.data)] <- c(duplicate.samples)


# class.data.long <- class.data %>% gather(Aceyl, AreaValue, -experiment.group)
# 
# class.data.long$Aceyl <- class.data.long$Aceyl %>% str_extract_all(., "\\(\\d*.*\\)") %>% str_remove_all(., "[\\(\\)]")

ggplot(class.data.long, aes(x=Aceyl, fill=experiment.group, y=AreaValue)) + geom_bar(position="dodge", stat="identity") +
  theme_bw() +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  labs(ylab="Main Area", fill="Groups", xlab="AT") +
  ggtitle(pick.class) +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle=45, hjust=1)) +
  xlab("Aceyl composition")
  





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
  refed_vs_fasted_in_WT=GPAT3_4i-DMSO,
  refed_vs_fasted_in_KO=DGAT2_GPAT3_4i-DGAT2i,
  ko_vs_wt_in_fasted=DGAT2i-DMSO,
  ko_vs_wt_in_refed=DGAT2_GPAT3_4i-GPAT3_4i,
  Diff=(DGAT2_GPAT3_4i-DGAT2i)-(GPAT3_4i-DMSO),levels=design)
##################################################################################################

fit2 <- contrasts.fit(fit, cont.matrix)

fit2 <- eBayes(fit2)

# Use topTable to obtain significantly altered lipids
Diff <- topTable(fit2, coef="Diff", adjust.method = "fdr",
                 lfc=0, number=nrow(filtered.lipids), sort.by="logFC")

refed_vs_fasted_in_WT <- topTable(fit2, coef="refed_vs_fasted_in_WT", adjust.method = "fdr",
                                  lfc=0, number=nrow(filtered.lipids), sort.by="logFC")

refed_vs_fasted_in_KO <- topTable(fit2, coef="refed_vs_fasted_in_KO", adjust.method = "fdr",
                                  lfc=0, number=nrow(filtered.lipids), sort.by="logFC")

ko_vs_wt_in_fasted <- topTable(fit2, coef="ko_vs_wt_in_fasted", adjust.method = "fdr",
                               lfc=0, number=nrow(filtered.lipids), sort.by="logFC")

ko_vs_wt_in_refed <- topTable(fit2, coef="ko_vs_wt_in_refed", adjust.method = "fdr",
                              lfc=0, number=nrow(filtered.lipids), sort.by="logFC")
# Write to .csv
write.csv(Diff, "data/Diff.csv")

write.csv(refed_vs_fasted_in_WT, "data/refed_vs_fasted_in_WT.csv")

write.csv(refed_vs_fasted_in_KO, "data/refed_vs_fasted_in_KO.csv")

write.csv(ko_vs_wt_in_fasted, "data/ko_vs_wt_in_fasted.csv")

write.csv(ko_vs_wt_in_refed, "data/ko_vs_wt_in_refed.csv")


#message("Please input the contrast groups from the toptable")
#input <- readline("Information retrieved from the toptable: ")
fold.change <- readline("Please input the times of Fold change: ")

input <- Diff

# setting treshold for the difference of experiment groups
#input$sig <- ifelse(-log10(input$adj.P.Val< 0.05 & abs(input$logFC) > fold.change), "significant", "nonsignificant")

input$sig <- ifelse(input$adj.P.Val< 0.05 & abs(input$logFC) > fold.change, "significant", "nonsignificant")
input <- input %>% group_by(sig)

points <- subset(input,abs(Diff$logFC) > 2 & -log10(Diff$adj.P.Val< 0.05) )
  volc = ggplot(input, aes(x=logFC, y=adj.P.Val)) + #volcanoplot with log2Foldchange versus pvalue
    geom_rect(xmin = -2, xmax = fold.change, ymin = -Inf, ymax = Inf,
            fill = "grey60")+
    geom_vline(xintercept=fold.change, linetype="solid", colour="grey60", size=0.4) +
    geom_vline(xintercept=-2, linetype="solid", colour="grey60", size=0.4) +
    #Add a horizontal line for P-value cut-off
    geom_hline(yintercept=0.05, linetype="solid", colour="grey60", size=0.4)+
    #scale_y_continuous(breaks=seq(,by=0.5))+  # Ticks from 0-10, every .25
    geom_point(aes(col=sig, size=AveExpr)) + #add points colored by significance, point size by average value
    #scale_size_continuous(range = c(0.25,2.5))+ #add size limit
    #scale_color_manual(values=c("dodgerblue3", "darkorange3", 'chartreuse4', 'grey45')) + 
    #ggtitle ('Total Lipid species') +
    #xlim(c(-2.7, 2.7)) + ylim(c(0, 2)) + #set the x,y axis limits.
    xlab("Log2 fold change ") + ylab("-Log10(q value)") + 
    theme_bw()+ # remove background
    theme(line=element_blank()) +
    geom_text_repel(data=subset(input, -log10(adj.P.Val)<0.05 & abs(logFC) > fold.change), size=2.5)+aes(label=rownames(input))#adding text for the FDR points.
  #Add a vertical line for fold change cut-offs
  
  





