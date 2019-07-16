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




## @knitr variablesXY
# Read in file, and store the data 
lipidomics <- read_csv(file.choose())

### REMOVE ODD NUMBER IN THE MOLEC
#standards <- c("LPA(17:1)", "PA(37:4)", "LPS(17:1)", "PS(17:0/20:4)")

message("If the file has extracted and unextracted samples, please input the Grades and APValue as guides.
        And if only one control sample need to delete, just type last sample one more time or 
        something unexist for the second option. e.g. Grade[s100], APValue[s99]")

message("delete the extracted and extracted sample grades, eg. Grade[s22]")
extracted.Grade <- as.character(readline("extracted.Grade -----> "))
unextracted.Grade <- as.character(readline("unextracted.Grade -----> "))


############################################### 
# interactive way
message("Please input the extrated and unextracted sample names, e.g. APValue[s17]")
extracted.P.No <- as.character(readline("extracted.P.No -----> "))
unextracted.P.No <- as.character(readline("unextracted.P.No -----> "))


### either -contains(colnumber of these samples) or -contains(name of the samples) // contains(pattern)
# eg. 
# extracted.P.No <- "APValue[s17]"   #### if you want interactive way, please comment this line with "#"
# 
# unextracted.P.No <- "APValue[s18]"  #### if you want interactive way, please comment this line with "#"

#standards <- c("TG(17:1/17:1/17:1)","PS(17:0/20:4)","PI(17:0/20:4)","PG(17:0/14:1)","PE(17:0/14:1)","LPS(17:1)","LPE(17:1)","LPE(17:1)","LPC(17:1)","DG(19:0/19:0)","CerG1(d18:1/12:0)","Cer(d18:1/17:0)","CL(14:1/14:1/15:1/14:1)")


# select and count the sample grades you want, filter the p value which less than 0.001
lipidCount <- lipidomics %>% 
  select(contains("Grade"), -contains("Grade[c]" ), -contains(extracted.Grade),-contains(unextracted.Grade), 
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

# check how many lipids passed filtering
describe(filtered.lipidomics$Class)

##########################################################################
# QC PLOT 1 - Abundance vs. retention time [requires ggplot2]
##############################################################################################
# Interactive way
##############################################################
message("Please input the MainArea of the sample name to check its abundance vs. retention time, 
        e.g. `MainArea[s1]`. 
        Please note that: ` ` is the backtick which is same position as tilde")

abundance.sample <- readline("abundance.sample -----> ")

ggplot(data=filtered.lipidomics, aes_string(x = sprintf("log10(%s)", abundance.sample), y = "BaseRt")) +
  geom_point() +
  theme_bw() +
  facet_grid(.~Class) +
  labs("Abundance VS. Retention time", x="lipid class (log10(MainArea))", y="Retention time")
##############################################################################################

################
# hard code way for ggplot
#################################################################################################
# ggplot(data=filtered.lipidomics, aes(x =log10(`MainArea[s1]`), y = BaseRt)) +
#   geom_point() +
#   theme_bw() +
#   facet_grid(.~Class) +
#   labs("Abundance VS. Retention time", x="lipid class (log10(MainArea))", y="Retention time")
##################################################################################################

message("Please input the name you want to store for the graph. e.g. retention.pdf")
plot.name <- readline("QC plot name: ")
ggsave(filename = plot.name, path = 'plot/', device = "pdf")

############################################################################


# making group
message("Please input the info of the experiment groups below and end with 'Enter'
        Please Only input the sample groups for analysis.
        e.g. Group number: 2   
        Name of Group 1: wild type
        Sample names of Group 1: s1 s2 s3 s4
        Name of Group 2: knockout
        Sample names of Group 2: s4 s5 s6 s7")

# info about group size
ngroups <- readline("Group number: ")

# store group info into list
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

# Type the group info and store it into a variable
total.info <- inputGroups(ngroups)
message("Please take a look at the group information below")
glimpse(total.info)

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
message("Take a look at the sample info and its position information below")
glimpse(sample.info)

# retrieve the group names
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

##########################################################################
# QC PLOT 2 - Pair-wise correlation between replicates

# retrieve the sample info position 
index <- 1:length(sample.info)
sample.index <- index[index %% 2 ==0]

# Storing pairwise plot in the plot directory
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

# making group repeats according to its position for making groups of later PCA
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

# making group info for PCA 
filtered.lipidomics.PCA$Group <- group.repeats
log2.filtered.lipidomics.PCA$Group <- group.repeats

# Perform PCA 
res.pca <-  PCA(filtered.lipidomics.PCA, scale.unit=TRUE, ncp=5, quali.sup=ncol(filtered.lipidomics.PCA), graph=T)
concat <-  cbind.data.frame(filtered.lipidomics.PCA[, ncol(filtered.lipidomics.PCA)], res.pca$ind$coord)
ellipse.coord = coord.ellipse(concat, bary=TRUE)

pdf(file="plot/sample.pca.pdf")
plot.PCA(res.pca, habillage=ncol(filtered.lipidomics.PCA), ellipse=ellipse.coord, cex=0.8, label="all")
dev.off()

pdf(file="plot/group.pca.pdf")
plot.PCA(res.pca, habillage=ncol(filtered.lipidomics.PCA), ellipse=ellipse.coord, cex=0.8, label="quali")
dev.off()

# Output the filtered lipidome
write.csv(filtered.lipidomics, "data/filtered.lipidomics.csv")
###################################################################################### 


#########################################################################
######## making barplot graphs
# For total classes
# select data contains sample info of the groups
filtered.groups <- filtered.lipidomics %>% select(Class, sample.list) 
# reform the data by aggregating data contains same class name
filtered.groups <- aggregate(.~Class, data=filtered.groups, FUN=sum)
# build new data frame contains the group info
total.class.data <- filtered.groups[,-1] %>% t %>% data.frame(group.repeats, .)
names(total.class.data) <- c("experiment.group",filtered.groups$Class) 
# calculate the standard deviation for each class in each group
total.class.sd <- aggregate(.~experiment.group, data=total.class.data, function(x) sd(x)) 
# reformat the data structure contains sd information into long data
total.class.sd <- total.class.sd %>% gather( class, sd, -experiment.group)

# aggregating data by group
total.class.wide <- total.class.data  %>% group_by(experiment.group) %>% summarise_all(.funs=sum) 
# reformat data into long data
total.class.long <- total.class.wide %>% gather(class, value, -experiment.group)
# combine data and its sd info for ggplot
total.class.long <- bind_cols(total.class.long, sd=total.class.sd$sd)

# plot the total class info by groups
# the y axis range can be editted by changing the parameter in "coord_cartesian(ylim=c(1e5,1e13))"
ggplot(total.class.long, aes(x=class, fill=experiment.group, y=value)) + 
  geom_bar(position="dodge", stat="identity") +
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

# input the plot name and store it
message("Please input the name you want to store for the graph. e.g. retention.pdf")
name.class  <- readline("Total lipid classes plot name: ")
ggsave(filename = name.class, path = 'plot/', device = "pdf")

############################
# For each class
message("How many classes you want to visualize for barplots?")
ntimes <- readline("Please input the numbers of classes you want to visualize: ")

# define the function for ploting 
eachClassPlots <- function(){
  message("which class you want to visualize for barplot?")
  pick.class <- as.character(readline("Please input the class for barplot visualization: "))
  # get information of the samples
  filtered.class <- filtered.lipidomics %>% 
    select(Class, LipidMolec, sample.list) %>% 
    filter(grepl(pick.class, Class))
  # aggregate the repeated data by its class and lipid molecular names
  class.data <- aggregate(.~Class+LipidMolec, data=filtered.class, FUN=sum) %>% select(-Class)
  # Making the group information for the data
  data.wide <- data.frame(t(class.data[,-1]), group.repeats)
  # retreive the numeric information of lipid molecular names for better visualization
  molec.names <- class.data[,1] %>% 
    unlist() %>% 
    str_extract_all(., "\\(\\d*.*\\)") %>% s
  tr_remove_all(., "[\\(\\)]") 
  # making names for the new data frame                
  names(data.wide) <- c(molec.names, "experiment.group")
  # calculate the sd information   
  class.sd <- aggregate(.~experiment.group, data=data.wide, function(x) sd(x)) 
  class.sd <- class.sd %>% gather( Aceyl, sd, -experiment.group)
  # reformat the data information
  class.long <- data.wide %>% group_by(experiment.group) %>% summarise_all(funs(sum)) 
  class.data.long <- class.long %>% gather(Aceyl, AreaValue, -experiment.group)
  # combine the standard deviation information with data information
  class.data.long <- bind_cols(class.data.long, sd=class.sd$sd)
  
  # plot for the individual class by groups
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
    coord_cartesian(ylim=c(1e5,1e10))
  
  message("Please input the name you want to store for the graph. e.g. TG.pdf")
  class.type <- readline("The lipid class plot name: ")
  ggsave(filename = class.type, path = 'plot/', device = "pdf")
}
# plot for each class
replicate(ntimes, eachClassPlots())
#################################################################################

source("a1.Rmd")


