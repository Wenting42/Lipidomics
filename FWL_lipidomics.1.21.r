####################################################################################
# Script: FWL_lipidomics.1.1.R
# Author: Niklas, Kenny, Wenting 
# Notes: This script helps generating the graph and data for the flowork of Lipidomics
#         1) interactive way. 
#         First please store the file name in the variable file.name for check 
#         Run the script by command "source("FWL_lipidomics.1.1.R") in console.
#         All the data are stored in the directory called "data", and plots in "plot"
#         2) You can change the parameters in the script and run it lines by lines.
#           Before you change any parameters, please copy the original script first.
#####################################################################################


rm(list=ls())  ##### empty the environment first. If you need to keep former variables in workspace, just comment this part

##### Please install the lacking packages ------> install.packages("package")
list.of.packages <- c("FactoMineR", "tidyverse", "ggplot2", "magrittr", "ggrepel", "reshape2",
                      "stargazer", "Hmisc", "limma", "factoextra", "scales", "RColorBrewer",
                      "stringr", "readxl")
need.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if(length(need.packages)) BiocManager::install(need.packages)

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

message('Please store the file name in the file.name for check. And you could go to script to add the file name in read_csv(file.choose("data_raw/")). 
        If you need to change some parameters or edit the script, please copy the original script first. ')
file.name <- readline("Please input the file name: ")

# Read in file, and store the data 
lipidomics <- read_csv(file.choose("data_raw/"))

### REMOVE ODD NUMBER IN THE MOLEC
#standards <- c("LPA(17:1)", "PA(37:4)", "LPS(17:1)", "PS(17:0/20:4)")

message("If your experiment includes standard control and solvent control, please input the information of the samples. And if only one control sample need to delete, just type this sample two times or plus something unexist for the second option. 
        e.g. Grade[s100], APValue[s99]")

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
# extracted.P.No <- "APValue[s17]"   #### if you want hard code way, please uncomment this part
# 
# unextracted.P.No <- "APValue[s18]"  #### if you want interactive way, please comment this line with "#"


# select and count the sample grades you want, filter the p value which less than 0.001
lipidCount <- lipidomics %>% 
  select(contains("Grade"), -contains("Grade[c]" ),
         -contains(extracted.Grade),-contains(unextracted.Grade), 
         contains("APValue"),-contains(extracted.P.No), 
         -contains(unextracted.P.No))  %>% 
  transmute(A = rowSums(. == "A"), B = rowSums(. == "B"), 
            C = rowSums(. == "C"), D = rowSums(. == "D"), 
            No.grade = rowSums(.=="-"), APvalue.001 = rowSums(. <= "0.001"))

# add this count into original strucutes
lipidSelect <- lipidomics %>% bind_cols(lipidCount)

# Filter the dataset based on your criteria  
# flexible parameter 4 which depends on the experiment for total number of A and B 
filtered.lipidomics <- lipidSelect %>% 
  rowwise() %>% 
  filter( Rej == 0 & sum(A, B) >= 4 & APvalue.001 >= 4)

### remove artifical lipids
# mark the artifical classes position
odd.index <- filtered.lipidomics$LipidMolec %>%
  str_locate(., "(\\d[13579]:)|([13579]:)")
percent.odd <- length(unique(which(!is.na(odd.index), arr.ind=TRUE)[,1]))/nrow(filtered.lipidomics)

a <- filtered.lipidomics %>% filter(Class=="TG") 
b <- a$LipidMolec %>%   str_locate(., "(17:)|(7:)")
b <- a$LipidMolec %>%   str_locate(., "(\\d[13579]:)|([13579]:)")
odd.chain <- length(unique(which(!is.na(b), arr.ind=TRUE)[,1]))/nrow(a)

# odd in total is 11%
# TG in total is 21.5%
nrow(a)/nrow(filtered.lipidomics)
# odd in TG is 13.8% while 17: odd chain is 7.3%

print(percent.odd)
print(odd.chain)

# #############    remove the artifical classes
# filtered.lipidomics <- filtered.lipidomics %>%
#                         slice(unique(which(is.na(odd.index), arr.ind=TRUE)[,1]))

# check how many lipids passed filtering
print(describe(filtered.lipidomics$Class))
print(summary(as.data.frame(filtered.lipidomics$Class), maxsum=nrow(filtered.lipidomics)))



##########################################################################
# QC PLOT 1 - Abundance vs. retention time [requires ggplot2]
##############################################################################################
# Interactive way
##############################################################
message("Please input the MainArea of the sample name to check its abundance vs. retention time, 
        e.g. `MainArea[s1]`. 
        Please note that: ` ` is the backtick which is same position as tilde")

abundance.sample <- readline("abundance.sample -----> ")

retention.plot <- ggplot(data=filtered.lipidomics, 
                         aes_string(x = sprintf("log10(%s)", abundance.sample), y = "BaseRt")) +
  geom_point() +
  theme_bw() +
  facet_grid(.~Class) +
  labs("Abundance VS. Retention time", 
       x="lipid class (log10(MainArea))", y="Retention time")
print(retention.plot)
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
ggsave(filename = plot.name, path = 'plot/', device = "pdf", width=15, height=15, dpi=300)

############################################################################


# making group
message("Please input the info of the experimental groups below and end with 'Enter'
        Please Only input the sample names (replicate) for analysis.
        e.g. Experimental Group number: 2   
        Name of Experimental Group 1: wild_type
        Sample names (replicate number from lipid search) : s1 s2 s3 s4 
        Name of Group 2: knockout
        Sample names of Group 2: s4 s5 s6")

# info about group size
ngroups <- readline("Group number: ")

# store group info into list
inputGroups   <- function(n){
  group.info  <- c()
  sample.info <- c()
  info <- list()
  for(i in 1:n){
    ms1                   <- paste("Name of Group ", i, ": ")
    group.info[i]         <- readline(prompt= ms1)
    ms2                   <- paste("Sample names of Group ", i, ": ")
    sample.info[i]        <- readline(prompt=ms2)
    info[[group.info[i]]] <- sample.info[i]
  }
  return(info)
}

# Type the group info and store it into a variable
total.info <- inputGroups(ngroups)
message("Please take a look at the group information below")
glimpse(total.info)

message("Do you need retype the group info? ")
group.condition <- readline("Y/N: ")
if(group.condition == "N"){
  message("Continue")
}else{
  total.info <- inputGroups(ngroups)
  glimpse(total.info)
}


# Grep index and sample names 
grepIndex <- function(sample, data){
  sample.names  <- list()
  group.names   <- c()
  col.index     <- list()
  for( i in 1:length(names(sample))){
    group.names[i]                  <- names(sample)[i]
    group.sample.names              <- sample[[i]] %>% 
      strsplit(., "\\s+") %>% 
      unlist() %>% 
      str_to_lower()
    sample.names[[group.names[i]]]  <- paste("MainArea[", group.sample.names, "]", sep="")
    col.index[[group.names[i]]]     <- which(colnames(data) %in% sample.names[[i]])
  }
  info      <- rbind(sample.names, col.index)
  return(info)
}

# store the sample and index in a list
sample.info <- grepIndex(total.info, filtered.lipidomics) 
message("Take a look at the sample info and its position information below")
glimpse(sample.info)

# retrieve the group names
group.names <- dimnames(sample.info)[[2]]

# Preparation for pair-wise correlations
inf2NA      <- function(x) { x[is.infinite(x)] <- NA; x }

panel.cor   <- function(x, y, digits = 2, cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  # correlation coefficient
  r   <- cor(inf2NA(x), inf2NA(y), use = "pairwise.complete.obs", method = "pearson")
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste("r= ", txt, sep = "")
  text(0.5, 0.6, txt)
  
  # p-value calculation
  p                       <- cor.test(x, y)$p.value
  txt2                    <- format(c(p, 0.123456789), digits = digits)[1]
  txt2                    <- paste("p= ", txt2, sep = "")
  if(is.na(p<0.01)) txt2  <- paste("p= ", "<0.01", sep = "")
  text(0.5, 0.4, txt2)
}

panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h   <- hist(x, plot = FALSE)
  breaks  <- h$breaks; nB <- length(breaks)
  y       <- h$counts; y  <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "#4393C3", ...)
}

##########################################################################
# QC PLOT 2 - Pair-wise correlation between replicates

# retrieve the sample info position 
index         <- 1:length(sample.info)
sample.index  <- index[index %% 2 ==0]

# Storing pairwise plot in the plot directory
for(i in 1:length(sample.index)){
  range     <- sample.index[i]
  plot.name <- paste("pairs.plot.", i, ".pdf",sep="")
  path      <- file.path("plot/", plot.name)
  pdf(file=path)
  pairs(log10(filtered.lipidomics[, sample.info[[range]]]), 
        lower.panel = panel.smooth, diag.panel=panel.hist, 
        upper.panel = panel.cor, xlim=c(4,12), ylim=c(4,12))
  #dev.off()
}

# making group repeats according to its position for making groups of later PCA
index.list  <- index[!index %% 2 ==0]
sample.list <- c()
for(i in length(index.list):1){
  j           <- index.list[i]
  sample.list <- c(sample.info[[j]], sample.list)
}

group.repeats <- c()
for(i in length(group.names):1){
  k             <- index.list[i]
  len           <- length(sample.info[[k]])
  repeats       <- rep(group.names[i], len)
  group.repeats <- c(repeats, group.repeats)
}


# Formatting the table for PCA
filtered.lipidomics.PCA <-  filtered.lipidomics %>% 
  select(sample.list) %>% 
  t() %>% 
  as.data.frame()

colnames(filtered.lipidomics.PCA) <- filtered.lipidomics$LipidMolec                           
log2.filtered.lipidomics.PCA      <- log((filtered.lipidomics.PCA+1), 2)

# making group info for PCA 
filtered.lipidomics.PCA$Group       <- group.repeats
log2.filtered.lipidomics.PCA$Group  <- group.repeats

# Perform PCA 
res.pca         <-  PCA(filtered.lipidomics.PCA, scale.unit=TRUE, ncp=5, 
                        quali.sup=ncol(filtered.lipidomics.PCA), graph=T)
concat          <-  cbind.data.frame(filtered.lipidomics.PCA[, ncol(filtered.lipidomics.PCA)], 
                                     res.pca$ind$coord)
ellipse.coord   <-  coord.ellipse(concat, bary=TRUE)

#pdf(file="plot/sample.pca.pdf")
plot.PCA(res.pca, habillage=ncol(filtered.lipidomics.PCA), 
         ellipse=ellipse.coord, cex=0.8, label="all")
dev.copy(pdf, "plot/sample.pca.pdf")
dev.off()


plot.PCA(res.pca, habillage=ncol(filtered.lipidomics.PCA), 
         ellipse=ellipse.coord, cex=0.8, label="quali")
dev.copy(pdf, "plot/group.pca.pdf")
dev.off()

# Output the filtered lipidome
write.csv(filtered.lipidomics, "data/filtered.lipidomics.csv")
###################################################################################### 


#########################################################################
######## making barplot graphs
# For total classes
totalClass <- function(filtered.lipidomics, sample.list, group.repeats){
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
  total.plot <- ggplot(total.class.long, aes(x=reorder(class, value), fill=experiment.group, y=value)) + 
    geom_bar(position="dodge", stat="identity") +
    geom_errorbar(aes(ymin=value, ymax=value+sd), width=.3,
                  position=position_dodge(.9)) +
    theme_bw() +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)), limits=c(0, 1e20)) +
    # scale_y_continuous(labels = scales::scientific) +
    #scale_y_continuous(trans="log10", limits=c( 1, NA)) +
    labs(ylab="Main Area", fill="Groups") +
    ggtitle("Total lipid classes") +
    theme(plot.title = element_text(hjust = 0.5), 
          axis.text.x = element_text(angle=45, hjust=1))+
    coord_cartesian(ylim=c(1e3,1e18)) ## these y axis range parameter can be changed
  print(total.plot)
  ggsave(filename="total.class.1.pdf",path = 'plot/', device = "pdf")
  
  message("Normalization for each class in each group.")
  control.group <- readline("Please input a group as control group for normalization, e.g. A : ")
  control.index <- which(total.class.wide$experiment.group %in% control.group)
  
  # normalization 
  # total.data <- log(total.class.data[,-1], 2)
  # total.d <- bind_cols(experiment.group=total.class.data$experiment.group, total.data)
  # 
  # this is just boxplot which is bizzare  
  # data.long <- total.class.data %>% gather(class, value, -experiment.group)
  # ggplot(data.long, aes(x=class, y=value, fill=experiment.group)) +
  #   geom_boxplot()
  #   
  
  
  # normalization of sum of main area
  normalization.factor <- apply(total.class.wide[-control.index,-1], 1, 
                                function(x)log2(x/total.class.wide[control.index,-1]+1)) %>% 
    unlist() %>%
    matrix(., nrow=nrow(total.class.wide)-1, byrow=TRUE)
  
  normal.log2.matrix <- total.class.wide
  normal.log2.matrix[-control.index, -1] <- normalization.factor 
  normal.log2.matrix[control.index, -1] <- rep(1, ncol(total.class.wide)-1)
  
  normal.log2.data <- normal.log2.matrix %>% gather(class, value, -experiment.group)
  
  
  total.normalization.plot <- ggplot(normal.log2.data, aes(x=reorder(class, value), fill=experiment.group, y=value)) +
    geom_col(position="dodge") +
    theme_bw() +
    coord_flip() +
    geom_hline(yintercept=1, linetype="dashed", colour="grey45", size=0.02) +
    labs(title="Total lipid classes",x="Lipid class", y="Log2 fold change") +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle=45, hjust=1))
  print(total.normalization.plot)
  ggsave(filename="total.class.foldchange.pdf",path = 'plot/', device = "pdf")
}

totalClass(filtered.lipidomics, sample.list, group.repeats)




# remove classes unwanted in the total class plot
message("Please check the plot total.class.1.pdf in the plot directory or look the Plots pannel in Rstudio first. Do you have any more classes to remove?")
judge <- readline("Y/N: ")
if(judge == "Y"){
  class.remove      <- readline("Please input the classes you want to remove, e.g(CerG3 ChE): ")
  class.remove      <- class.remove %>% str_split(., "\\s+") %>% unlist()
  filtered.lipidomics  <- filtered.lipidomics %>% filter(!Class %in% class.remove)
}






# # input the plot name and store it
# message("Please input the name you want to store for the graph. e.g. total.class.pdf")
# name.class  <- readline("Total lipid classes plot name: ")
# ggsave(filename = name.class, path = 'plot/', device = "pdf")

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
    filter(Class==pick.class)
  # aggregate the repeated data by its class and lipid molecular names
  class.data <- aggregate(.~Class+LipidMolec, data=filtered.class, FUN=sum) %>%  select(-Class)
  # rename the class
  molec.names <- class.data[,1] %>% 
    unlist() %>% 
    str_extract_all(., "\\(\\d*.*\\)") %>% 
    str_remove_all(., "[\\(\\)]") 
  
  # Making the group information for the data
  data.wide <- data.frame(group.repeats, t(class.data[, -1]))
  # making names for the new data frame                
  names(data.wide) <- c("experiment.group", molec.names)
  
  # calculate the sd information   
  class.sd <- aggregate(.~experiment.group, data=data.wide, function(x) sd(x)) %>%      
    gather( Acyl, sd, -experiment.group)
  
  # reformat the data information
  class.long <- aggregate(.~experiment.group, data=data.wide, FUN=sum) %>% 
    gather(Acyl, AreaValue, -experiment.group) %>% 
    merge.data.frame(., class.sd)
  
  # plot for the individual class by groups
  class.plot <- ggplot(class.long, aes(x=reorder(Acyl, AreaValue), fill=experiment.group, y=AreaValue)) + 
    geom_bar(position="dodge", stat="identity") +
    geom_errorbar(aes(ymin=AreaValue, ymax=AreaValue+sd), 
                  width=.8, position=position_dodge(.9)) +
    theme_bw() +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), # base of 10
                  labels = trans_format("log10", math_format(10^.x)), limits=c(0, 1e20)) +
    #scale_y_continuous(labels = scales::scientific)+    # scientific notation
    #scale_y_continuous(trans="log10", limits=c(NA, 1))
    #labs(ylab="Main Area", fill="Groups", xlab="AT") +
    ggtitle(pick.class) +
    theme(plot.title = element_text(hjust = 0.5), 
          axis.text.x = element_text(angle=45, hjust=1)) +
    xlab("Acyl composition") +
    ylab("Main Area") +
    coord_cartesian(ylim=c(1e3,1e18)) 
  
  print(class.plot)
  
  message("Please input the name you want to store for the graph. e.g. TG.pdf")
  class.type <- readline("The lipid class plot name: ")
  ggsave(filename = class.type, path = 'plot/', device = "pdf", width=15, height=15, dpi=300)
}
# plot for each class
replicate(ntimes, eachClassPlots())
###################################################################################


# volcano plot

# Create a design matrix 
samples <- factor(group.repeats)
design <- model.matrix(~0+samples)
colnames(design) <- levels(samples)

# Reformat the filtered_lipidomics_PCA, consolidate duplicated lipids and log2+1 transform 
filtered.lipids <- filtered.lipidomics %>% select(Class, LipidMolec, sample.list)
# aggregate the same classes
filtered.lipids <- aggregate(. ~LipidMolec+Class, data=filtered.lipids, FUN=sum)
rownames(filtered.lipids) <- filtered.lipids$LipidMolec

# filtered data
filtered.lipids <- filtered.lipids %>% select(-LipidMolec, -Class)
log2.filtered.lipids <- log((filtered.lipids+1), 2)

# write the filtered data information into csv
filtered.lipids %>% write_csv('data/filtered.lipids.csv')

# Fit model and extract contrasts 
fit <- lmFit(log2.filtered.lipids, design)

###################################################################################################
message("Please note that you NEED make contrast groups manually to Compare the difference between/among groups.
        e.g. compare B+C against A: A-(B+C)/2; A against B:  A-B; A-B against C-D: (A-B)-(C-D), etc.")
message("How many groups you designed for making comparison?")
n.comparisons <- readline("Number of comparison groups-----> ")

replicate( n.comparisons, {
  cont.matrix <- makeContrasts(
    readline("Please input the comparison, e.g A-B: "),
    levels = design
  )
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  comparison <- toptable_OUTPUT1 <- topTable(fit2, coef=1, adjust.method = 'fdr',
                                             lfc=0, number=nrow(filtered.lipids)) %>% 
    as.data.frame()
  # store the result into csv
  name1 <- readline("Please input the name for comparison data file, e.g. A.B.csv: ")
  write_csv(comparison, file.path("data/", name1))
  # store significant result into csv
  toptable_OUTPUT1.1 <- topTable(fit2, coef=1, adjust.method = 'fdr',
                                 p.value=0.05, lfc=log2(1.5), number=nrow(filtered.lipids)) %>% 
    as.data.frame()
  name1.1 <- readline("Please input the name for significant comparison data file, e.g. A.B.sig.csv: ")
  
  write_csv(x=toptable_OUTPUT1.1, file.path("data", name1.1))
  
  # volcano plot
  fold.change <- as.numeric(readline("Please input the fold change treshold for the volcano plot, e.g 1: "))
  
  # volcano input info. 
  input <- comparison  # input of comparison group info which can be changed
  # significant data
  input$sig <- factor(input$adj.P.Val< 0.05 & abs(input$logFC) > fold.change) 
  #  input <- input %>% group_by(sig) 
  
  # size of significant data
  sig.lipids <- sum(input$adj.P.Val< 0.05 & abs(input$logFC) > fold.change)
  message("When the tests' q value treshold is 0.05 and the fold change threshold is ",fold.change,". The number of lipids which are statistically significant are: ", sig.lipids)
  
  # Define significant data for volcano plot graphing
  points <- input %>% 
    rownames_to_column('lipid') %>% 
    filter(abs(logFC) > fold.change & adj.P.Val< 0.05 ) %>% 
    column_to_rownames('lipid')
  
  # making volcano plot        
  volc.plot<- ggplot(input, aes(logFC, -log10(adj.P.Val))) + #volcanoplot with log2Foldchange versus pvalue
    geom_rect(aes(xmin = -fold.change, xmax = fold.change, ymin = -Inf, ymax = Inf),
              fill = "grey90", linetype="dashed", color="red", size=0.2)+
    geom_point() +
    geom_vline(xintercept=0, linetype="solid", colour="grey60", size=0.2) +
    #Add a horizontal line for P-value cut-off
    geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="red", size=0.2)+
    #geom_hline(yintercept=-log10(0.01), linetype="dashed", colour="red", size=0.4)+
    scale_y_continuous()+  # Ticks from 0-10, every .25
    geom_point(aes(col=sig, size=AveExpr)) + #add points colored by significance, point size by average value
    scale_size_continuous(range = c(0.25,2.5), breaks=c(10, 20, 30))+ #add size limit
    scale_color_manual(values=c("#bdbdbd", "#3182bd", 
                                'chartreuse4', 'grey45',"darkorange3"))+
    #ggtitle ('Total Lipid species') +
    #xlim(c(-2.7, 2.7)) + ylim(c(0, 2)) + #set the x,y axis limits.
    xlab("Log2 fold change ") + ylab("-Log10(q value)") + 
    theme_bw()+ # remove background
    theme(line=element_blank())+
    #scale_y_continuous(expand = c(0, 0), limits=c(0,15)) +       ###### start at 0
    #geom_text_repel(data=points, aes(label=rownames(points)), size=2)#adding text for the FDR points.
    geom_text_repel(data=points, aes(label=rownames(points)), size=2)
  
  
  print(volc.plot)
  plot.name <- readline("Please input the volcano plot name: ")
  ggsave(filename = plot.name, path = 'plot/', device = "pdf", width=15, height=15, dpi=300)
  #dev.off()
  
  
  
  # build mutated data frame
  class.names <- rownames(input) %>% str_extract_all(., "(.+)\\(") %>% str_remove_all(., "\\(")
  
  #Add a new column to specify lipids
  input.lip <- input
  input.lip$lipclass <- class.names
  input.lip$lipclass <- ifelse(input.lip$sig=="TRUE", 
                               case_when( (input.lip$lipclass %in% c("PC","PE","PG","PI","PS","LPC", "LPE",
                                                                     "LPI","dMePE","CL"))~ "Glycerophospholipids",
                                          input.lip$lipclass %in% c("TG", "DG") ~ "Neutral lipids",
                                          input.lip$lipclass %in% c("SM","So","Cer","CerG1","CerG2","GD1a",
                                                                    "GD3","GM1","GM2","GM3","GT3") ~ "Sphingolipids",
                                          input.lip$lipclass %in% c("ChE","Cholestoral") ~ "Sterols",
                                          TRUE ~ "Other lipids"), 
                               "n.s")
  # making levels for the class category
  class.levels <- c("Glycerophospholipids", "Neutral lipids", "Sphingolipids", "Sterols", "Other lipids", "n.s")
  # make the class category as factors
  input.lip$lipclass <- factor(input.lip$lipclass, levels=class.levels)
  # rearrange the data by its levels 
  input.lip <-  input.lip %>% rownames_to_column('lipid') %>% arrange(lipclass) %>% column_to_rownames('lipid')
  # subset the significant points for text information
  sig.points <- input.lip %>% 
    rownames_to_column('lipid') %>% 
    filter(!lipclass=="n.s") %>% 
    column_to_rownames('lipid')
  
  # Recolor volcano plot 
  volc.plot.color <- ggplot(input.lip, aes(logFC, -log10(adj.P.Val)), fill=lipclass) + #volcanoplot with log2Foldchange versus pvalue
    geom_rect(aes(xmin = -fold.change, xmax = fold.change, ymin = -Inf, ymax = Inf),
              fill = "grey90", linetype="dashed", color="red", size=0.2)+
    geom_point() +
    geom_vline(xintercept=0, linetype="solid", colour="grey60", size=0.2) +
    #Add a horizontal line for P-value cut-off
    geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="red", size=0.2)+
    #geom_hline(yintercept=-log10(0.01), linetype="dashed", colour="red", size=0.4)+
    scale_y_continuous()+  # Ticks from 0-10, every .25
    geom_point(aes(col=lipclass, size=AveExpr)) + #add points colored by significance, point size by average value
    scale_size_continuous(range = c(0.25,2.5), breaks=c(10, 20, 30))+ #add size limit
    scale_color_manual(values=c("darkorange3","dodgerblue3",  'chartreuse4', '#756bb1','#c994c7', '#636363'), drop=FALSE)+
    #ggtitle ('Total Lipid species') +
    #xlim(c(-2.7, 2.7)) + ylim(c(0, 2)) + #set the x,y axis limits.
    xlab("Log2 fold change ") + ylab("-Log10 (q value)") + 
    theme_bw()+ # remove background
    theme(line=element_blank())+
    #scale_y_continuous(expand = c(0, 0), limits=c(0,15)) +       ###### start at 0
    #geom_text_repel(data=points, aes(label=rownames(points)), size=2)#adding text for the FDR points.
    geom_text_repel(data=points, aes(label=rownames(points)), size=2)
  
  print(volc.plot.color)
  plot.name <- readline("Please input the volcano plot name: ")
  ggsave(filename = plot.name, path = 'plot/', device = "pdf", width=15, height=15, dpi=300)
  
})






