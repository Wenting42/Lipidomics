####################################################################################
# Script: FWL_lipidomics.1.1.R
# Author: Wenting 
# Notes: 
#         This script is based on Niklas and Kenny's previous massive work. 
#         It helps generating the graph and data for the flowork of Lipidomics. 
#
#         To start, typing command in the console-----> source("FWL_lipidomics2.0.R") 
#
#         1) interactive way.  
#         First please store the file name in the variable file.name for check 
#         Run the script by command "source("FWL_lipidomics2.0.R") in console.
#         All the data are stored in the directory called "data", and plots in "plot"
#         2) You can change the parameters in the script and run it lines by lines.
#           Before you change any parameters, please copy the original script first.
#           
# Warning: Please DO NOT change png format to pdf which 
#           might induce troubles for displaying plots! 
#         ctrol+c is the function to quit the script
#####################################################################################


rm(list=ls())  ##### empty the environment first. If you need to keep former variables in workspace, just comment this part
setSessionTimeLimit(cpu = Inf, elapsed = Inf)
##### Please install the lacking packages ------> install.packages("package")
list.of.packages <- c("FactoMineR", "tidyverse", "ggplot2", "magrittr", "ggrepel", "reshape2",
                      "stargazer", "Hmisc", "limma", "factoextra", "scales", "RColorBrewer",
                      "stringr", "readxl")
need.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if(length(need.packages)) BiocManager::install(need.packages)


lapply(list.of.packages, function(x) 
  suppressMessages(require(x, character.only = TRUE, quietly = TRUE)))

mkdirs <- function(x){
  for(i in 1:length(x))
  if(!file.exists(x[i])){
    mkdirs(dirname(x[i]))
    dir.create(x[i])
  }
}

dirs <- c("plot", "data", "plot/classes")
mkdirs(dirs)

message('Please store the file name in the file.name for check. And you could go to script to add the file name in read_csv(file.choose("data_raw/")). 
        If you need to change some parameters or edit the script, please copy the original script first. ')
file.name <- readline("Please input the file name: ")

# Read in file, and store the data 
lipidomics <- read_csv(file.choose("data_raw/"), col_types = cols())

### REMOVE ODD NUMBER IN THE MOLEC
#standards <- c("LPA(17:1)", "PA(37:4)", "LPS(17:1)", "PS(17:0/20:4)")


message("Do you have internal controls like the standards or solvent?")
internal.check <- readline("Please type Y/N: ") %>% str_to_lower(.)
if(internal.check!="n"){
  
  message("If your experiment includes standard control and solvent control, please input the number of the sample. And if only one control sample need to delete, just type the same sample two times or some unexist sample for the other option. 
          e.g. sxx")
  
  message("Sample ID for extracted and unextracted internal control, eg. s22")
  extracted.sample <- as.character(readline("extracted.sample -----> ")) %>% str_to_lower(.)
  extracted.Grade <- paste("Grade[", extracted.sample, "]", sep="")
  extracted.P.No <- paste("APValue[", extracted.sample, "]", sep="")
  unextracted.sample <- as.character(readline("unextracted.sample -----> ")) %>% str_to_lower(.)
  unextracted.Grade <- paste("Grade[", unextracted.sample, "]", sep="")
  unextracted.P.No <- paste("APValue[", unextracted.sample, "]", sep="")
  
  
  
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
  
}else{
  lipidCount <- lipidomics %>% 
    select(contains("Grade"), -contains("Grade[c]" ),
           contains("APValue"))  %>% 
    transmute(A = rowSums(. == "A"), B = rowSums(. == "B"), 
              C = rowSums(. == "C"), D = rowSums(. == "D"), 
              No.grade = rowSums(.=="-"), APvalue.001 = rowSums(. <= "0.001"))
  
}
# add this count into original strucutes
lipidSelect <- lipidomics %>% bind_cols(lipidCount)

# Filter the dataset based on your criteria  
# flexible parameter 4 which depends on the experiment for total number of A and B 
message("Please note the file is filtered by three standards. For same lipidmolecule: \n 1. Rej (rejection identification) =0; 2. at least k number of Grade B; 3. there are at least 4 p value less or equal than 0.001.\n The filtered data is stored in the filtered.data.csv.")
k <- readline("Please type the number of at least B grade you want, e.g. 4: ") %>% as.numeric(.)
j <- readline("Please type the number of at least p values less or equal than 0.001 you want, e.g. 4: ") %>% as.numeric(.)
filtered.lipidomics <- lipidSelect %>% 
  rowwise() %>% 
  filter( Rej == 0 & sum(A, B) >= k & APvalue.001 >= j)

### remove artifical lipids
# mark the artifical classes position
# odd.index <- filtered.lipidomics$LipidMolec %>%
#   str_locate(., "(\\d[13579]:)|([13579]:)")
# percent.odd <- length(unique(which(!is.na(odd.index), arr.ind=TRUE)[,1]))/nrow(filtered.lipidomics)
# 
# a <- filtered.lipidomics %>% filter(Class=="TG") 
# b <- a$LipidMolec %>%   str_locate(., "(17:)|(7:)")
# b <- a$LipidMolec %>%   str_locate(., "(\\d[13579]:)|([13579]:)")
# odd.chain <- length(unique(which(!is.na(b), arr.ind=TRUE)[,1]))/nrow(a)
# 
# # odd in total is 11%
# # TG in total is 21.5%
# nrow(a)/nrow(filtered.lipidomics)
# # odd in TG is 13.8% while 17: odd chain is 7.3%

# print(percent.odd)
# print(odd.chain)

# #############    remove the artifical classes
# filtered.lipidomics <- filtered.lipidomics %>%
#                         slice(unique(which(is.na(odd.index), arr.ind=TRUE)[,1]))
message("The info below and summary plot will show the summary information of classes after filtering the data")
# check how many lipids passed filtering
print(describe(filtered.lipidomics$Class))
print(summary(as.data.frame(filtered.lipidomics$Class), maxsum=nrow(filtered.lipidomics)))


prop.bw.plot <- ggplot(filtered.lipidomics, aes(x=reorder(filtered.lipidomics$Class, filtered.lipidomics$Class, function(x)length(x)))) + 
  geom_bar(aes(y = (..count..)/sum(..count..)), stat="count")+ 
  scale_y_continuous(labels=scales::percent) +
  theme_bw() +
  ylab("relative frequencies") +
  scale_colour_grey(end=0) +
  ggtitle("Class Summary plot: \n Proportion of different class numbers among all the samples")+
  xlab("classes")+ ylab("Relative frequencies")  +
  coord_flip() 

print(prop.bw.plot) 

ggsave(filename = "Prop.bw.summary.pdf", path = 'plot/', device = "pdf")



prop.plot <- ggplot(filtered.lipidomics, aes(x=reorder(Class, Class, function(x)length(x)), fill = Class)) + 
  geom_bar(aes(y = (..count..)/sum(..count..)), stat="count")+ 
  scale_y_continuous(labels=scales::percent) +
  theme_bw() +
  ylab("relative frequencies") +
  scale_color_brewer() +
  ggtitle("Class Summary plot: \n Proportion of different class numbers among all the samples")+
  xlab("lipid classes")+ ylab("Relative frequencies (%)")  +
  coord_flip() 

print(prop.plot) 

ggsave(filename = "Prop.summary.pdf", path = 'plot/', device = "pdf")


##########################################################################
# QC PLOT 1 - Abundance vs. retention time [requires ggplot2]
##############################################################################################
# Interactive way
##############################################################
# individual sample retention plot

# message("Please input the MainArea of the sample name to check its abundance vs. retention time, 
#         e.g. `MainArea[s1]`. 
#         Please note that: ` ` is the backtick which is same position as tilde")
# 
# abundance.sample <- readline("Please type the sample name as this format. e.g. `MainArea[s1]` -----> ")
# 
# retention.plot <- ggplot(data=filtered.lipidomics, 
#                          aes_string(x = sprintf("log10(%s)", abundance.sample), y = "BaseRt")) +
#   geom_point() +
#   theme_bw() +
#   facet_grid(.~Class) +
#   labs("Abundance VS. Retention time", 
#        x="lipid class (log10(MainArea))", y="Retention time")
# print(retention.plot)
# message("Please input the name you want to store for the graph. e.g. retention.pdf")
# plot.name <- readline("QC plot name: ")
# ggsave(filename = plot.name, path = 'plot/', device = "pdf", dpi=300)



retention.data <- filtered.lipidomics %>% select(contains("MainArea[s"), BaseRt, Class) %>% gather(sample, MainArea, -c(BaseRt, Class))
retention.all.plot <- ggplot(data=retention.data, 
                             aes(x=log10(MainArea), y = BaseRt)) +
  geom_point() +
  theme_bw() +
  facet_grid(.~Class) +
  labs("Abundance VS. Retention time", 
       x="All samples (log10(MainArea))", y="Retention time (mins)")
print(retention.all.plot)
ggsave(filename = "all.retention.png", path = 'plot/', device = "png", dpi=300)
#ggsave(filename = "all.retention.png", path = 'plot/', device = "png", width = 10, height = 8, dpi = 150, units = "in")


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



############################################################################
options(scipen=999)

# making group
message("Please input the info of the experimental groups below and end with 'Enter'
        Please Only input the sample names (replicate) for analysis.
        e.g. Experimental Group number: 2   
        Name of Experimental Group 1: wild_type
        Sample names (replicate number from lipid search) : s1 s2 s3 s4 
        Name of Group 2: knockout
        Sample names of Group 2: s4 s5 s6")


Input <- function(filtered.lipidomics){
# info about group size
ngroups <- readline("Group number: ") %>% as.numeric(.)
CheckType <- function(x){
  if((x%%1 != 0) | (is.na(x))){
    print("Group number must be numeric!")
    x <- readline("Group number: ") %>% as.numeric(.)
    return(CheckType(x))
  }else{
    if(x==0){
      print("You input 0 groups")
      x <- readline("Group number: ") %>% as.numeric(.)
      return(CheckType(x))
    }else{
      return(x)
    }
  }
}
ngroups <- CheckType(ngroups)

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
while(str_to_lower(group.condition) == "y"){
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
}

# store the sample and index in a list
sample.info <- grepIndex(total.info, filtered.lipidomics) 
message("Take a look at the sample info and its column position information in the file below")
glimpse(sample.info)

# retrieve the group names
group.names <- dimnames(sample.info)[[2]]

return(list(sample.info, group.names, ngroups))

}

samples <- Input(filtered.lipidomics)
sample.info <- samples[[1]]
group.names <- samples[[2]]
ngroups <- samples[[3]]









pca.pairs.plots <- function(sample.info, group.names, filtered.lipidomics, mark){

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



#############################################################################################################
# this part reformats the axis label positions on same sides. Source: stackflow
pairs2 <- function (x, labels, panel = points, ..., 
                    lower.panel = panel.smooth, diag.panel=panel.hist, 
                    upper.panel = panel.cor, text.panel = textPanel, 
                    label.pos = 0.5 + has.diag/3, cex.labels = NULL, font.labels = 1, 
                    row1attop = TRUE, gap = 1) {
  textPanel <- function(x = 0.5, y = 0.5, txt, cex, font) text(x, 
                                                               y, txt, cex = cex, font = font)
  localAxis <- function(side, x, y, xpd, bg, col = NULL, main, 
                        oma, ...) {
    if (side%%2 == 1) 
      Axis(x, side = side, xpd = NA, ...)
    else Axis(y, side = side, xpd = NA, ...)
  }
  localPlot <- function(..., main, oma, font.main, cex.main) plot(...)
  localLowerPanel <- function(..., main, oma, font.main, cex.main) lower.panel(...)
  localUpperPanel <- function(..., main, oma, font.main, cex.main) upper.panel(...)
  localDiagPanel <- function(..., main, oma, font.main, cex.main) diag.panel(...)
  dots <- list(...)
  nmdots <- names(dots)
  if (!is.matrix(x)) {
    x <- as.data.frame(x)
    for (i in seq_along(names(x))) {
      if (is.factor(x[[i]]) || is.logical(x[[i]]))
        x[[i]] <- as.numeric(x[[i]])
      if (!is.numeric(unclass(x[[i]])))
        stop("non-numeric argument to 'pairs'")
    }
  }else if(!is.numeric(x)){
    stop("non-numeric argument to 'pairs'")
  }
  
  panel <- match.fun(panel)
  if ((has.lower <- !is.null(lower.panel)) && !missing(lower.panel))
    lower.panel <- match.fun(lower.panel)
  if ((has.upper <- !is.null(upper.panel)) && !missing(upper.panel))
    upper.panel <- match.fun(upper.panel)
  if ((has.diag <- !is.null(diag.panel)) && !missing(diag.panel))
    diag.panel <- match.fun(diag.panel)
  if (row1attop) {
    tmp <- lower.panel
    lower.panel <- upper.panel
    upper.panel <- tmp
    tmp <- has.lower
    has.lower <- has.upper
    has.upper <- tmp
  }
  nc <- ncol(x)
  if (nc < 2)
    stop("only one column in the argument to 'pairs'")
  has.labs <- TRUE
  if (missing(labels)) {
    labels <- colnames(x)
    if (is.null(labels))
      labels <- paste("var", 1L:nc)
  }else if (is.null(labels))
    has.labs <- FALSE
  oma <- if ("oma" %in% nmdots){
    dots$oma
  }else NULL
  main <- if ("main" %in% nmdots){
    dots$main
  }else NULL
  if (is.null(oma)) {
    oma <- c(4, 4, 4, 4)
    if (!is.null(main))
      oma[3L] <- 6
  }
  opar <- par(mfrow = c(nc, nc), mar = rep.int(gap/2, 4), oma = oma)
  on.exit(par(opar))
  dev.hold()
  on.exit(dev.flush(), add = TRUE)
  for (i in if (row1attop)
    1L:nc
    else nc:1L) for (j in 1L:nc) {
      localPlot(x[, j], x[, i], xlab = "", ylab = "", axes = FALSE,
                type = "n", ...)
      if (i == j || (i < j && has.lower) || (i > j && has.upper)) {
        box()
        # edited here...
        #           if (i == 1 && (!(j%%2) || !has.upper || !has.lower))
        #           localAxis(1 + 2 * row1attop, x[, j], x[, i],
        #                       ...)
        # draw x-axis
        if (i == nc & j != nc)
          localAxis(1, x[, j], x[, i],
                    ...)
        # draw y-axis
        if (j == 1 & i != 1)
          localAxis(2, x[, j], x[, i], ...)
        #           if (j == nc && (i%%2 || !has.upper || !has.lower))
        #             localAxis(4, x[, j], x[, i], ...)
        mfg <- par("mfg")
        if (i == j) {
          if (has.diag)
            localDiagPanel(as.vector(x[, i]), ...)
          if (has.labs) {
            par(usr = c(0, 1, 0, 1))
            if (is.null(cex.labels)) {
              l.wid <- strwidth(labels, "user")
              cex.labels <- max(0.8, min(2, 0.9/max(l.wid)))
            }
            text.panel(0.5, label.pos, labels[i], cex = cex.labels,
                       font = font.labels)
          }
        }
        else if (i < j)
          localLowerPanel(as.vector(x[, j]), as.vector(x[,
                                                         i]), ...)
        else localUpperPanel(as.vector(x[, j]), as.vector(x[,
                                                            i]), ...)
        if (any(par("mfg") != mfg))
          stop("the 'panel' function made a new plot")
      }
      else par(new = FALSE)
    }
  if (!is.null(main)) {
    font.main <- if ("font.main" %in% nmdots)
      dots$font.main
    else par("font.main")
    cex.main <- if ("cex.main" %in% nmdots)
      dots$cex.main
    else par("cex.main")
    mtext(main, 3, 3, TRUE, 0.5, cex = cex.main, font = font.main)
  }
  invisible(NULL)
}
######################################################################


# retrieve the sample info position 
index         <- 1:length(sample.info)
sample.index  <- index[index %% 2 ==0]

######################################################################
# QC PLOT 2 - Pair-wise correlation between replicates, reformatted axis
for(i in 1:length(sample.index)){
  range     <- sample.index[i]
  plot.name <- paste("pairs.plot.", i, ".", mark, ".pdf",sep="")
  path      <- file.path("plot/", plot.name)
  pairs2(log10(filtered.lipidomics[, sample.info[[range]]]), 
         lower.panel = panel.smooth, diag.panel=panel.hist, 
         upper.panel = panel.cor, xlim=c(4,12), ylim=c(4,12))
  dev.copy(pdf, path)
  dev.off()
}

#################################################################


# Storing pairwise plot in the plot directory
for(i in 1:length(sample.index)){
  range     <- sample.index[i]
  plot.name <- paste("pairs.plot.", i, ".", mark, ".oldversion.pdf",sep="")
  path      <- file.path("plot/", plot.name)
  pairs(log10(filtered.lipidomics[, sample.info[[range]]]), 
        lower.panel = panel.smooth, diag.panel=panel.hist, 
        upper.panel = panel.cor, xlim=c(4,12), ylim=c(4,12))
  dev.copy(pdf, path)
  dev.off()
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

# if just show the sample names by removing the MainArea prefix
rownames(filtered.lipidomics.PCA) <- filtered.lipidomics.PCA %>% 
  rownames(.) %>% 
  str_remove_all(., "MainArea\\[") %>% str_remove_all(., "\\]")

log2.filtered.lipidomics.PCA      <- log((filtered.lipidomics.PCA+1), 2)

# making group info for PCA 
filtered.lipidomics.PCA$Group       <- group.repeats
log2.filtered.lipidomics.PCA$Group  <- group.repeats

# Perform PCA 
res.pca         <-  PCA(filtered.lipidomics.PCA, scale.unit=TRUE, 
                        quali.sup=ncol(filtered.lipidomics.PCA), graph=T)



print(get_eigenvalue(res.pca))





concat          <-  cbind.data.frame(filtered.lipidomics.PCA[, ncol(filtered.lipidomics.PCA)], 
                                     res.pca$ind$coord)
ellipse.coord   <-  coord.ellipse(concat, bary=TRUE)




#pdf(file="plot/sample.pca.pdf")
plot.PCA(res.pca, habillage=ncol(filtered.lipidomics.PCA), 
         ellipse=ellipse.coord, cex=0.8, label="all", ncp=5)
dev.copy(pdf, "plot/sample.pca.pdf")
dev.off()


plot.PCA(res.pca, habillage=ncol(filtered.lipidomics.PCA), 
         ellipse=ellipse.coord, cex=0.8, label="quali")
dev.copy(pdf, "plot/group.pca.pdf")
dev.off()

# Output the filtered lipidome
write.csv(filtered.lipidomics, "data/filtered.lipidomics.csv")
###################################################################################### 
return(list(sample.list, group.repeats))
}
label <- "first"
info.list <-  pca.pairs.plots(sample.info, group.names, filtered.lipidomics, label)





message("Please check the plots in the plot directory or r studio plots pannel.\nDo you need to re-enter the sample info for analysis?")
pca.check <- readline("Please type Y/N: ")
if(str_to_lower(pca.check)=="y"){
  samples <- Input(filtered.lipidomics)
  sample.info <- samples[[1]]
  group.names <- samples[[2]]
  ngroups <- samples[[3]]
  label <- "new"
  info.list <-  pca.pairs.plots(sample.info, group.names, filtered.lipidomics, label)
}

sample.list <- info.list[[1]]
group.repeats <- info.list[[2]]
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
  
  total.class.var <- total.class.sd %>% transmute(experiment.group=experiment.group, class=class, var=sd^2)
  
  
 
  
  # aggregating data by group
  total.class.wide <- total.class.data  %>% group_by(experiment.group) %>% summarise_all(.funs=sum)  %>% view()
  # reformat data into long data
  total.class.long <- total.class.wide %>% gather(class, value, -experiment.group)
  # combine data and its sd info for ggplot
  total.class.long <- bind_cols(total.class.long, sd=total.class.sd$sd)
  
  
  # plot the total class info by groups
  # the y axis range can be editted by changing the parameter in "coord_cartesian(ylim=c(1e5,1e13))"
  #quartz()
  total.plot <- ggplot(total.class.long, aes(x=reorder(class, value), fill=experiment.group, y=value)) + 
    geom_bar(position="dodge", stat="identity") +
    geom_errorbar(aes(ymin=value, ymax=value+sd), width=.3,
                  position=position_dodge(.9)) +
    scale_y_continuous(trans=log10_trans(),
                       labels = trans_format("log10", math_format(10^.x)), 
                       limits=c(1, 1e20)) +
    theme_bw() +
    xlab("Lipid Classes") +
    ylab("AUC (Area under curve)") +
    ggtitle("Total lipid classes") +
    theme(plot.title = element_text(hjust = 0.5), 
          axis.text.x = element_text(angle=45, hjust=1)) +
    coord_cartesian(ylim=c(1e3, 1e14)) ## these y axis range parameter can be changed
  print(total.plot)
  ggsave(filename="total.class.png", path = 'plot/', device = "png", dpi = 300)
  
  
  FoldPlot <- function(data){
    ggplot(data, aes(x=class, fill=experiment.group, y=value)) +
      geom_col(position="dodge") +
      geom_hline(yintercept=1, linetype="dashed", colour="red", size=0.2) +
      theme_bw() +
      coord_flip() +
      geom_hline(yintercept=1, linetype="dashed", colour="grey45", size=0.02) +
      theme(plot.title = element_text(hjust = 0.5))
  }
  
  # design control groups:
  n.control <- readline("Please input how many controls for normalizations you want, e.g. 1: ")
  control.group <- c()
  for(i in 1:n.control){
    info <- paste("Please input name of control Group", i," e.g. A: ")
    control.group[i] <- readline(prompt = info)
    contrast.groups <- readline("Please input the groups which are split by space for contrast, e.g. B C D: ")
    contrast.groups <- contrast.groups %>% str_split(., "\\s+") %>% unlist()
    control.index <- which(total.class.wide$experiment.group %in% control.group[i])
    contrast.index <- which(total.class.wide$experiment.group %in% contrast.groups)
    normalization.factor <- apply(total.class.wide[contrast.index, -1], 1, 
                                  function(x)log2(x/total.class.wide[control.index,-1]+1)) %>% 
      unlist() %>%
      matrix(., nrow=length(contrast.index), byrow=TRUE) %>% 
      cbind(total.class.wide[contrast.index, 1], .)  
    
    dd <- total.class.wide %>% slice(control.index, contrast.index)
    log.dd <- dd
    log.dd[, -1] <- log2(dd[,-1]) %>% view()
    
    subc <- log.dd[1, -1] %>% unlist(.)
    log.dd[ ,-1] <- sweep(log.dd[,-1], 2, subc, "-")
    log.dd[1, -1] <- apply(log.dd[1,-1], 1, function(x) ifelse(x==0, 1, 1))
    
    od <- total.class.wide %>% filter(experiment.group==control.group) %>% 
      select(-experiment.group) %>% apply(., 1, function(x) x <- order(x)) 
    odf <- total.class.wide %>% select(1+od) %>% colnames(.)
    
    d.long <- log.dd %>% gather(class, value, -experiment.group)
    o.long <- d.long %>% mutate(class= factor(class, levels = odf))  %>% arrange(class)
    
    ordered.foldchange.plot <- FoldPlot(o.long) + 
      labs(title="Total lipid classes",x="Lipid class", y="Log2 fold change", 
           caption = "The axis of classes are ordered decreasingly by control group") 
    
    print(ordered.foldchange.plot)
   
    filenames <- paste("orderedControl.class.log2.oo", ".png")
    ggsave(filename=filenames, path = 'plot/', device = "png", dpi = 300)
    
    
    
    colnames(normalization.factor) <- colnames(total.class.wide)
    
    normal.log2.matrix <- bind_rows(total.class.wide[control.index, ], normalization.factor)            
    normal.log2.matrix[1, -1] <- rep(1, ncol(normal.log2.matrix)-1)
    
    normal.log2.data <- normal.log2.matrix %>% gather(class, value, -experiment.group) 
   
    
    #rank of variance
    control.var <- total.class.var %>%
                    filter(experiment.group==control.group[i]) %>%
                    arrange(var) %>%
                    transmute(experiment.group=experiment.group, class=class, var=var, rank=log2(var/var[1]+1))


    rank.control.class.var <- ggplot(data=control.var)+ aes(x= reorder(class, rank), y=rank) +
                              geom_bar(position="dodge", stat="identity") +
                              labs(title = "Relative variances of classes in the control group",
                                   x = "rank of variance",
                                   y = "variance of classes")

    print(rank.control.class.var)

    
    total.normalization.plot <-  ggplot(normal.log2.data, aes(x=reorder(class, value), fill=experiment.group, y=value)) +
                                  geom_col(position="dodge") +
                                  geom_hline(yintercept=1, linetype="dashed", colour="red", size=0.2) +
                                  theme_bw() +
                                  coord_flip() +
                                  geom_hline(yintercept=1, linetype="dashed", colour="grey45", size=0.02) +
                                  theme(plot.title = element_text(hjust = 0.5)) +
                                labs(title="Total lipid classes",x="Lipid class", y="Log2 fold change")
    
    print(total.normalization.plot)
    filenames <- paste("class.log2.", i, ".png")
    ggsave(filename=filenames, path = 'plot/', device = "png", dpi = 300)
   
    # order the class axis by control group of decreasing order 
    order.factor <- total.class.wide %>% 
      filter(experiment.group==control.group[i]) %>% 
      select(-1) %>% 
      apply(., 1, function(x)sort(x)) %>% rownames(.)
    
    order.dt <- normal.log2.data %>% mutate(class=factor(class, levels = order.factor)) %>% arrange(class)
    ordered.foldchange.plot <- FoldPlot(order.dt) + 
                                labs(title="Total lipid classes",x="Lipid class", y="Log2 fold change", 
                                     caption = "The axis of classes are ordered decreasingly by control group") 
                                
    print(ordered.foldchange.plot)
    filenames <- paste("orderedControl.class.log2.", i, ".png")
    ggsave(filename=filenames, path = 'plot/', device = "png", dpi = 300)
    
  }
  
}

totalClass(filtered.lipidomics, sample.list, group.repeats)

































# # input the plot name and store it
# message("Please input the name you want to store for the graph. e.g. total.class.pdf")
# name.class  <- readline("Total lipid classes plot name: ")
# ggsave(filename = name.class, path = 'plot/', device = "pdf")

############################
# For each class


###########################

# plot every class
classes <- unique(filtered.lipidomics$Class)
for(i in 1:length(classes)){
  options(warn=-1)
  pick.class <- classes[i]
  print(pick.class)
  
  
  # get information of the samples
  filtered.class <- filtered.lipidomics %>% 
    select(Class, LipidMolec, sample.list) %>% 
    filter(Class==pick.class)
  # aggregate the repeated data by its class and lipid molecular names
  class.data <- aggregate(.~Class+LipidMolec, data=filtered.class, FUN=sum) %>%  select(-Class)
  
  # rename the class
  molec.names <- class.data[, 1] %>% 
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
    merge.data.frame(., class.sd) %>% 
    arrange(Acyl)
  
  
  
  
  # setting the plot limits when the bar numbers exceed the threshold n.bar
  n.bar <- 60     # estimation of threshold which can be modified and at least bigger than group number
  n.observations <- nrow(class.long) 
  n.groups <- as.numeric(ngroups) 
  
  quartz()
  classPlot <- function(data){ggplot(data, aes(x=reorder(Acyl, AreaValue), fill=experiment.group, y=AreaValue)) + 
      geom_bar(position="dodge", stat = "identity") +
      geom_errorbar(aes(ymin=AreaValue, ymax=AreaValue+sd), 
                    width=.8, position=position_dodge(.9)) +
      theme_bw() +
      scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), # base of 10
                    labels = trans_format("log10", math_format(10^.x)), limits=c(NA, 1e20)) +
      theme(plot.title = element_text(hjust = 0.5), 
            axis.text.x = element_text(angle=45, hjust=1)) +
      labs(title= pick.class, x="Acyl composition", y="Main Area", caption="Error bar is the standard deviation for each class in each group") +
      coord_cartesian(ylim=c(1e3,1e14)) 
  }
  
  
  
  if(n.observations <= n.bar){
    p1 <- classPlot(class.long)
    print(p1)
    #ggsave(filename = paste(pick.class, ".png", sep=""), path = 'plot/classes', device = "png", width=15, height=15, dpi=300)
    ggsave(filename = paste(pick.class, ".png", sep=""), path = 'plot/classes', device = "png")
    
  }else{
    #class.long$bins <- 1:n.observations
    if(n.bar %% n.groups != 0){
      n.bar <- (n.bar%/% n.groups)*n.groups
    }
    n.facet <- n.observations %/% n.bar
    for(k in 0:n.facet){
      #m <- paste("a",i, sep="")
      if(k < n.facet){
        #class.long$bins[((n.bar*i+1):(n.bar*(i+1)))] <-  rep(m, n.bar)
        row.range <- (n.bar*k+1):(n.bar*(k+1))
        data <- class.long %>% slice(row.range) 
        p2 <- classPlot(data)
        print(p2)
        ggsave(filename = paste(pick.class, ".", k+1, ".png", sep=""), path = 'plot/classes/', device = "png")
      }else{
        #class.long$bins[(n.bar*i):nrow(class.long)] <- rep(m, nrow(class.long)%%n.bar)
        row.range <- (n.bar*k+1):nrow(class.long)
        data <- class.long %>% slice(row.range) 
        p2 <- classPlot(data)
        print(p2)
        ggsave(filename = paste(pick.class, ".", k+1, ".png", sep=""), path = 'plot/classes/', device = "png")
      }
    }
  }
}


# 
# # message("How many more classes you want to visualize for barplots?")
# ntimes <- readline("Please input the numbers of classes you want to visualize: ")

# define the function for ploting
# eachClassPlots <- function(){
#   message("which class you want to visualize for barplot?")
#   pick.class <- as.character(readline("Please input the class for barplot visualization: "))
#   # get information of the samples
#   filtered.class <- filtered.lipidomics %>%
#     select(Class, LipidMolec, sample.list) %>%
#     filter(Class==pick.class)
#   # aggregate the repeated data by its class and lipid molecular names
#   class.data <- aggregate(.~Class+LipidMolec, data=filtered.class, FUN=sum) %>%  select(-Class)
#   # rename the class
#   molec.names <- class.data[,1] %>%
#     unlist() %>%
#     str_extract_all(., "\\(\\d*.*\\)") %>%
#     str_remove_all(., "[\\(\\)]")
# 
#   # Making the group information for the data
#   data.wide <- data.frame(group.repeats, t(class.data[, -1]))
#   # making names for the new data frame
#   names(data.wide) <- c("experiment.group", molec.names)
# 
#   # calculate the sd information
#   class.sd <- aggregate(.~experiment.group, data=data.wide, function(x) sd(x)) %>%
#     gather( Acyl, sd, -experiment.group)
# 
#   # reformat the data information
#   class.long <- aggregate(.~experiment.group, data=data.wide, FUN=sum) %>%
#                 gather(Acyl, AreaValue, -experiment.group) %>%
#                 merge.data.frame(., class.sd)
# 
#   # plot for the individual class by groups
#   class.plot <- ggplot(class.long, aes(x=reorder(Acyl, AreaValue), fill=experiment.group, y=AreaValue)) +
#     geom_bar(position="dodge", stat="identity") +
#     geom_errorbar(aes(ymin=AreaValue, ymax=AreaValue+sd),
#                   width=.8, position=position_dodge(.9)) +
#     theme_bw() +
#     scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), # base of 10
#                   labels = trans_format("log10", math_format(10^.x)), limits=c(NA, 1e20)) +
#     #scale_y_continuous(labels = scales::scientific)+    # scientific notation
#     #scale_y_continuous(trans="log10", limits=c(NA, 1))
#     #labs(ylab="Main Area", fill="Groups", xlab="AT") +
#     ggtitle(pick.class) +
#     theme(plot.title = element_text(hjust = 0.5),
#           axis.text.x = element_text(angle=45, hjust=1)) +
#     xlab("Acyl composition") +
#     ylab("Main Area") +
#     coord_cartesian(ylim=c(1e3,1e14))
# 
#   print(class.plot)
# 
#    message("Please input the name you want to store for the graph. e.g. TG.pdf")
#    class.type <- readline("The lipid class plot name: ")
#    ggsave(filename = class.type, path = 'plot/', device = "pdf", width=15, height=15, dpi=300)
#  }
# # # plot for missed classes
#  if(ntimes != "0"){
#   replicate(ntimes, eachClassPlots())
# }
# #
# # 




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
message("Begin to plot volcano \n Please note that you NEED make contrast groups manually to Compare the difference between/among groups.
        e.g. compare B+C against A: A-(B+C)/2; A against B:  A-B; A-B against C-D: (A-B)-(C-D), etc.")
message("How many groups you designed for making comparison?")
n.comparisons <- readline("Number of comparison groups-----> ")

replicate( n.comparisons, {
  cont.matrix <- makeContrasts(
    readline("Please input the groups names for comparison, '-' means the contrast, e.g A-B: "),
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
                                 p.value =0.05, lfc=log2(1.5), number=nrow(filtered.lipids)) %>% 
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
  
  #contrast.info <- dimnames(cont.matrix)[[2]] %>% str_split(., "-") %>% unlist()
  #title.info <- paste(contrast.info)
  # making volcano plot        
  volc.plot<- ggplot(input, aes(logFC, -log10(adj.P.Val))) + #volcanoplot with log2Foldchange versus pvalue
    geom_rect(aes(xmin = -fold.change, xmax = fold.change, ymin = -Inf, ymax = Inf),
              fill = "#f0f0f0", linetype="dashed", color="red", size=0.2)+
    geom_point() +
    geom_vline(xintercept=0, linetype="solid", colour="grey60", size=0.2) +
    #Add a horizontal line for P-value cut-off
    geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="red", size=0.2)+
    #geom_hline(yintercept=-log10(0.01), linetype="dashed", colour="red", size=0.4)+
    scale_y_continuous()+  # Ticks from 0-10, every .25
    geom_point(aes(col=sig, size=AveExpr)) + #add points colored by significance, point size by average value
    scale_size_continuous(range = c(0.25,2.5), breaks=c(10, 20, 30))+ #add size limit
    scale_color_manual(values=c("#bdbdbd", "#3182bd", 
                                'chartreuse4', '#deebf7',"darkorange3"))+
    xlab(bquote("Fold change, " ~ log[2]~"("~textstyle(frac({x}, {y}))~scriptstyle(, AUC)~")")) + ylab("-Log10(q value)") + 
    ggtitle("volcano plot") +
    theme_bw()+ # remove background
    theme(line=element_blank()) +
    #scale_y_continuous(expand = c(0, 0), limits=c(0,15)) +       ###### start at 0
    #geom_text_repel(data=points, aes(label=rownames(points)), size=2)#adding text for the FDR points.
    geom_text_repel(data=points, aes(label=rownames(points)), size=2, show.legend =FALSE) 
    
  
  
  
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
  
  
  fold.points <- input %>% 
    rownames_to_column('lipid') %>% 
    filter(abs(logFC) > 5 & adj.P.Val>= 0.05 ) %>% 
    column_to_rownames('lipid')
  
  
  # Recolor volcano plot 
  volc.plot.color <- ggplot(input.lip, aes(logFC, -log10(adj.P.Val))) + #volcanoplot with log2Foldchange versus pvalue
    geom_rect(aes(xmin = -fold.change, xmax = fold.change, ymin = -Inf, ymax = Inf),
              fill = "#f0f0f0", linetype="dashed", color="red", size=0.2)+
    geom_point()+
    geom_vline(xintercept=0, linetype="solid", colour="grey60", size=0.2) +
    #Add a horizontal line for P-value cut-off
    geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="red", size=0.2)+
    #geom_hline(yintercept=-log10(0.01), linetype="dashed", colour="red", size=0.4)+
    scale_y_continuous()+  # Ticks from 0-10, every .25
    #geom_point(aes(col=lipclass, size=AveExpr)) + #add points colored by significance, point size by average value
    scale_size_continuous(range = c(0.25,2.5), breaks=c(10, 20, 30))+ #add size limit
    
    #ggtitle ('Total Lipid species') +
    #xlim(c(-2.7, 2.7)) + ylim(c(0, 2)) + #set the x,y axis limits.
    xlab(bquote("Fold change, " ~ log[2]~"("~textstyle(frac({x}, {y}))~scriptstyle(, AUC)~")")) + ylab("-Log10 (q value)") + 
    theme_bw()+ # remove background
    theme(line=element_blank())+
    #scale_y_continuous(expand = c(0, 0), limits=c(0,15)) +       ###### start at 0
    #geom_text_repel(data=points, aes(label=rownames(points)), size=2)#adding text for the FDR points.
    geom_text_repel(data=sig.points, aes(label=rownames(sig.points)), size=2)  
  
    volc.plot.color <- volc.plot.color+ 
      geom_text_repel(data=fold.points, aes(label=rownames(fold.points)), size=2) + scale_color_manual(values=c("darkorange3","dodgerblue3",  'chartreuse4', '#c51b8a', '#756bb1','#969696'), drop=FALSE) + 
    geom_point(aes(col=lipclass, size=AveExpr)) 
  #geom_rug()
  
  print(volc.plot.color)
  plot.name <- readline("Please input the volcano plot name: ")
  ggsave(filename = plot.name, path = 'plot/', device = "pdf", width=15, height=15, dpi=300)
 # dev.off()

})





