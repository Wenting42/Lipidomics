
####################################################################################
# Script: FWL_lipidomics.2.2.functions.R
# Author: Niklas, Kenny, Wenting
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

##### Please install the lacking packages ------> install.packages("package")
source("https://install-github.me/tidyverse/rlang")
if(!length("BiocManager")) install.packages("BiocManager")
list.of.packages <- c("rlang","FactoMineR", "scales", "magrittr", "ggrepel", "reshape2",
                      "stargazer", "Hmisc", "limma", "factoextra", "scales", "RColorBrewer",
                      "stringr", "readxl", "RCy3", "igraph", "tidyverse")
need.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if(length(need.packages)) BiocManager::install(need.packages)

# Suppress the package messages
lapply(list.of.packages, function(x) 
  suppressMessages(require(x, character.only = TRUE, quietly = TRUE)))


########################################################################
### Functions
########################################################################
# function name: mkdirs
# utility: Check directory existence, if not then make one
########################################################################
mkdirs <- function(x){
  for(i in 1:length(x))
    if(!file.exists(x[i])){
      mkdirs(dirname(x[i]))
      dir.create(x[i])
    }
}


########################################################################
# function name: InternalCheck
# parameter: internals (y/n)
# utility: This function will check if there are internal control samples exist. 
#          It will delete the internal controls for later PCA and correlations analysis. 
#         It will also make new columns which store the information of sums of Grade A, B, 
#         C and D, plus the p values which are less or equal than 0.001 for each row. It will return the 
#         the extra information
########################################################################

InternalCheck <- function(internals){
  if(internals=="y"){
    message("If your experiment includes standard control and solvent control, please input the number of the sample. And if only one control sample need to delete, just type the same sample two times or some unexist sample for the other option.
        e.g. sxx")
    message("Sample ID for extracted and unextracted internal control, eg. s22")
    extracted.sample <- as.character(readline("extracted.sample -----> ")) %>% str_to_lower(.)
    extracted.Grade <- paste("Grade[", extracted.sample, "]", sep="")
    extracted.P.No <- paste("APValue[", extracted.sample, "]", sep="")
    unextracted.sample <- as.character(readline("unextracted.sample -----> ")) %>% str_to_lower(.)
    unextracted.Grade <- paste("Grade[", unextracted.sample, "]", sep="")
    unextracted.P.No <- paste("APValue[", unextracted.sample, "]", sep="")
   
    # select and count the sample grades you want, filter the p value which less than 0.001
    lipidCount <- lipidomics %>%
      select(contains("Grade"), -contains(extracted.Grade),-contains(unextracted.Grade),
             contains("APValue"), -contains(extracted.P.No), -contains(unextracted.P.No))  %>%
      transmute(A = rowSums(. == "A"), B = rowSums(. == "B"),
                C = rowSums(. == "C"), D = rowSums(. == "D"),
                No.grade = rowSums(.=="-"), APvalue.001 = rowSums(. <= "0.001"))
    
  }else if(internals=="n"){
    lipidCount <- lipidomics %>%
      select(contains("Grade"), contains("APValue"))  %>%
      transmute(A = rowSums(. == "A"), B = rowSums(. == "B"),
                C = rowSums(. == "C"), D = rowSums(. == "D"),
                No.grade = rowSums(.=="-"), APvalue.001 = rowSums(. <= "0.001"))
  }
  else{
    print("You are typing wrong information, please type again") 
    internal <- readline("Please type Y/N: ") %>% str_to_lower(.)
    return(InternalCheck(internal))
  }
  return(lipidCount)
}


########################################################################
# function name: PropPlot
# parameter: data (filtered_lipidomics)
# utility: plot the classes proportion summary for filtered data
########################################################################
PropPlot <- function(data){
  prop.plot <- ggplot(data, aes_string(x=reorder(data$Class, data$Class, function(x)length(x)))) +
    geom_bar(aes(y = (..count..)/sum(..count..)), stat="count")+
    scale_y_continuous(labels=scales::percent) +
    theme_bw() +
    ylab("relative frequencies") +
    ggtitle("Class Summary plot: \n Proportion of different class numbers among all the samples")+
    xlab("classes")+ ylab("Relative frequencies (%)")  +
    coord_flip()
  print(prop.plot)
}


########################################################################
# function name: IndividuleRetentionPlot
# parameter: lipid.data, sample 
# utility: This function is unused. If you want to use it, please uncomment the 
#          correponding part in the main script. it will 
#           plot retention time for individule sample and save the plot 
########################################################################
IndividuleRetentionPlot <- function(lipid.data, sample){
  retention.plot <-   ggplot(data=lipid.data,
                             aes_string(x = sprintf("log10(%s)", sample), y = "BaseRt")) +
    geom_point() +
    theme_bw() +
    facet_grid(.~Class) +
    labs("Abundance VS. Retention time",
         x="lipid class (log10(MainArea))", y="Retention time")
  print(retention.plot)
  message("Please input the name you want to store for the graph. e.g. retention.pdf")
  plot_name <- readline("QC plot name: ")
  ggsave(filename = plot_name, path = 'plot/', device = "pdf", dpi=300)
}


########################################################################
# function name: AllRetentionPlot
# parameter: new.data (retention_data)
# utility: plot the retention time for all sample abundances and save the plot as 
#           all.retention.png in the plot directory
########################################################################
AllRetentionPlot <- function(new.data){
  retention.all.plot <- ggplot(data=new.data, aes(x=log10(MainArea), y = BaseRt)) +
    geom_point() +
    theme_bw() +
    facet_grid(.~Class) +
    labs("Abundance VS. Retention time",
         x="All samples (log10(MainArea))", y="Retention time (mins)")
  print(retention.all.plot)
  ggsave(filename = "all.retention.png", path = 'plot/', device = "png", dpi=300)
  #ggsave(filename = "all.retention.png", path = 'plot/', device = "png", width = 10, height = 8, dpi = 150, units = "in")
}


########################################################################
# function name: GrepIndex
# parameter: sample, data
# utility: get sample names and its column position information from filtered data csv file
########################################################################
GrepIndex <- function(sample, data){
  sample.names  <- list()
  group_names   <- c()
  col.index     <- list()
  for( i in 1:length(names(sample))){
    group_names[i]                  <- names(sample)[i]
    group_sample_names              <- sample[[i]] %>%
      strsplit(., "\\s+") %>%
      unlist() %>%
      str_to_lower()
    sample.names[[group_names[i]]]  <- paste("MainArea[", group_sample_names, "]", sep="")
    col.index[[group_names[i]]]     <- which(colnames(data) %in% sample.names[[i]])
  }
  info      <- rbind(sample.names, col.index)
  return(info)
}


########################################################################
# function name: Input
# parameter: data (filtered_lipidomics)
# utility: input group information of samples from the console and reterieve 
#           corresponding sample information from the csv file
########################################################################
Input <- function(data){
  # info about group size
  ngroups <- readline("Group number: ") %>% as.numeric(.)
  ngroups <- CheckType(ngroups)
  
  # Type the group info and store it into a variable
  total_info <- InputGroups(ngroups)
  
  # store the sample and index in a list
  sample_info <- GrepIndex(total_info, data)
  message("Take a look at the sample info and its column position information in the file below")
  glimpse(sample_info)
  
  # retrieve the group names
  group_names <- dimnames(sample_info)[[2]]
  
  return(list(sample_info, group_names, ngroups))
  
}


########################################################################
# function name: CheckType
# parameter: x (ngroups from Input function/group number)
# utility: check if the group number is numeric and store the correct form
#         of group number information
########################################################################
CheckType <- function(x){
  if((x%%1 != 0) | (is.na(x))){
    print("Group number must be numeric!")
    x <- readline("Group number: ") %>% as.numeric(.)
    return(CheckType(x))
  }else{
    if(x<=0){
      print("You input wrong group number, please retype")
      x <- readline("Group number: ") %>% as.numeric(.)
      return(CheckType(x))
    }else{
      return(x)
    }
  }
}


########################################################################
# function name: InputGroups
# parameter: n (ngroups from input function/group number)
# utility: read standard input from console and check if it's correct and
#           store group info into list
########################################################################
InputGroups <- function(n){
  group_info  <- c()
  sample_info <- c()
  info <- list()
  for(i in 1:n){
    ms1                   <- paste("Name of Group ", i, ": ")
    group_info[i]         <- readline(prompt= ms1)
    ms2                   <- paste("Sample names of Group ", i, ": ")
    sample_info[i]        <- readline(prompt=ms2)
    info[[group_info[i]]] <- sample_info[i]
  }
  message("Please take a look at the group information below")
  glimpse(info)
  message("Do you need retype the group info? ")
  group_condition <- readline("Y/N: ")
  if(str_to_lower(group_condition) !="n"){
    if(str_to_lower(group_condition) =="y"){
      return(InputGroups(n))
    }else{
      print("You typed wrong condition info, and pipeline will assume you want to type again.")
      return(InputGroups(n))
    }
  }
  return(info)
}


########################################################################
# function name: PCA_pairs_Plot
# parameter: sample_info, group_names, filtered_lipidomics, mark
# utility: this function will plot PCA and pair-wise correlations among the samples.
#           This will also produce two kinds of correlation plots based on different
#           axis display style. 
#           Please notice that this function is really long format. 
########################################################################

PCA_pairs_Plot <- function(sample_info, group_names, filtered_lipidomics, mark){
  
  # Preparation for pair-wise correlations
  inf2NA      <- function(x) { x[is.infinite(x)] <-  NA; x }
  
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
  index         <- 1:length(sample_info)
  sample_index  <- index[index %% 2 ==0]
  
  ######################################################################
  # QC PLOT 2 - Pair-wise correlation between replicates, reformatted axis
  for(i in 1:length(sample_index)){
    range     <- sample_index[i]
    plot_name <- paste("pairs.plot.", i, ".", mark, ".pdf",sep="")
    path      <- file.path("plot/", plot_name)
    pairs2(log10(filtered_lipidomics[, sample_info[[range]]]), 
           lower.panel = panel.smooth, diag.panel=panel.hist, 
           upper.panel = panel.cor, xlim=c(4,12), ylim=c(4,12))
    dev.copy(pdf, path)
    dev.off()
  }
  
  #################################################################
  
  
  # Storing pairwise plot in the plot directory
  for(i in 1:length(sample_index)){
    range     <- sample_index[i]
    plot_name <- paste("pairs.plot.", i, ".", mark, ".oldversion.pdf",sep="")
    path      <- file.path("plot/", plot_name)
    pairs(log10(filtered_lipidomics[, sample_info[[range]]]), 
          lower.panel = panel.smooth, diag.panel=panel.hist, 
          upper.panel = panel.cor, xlim=c(4,12), ylim=c(4,12))
    dev.copy(pdf, path)
    dev.off()
  }
  
  # making group repeats according to its position for making groups of later PCA
  index_list  <- index[!index %% 2 ==0]
  sample_list <- c()
  for(i in length(index_list):1){
    j           <- index_list[i]
    sample_list <- c(sample_info[[j]], sample_list)
  }
  
  group_repeats <- c()
  for(i in length(group_names):1){
    k             <- index_list[i]
    len           <- length(sample_info[[k]])
    repeats       <- rep(group_names[i], len)
    group_repeats <- c(repeats, group_repeats)
  }
  
  
  # Formatting the table for PCA
  filtered_lipidomics_PCA <-  filtered_lipidomics %>% 
    select(sample_list) %>% 
    t() %>% 
    as.data.frame()
  
  colnames(filtered_lipidomics_PCA) <- filtered_lipidomics$LipidMolec    
  
  # if just show the sample names by removing the MainArea prefix
  rownames(filtered_lipidomics_PCA) <- filtered_lipidomics_PCA %>% 
    rownames(.) %>% 
    str_remove_all(., "MainArea\\[") %>% str_remove_all(., "\\]")
  
  log2_filtered_lipidomics_PCA      <- log((filtered_lipidomics_PCA+1), 2)
  
  # making group info for PCA 
  filtered_lipidomics_PCA$Group       <- group_repeats
  log2_filtered_lipidomics_PCA$Group  <- group_repeats
  
  # Perform PCA 
  res.pca         <-  PCA(filtered_lipidomics_PCA, scale.unit=TRUE, 
                          quali.sup=ncol(filtered_lipidomics_PCA), graph=T)
  
  #print(get_eigenvalue(res.pca))
  
  fviz_screeplot(res.pca, ncp=10)
  
  concat          <-  cbind.data.frame(filtered_lipidomics_PCA[, ncol(filtered_lipidomics_PCA)], 
                                       res.pca$ind$coord)
  ellipse.coord   <-  coord.ellipse(concat, bary=TRUE)
  
  
  
  #pdf(file="plot/sample.pca.pdf")
  plot.PCA(res.pca, habillage=ncol(filtered_lipidomics_PCA), 
           ellipse=ellipse.coord, cex=0.8, label="all", ncp=5)
  dev.copy(pdf, "plot/sample.pca.pdf")
  dev.off()
  
  
  plot.PCA(res.pca, habillage=ncol(filtered_lipidomics_PCA), 
           ellipse=ellipse.coord, cex=0.8, label="quali")
  dev.copy(pdf, "plot/group.pca.pdf")
  dev.off()
  
  # # Output the filtered lipidome
  # write.csv(filtered_lipidomics, "data/filtered_lipidomics.csv")
  ###################################################################################### 
  return(list(sample_list, group_repeats))
}


########################################################################
# function name: PCAcheck
# parameter: pca_check
# utility: check if deleting samples needed and replot PCA
########################################################################
PCAcheck <- function(pca_check, data){
  if(str_to_lower(pca_check)!="n"){
    if(str_to_lower(pca_check)=="y"){
      samples <- Input(data)
      sample_info <- samples[[1]]
      group_names <- samples[[2]]
      ngroups <- samples[[3]]
      label <- "new"
      info_list <-  PCA_pairs_Plot(sample_info, group_names, data, label)
    }else{
      print("You typed wrong, please retype.")
      pca_check <- readline("Please type Y/N: ")
      return(PCAcheck(pca_check, data))
    }
  }
  return(info_list)
}  


#########################################################################
# function name: SD
# parameter: dt
# utility: calculate sd for each class in each group
########################################################################
SD <- function(dt){
  data.sd <- aggregate(.~experiment.group, data=dt, function(x) sd(x)) %>%
    gather(class, sd, -experiment.group)
  return(data.sd)
}


########################################################################
# function name: WideToLong
# parameter: wide_data, sd
# utility: transform wide data to long data format
########################################################################
WideToLong <- function(wide.data, sd){
  data.long <- aggregate(.~experiment.group, data=wide.data, FUN=sum) %>%
    gather(class, value, -experiment.group) %>%
    merge.data.frame(., sd) %>%
    arrange(class)
  return(data.long)
}


########################################################################
# function name: Subtra
# parameter: x
# utility: subtracted background (blank sample c) from all the samples
########################################################################
Subtra <- function(x, na.rm=FALSE){x-filtered_lipidomics$`MainArea[c]`}


########################################################################
# function name: ReplaceInf
# parameter: x
# utility: replace all the zero values to negative values, such that the 
#         log transformed data of zero value will produce NaN other than -Inf
########################################################################
#ReplaceInf <- function(x, na.rm=FALSE){x=replace(x, x == 0, -1)}
ReplaceInf <- function(x, na.rm=FALSE){x=replace(x, x == -Inf, NaN)}


########################################################################
# function name: log2trans
# parameter: x
# utility: take log transformation of the data
########################################################################
log2trans <-  function(x, na.rm=FALSE){log2(x)}


########################################################################
# function name: ClassPlot
# parameter: data
# utility: make bar plots
########################################################################
ClassPlot <- function(data){
  classplot <- ggplot(data, aes(x=reorder(class, value), fill=experiment.group, y=value)) +
    geom_bar(position="dodge", stat="identity") +
    geom_errorbar(aes(ymin=value, ymax=value+sd), width=.3,
                  position=position_dodge(.9)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle=45, hjust=1)) + 
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    )
  
  return(classplot)
  #+ coord_cartesian(ylim=c(1e3, 1e14)) ## these y axis range parameter can be changed
}


########################################################################
# function name: FoldPlot
# parameter: data
# utility: make fold change bar plots 
########################################################################
FoldPlot <- function(data){
  ggplot(data, aes(x=class, fill=experiment.group, y=value)) +
    geom_bar(position="dodge", stat="identity") +
    geom_hline(yintercept=1, linetype="dashed", colour="red", size=0.2) +
    theme_bw() +
    coord_flip() +
    #    geom_hline(yintercept=1, linetype="dashed", colour="grey45", size=0.02) +
    theme(plot.title = element_text(hjust = 0.5))
}



########################################################################
# function name: TotalClass
# parameter: data, sample, group_info
# utility: get aggregated main area data and standard deviation for samples. 
########################################################################
# For total classes
TotalClass <- function(data, sample, group_info){
  # select data contains sample info of the groups
  filtered_groups <- data %>% select(Class, sample)
  # reform the data by aggregating data contains same class name
  filtered_groups <- aggregate(.~Class, data=filtered_groups, FUN=sum)
  # build new data frame contains the group info
  total_class_data <- filtered_groups[,-1] %>% t %>% data.frame(group_info, .)
  names(total_class_data) <- c("experiment.group",filtered_groups$Class)
  # calculate the standard deviation for each class in each group
  total_class_sd <- SD(total_class_data)
  
  # # get median sample main area value for each group 
  # total_class_median <- total_class_data %>% group_by(experiment.group) %>% summarise_all(.funs=median)
  
  # aggregating data by group
  total_class_wide <- total_class_data  %>% group_by(experiment.group) %>% summarise_all(.funs=sum)
  # reformat data into long data
  total_class_long <- total_class_wide %>% gather(class, value, -experiment.group)
  
  # combine data and its sd info for ggplot
  total_class_long <- bind_cols(total_class_long, sd=total_class_sd$sd)
  #total_class_long <- WideToLong(total_class_wide, total_class_sd)
  return(list(total_class_long, total_class_wide))  
}


########################################################################
# function name: PassToNormal
# parameter: data, n_control
# utility: pick several control groups and plot fold changes bar plots among all the groups 
#           The function will produce two kind of display of plots. One will be particularly
#           order the class axis based on its value in control group
########################################################################
PassToNormal <- function(data, n_control){
  total_log_dt<- TotalNormalClass(data, n_control)
  total_log_data <- total_log_dt[[1]]
  order_total <- total_log_dt[[2]]
  for(i in 1:n_control){ 
    total_log <- total_log_data[[i]]
    total_log2_plot <- FoldPlot(total_log) +
      labs(title="Total lipid classes", x="Lipid class (background subtracted)", y="Log2 fold change of imputated data")
    print(total_log2_plot)
    filenames <- paste("class.log2.", i, ".png", sep = "")
    ggsave(filename=filenames, path = 'plot/', device = "png", dpi = 300)
    
    # log2 fold changes for each class among different groups AND ordered the class axis by the control group
    ordered_data<- order_total[[i]]
    ordered_fc_plot <- FoldPlot(ordered_data) +
      labs(title="Total lipid classes", x="Lipid class (background subtracted)", y="Log2 fold change of imputated data",
           caption = "The axis of classes are ordered decreasingly by control group")
    
    print(ordered_fc_plot)
    filenames <- paste("orderedControl.log2.", i, ".png", sep = "")
    ggsave(filename=filenames, path = 'plot/', device = "png", dpi = 300)
  }
}



getMedian <- function(data){
  sample_median <- data %>% group_by(experiment.group) %>% summarise_all(median) 
  }
getMean <- function(data) {
  sample_mean <- data %>% group_by(experiment.group) %>% summarise_all(mean) 
  }


########################################################################
# function name: TotalNormalClass
# parameter: data, n_control
# utility: calculate the fold changes among median samples in groups
########################################################################
TotalNormalClass <- function(data, n_control){
  fc_long <- list()
  order_fc_long <- list()
  for(i in 1:n_control){
    #control_group <- c()
    print("Yes")
    # design control groups:
    info <- paste("Please input name of control Group", i," e.g. A: ")
    control_group <- readline(prompt = info)
    contrast_groups <- readline("Please input the groups which are split by space for contrast, e.g. B C D: ")
    contrast_groups <- contrast_groups %>% str_split(., "\\s+") %>% unlist()
    print(contrast_groups)
    # get the median sample for each group
    median_sample <- getMedian(data)     ##### pass get Median function
    control_index <- which(median_sample$experiment.group %in% control_group)
    contrast_index <- which(median_sample$experiment.group %in% contrast_groups)
    # calculate the fold changes among the groups
    subtraction <- median_sample %>% select(-experiment.group) %>% slice(control_index) %>% unlist(.)
    fold_samples <- median_sample %>% slice(control_index, contrast_index)
    fold_samples[, -1] <- sweep(fold_samples[, -1], 2, subtraction, "-" )
    fold_samples[1, -1] <- apply(fold_samples[1, -1], 1, function(x)ifelse(x==0, 1, 1) )
    
    # order the class axis by control group of decreasing order in log transformation data
    order_class <- median_sample %>%
      filter(experiment.group==control_group) %>%
      select(-experiment.group) %>% apply(., 1, function(x) x <- order(x))  
    order_factor <- median_sample %>% select(1+ order_class) %>% colnames(.)
    # reformat the data form into long data
    fc_long[[i]] <- fold_samples %>% gather(class, value, -experiment.group)
    order_fc_long[[i]] <- fc_long[[i]] %>% mutate(class= factor(class, levels = order_factor))  %>% arrange(class) 
  }
  return(list(fc_long, order_fc_long))
}


########################################################################
# function name: EachClass
# parameter: data
# utility: plot for each class
########################################################################
EachClass <- function(data){
  classes <- unique(data$Class)
  class_long_list <- list()
  class_wide_list <- list()
  for(i in 1:length(classes)){
    options(warn=-1)
    pick_class <- classes[i]
    print(pick_class)
    # get information of the samples
    filtered_class <- data %>%
      select(Class, LipidMolec, sample_list) %>%
      filter(Class==pick_class)
    # aggregate the repeated data by its class and lipid molecular names
    class_data <- aggregate(.~Class+LipidMolec, data=filtered_class, FUN=sum) %>%  select(-Class)
    
    # rename the class
    molec_names <- class_data[, 1] %>%
      unlist() %>%
      str_extract_all(., "\\(\\d*.*\\)") %>%
      str_remove_all(., "[\\(\\)]")
    
    mark <- paste("(", pick_class, ")", molec_names, sep="")
    # Making the group information for the data
    data_wide <- data.frame(group_repeats, t(class_data[, -1])) 
    # making names for the new data frame
    names(data_wide) <- c("experiment.group", molec_names)
    class_wide_list[[i]] <- data_wide
    names(class_wide_list)[i] <- pick_class
    colnames(class_wide_list[[i]]) <- c("group", mark)
    # calculate the sd information
    class_sd <- SD(data_wide)
    
    # reformat the data information
    class_long <- WideToLong(data_wide, class_sd)
    # store the data
    class_long_list[[i]] <- class_long
    names(class_long_list)[i] <- pick_class
  }
  return(list(class_long_list, class_wide_list))
}


########################################################################
# function name: EachClassPlot
# parameter: data, n_control
# utility: calculate the fold changes among median samples in groups
#         if you are a mac user, please install quartz via 
########################################################################
EachClassPlot <- function(n_bar, n_groups, long_data){
  quartz()
  classes <- names(long_data)
  for(i in 1:length(classes)){
    options(warn=-1)
    pick_class <- classes[i]
    print(pick_class)
    observations <- nrow(long_data[[i]])
    if(observations <= n_bar){
      p1 <- ClassPlot(long_data[[i]]) 
      p1 <- p1 +labs(title= pick_class, x="Acyl composition", y="Main Area", caption="Error bar is the standard deviation for each class in each group") 
      print(p1)
      #ggsave(filename = paste(pick_class, ".png", sep=""), path = 'plot/classes', device = "png", width=15, height=15, dpi=300)
      ggsave(filename = paste(pick_class, ".png", sep=""), path = 'plot/classes', device = "png")
      
    }else{
      if(n_bar %% n_groups != 0){
        n_bar <- (n_bar%/% n_groups)*n_groups
      }
      nfacet <- observations %/% n_bar
      for(k in 0:nfacet){
        if(k < nfacet){
          ranges <- (n_bar*k+1):(n_bar*(k+1))
          data <- long_data[[i]] %>% slice(ranges)
          p2 <- ClassPlot(data)
          p2 <- p2 + labs(title= pick_class, x="Acyl composition", y="Main Area", caption="Error bar is the standard deviation for each class in each group") 
          print(p2)
          ggsave(filename = paste(pick_class, ".", k+1, ".png", sep=""), path = 'plot/classes/', device = "png")
        }else{
          ranges <- (n_bar*k+1):nrow(long_data[[i]])
          data <- long_data[[i]] %>% slice(ranges)
          p2 <- ClassPlot(data)
          p2 <- p2 + labs(title= pick_class, x="Acyl composition", y="Main Area", caption="Error bar is the standard deviation for each class in each group") 
          print(p2)
          ggsave(filename = paste(pick_class, ".", k+1, ".png", sep=""), path = 'plot/classes/', device = "png")
        }
      }
    } 
  }
  dev.off()
}


########################################################################
# function name: ImputeMinProb
# parameter: data, q, tune.sigma
# utility: this function applys the imputation of left-censored missing data 
#           by random draws from a Gaussian distribution centered in a minimal 
#           value. The minimal value observed is estimated as being the q-th 
#           quantile (e.g. q = 0.01) of the observed values in that sample. 
#           q is a scalar used to determine a low expression value to be used 
#           for missing data imputation. 0 < q < 1, in this case q should be 
#           set to a low value. The default value is q = 0.01. tune.sigma is 
#           a scalar used to control the standard deviation of the Gaussian 
#           distribution used for random draws. If the sd is overestimated, 
#           than 0 < sigma.coef < 1. The default value is tune.sigma = 1.
#           The function is obtained from imputeLCMD package
########################################################################
ImputeMinProb <- function (data, q = 0.01, tune.sigma = 1) 
{
  n_samples <-  dim(data)[2]
  n_observations <-  dim(data)[1]
  imputated_data <-  data
  min_samples <-  apply(imputated_data, 2, quantile, prob = q, na.rm = T)
  count_NAs <-  apply(!is.na(data), 1, sum)
  count_NAs <-  count_NAs/n_samples
  filtered_data <-  data[which(count_NAs > 0.5), ]
  data_sd <-  apply(filtered_data, 1, sd)
  sd_temp <-  median(data_sd, na.rm = T) * tune.sigma
  print(sd_temp)
  for (i in 1:(n_samples)) {
    input <-  rnorm(n_observations, 
                    mean = min_samples[i], 
                    sd = sd_temp)
    imputated_data[which(is.na(data[, i])), i] <- input[which(is.na(data[,i]))]
  }
  return(imputated_data)
}


########################################################################
# function name: BuildContrast
# parameter: fit, filtered_data
# utility: Build contrast matrix for limma 
########################################################################
BuildContrast <- function(fit, filtered_data){
  contmatrix <- makeContrasts(readline("Please input the groups names for comparison, '-' means the contrast, e.g A-B: "),
                               levels = design
  )
  fit2 <- contrasts.fit(fit, contmatrix)
  fit2 <- eBayes(fit2)
  comparison <- toptable_OUTPUT1 <- topTable(fit2, coef=1, adjust.method = 'fdr',
                                             lfc=0, number=nrow(filtered_data)) %>% as.data.frame()
  # store the result into csv
  name1 <- readline("Please input the name for comparison data file, e.g. A.B.csv: ")
  write_csv(comparison, file.path("data/", name1))
  # store significant result into csv
  toptable_OUTPUT1.1 <- topTable(fit2, coef=1, adjust.method = 'fdr',
                                 p.value =0.05, lfc=log2(1.5), number=nrow(filtered_data)) %>% as.data.frame()
  name1.1 <- readline("Please input the name for significant data (p value less than 0.05) file, e.g. A.B.sig.csv: ")
  write_csv(x=toptable_OUTPUT1.1, file.path("data", name1.1))
  # volcano plot
  fold_change <- as.numeric(readline("Please input the fold change treshold for the volcano plot, e.g 1: "))
  # volcano input info.
  input <- comparison  # input of comparison group info which can be changed
  # significant data
  input$sig <- factor(input$adj.P.Val< 0.05 & abs(input$logFC) > fold_change)
  # size of significant data
  significant_lipids <- sum(input$adj.P.Val< 0.05 & abs(input$logFC) > fold_change)
  message("When the tests' q value treshold is 0.05 and the fold change threshold is ",fold_change,". The number of lipids which are statistically significant are: ", significant_lipids)
  
  # Define significant data for volcano plot graphing
  points <- input %>%
    rownames_to_column('lipid') %>%
    filter(abs(logFC) > fold_change & adj.P.Val< 0.05 ) %>%
    column_to_rownames('lipid')
  # making volcano plot
  volc1 <- PlotVolc(input, points, fold_change) 
  volc1 <- volc1 + geom_point(aes(col=sig, size=AveExpr)) + 
                    scale_size_continuous(range = c(0.25,2.5), breaks=c(10, 20, 30))+
                    scale_color_manual(values=c("#bdbdbd", "#3182bd"))
  print(volc1)
  plot_name <- readline("Please input the volcano plot name: ")
  ggsave(filename = plot_name, path = 'plot/', device = "pdf", width=15, height=15, dpi=300)
  # build mutated data frame
  class_names <- rownames(input) %>% str_extract_all(., "(.+)\\(") %>% str_remove_all(., "\\(")
  #Add a new column to specify lipids
  input_lipids <- input
  input_lipids$lipclass <- class_names
  input_lipids$lipclass <- ifelse(input_lipids$sig=="TRUE", 
                               case_when( (input_lipids$lipclass %in% c("PC", "PE", "PG", "PI", "PS", "LPC", "LPE",
                                                                     "LPI", "dMePE", "CL"))~ "Glycerophospholipids",
                                          input_lipids$lipclass %in% c("TG", "DG", "MG") ~ "Neutral lipids",
                                          input_lipids$lipclass %in% c("SM", "So", "SoG1", "Cer", "CerG1","CerG2",
                                                                    "GM2","GM3", "CerG2GNAc1", "CerG3GNAc1", "CerG3NAc2") ~ "Sphingolipids",
                                          input_lipids$lipclass %in% c("ChE","Cholestoral") ~ "Sterols",
                                          TRUE ~ "Other lipids"), 
                               "n.s")
  # making levels for the class category
  class_levels <- c("Glycerophospholipids", "Neutral lipids", "Sphingolipids", "Sterols", "Other lipids", "n.s")
  # make the class category as factors
  input_lipids$lipclass <- factor(input_lipids$lipclass, levels=class_levels)
  # rearrange the data by its levels 
  input_lipids <-  input_lipids %>% 
                  rownames_to_column('lipid') %>% 
                  arrange(lipclass) %>% 
                  column_to_rownames('lipid')
  # subset the significant points for text information
  significant_points <- input_lipids %>% 
                        rownames_to_column('lipid') %>% 
                        filter(!lipclass=="n.s") %>% 
                        column_to_rownames('lipid')
  # subset the foldchange larger than 5 for text information                    
  fc_points <- input %>% 
    rownames_to_column('lipid') %>% 
    filter(abs(logFC) > 5 & adj.P.Val>= 0.05 ) %>% 
    column_to_rownames('lipid')
  # plot the color version 
  volc2 <- PlotVolc(input_lipids, significant_points, fold_change)
  volc2 <- volc2 + geom_text_repel(data=fc_points, aes(label=rownames(fc_points)), size=2) + 
    geom_point(aes(col=lipclass, size=AveExpr))  +
    scale_size_continuous(range = c(0.25,2.5), breaks=c(10, 20, 30)) +
    scale_color_manual(values=c( "darkorange3","dodgerblue3", 'chartreuse4', 
                                 '#c51b8a', '#756bb1',"#bdbdbd"), drop=FALSE)
  print(volc2)
  plot_name <- readline("Please input the volcano plot name: ")
  ggsave(filename = plot_name, path = 'plot/', device = "pdf", width=15, height=15, dpi=300)
}


########################################################################
# function name: PlotVolc
# parameter: input, points, fold_change
# utility: plot volcano graph
########################################################################
PlotVolc <- function(input, points, fold_change){
  plot <- ggplot(input, aes(logFC, -log10(adj.P.Val))) + 
    geom_rect(aes(xmin = -fold_change, xmax = fold_change, ymin = -Inf, ymax = Inf),
              fill = "#f0f0f0", linetype="dashed", color="red", size=0.2)+
    geom_point() +
    geom_vline(xintercept=0, linetype="solid", colour="grey60", size=0.2) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="red", size=0.2)+
    xlab(bquote("Fold change, " ~ log[2]~"("~textstyle(frac({x}, {y}))~scriptstyle(, AUC)~")")) + 
    ylab("-Log10(q value)") +
    ggtitle("volcano plot") +
    theme_bw()+ 
    theme(line=element_blank()) +
    geom_text_repel(data=points, aes(label=rownames(points)), size=2) 
  return(plot)
}


########################################################################
# function name: VolPlot
# parameter: n_comparisons, fit, filtered_data
# utility: passing parameters to BuildContrast function to plot volcano graph
#         for number of n_comparisons
########################################################################
VolPlot <- function(n_comparisons, fit, filtered_data){
  replicate( n_comparisons, BuildContrast(fit, filtered_data))
}


########################################################################
# function name: checkDownload
# parameter: x
# utility: Check if Cytoscape is download and open
########################################################################
checkDownload <- function(x){
  if(x != "y"){
    if(x == "n"){
      message("Please download Cytoscape from the link above and open it")
      x <- readline("Do you download Cytoscape 3.6.1 or higher and opened it? Y/N: ") %>% str_to_lower(.)
      return(checkDownload(x))
    }else{
      message("Please type again.")
      x <- readline("Do you download Cytoscape 3.6.1 or higher and opened it? Y/N: ") %>% str_to_lower(.)
      return(checkDownload(x))
    }
  }
}



########################################################################
# function name: buildPathWay
# parameter: none
# utility: build initial grapha
########################################################################
buildPathWay <- function(){
    g = new ('graphNEL', edgemode='directed')
    g = graph::addNode ('fatty acids', g)
    g = graph::addNode ('G3P', g)
    g = graph::addNode ('LCBs', g, edges = list('fatty acids'))
    g = graph::addNode ('Cer', g, edges = list("LCBs"))
    g = graph::addNode ('SM', g, edges = list("Cer"))
    g = graph::addNode ('LPA', g, edges = list(c("G3P", "fatty acids")))
    g = graph::addNode ('PA', g, edges = list('LPA'))
    g = graph::addNode ('DAG', g, edges = list('PA'))
    g = graph::addNode ('TAG', g, edges = list('DAG'))
    g = graph::addNode ('PC', g, edges = list('DAG'))
    g = graph::addNode ('LPC', g, edges = list('PC'))
  #  g = graph::addEdge ('LPC', 'PC', g)
    g = graph::addNode ('PE', g, edges = list('DAG'))
    g = graph::addNode ('LPE', g, edges = list('PE'))
 #   g = graph::addEdge ('LPE', 'PE', g)
  #  g = graph::addEdge ('PE', 'PC', g)
    g = graph::addNode ('CDP-DAG', g, edges = list('PA'))
    g = graph::addNode ('PI', g, edges = list('CDP-DAG'))
    g = graph::addNode ('LPI', g, edges = list("PI"))
  #  g = graph::addEdge ('LPI', 'PI', g)
    g = graph::addNode ('PG', g, edges = list('CDP-DAG'))
    g = graph::addNode ('LPG', g, edges = list('PG'))
   # g = graph::addEdge ('LPG', 'PG', g)
    g = graph::addNode ('PS', g, edges = list(c('CDP-DAG', 'PE')))
    g = graph::addNode ('LPS', g, edges = list('PS'))
    g <- createNetworkFromGraph (g, title='simple network', collection='GraphNEL Example')
    addCyEdges(list("LPC", "PC"), edgeType = "interacts with", directed = FALSE, network = g)
    addCyEdges(list("LPE", "PE"), edgeType = "interacts with", directed = FALSE, network = g)
    addCyEdges(list('PE', 'PC'), edgeType = "interacts with", directed = FALSE, network = g)
    addCyEdges(list('LPI', 'PI'), edgeType = "interacts with", directed = FALSE, network = g)
    addCyEdges(list('LPG', 'PG'), edgeType = "interacts with", directed = FALSE, network = g)
    addCyEdges(list('LPS', 'PS'), edgeType = "interacts with", directed = FALSE, network = g)
    
    #g = graph::addEdge ('LPS', 'PS', g)
    return(g)
}

