####################################################################################
# Script: fattyAcidsSaturation_analysis.FUNCTIONS.2.2.3.r
# Author: Wenting
# Notes: This script is the assistant script of fattyAcidsSaturation_analysis2.1.r
#         which helps calculate the SFA, MUFA, PUFA and the aggregated main area for each sample.
#       First, Make sure that your R and Rstudio are newly enough for installing the 
#             packages needed in the script. Otherwise the script will pop out warnings 
#             for packages and won't run.
#       Second, please set the working environment first. Users can copy the project and
#               put it in the directory with the scripts. After that please click the
#               project first and then run the script.
#       Third, to run the script, please type the command below in the console:
#                 source("fattyAcidsSaturation_analysis2.2.3.r") #
#####################################################################################

# load packages needed
initial.package <- "BiocManager"
need.package <- initial.package[!(initial.package %in% installed.packages()[, "Package"])]
if(length(need.package) > 0) install.packages(need.package)


list.of.packages <- c( "FactoMineR", "factoextra", "scales", "magrittr", "ggrepel", "reshape2",
                       "stargazer", "Hmisc", "limma", "factoextra", "scales", "RColorBrewer",
                       "stringr", "readxl", "RCy3", "igraph", "tidyverse", "dplyr", "rlang", "ggforce")
need.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if(length(need.packages) > 0) BiocManager::install(need.packages)


# Suppress the package messages
lapply(list.of.packages, function(x) 
  suppressMessages(require(x, character.only = TRUE, quietly = TRUE)))


# color which are kind for color blinded people, source from Internet
clPalette1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

clPalette2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


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
# function name: getMedians
# utility: find median
########################################################################
getMedians <- function(x){
  median <- median(x)
  return(median)
}




#########################################################################
# function name: CountPatternByRow
# parameter: x (a list or vector of lipid molecules)
# utility: find the saturation pattern and counts
#           SFA pattern is defined as ":0", MUFA is defined as ":1" 
#           and PUFA is defined as ":[2-9]" (colon followed by any number from 2 to 9)
########################################################################
CountPatternByRow <- function(x){
  pattern_name <- i <- j <- k <- c()
  i <- str_count(x, ":0")
  j <- str_count(x, ":1")
  k <- str_count(x, ":[2-9]")
  SFA <- paste(rep("SFA", i), collapse = "/")
  MUFA <- paste(rep("MUFA", j), collapse = "/")
  PUFA <- paste(rep("PUFA", k), collapse = "/")
  pattern_name <- paste( SFA, MUFA, PUFA, sep = "/")
  patterns <- data.frame(pattern = pattern_name, SFA = i, MUFA = j, PUFA = k, 
                         stringsAsFactors = FALSE) %>% as.matrix(.)
  return(patterns)
}



##############################################
# function name: addquotes
# parameter: ...
# utility: add quotes for args
#################################################
addquotes <- function(...){
  args <- ensyms(...)
  paste(purrr::map(args, rlang::as_string), collapse = "")
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
# parameter: data (filtered_lipids)
# utility: input group information of samples from the console and reterieve 
#           corresponding sample information from the csv file
########################################################################
Input <- function(data){
  ###  making group
  message("Please input the info of the experimental groups below and end with 'Enter'
        Please Only input the sample names (replicate) for analysis.
        e.g. Experimental Group number: 2   
        Name of Experimental Group 1: wild_type
        Sample names (replicate number from lipid search) : s1 s2 s3 s4 
        Name of Group 2: knockout
        Sample names of Group 2: s4 s5 s6")
  # info about group size
  ngroups <- readline("Group number: ") %>% as.numeric(.)
  ngroups <- CheckType(ngroups)
  
  # Type the group info and store it into a variable
  total_info <- InputGroups(ngroups)
  
  # store the sample and index in a list
  info <- GrepIndex(total_info, data)
  message("Take a look at the sample info and its column position information in the file below")
  glimpse(info)
  
  # retrieve the group names
  group_names <- dimnames(info)[[2]]
  
  return(list(info, group_names, ngroups))
  
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
  info <- c()
  total_info <- list()
  for(i in 1:n){
    ms1                   <- paste("Name of Group ", i, ": ")
    group_info[i]         <- readline(prompt= ms1)
    ms2                   <- paste("Sample names of Group ", i, ": ")
    info[i]        <- readline(prompt=ms2)
    total_info[[group_info[i]]] <- info[i]
  }
  message("Please take a look at the group information below")
  glimpse(total_info)
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
  return(total_info)
}


########################################################################
# function name: RetrieveInfo
# parameter: info
# utility: retrieve group and its sample information
########################################################################
RetrieveInfo <- function(info){
  # retrieve the sample info position 
  index         <- 1:length(info)
  sample_index  <- index[index %% 2 ==0]
  # making group repeats according to its position for making groups of later PCA
  index_list  <- index[!index %% 2 ==0]
  sample_list <- c()
  for(i in length(index_list):1){
    j           <- index_list[i]
    sample_list <- c(info[[j]], sample_list)
  }
  
  group_repeats <- c()
  for(i in length(group_names):1){
    k             <- index_list[i]
    len           <- length(info[[k]])
    repeats       <- rep(group_names[i], len)
    group_repeats <- c(repeats, group_repeats)
  }
  return(list(group_repeats, sample_list, sample_index))
}



########################################################################
# function name: getSameple
# parameter: fa_percent, group_info
# utility: find mean or median of samples in each group
########################################################################
CalcGroup <- function(percent_info, groups){
  group_list <- unique(groups[, 2]) %>% unlist()
  all_samples <- unique(groups[, 1]) %>% unlist()
  names(all_samples) <- all_samples
  dt <- percent_info %>% 
    mutate_at(all_samples, list(SFA = ~.*`%SFA`, MUFA = ~.*`%MUFA`, PUFA = ~.*`%PUFA`))
  write_csv(dt, "data/individual_saturations.csv")
  dd <- dt %>% 
    select(Class, contains("MainArea"), -ends_with("]")) %>% 
    group_by(Class) %>% 
    summarise_all(sum)
  write_csv(dd, "data/class_saturations.csv")
  for (i in group_list){
    samples<- filter(groups, groups == i) %>% ungroup() %>% select(samples)
    sample_list <- samples$samples
    # bulid names
    nameMedian <- addquotes(!!i, "_median")
    nameMean <- addquotes(!!i, "_mean")
    nameSD <- addquotes(!!i, "_sd")
    k <- c("_SFA", "_MUFA", "_PUFA")
    medians <- sapply(k, function(x)addquotes(!!nameMedian, !!x))
    means <- sapply(k, function(x)addquotes(!!nameMean, !!x))
    sds <- sapply(k, function(x)addquotes(!!nameSD, !!x))
    samples_SFA <- sapply(sample_list, function(x)addquotes(!!x, "_SFA"))
    samples_MUFA <- sapply(sample_list, function(x)addquotes(!!x, "_MUFA"))
    samples_PUFA <- sapply(sample_list, function(x)addquotes(!!x, "_PUFA"))
    # get median, mean and sd for samples in each group 
    dt1 <-  dd %>% rowwise() %>% transmute(!!(medians[1]) := median(c(!!!syms(samples_SFA))),
                                           !!(medians[2]) := median(c(!!!syms(samples_MUFA))),
                                           !!(medians[3]) := median(c(!!!syms(samples_PUFA))),
                                           !!(means[1]) := mean(c(!!!syms(samples_SFA))),
                                           !!(means[2]) := mean(c(!!!syms(samples_MUFA))),
                                           !!(means[3]) := mean(c(!!!syms(samples_PUFA))),
                                           !!(sds[1]) := sd(c(!!!syms(samples_SFA))),
                                           !!(sds[2]) := sd(c(!!!syms(samples_MUFA))),
                                           !!(sds[3]) := sd(c(!!!syms(samples_PUFA))))
    dd <- cbind(dd, dt1)
  }
  return(dd)
}





########################################################################
# function name: FormatData
# parameter: data, pick
# utility: reformat the mean or median based data for plots
########################################################################
FormatData <- function(data, pick){
  dt <- data %>% select(Class, contains(pick))
  pattern1 <- addquotes("_", !!pick, "_.*")
  pattern2 <- addquotes(".*", !!pick, "_")
  dt <- dt %>% gather("TYPE", !!pick, -Class)
  dt$Groups <- str_replace_all(dt$TYPE, pattern1, "")
  dt$TYPE <- str_replace_all(dt$TYPE, pattern2, "")
  return(dt)
}




########################################################################
# function name: PlotTypes
# parameter: data
# utility: plot for SFA, MUFA, PUFA for each lipid class
########################################################################
PlotTypes <- function(data, m){
  m <- syms(m)
  p1 <- m[[1]]
  p2 <- m[[2]]
  p3 <- m[[3]]
  plot_types <-   ggplot(data, aes(x=eval(p1), y=eval(p2), fill = eval(p3))) +
    # facet_wrap(~Class, scales = "free") +
    scale_fill_manual(values = clPalette1) +
    theme_bw()+
    theme(
      axis.text.x = element_text(size = 6, angle = 60, hjust = 1),
      panel.grid.major = element_blank(),
      #panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      #strip.background = element_rect(colour="black", fill="white", 
      #                               size= 0.5, linetype="solid")
      strip.background = element_blank(),
      axis.line = element_line(colour = "black", 
                               size = .2, 
                               arrow = arrow(length = unit(0.06, "cm"), type = "closed")),
      panel.border = element_blank()) 
      # scale_y_log10(
      #   breaks = scales::trans_breaks("log10", function(x) 10^x),
      #   labels = scales::trans_format("log10", scales::math_format(10^.x))
    
  return(plot_types) 
} 




########################################################################
# function name: check_group
# parameter: groups
# utility: check if the user type correct group name
########################################################################
check_group <- function(groups){
  control <- readline("Please input the control group for normalization: ")
  if(!control %in% groups){
    message("The input is wrong, please type again.")
    control <- check_group(groups)
  } 
  return(control)
}



########################################################################
# function name: PassToNormal
# parameter: data
# utility: transform the data to get fold change information based on control group
########################################################################
FoldChange <- function(data){
  groups <- unique(data$Groups)
  control <- check_group(groups)
  n <- length(groups)
  # get fold change 
  control_values <- data %>% filter(Groups == control) %>% select(mean) %>% unlist()
  values <- rep(control_values, n)
  data$mean <- data$mean/values
  # replace NaN with 0
  data$mean <- lapply(data$mean, function(x)ifelse(x == "NaN", 0, x)) %>% unlist()
  return(data)
}







