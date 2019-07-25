####################################################################################
# Script: fattyAcidsSaturation_analysis.FUNCTIONS.2.2.r
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
#                 source("fattyAcidsSaturation_analysis2.2.r") #
#  fixed some calculation bugs
#####################################################################################

# load packages needed
initial.package <- "BiocManager"
need.package <- initial.package[!(initial.package %in% installed.packages()[, "Package"])]
if(length(need.package) > 0) install.packages(need.package)


list.of.packages <- c( "FactoMineR", "factoextra", "scales", "magrittr", "ggrepel", "reshape2",
                       "stargazer", "Hmisc", "limma", "factoextra", "scales", "RColorBrewer",
                       "stringr", "readxl", "RCy3", "igraph", "tidyverse", "dplyr", "rlang")
need.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if(length(need.packages) > 0) BiocManager::install(need.packages)


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
  for (i in group_list){
    samples<- filter(groups, groups == i) %>% ungroup() %>% select(samples)
    sample_list <- samples$samples
    nameMedian <- addquotes(!!i, "_median")
    nameMean <- addquotes(!!i, "_mean")
  
    # get median and mean for samples in each group 
   dt1 <-  percent_info %>% 
           rowwise() %>% 
           transmute(!!nameMedian := median(c(!!!syms(sample_list))), 
                     !!nameMean := mean(c(!!!syms(sample_list))))
   percent_info <- cbind(percent_info, dt1)
    
   k1 <- c("%SFA", "%MUFA", "%PUFA")
   exprMedian <- sapply(k1, function(x)expr(!!sym(x) * !!sym(nameMedian)))
   exprMean <- sapply(k1, function(x)expr(!!sym(x) * !!sym(nameMean)))
   
   k2 <- c("_SFA", "_MUFA", "_PUFA")  
   medians <- sapply(k2, function(x)addquotes(!!nameMedian, !!x))
   means <- sapply(k2, function(x)addquotes(!!nameMean, !!x))

   
   # do saturation analysis for each group 
   dt2 <- percent_info %>% 
          rowwise() %>% 
          transmute(!!(medians[1]) := eval(exprMedian[[1]]),
                    !!(medians[2]) := eval(exprMedian[[2]]), 
                    !!(medians[3]) := eval(exprMedian[[3]]),
                    !!(means[1]) := eval(exprMean[[1]]),
                    !!(means[2]) := eval(exprMean[[2]]),
                    !!(means[3]) := eval(exprMean[[3]]))
   percent_info <- cbind(percent_info, dt2)
  }
  return(percent_info)
}



########################################################################
# function name: ReformatData
# parameter: data, pick
# utility: Transforme data by pass mean or median based sample in 
#         each group 
########################################################################
ReformatData <- function(data, pick){
  methods <- c("mean", "median")
  exclude <- methods[! methods %in% pick]  
  info <- data %>% select(Class, FA.name, contains(!!pick), -contains(!!exclude))
  data1 <- info %>% group_by(Class) %>% select(-FA.name, -ends_with(!!pick)) %>% summarise_all(sum)
  pattern <- addquotes(".*", !!pick, "_")
  # reformat the data
  dt <- data1 %>% 
    gather("Pattern", "Value", - Class) %>% 
    mutate(Group = str_remove(Pattern, "_.*")) %>% 
    mutate_at(vars(Pattern), list(~str_remove(Pattern, pattern)))
  return(list(info, dt))
}



########################################################################
# function name: PlotStack
# parameter: data
# utility: plot stack plots of SFA, MUFA, PUFA for each lipid class
########################################################################
PlotStack <- function(data){
  stackplots <- ggplot(data, aes(x=Group, y=Value, fill = Pattern)) +
                geom_bar(stat = "identity") +
                facet_wrap(~Class, scales = "free") +
                scale_fill_brewer(palette = "Set2") +
                theme_bw()+
                theme(axis.text.x = element_text(size = 8),
                      panel.grid.major = element_blank(),
                      #panel.grid.minor = element_blank(),
                      legend.title = element_blank(),
                      #strip.background = element_rect(colour="black", fill="white", 
                       #                               size= 0.5, linetype="solid")
                      strip.background = element_blank(),
                      #axis.line = element_line(colour = "black"),
                      panel.border = element_blank(),)
}



########################################################################
# function name: PlotGroups
# parameter: data
# utility: plot SFA, MUFA and PUFA for each class
########################################################################
PlotGroups <- function(data){
fattyplots <- ggplot(data, aes(x=Class, y=Value, fill=Group)) +
  geom_bar(stat="identity", position="dodge") +
  scale_fill_brewer(palette = "Set2")
}


check_group <- function(groups){
  control <- readline("Please input the control group for normalization: ")
  if(!control %in% groups){
    message("The input is wrong, please type again.")
    control <- check_group(groups)
  } 
  return(control)
}



########################################################################
# function name: PlotGroups
# parameter: data
# utility: plot SFA, MUFA and PUFA for each class
########################################################################
PassToNormal <- function(data){
  groups <- unique(data$Group)
  control <- check_group(groups)
  data
}







