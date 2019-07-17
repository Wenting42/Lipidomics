####################################################################################
# Script: fattyAcidsSaturation_analysis.FUNCTIONS.2.1.r
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
#                 source("fattyAcidsSaturation_analysis2.1.r")   
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
# function name: CalculateMainArea
# parameter: ...
# utility: calculate the main area for SFA, MUFA, PUFA in each group
########################################################################
CalculateMainArea <- function(...){
  vars <- ensyms(...)
  group_names <- c()
  for( i in seq_along(vars)){
    res <- vars[[i]]
    sfa <- expr(SF * (!!res))
    mufa <- expr(MUF * (!!res))
    pufa <- expr(PUF * (!!res))
    percent_info <- fa_percent %>% transmute(!!sfa, !!mufa, !!pufa)
    
    each_fa <- cbind(each_fa, percent_info)
    sfs_name <- addquotes(!!res, "_SFA")
    mufa_name <- addquotes(!!res, "_MUFA")
    pufa_name <- addquotes(!!res, "_PUFA")
    groups <- c(sfs_name, mufa_name, pufa_name)
    group_names <- c(group_names, groups)
    
  }
  colnames(each_fa) <- c(group_names)
  each_fa <- bind_cols(each_fa, fa_percent)
  return(each_fa)
}

