####################################################################################
# Script: fattyAcids.02.r:w
# Author: Wenting Lyu
# Notes: This script helps calculate the 
#        SFA, MUFA, PUFA and the aggregated main area for each sample.
#       First, Make sure that your R and Rstudio are newly enough for installing the 
#             packages needed in the script. Otherwise the script will pop out warnings 
#             for packages and won't run.
#       Second, click the "Laura.project.Rproj" 
#       Third, to run the script, please type the command below in the console:
#                 source("Laura.fa.2.0.r")   
#              or, click the button "source"
#         
#####################################################################################
rm(list=ls())
# 
# list.of.packages <- c("BiocManager","tidyverse", "readr", "magrittr")
# need.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
# if(length(need.packages)) install.packages(need.packages)
# 
# #######################
# # Required library
# library("tidyverse")
# library("readr")
# library("magrittr")
#######################
#########################################################################################
message("Please run script txt_to_csv.R for converting original file from LipidSearch to get standard file for analysis. \n If you already have file generated by this script, you can continue. ")




source("FWL_lipidomics.2.5.functions.R")
############################################################################################ 
# filter data 
# Read in file, and store the data 
message("Please open the standand file generated from txt_to_csv.R")
lipidomics <- read_csv(file.choose("Please open the standard file"), col_types = cols())

# select and count the sample grades you want, filter the p value which less than 0.001
# this step extract the needed sample information from the file automatically
lipidCount <- lipidomics %>% 
  select(contains("Grade"), -contains("Grade[c]" ),
         contains("APValue"))  %>% 
  transmute(A = rowSums(. == "A"), B = rowSums(. == "B"), 
            C = rowSums(. == "C"), D = rowSums(. == "D"), 
            No.grade = rowSums(.=="-"), APvalue.001 = rowSums(. <= "0.001"))

# add this count into original strucutes
lipidSelect <- lipidomics %>% bind_cols(lipidCount)

# Filter the dataset based on your criteria  
message("Please note the file is filtered by three standards. For same lipidmolecule: \n 1. Rej (rejection identification) =0; 2. at least group_info number of Grade B; 3. there are at least 4 p value less or equal than 0.001.\n The filtered data is stored in the filtered.data.csv.")
k <- readline("Please type the number of at least B grade you want, e.g. 4: ") %>% as.numeric(.)
j <- readline("Please type the number of at least p values less or equal than 0.001 you want, e.g. 4: ") %>% as.numeric(.)
filtered_lipidomics <- lipidSelect %>% 
  rowwise() %>% 
  filter( Rej == 0 & sum(A, B) >= k & APvalue.001 >= j)

# store the filtered data into the filtered.data.csv file  
message("Please note that the filtered data are stored in filtered.data.csv.")
write.csv(filtered_lipidomics, "filtered.data.csv")



# input the group information and reterieve the data from the csv file.
samples <- Input(filtered_lipidomics)
# sample info
sample_info <- samples[[1]]
# group name info
group_names <- samples[[2]]
# group number info
ngroups <- samples[[3]]

# retrieve the group and its samples information
info_list <- RetrieveInfo(sample_info)
repeats <- info_list[[1]]
sample_list <- info_list[[2]]



# function to help finding the saturation pattern and counts. 
# SFA pattern is defined as ":0", MUFA is defined as ":1" 
# and PUFA is defined as ":[2-9]" (colon followed by any number from 2 to 9)
message("Please note that in the LipidMolec column, \n Pattern of SFA is defined as *:0 (colon followed by 0)\n MUFA is defined as *:1 (colon followed by) \n PUFA is defined as *:[2-9] (colon followed by any number from ranges of 2-9)")
CountPatternByRow <- function(x){
  pattern_name <- i <- j <- group_info <- c()
  i <- str_count(x, ":0")
  j <- str_count(x, ":1")
  group_info <- str_count(x, ":[2-9]")
  SFA <- paste(rep("SFA", i), collapse = "/")
  MUFA <- paste(rep("MUFA", j), collapse = "/")
  PUFA <- paste(rep("PUFA", group_info), collapse = "/")
  pattern_name <- paste( SFA, MUFA, PUFA, sep = "/")
  patterns <- data.frame(pattern= pattern_name, SFA=i, MUFA=j, PUFA=group_info, stringsAsFactors = FALSE) %>% as.matrix(.)
  return(patterns)
}

# make a data store the pattern and count information
lipid_list <- sapply(filtered_lipidomics$LipidMolec, function(x) CountPatternByRow(x)) 
lipid_list <- t(lipid_list) %>% data.frame(., stringsAsFactors = FALSE)
# name the pattern data
colnames(lipid_list) <- c("FA.name", "SFA", "MUFA", "PUFA")

# reformat the FA.name style by deleting extra speration mark "/"
message("Please note that any lipid molecule pattern which are not from the SFA, MUFA, PUFA are defined as None, e.g Q10")
lipid_list$FA.name <- apply(lipid_list, 1, function(x)x[1] <- x[1] %>% str_replace(., "//", "/") %>% str_remove(., "^/") %>% str_remove(., "/$") )
# assign the other pattern as None, e.g. Q10
lipid_list$FA.name[lipid_list$FA.name==""] <- "None"

# reorder the column position and put the counting patterns in between FA and FA Group key columns
count_lipid <- filtered_lipidomics %>% select(1:4) %>% cbind(., lipid_list)
count_lipid <- count_lipid %>% cbind(., filtered_lipidomics[5:ncol(filtered_lipidomics)]) 
rownames(count_lipid) <- NULL
# write the count of pattern into count_lipid.csv file. 
message("Please note the file contains the count for pattern SFA, MUFA and PUFA are stored in the count_lipid.csv")
write.csv(count_lipid, "data/count_lipid.csv")

# select columns of information for calculating the aggregation of MainArea of samples
selected_lipids <- count_lipid %>% select(Class,  FA.name, SFA, MUFA, PUFA, contains("MainArea"))
# transform attributes of SFA, MUFA, PUFA for calculations
names <- c("SFA", "MUFA", "PUFA")
DataAttribute <- function(x){class(x) <- as.numeric(x)}
transformed_lipid <- selected_lipids %>% mutate_at(names, DataAttribute)

#calculate the aggregation of lipid by same class and FA.name
aggregate_lipids <- transformed_lipid %>% group_by(Class, FA.name) %>% summarise_all(funs(sum)) 
observation_count <- transformed_lipid %>% group_by(Class, FA.name) %>% tally()
aggregate_lipids <- left_join(aggregate_lipids, observation_count) 
aggregate_lipids <- aggregate_lipids %>% select(Class, n, FA.name, SFA, MUFA, PUFA, sample_list)

# write the aggregation information into aggregated.csv file
message("Please note that the aggregated csv file will be stored in aggregated.csv")
aggregate_lipids %>% arrange(Class, FA.name) %>% write.csv(., "data/aggregated.csv")

# select three variables Class, n, FA.name for analysis
group_lipids <- aggregate_lipids %>% select(Class, n, FA.name)

# make a data frame contains sample information and group information
group_info <- data.frame(samples=sample_list, groups=repeats, stringsAsFactors = FALSE) %>% group_by(groups) 

# get the MainArea information for the groups
extracted_group<- list()
group_lists <- c()
group_area <- data.frame(row.names = 1:nrow(aggregate_lipids))
# extract sample information by its group
for(i in 1:ngroups){
  extracted_group[[i]] <- subset(group_info, groups==group_names[i], select = samples)
  assign(group_names[i], unlist(extracted_group[[i]]))
  group_lists[[i]] <- get(group_names[i])
  sum_samples <-  aggregate_lipids %>% ungroup() %>% transmute(area=rowSums(.[group_lists[[i]]]))
  group_area <- cbind(group_area, sum_samples)
}
colnames(group_area) <- group_names
  
# calculate the SFA, MUFA and PUFA's percentages 
fa_percent <- group_lipids %>% rowwise() %>% 
  # calculate how many SFA, MUFA and PUFA by row (by different SFA, MUFA and PUFA combination patterns)
  mutate(SFA = str_count(FA.name, "SFA"),
         MUFA= str_count(FA.name, "MUFA"), 
         PUFA= str_count(FA.name, "PUFA"), ) %>%   
  # calcutate SFA, MUFA and PUFA's percentage by row
  mutate(SF = SFA/(SFA+ MUFA + PUFA), MUF = MUFA/(SFA+ MUFA + PUFA), PUF= PUFA/(SFA+ MUFA + PUFA))  %>% 
  cbind(group_area, .) 

# calculate different groups' MainArea by its SFA, MUFA and PUFA percentage
each_fa <- data.frame(row.names = 1:nrow(aggregate_lipids))
f <- function(...){
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

names <- syms(group_names) 
# calculate the main area for SFA, MUFA, PUFA in each group
data <- f(!!!names)

# calculate the MainArea for SFA, MUFA, PUFA by Class in each group
area_groups <- data %>% ungroup() %>% group_by(Class) %>% summarise_at(.vars = colnames(.)[1:(4*ngroups)], sum)

write_csv(area_groups, "data/out.csv")
# #%>%   
#   mutate(SFA = SF*n, MUFA= MUF*n, PUFA=PUF*n) %>% ungroup() %>% group_by(Class) %>% 
#   summarise(SFA= sum(SFA), MUFA = sum(MUFA), PUFA = sum(PUFA))
# 
