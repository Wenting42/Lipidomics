####################################################################################
# Script: FWL_FattyAcids.functions.2.9.R
# Author: Wenting Lyu
# Notes: This script assist executing for main script FWL_lipidomics2.9.R which
#         helps generating the graph and data for the workflow of lipidomics.
#         it also calculate the SFA, MUFA, PUFA and the aggregated main area 
#         for each sample.
#         First, Make sure that your R and Rstudio are newly enough for installing the 
#         packages needed in the script. Otherwise the script will pop out warnings 
#         for packages and won't run.
#         Second, typing command in the console-----> source("FWL_lipidomics2.9.R")
#         or press the source button.
#         Third, users can independently implement this analysis by running 
#         "fattyAcidsaturation_analysis2.1.r" in directory fattyAcids_saturation_analysis.
#         This script is derived from Laura's project
#####################################################################################


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
aggregate_lipids <- transformed_lipid %>% group_by(Class, FA.name) %>% summarise_all(list(sum))
observation_count <- transformed_lipid %>% group_by(Class, FA.name) %>% tally()
aggregate_lipids <- left_join(aggregate_lipids, observation_count) 
aggregate_lipids <- aggregate_lipids %>% select(Class, n, FA.name, SFA, MUFA, PUFA, sample_list)

# write the aggregation information into aggregated.csv file
message("Please note that the aggregated csv file will be stored in aggregated.csv")
aggregate_lipids %>% arrange(Class, FA.name) %>% write.csv(., "data/aggregated.csv")

# select three variables Class, n, FA.name for analysis
group_lipids <- aggregate_lipids %>% select(Class, n, FA.name)

# make a data frame contains sample information and group information
group_info <- data.frame(samples=sample_list, groups=group_repeats, 
                         stringsAsFactors = FALSE) %>% 
  group_by(groups) 

# calculate the SFA, MUFA and PUFA's percentages 
fa_percent <- group_lipids %>% rowwise() %>% 
  # calculate how many SFA, MUFA and PUFA by row (by different SFA, MUFA and PUFA combination patterns)
  mutate(SFA = str_count(FA.name, "SFA"),
         MUFA= str_count(FA.name, "MUFA"), 
         PUFA= str_count(FA.name, "PUFA"), ) %>%   
  # calcutate SFA, MUFA and PUFA's percentage by row
  mutate("%SFA" = SFA/(SFA+ MUFA + PUFA), "%MUFA" = MUFA/(SFA+ MUFA + PUFA), "%PUFA"= PUFA/(SFA+ MUFA + PUFA)) 

# combine the percentage information with aggregated lipid info
fa_percent <- aggregate_lipids %>% ungroup() %>% select(sample_list) %>% cbind(fa_percent, .)


# find median and mean for each sample in each group and calculate its value
data <- CalcGroup( fa_percent, group_info)
write.csv(data, "data/mean_median.csv")

# reformat the data for later visualization
mean_base <- ReformatData(data, "mean")
mean1 <- mean_base[[1]]
write.csv(mean1, "data/mean_data.csv")
mean2 <- mean_base[[2]]
write.csv(mean2, "data/eachClass.mean.csv")

# visualization mean based sample in each group
p1 <- PlotStack(mean2) + labs(title = "mean based stack plots") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
print(p1)
ggsave(filename = "meanBased.stackplots.pdf", path="plot/", device="pdf")

# # visualiaztion based 
# p2 <- PlotGroups(mean2)
# print(p2)
# ggsave(filename = "meanBased.fattyAcids.pdf", path="plot/", device="pdf")
# 
# # normalized the data for visualization
# dt <- PassToNormal(mean2)
# p3 <- PlotGroups(dt)
# print(p3)