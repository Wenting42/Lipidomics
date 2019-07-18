####################################################################################
# Script: FWL_FattyAcids.functions.R
# Author: Wenting Lyu
# Notes: This script assist executing for main script FWL_lipidomics2.8.R which
#         helps generating the graph and data for the workflow of lipidomics.
#         it also calculate the SFA, MUFA, PUFA and the aggregated main area 
#         for each sample.
#         First, Make sure that your R and Rstudio are newly enough for installing the 
#         packages needed in the script. Otherwise the script will pop out warnings 
#         for packages and won't run.
#         Second, typing command in the console-----> source("FWL_lipidomics2.8.R")
#         or press the source button.
#         Third, users can independently implement this analysis by running 
#         "fattyAcidsaturation_analysis2.1.r" in directory fattyAcids_saturation_analysis.
#         This project is derived from Laura's idea
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
group_info <- data.frame(samples=sample_list, groups=group_repeats, stringsAsFactors = FALSE) %>% group_by(groups) 

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


names <- syms(group_names) 
# calculate the main area for SFA, MUFA, PUFA in each group
data <- CalculateMainArea(!!!names)
message("Please note that the main area info for FA in each group is stored in eachFA.csv")
write_csv(data, "data/eachFA.csv")

# calculate the MainArea for SFA, MUFA, PUFA by Class in each group
area_groups <- data %>% ungroup() %>% group_by(Class) %>% summarise_at(.vars = colnames(.)[1:(4*ngroups)], sum)
message("Please note that MainArea for SFA, MUFA, PUFA by class in each group information is stored in classFA.csv")
write_csv(area_groups, "data/classFA.csv")


# reformat data structure for visualization
rdata <- area_groups %>% 
  select(-group_names) %>% 
  gather(., key="saturation_of_group", value="value", -Class) %>% 
  mutate(experiment_group = str_remove(saturation_of_group, "_.*")) %>% 
  mutate_at(vars(saturation_of_group),  list(~str_remove(saturation_of_group, ".*_")))


# plot saturantion status in different groups of different class
stackplots <- ggplot(data = rdata, aes(x=experiment_group, y=value, fill = saturation_of_group)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Class, scales = "free") +
  scale_fill_brewer(palette = "Pastel2")
print(stackplots)
ggsave(filename = "stackplots.pdf", path="plot/", device="pdf")

# plot saturation status in different groups among all class
fattyplots <- ggplot(data = rdata, aes(x=Class, y=value, fill=experiment_group)) + 
  geom_bar(stat="identity", position="dodge") + 
  scale_fill_manual(values = clPalette1)
print(fattyplots)
ggsave(filename = "fatty_acids.pdf", path="plot/", device="pdf")
