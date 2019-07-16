####################################################################################
# Script: FWL_lipidomics.2.2.R
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
source("FWL_lipidomics.2.4.functions.R")

# Check directory existence, if not then make one
dirs <- c("plot", "data", "plot/classes")
mkdirs(dirs)

# store the file name in the variable file_name
message('Please store the file name in the file.name for check. And you could go to script to add the file name in read_csv(file.choose("data_raw/")). 
        If you need to change some parameters or edit the script, please copy the original script first. ')
file_name <- readline("Please input the file name: ")

# Read in file, and store the data 
lipidomics <- read_csv(file.choose("data_raw/"), col_types = cols())

# check if there is internal control samples, and extract information of different grades and p value 
# which are less or equal than 0.001 for each molecule and store this info in the variable lipid_count
message("Do you have internal controls like the standards or solvent?")
internal <- readline("Please type Y/N: ") %>% str_to_lower(.)
lipid_count <- InternalCheck(internal)

# add this count into original strucutes
lipid_select <- lipidomics %>% bind_cols(lipid_count)

# Filter the dataset based on your criteria  
# flexible parameter 4 which depends on the experiment for total number of A and B 
message("Please note the file is filtered by three standards. For same lipidmolecule: \n 1. Rej (rejection identification) =0; 2. at least k number of Grade B; 3. there are at least 4 p value less or equal than 0.001.\n The filtered data is stored in the filtered.data.csv.")
k <- readline("Please type the number of at least B grade you want, e.g. 4: ") %>% as.numeric(.)
j <- readline("Please type the number of at least p values less or equal than 0.001 you want, e.g. 4: ") %>% as.numeric(.)
filtered_lipidomics <- lipid_select %>% 
  rowwise() %>% 
  filter( Rej == 0 & sum(A, B) >= k & APvalue.001 >= j) %>% 
  as.data.frame(.)

write.csv(filtered_lipidomics, "data/filtered.raw.data.csv")


message("The info below and summary plot will show the summary information of classes after filtering the data")
# two methods check how many lipids passed filtering
print(describe(filtered_lipidomics$Class))
print(summary(as.data.frame(filtered_lipidomics$Class), maxsum=nrow(filtered_lipidomics)))


############################################################################
###  making group
message("Please input the info of the experimental groups below and end with 'Enter'
        Please Only input the sample names (replicate) for analysis.
        e.g. Experimental Group number: 2   
        Name of Experimental Group 1: wild_type
        Sample names (replicate number from lipid search) : s1 s2 s3 s4 
        Name of Group 2: knockout
        Sample names of Group 2: s4 s5 s6")

# input the group information and reterieve the data from the csv file.
samples <- Input(filtered_lipidomics)
# sample info
sample_info <- samples[[1]]
# group name info
group_names <- samples[[2]]
# group number info
ngroups <- samples[[3]]

label <- "first"
info_list <- PCA_pairs_Plot(sample_info, group_names, filtered_lipidomics, label)

# check if deleting outlier samples needed and plot new PCA
message("Please check the plots in the plot directory or r studio plots pannel.\nDo you need to re-enter the sample info for analysis?")
pca_check <- readline("Please type Y/N: ")
info_list <- PCAcheck(pca_check, filtered_lipidomics)

sample_list <- info_list[[1]]
group_repeats <- info_list[[2]]



# total class lipids of raw data by subtracting background (blank sample c) from all samples
subtracted_lipids <- filtered_lipidomics %>% select(Class, LipidMolec, contains("MainArea[c]"), sample_list) %>% mutate_at(sample_list, Subtra)
message("The subtracted samples from blank sample c are stored in the file subtracted_lipids")
write.csv(subtracted_lipids, "data/subtracted_lipids.csv")
# aggregat the lipid molecule that are same but have different retention time
subtracted_lipids <- aggregate(.~Class + LipidMolec, data=subtracted_lipids, FUN=sum)


# imputate data for each molecule which is same process as above 
log2_lipids <- subtracted_lipids %>% mutate_at(sample_list, log2trans) 
message("Please note that the log 2 transformed data for all molecules are stored in log.molec.csv")
write.csv(log2_lipids, "data/UnimputMolec.csv")
replace_inf_lipids <- log2_lipids %>% mutate_at(sample_list, ReplaceInf)
names <- replace_inf_lipids$LipidMolec
preImputated_lipids <- replace_inf_lipids %>% select(sample_list) 
imputated_lipids <- ImputeMinProb(preImputated_lipids, 0.01, 1)
rownames(imputated_lipids) <- names
write.csv(imputated_lipids, "data/imputeMolec.csv")

