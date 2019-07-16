####################################################################################
# Script: FWL_lipidomics.2.7.R
# Author: Wenting 
# Notes: 
#         This script is based on Niklas and Kenny's previous massive work. 
#         It helps generating the graph and data for the flowork of Lipidomics. 
#
#         To start, typing command in the console-----> source("FWL_lipidomics2.0.R") 
#
#         1) interactive way.  
#         First please store the file name in the variable file.name for check 
#         Run the script by command "source("FWL_lipidomics2.7.R") in console.
#         All the data are stored in the directory called "data", and plots in "plot"
#         2) You can change the parameters in the script and run it lines by lines.
#           Before you change any parameters, please copy the original script first.
#           

#####################################################################################


rm(list=ls())  ##### empty the environment first. If you need to keep former variables in workspace, just comment this part
setSessionTimeLimit(cpu = Inf, elapsed = Inf)
# source function from FWL_lipidomics.**.functions.R script
source("FWL_lipidomics.2.7.functions.R")
# Check directory existence, if not then make one
# all the plots will stored in plot directory and data in the data directory
dirs <- c("plot", "data", "plot/classes")
mkdirs(dirs)

# store the file name in the variable file_name for check
message('Please store the file name in the file.name for check. And you could go to script to add the file name in read_csv(file.choose("data_raw/")). 
        If you need to change some parameters or edit the script, please copy the original script first. ')
file_name <- readline("Please input the file name: ")

# Read in file, and store the data 
lipidomics <- read_csv(file.choose("data_raw/"), col_types = cols())

# check if there is internal control samples, and extract information of different grades and p value 
# which are less or equal than 0.001 for each molecule and store this info in the variable lipid_count
message("Do you have internal controls like the standards or solvent?")
internal <- readline("Please type Y/N: ") %>% str_to_lower(.)
lipid_check <- InternalCheck(internal, lipidomics)
lipid_count <- lipid_check[[1]]

# remove standards
lipid_remove <- lipid_check[[2]]

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


# mark the artifical classes position
odd.index <- filtered_lipidomics$LipidMolec %>%
  str_locate(., "(\\d[13579]:)|([13579]:)")
percent.odd <- length(unique(which(!is.na(odd.index), arr.ind=TRUE)[,1]))/nrow(filtered_lipidomics) 
percent.odd <- percent(percent.odd*100/100)
message("The odd chain of fatty acids percent is ", percent.odd, " in total.")



########################3
message("The info below and summary plot will show the summary information of classes after filtering the data")
# two methods check how many lipids passed filtering
print(describe(filtered_lipidomics$Class))
print(summary(as.data.frame(filtered_lipidomics$Class), maxsum=nrow(filtered_lipidomics)))

# plot class proportion black white plot
prop.bw.plot <- PropPlot(filtered_lipidomics) + scale_fill_grey()
ggsave(filename = "Prop.bw.summary.pdf", path = 'plot/', device = "pdf")

##########################################################################
# QC PLOT 1 - Abundance vs. retention time [requires ggplot2]
##############################################################################################
### individual sample retention plot
## if you want individule sample retention time vs. abundance plot, just uncomment the part below
# message("Please input the MainArea of the sample name to check its abundance vs. retention time, 
#         e.g. `MainArea[s1]`. 
#         Please note that: ` ` is the backtick which is same position as tilde")
# 
# abundance.sample <- readline("Please type the sample name as this format. e.g. `MainArea[s1]` -----> ")
# IndividuleRetentionPlot(filtered_lipidomics, abundance.sample) 

### hard code way for ggplot
# ggplot(data=filtered_lipidomics, aes(x =log10(`MainArea[s1]`), y = BaseRt)) +
#   geom_point() +
#   theme_bw() +
#   facet_grid(.~Class) +
#   labs("Abundance VS. Retention time", x="lipid class (log10(MainArea))", y="Retention time")


###  Abundance vs. retention time for all samples
retention_data <- filtered_lipidomics %>% 
  select(contains("MainArea[s"), BaseRt, Class) %>% 
  gather(sample, MainArea, -c(BaseRt, Class))
AllRetentionPlot(retention_data)
##############################################################################################


################################################################################
# detect same lipid molecules with different retention time
### Filter the lipid molecule contains same name but different retention time based on your criteria  
duplicate_molecs <- Detect_duplicates(filtered_lipidomics)
# diff_baseRT will store the retention time differences for the same lipid molecule
diff_baseRT <- duplicate_molecs %>% select(LipidMolec, BaseRt) %>% group_by(LipidMolec) %>% summarise(gap=max(BaseRt)-min(BaseRt)) %>% view()
message("Please note that the difference of retention time for same lipid molecule is stored in diffRT.csv")
write_csv(diff_baseRT, "data/diffRT.csv")

# if move on, fix method for duplicated molecules
if(nrow(duplicate_molecs) != 0){
  message("!!!warning: There exist same lipid molecules but different retention time.\n !!!!!! There might be contaminations in the samples. \n If you check and still wanna move on, there are two methods for handling the data. \n\nPlease note that the duplicated molecule with different retention time is stored in data/duplicated.molecules.csv")
  write_csv(duplicate_molecs, "data/duplicated.molecules.csv")
  # get all names of sample
  all_names <- colnames(duplicate_molecs)[-c(1:3)]
  # two methods to fix duplicate lipid molecules
  message("Here provides two methods to deal with the duplicated molecules. \n Options:\n A: get the lipid molecule which summation of main area is bigger. \n B: aggregate the two different retention time molecule \n Please note that the aggregation method will lose BaseRt information")
  fix_method <- readline("Please input the method you like to handling the duplicated molecules, please type A or B: ")
  fixed_dt<- fix_duplicate(filtered_lipidomics, duplicate_molecs, fix_method)
  # transform duplicate lipid molecules
  duplicate_output <- fixed_dt[[1]]
  # store reserved duplicate lipid molecules
  message("Please note that the reserved duplicate lipid molecules are stored in reserved_duplicates.csv")
  write_csv(duplicate_output, "data/reserved_duplicates.csv")
  # get filtered lipid information within processed duplicate lipid molecules
  cleaned_dt <- fixed_dt[[2]]
  filtered_lipidomics <- cleaned_dt
  message("Please note that the filtered lipid molecules after removing the duplicates are stored in removeduplicates.csv")
  write_csv(filtered_lipidomics, "data/removeduplicates.csv")
}


# plot the background information of the blank sample c
blank_sample <-filtered_lipidomics %>% select(Class, contains("MainArea[c]")) %>% 
  group_by(Class) %>% summarise_all(.funs=sum)
colnames(blank_sample) <- c("class", "value")

# since the value of blank sample is very large and by adding 1 to fix the log transformation of 0 value 
# for visualization
blank_sample$value <- blank_sample$value + 1
# plot the blank sample
ggplot(blank_sample, aes(x=reorder(class, value), y=value)) + geom_bar(position="dodge", stat="identity") +
  scale_y_log10( breaks = scales::trans_breaks("log2", function(x) 10^x),
                 labels = scales::trans_format("log2", scales::math_format(10^.x)))  +
  labs(title = "Area under Curve for blank sample (background) in each class", 
       x = "Lipid Class", y="Main Area") + 
  theme_bw() +
  coord_cartesian(ylim=c(1, 1e10))

ggsave(filename = "background.png", path="plot/", device = "png")

############################################################################
# Make groups
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


# detect empty and negative lipids
empty_data <- filtered_lipidomics %>% filter_at(sample_list, any_vars(.<= 0))
counts <- DetectEmpty(ngroups, sample_info, empty_data)
write_csv(counts, "data/checkInvalid.csv")
message("Please mannualy check the data contains potential invalid values.\n And the delete the row of molecules which is invalid to you. And save the changes")


# deleting the lipid molecules you select 
message("Now we need to open the changed file checkInvalid.csv which delete the molecules you don't trust")
filter_empty <- read_csv(file.choose("data/checkInvalid.csv"))
delete_empty <- anti_join(empty_data, filter_empty, by = "LipidMolec")
filtered_lipidomics <- anti_join(filtered_lipidomics, delete_empty, by = "LipidMolec")
write_csv(filtered_lipidomics, "data/post.filtered.lipidomics.csv")

# calculate the saturation for different lipid class
source("FattyAcids.R")
message("The SFA, MUFA, PUFA information will be stored in the count_lipid.csv and aggregated.csv")

#########################################################################
######## making barplot graphs
# raw data 
### For total classes, reformatted the data
total_data1 <- TotalClass(filtered_lipidomics, sample_list, group_repeats) 

# total class lipids of unprocessed raw data visualization
total_plot <- ClassPlot(total_data1[[1]]) + 
  labs(title = "Total lipid classes", x = "Lipid Classes", y= "AUC (Area under curve)") + 
  geom_errorbar(aes(ymin=value, ymax=value+sd), width=.3, 
                position=position_dodge(.9)) +
  scale_fill_manual(values = clPalette1)
print(total_plot)                                                                                                                 
ggsave(filename="total.class.png", path = 'plot/', device = "png", dpi = 300)



# total class lipids of raw data by subtracting background (blank sample c) from all samples
subtracted_lipids <- filtered_lipidomics %>% select(Class, LipidMolec, contains("MainArea[c]"), sample_list) %>% mutate_at(c("MainArea[c]", sample_list), list(~.-`MainArea[c]`))
message("The subtracted samples from blank sample c are stored in the file subtracted_lipids")






# Make bar plots for the subtracted data
total_data2 <- TotalClass(subtracted_lipids, sample_list, group_repeats) 
total_data_long <- total_data2[[1]]
total_data_wide <- total_data2[[2]]
#total_data_long$value <- ifelse(total_data_long$value<0, 0, total_data_long$value)
total_plot2 <- ClassPlot(total_data_long) + 
  labs(title = "Total lipid classes", x = "Lipid Classes", y= "Main Area (background subtracted)") + 
  geom_errorbar(aes(ymin=value, ymax=value+sd), width=.3, 
                position=position_dodge(.9)) 
print(total_plot2)                                                                                                                 
ggsave(filename="subtracted.total.class.png", path = 'plot/', device = "png")

# log transformation and visualize approximately normal distribution of the transformed data
log2_lipids <- subtracted_lipids %>% mutate_at(sample_list, log2trans)
ggplot(data=log2_lipids, aes(x=`MainArea[s3]`)) + geom_density() 
ggplot(data=log2_lipids, aes(x=`MainArea[s3]`)) + geom_histogram(bins=30) 



## Fold change plot by taking median sample in the control group
control <- readline("Please input the number of controls for visualization fold change, e.g. 1: ")
# aggregate same class lipid
t_data <- subtracted_lipids %>% select(-LipidMolec)
f_input <- aggregate(.~Class, data=t_data, FUN = sum)
# take log transformation for subtracted data 
log2_input <- f_input %>% mutate_at(sample_list, log2trans)
# replace negative infinity value with NA
inf_input <- log2_input %>% mutate_at(sample_list, ReplaceInf)
# imputate data for each class
fold_names <- inf_input$Class
imputated_class <- inf_input %>% select(sample_list) %>% ImputeMinProb(., 0.01, 1)
rownames(imputated_class) <- fold_names
message("The total class log transformed data will stored in log2.data.csv")
write.csv(imputated_class, "data/imputClass.csv")
# get fold change information for each class in each group
foldchange_input <- t(imputated_class) %>% cbind(experiment.group = group_repeats, .) %>% as.data.frame(.)
foldchange_input[, -1]<- sapply(foldchange_input[, -1], as.numeric) %>% as.data.frame(.)
message("The total class plot reflect fold changes among the groups by using the imputed log transformed median sample")
# plot for imputated class data 
PassToNormal(foldchange_input, control)


############################
### plot for each class
class_data <- EachClass(subtracted_lipids)
class_long <- class_data[[1]]
class_wide <- class_data[[2]]
class_wide_list <- sapply(class_wide, function(x) x <- t(x[,-1]))
lipidmolecNO_max <- sapply(class_wide_list, function(x) nrow(x)) %>% sum()

# setting the plot limits when the bar numbers exceed the threshold nbar
nbar <- 64     # estimation of threshold which can be modified and at least bigger than group number
EachClassPlot(nbar, ngroups, class_long)

# overview of each class plot
nbar <- lipidmolecNO_max
EachClassPlot(nbar, ngroups, class_long)

###################################################################################
# volcano plot
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



# Create a design matrix 
samples <- factor(group_repeats)
design <- model.matrix(~0+samples)
colnames(design) <- levels(samples)
message("Begin to plot volcano \nPlease note that you NEED make contrast groups manually to Compare the difference between/among groups.\ne.g. compare B+C against A: A-(B+C)/2; A against B:  A-B; A-B against C-D: (A-B)-(C-D), etc.")
ncomparisons <- readline("Please input the Number of comparison groups, e.g. 2: ")


# Fit model and extract contrasts 
fit <- lmFit(imputated_lipids, design)
# plot volcano graph
VolPlot(ncomparisons, fit, imputated_lipids)
###################################################################################################



### pathway
message("Please install Cytoscape 3.6.1 or greater, which can be downloaded from \nhttp://www.cytoscape.org/download.php.")
message("Once download, please make sure that you open the cytoscape while running for pathway part")
check_cytoscape<- readline("Do you download Cytoscape 3.6.1 or higher and opened it? Y/N: ") %>% str_to_lower(.)
checkDownload(check_cytoscape)

# confirmed that cytoscape is connected
cytoscapePing ()
cytoscapeVersionInfo ()

# if you install cytoscape 3.7 or higher then apps can be installed directly from R.
installation_responses <- c()
#list of app to install
cyto_app_toinstall <- c("clustermaker2", "enrichmentmap", "autoannotate", "wordcloud", "stringapp", "aMatReader")
cytoscape_version <- unlist(strsplit( cytoscapeVersionInfo()['cytoscapeVersion'],split = "\\."))
if(length(cytoscape_version) == 3 && as.numeric(cytoscape_version[1]>=3) && as.numeric(cytoscape_version[2]>=7)){
  for(i in 1:length(cyto_app_toinstall)){
    #check to see if the app is installed.  Only install it if it hasn't been installed
    if(!grep(commandsGET(paste("apps status app=\"", cyto_app_toinstall[i],"\"", sep="")), 
             pattern = "status: Installed")){
      installation_response <-commandsGET(paste("apps install app=\"", 
                                                cyto_app_toinstall[i],"\"", sep=""))
      installation_responses <- c(installation_responses,installation_response)
    } else{
      installation_responses <- c(installation_responses,"already installed")
    }
  }
  installation_summary <- data.frame(name = cyto_app_toinstall, 
                                     status = installation_responses)
  list(installation_summary)
}

# help(package=RCy3)

setNodeShapeDefault("ellipse")
setNodeSizeDefault("30")
setNodeFontSizeDefault("20")
initial_graph <- buildPathWay()
# make duplicate graph for contrast 
h <- cloneNetwork(initial_graph)

message("Please type the control group and contrast groups for making pathway plot")
fc_data <- TotalNormalClass(foldchange_input, 1)
fc_data_long <- fc_data[[2]][[1]]
fc_class <- fc_data_long %>% spread(., key = experiment.group, value = value)


control.group <- fc_data_long$experiment.group[1]
path.class <- c('fatty acids', 'G3P', 'LCBs', 'Cer', 'SM', 'LPA', 'PA', 'DAG', 'TAG', 'PC',
                'LPC', 'PE', 'LPE', 'CDP-DAG', 'PI', 'LPI', 'PG', 'LPG', 'PS', 'LPS')
exist.class <- intersect(fc_class$class, path.class)
detect.lipid <- fc_class %>% filter(class %in% exist.class)
empty.lipid <- path.class[!path.class %in% exist.class]
empty.class <- sapply(group_names, function(x) ifelse(x == control.group, return(rep(1, length(empty.lipid))), return(rep(0, length(empty.lipid)))))
empty.data <- data.frame(class=empty.lipid, empty.class)
path.data <- rbind(detect.lipid, empty.data)
nodes <- path.data$class %>% as.character(.)
foldpath <- data.frame(log2fc = path.data$G2, row.names = nodes, stringsAsFactors = FALSE, key.column = "class")
write_csv(foldpath, "data/path.csv")
loadTableData(foldpath)





#dev.off()
