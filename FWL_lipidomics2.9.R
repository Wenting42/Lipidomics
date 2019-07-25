# Script: FWL_lipidomics2.8.R
# Author: Wenting 
# Notes:  To start, typing command in the console-----> source("FWL_lipidomics2.9.R") 
#         or press the source butthon. 
#         Please open the cytoscape on your computer and make sure Mac users installed XQuartz.
#         This script is designed as a reference for lipidomics experiment. 
#         It is based on Niklas and Kenny's previous work (Their work files can be found in folders 
#         quality_control and statistics_quantification). Acknowledge to Weng's technical 
#         guidance, Laura's fatty acid saturation analysis project and Sebstian's shadow experiment
#
# Usage:  1) Type the name of transformed data file you want for analysis and choose it.
#         2) Input the sample name of internal standandards seperated by space, the samples are
#             reformated and case insensitive, e.g. S18 S19 S20. 
#         3) Input the quality control standandards for Lipid Search software, the usual 
#            standards are 4 for numbers of grades and numbers of APValue less than 0.001.
#         4) Please check the proportion difference of retention time plots for quality control.
#            And check the retention time difference for same lipid molecule in data/diff_RT.csv
#            If the time gaps are too large. Users will need consider redo the experiment. 
#            If continue, users will need choose from 2 methods for dealing with duplicated lipid 
#             molecules by typing A or B. 
#            # Method A will pick the main lipid molecule which area under curve is largest.
#            # Method B will aggregate the duplicated molecules.
#         5) Please read background.png for estimating background signal. The background will
#             be removed from the main area under curve for later analysis.
#         6) After User input the group information, please check the file data/checkInvalid.csv 
#             which store the lipid molecule information after removing background. 
#             Users need delete lipid molecule based on the experiment design in this file or not,
#             then choose this file for open. 
#            # Please note that the rest invalid lipid molecules will be imputated later.
#         7) saturation analysis will combined into this step for the filtered data.
#         8) The data will be log transformed (data/log2.data.csv) and imputated (data/imputClass.csv)
#            for class summary plots.
#         9) For vocano plots, the lipid molecules data will also take log transformed (data/log.molec.csv) 
#            and imputated (data/imputeMolec.csv). Users need input contrast groups, e.g, FLOX - LKO.
#            The other outputs will be stored in files named by users.
#         10) The pathway visualization part is connected to cytoscape. Please be sure that the cytoscape is
#              open when running the pipeline.
#      
# Status:  In progress             

#####################################################################################


rm(list=ls())  ##### empty the environment first. If you need to keep former variables in workspace, just comment this part
#setSessionTimeLimit(cpu = Inf, elapsed = Inf)
# source function from FWL_lipidomics.**.functions.R script
source("FWL_lipidomics.FUNCTIONS2.9.r")
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

# add this count into original strucutes and remove the internal standards
lipid_select <- lipid_remove %>% bind_cols(lipid_count)


## or without remove the internal standards
#lipid_select <- lipidomics %>% bind_cols(lipid_count)


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
              str_detect(., "(\\d[13579]:.*)|([13579]:.*)")
percent.odd <- length(which(odd.index == TRUE))/nrow(filtered_lipidomics) 
percent.odd <- percent(percent.odd*100/100)
message("The odd chain of fatty acids percent is ", percent.odd, " in total.")

# mark the odd chains in TG class
message("The odd chain of TG information is stored in TG.odd.csv. Please have a look and 
        decide if you wanna delete any lipid molecules. If you do, just delete the corresponding
        lipid in TG.odd.csv. The pipeline will take care of the rest")
TGs <- filtered_lipidomics[which(odd.index == TRUE), ] %>% filter(Class=="TG") 
write.csv(TGs, "data/TG.odd.csv")


TTG17 <- TGs[which(TGs$LipidMolec == "TG(17:1/17:1/17:1)"), ] %>% select(contains("MainArea"))
all_lipids <- filtered_lipidomics[, which(str_detect(colnames(TGs), "MainArea"))] %>% summarise_all(sum)
other_lipids <- all_lipids - TG17
lipids <- rbind(TG17, all_lipids)
rownames(lipids) <- c("TG(17:1/17:1/17:1)", "Other_LipidMolecs")
colnames(lipids) <- colnames(lipids) %>% str_replace_all(., "MainArea\\[", "") %>% str_remove_all(., "\\]")
lipids$Pattern <- rownames(lipids)
transformed_lipids <- gather(lipids, Sample, Value, -Pattern)

# plot standard TG(17:1/17:1/17:1) in all samples
stackplots <- ggplot(transformed_lipids, aes(x=reorder(Sample, Value), y = Value, fill = Pattern)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Set2") +
  theme_bw()+
  theme(axis.text.x = element_text(size = 8),
        panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        #strip.background = element_rect(colour="black", fill="white", 
        #                               size= 0.5, linetype="solid")
        #strip.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        ) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  labs(title = "TG(17:1/17:1/17:1) lipid molecules in all samples", xlab= "samples",
       ylab = "aggregated AUC (Main Area under curve)")
print(stackplots)
ggsave(filename = "TG17.pdf", path="plot/", device="pdf")

# after manually deleting some odd chains
filter_odd <- read_csv(file.choose("data/TG.odd.csv"))
delete_odd <- anti_join(TGs, filter_odd, by = "LipidMolec")
filtered_lipidomics <- anti_join(filtered_lipidomics, delete_odd, by = "LipidMolec")
write_csv(filtered_lipidomics, "data/deleteODDTG.lipidomics.csv")

########################3
message("The info below and summary plot will show the summary information of classes after filtering the data")
# two methods check how many lipids passed filtering
print(describe(filtered_lipidomics$Class))
print(summary(as.data.frame(filtered_lipidomics$Class), maxsum=nrow(filtered_lipidomics)))

# plot class proportion black white plot
prop.bw.plot <- PropPlot(filtered_lipidomics) + scale_fill_grey()
print(prop.bw.plot)
ggsave(filename = "Prop.bw.summary.pdf", path = 'plot/', device = "pdf")


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
diff_baseRT <- duplicate_molecs %>% select(LipidMolec, BaseRt) %>% group_by(LipidMolec) %>% summarise(gap=max(BaseRt)-min(BaseRt)) 
message("Please note that the difference of retention time for same lipid molecule is stored in diff_RT.csv")
write_csv(diff_baseRT, "data/diff_RT.csv")

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
  theme( panel.grid.major = element_blank(),
         )
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

# detect invalid lipids (AUC of lipid molecule is 0 or negative) in raw data
empty_data <- filtered_lipidomics %>% filter_at(sample_list, any_vars(.<= 0))
counts <- DetectEmpty(ngroups, sample_info, empty_data)
write_csv(counts, "data/checkInvalid.raw.csv")


# deleting the lipid molecules you select 
message("Now we need to open the changed file checkInvalidRaw.csv which delete the molecules you don't trust")

filter_empty <- read_csv(file.choose("data/checkInvalidRaw.csv"))
delete_empty <- anti_join(empty_data, filter_empty, by = "LipidMolec")
filtered_lipidomics <- anti_join(filtered_lipidomics, delete_empty, by = "LipidMolec")
write_csv(filtered_lipidomics, "data/post.filtered.lipidomics.csv")


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
  theme_bw()+
  theme(axis.text.x = element_text(size = 8),
         panel.grid.major = element_blank(),
         #panel.grid.minor = element_blank(),
         legend.title = element_blank(),
         #strip.background = element_rect(colour="black", fill="white", 
         #                               size= 0.5, linetype="solid")
         strip.background = element_blank(),
         axis.line = element_line(colour = "black"),
         panel.border = element_blank(),
         )
print(total_plot)                                                                                                                 
ggsave(filename="total.class.png", path = 'plot/', device = "png", dpi = 300)



# total class lipids of raw data by subtracting background (blank sample c) from all samples
subtracted_lipids <- filtered_lipidomics %>% 
  mutate_at(c("MainArea[c]", sample_list), list(~.-`MainArea[c]`))
message("The subtracted samples from blank sample c are stored in the file subtracted_lipids.csv")
write.csv(subtracted_lipids, "data/subtracted_lipids.csv")


# Make bar plots for the subtracted data
total_data2 <- TotalClass(subtracted_lipids, sample_list, group_repeats) 
total_data_long <- total_data2[[1]]
total_data_wide <- total_data2[[2]]
#total_data_long$value <- ifelse(total_data_long$value<0, 0, total_data_long$value)
total_plot2 <- ClassPlot(total_data_long) + 
  labs(title = "Total lipid classes", x = "Lipid Classes", y= "Main Area (background subtracted)") + 
  geom_errorbar(aes(ymin=value, ymax=value+sd), width=.3, 
                position=position_dodge(.9)) +
  theme_bw()+
  theme(axis.text.x = element_text(size = 8),
        panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        #strip.background = element_rect(colour="black", fill="white", 
        #                               size= 0.5, linetype="solid")
        strip.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
  )
print(total_plot2)                                                                                                                 
ggsave(filename="subtracted.total.class.png", path = 'plot/', device = "png")

# log transformation and visualize approximately normal distribution of the transformed data
log2_lipids <- subtracted_lipids %>% mutate_at(sample_list, log2trans)

# randomly choose a sample to see its distribution approximately normal distribution
i <- sample(length(unique(sample_list)), 1)
ggplot(data=log2_lipids, aes(x = !!sym(sample_list[i]))) + geom_density() 
ggplot(data=log2_lipids,  aes(x = !!sym(sample_list[i]))) + geom_histogram(bins=35) + 
  geom_bar(aes(y = (..count..)/sum(..count..))) +
  theme_bw() +
  theme(text = element_text(size = 8),
        axis.text = element_text(size = 6), 
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        axis.title = element_text(size = 10),
        panel.background = element_rect(fill = "white")) + 
  labs(x = "sample distribution", y = "lipid molecule count") 


# detect invalid lipids (AUC of lipid molecule is 0 or negative) in subtracted data
empty_data <- subtracted_lipids %>% filter_at(sample_list, any_vars(.<= 0))
counts <- DetectEmpty(ngroups, sample_info, empty_data)
write_csv(counts, "data/checkInvalid.csv")


# deleting the lipid molecules you select 
message("Now we need to open the changed file checkInvalid.csv which delete the molecules you don't trust")

filter_empty <- read_csv(file.choose("data/checkInvalid.csv"))
delete_empty <- anti_join(empty_data, filter_empty, by = "LipidMolec")
subtracted_lipids <- anti_join(subtracted_lipids, delete_empty, by = "LipidMolec")
write_csv(subtracted_lipids, "data/post.subtracted.lipidomics.csv")

# calculate the saturation for different lipid class
source("FWL_FattyAcids.functions.2.9.R")
message("The SFA, MUFA, PUFA information will be stored in the count_lipid.csv and aggregated.csv")



## Fold change plot by taking median sample in the control group
control <- readline("Please input the number of controls for visualization fold change, e.g. 1: ")
# aggregate same class lipid
t_data <- subtracted_lipids %>% select(Class, contains("MainArea[c]"), sample_list)
f_input <- aggregate(.~Class, data=t_data, FUN = sum)
# take log transformation for subtracted data 
log2_input <- f_input %>% mutate_at(sample_list, log2trans)
# replace negative infinity value with NA
inf_input <- log2_input %>% mutate_at(sample_list, ReplaceInf)


# imputate data for each class
fold_names <- inf_input$Class

# impuate with random minimal value from normal distribution
imputated_class <- inf_input %>% select(sample_list) %>% ImputeMinProb(., 0.01, 1)
# imputate with 1
imputated_class <- inf_input %>% mutate_at(sample_list, ReplaceNa) %>% select(-Class, contains("[c]"))

rownames(imputated_class) <- fold_names
message("The total class log transformed data will stored in log2.data.csv")
write.csv(imputated_class, "data/imputClass.csv")
# get fold change information for each class in each group
foldchange_input <- t(imputated_class) %>% 
  cbind(experiment.group = group_repeats, .) %>% 
  as.data.frame(.)
foldchange_input[, -1]<- sapply(foldchange_input[, -1], as.numeric) %>% 
  as.data.frame(.)
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

dev.off()
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
write.csv(log2_lipids, "data/log.molec.csv")
replace_inf_lipids <- log2_lipids %>% mutate_at(sample_list, ReplaceInf)
names <- replace_inf_lipids$LipidMolec
preImputated_lipids <- replace_inf_lipids %>% select(sample_list) 
imputated_lipids <- ImputeMinProb(preImputated_lipids, 0.01, 1)
rownames(imputated_lipids) <- names
message("Please note that the imputated log transformed data will stored in imputeMolec.csv")
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


contmatrix <- makeContrasts(FLOX-LKO, levels=design)
fit2 <- contrasts.fit(fit, contmatrix)
fit3 <- eBayes(fit2)
comparison <- output1 <- topTable(fit3, coef=1, adjust.method = 'fdr',
                                  lfc=0, number=nrow(imputated_lipids)) %>% as.data.frame()


comparison$sig <- factor(comparison$adj.P.Val< 0.05 & abs(comparison$logFC) > 1)

points <- comparison %>%  rownames_to_column('lipid') %>%
  filter(abs(logFC) > 1 & adj.P.Val< 0.05 ) %>%
  column_to_rownames('lipid')

foldpoints <-  comparison %>%  rownames_to_column('lipid') %>%
  filter(abs(logFC) > 1 & adj.P.Val > 0.05 ) %>%
  column_to_rownames('lipid')

ppoints <-  comparison %>%  rownames_to_column('lipid') %>%
  filter(abs(logFC) < 1 & adj.P.Val < 0.05 ) %>%
  column_to_rownames('lipid')

npoints <-  comparison %>%  rownames_to_column('lipid') %>%
  filter(abs(logFC) < 1 & adj.P.Val > 0.05 ) %>%
  column_to_rownames('lipid')



comparison$sig <- ifelse(comparison$sig=="FALSE", 
                                case_when( 
                                            (abs(comparison$logFC) < 1 & comparison$adj.P.Val > 0.05) ~ "NS", 
                                            (abs(comparison$logFC) > 1 & comparison$adj.P.Val > 0.05) ~ "foldPoints", 
                                            (abs(comparison$logFC) < 1 & comparison$adj.P.Val < 0.05) ~ "pPoints",
                                            (abs(comparison$logFC) > 5 & comparison$adj.P.Val > 0.05) ~ "5fold"
                                          ), 
                                "significant")



class_names <- rownames(comparison) %>% str_extract_all(., "(.+)\\(") %>% str_remove_all(., "\\(")
input_lipids <- comparison

input_lipids$lipclass <- class_names
input_lipids$lipclass <- case_when( (input_lipids$lipclass %in% c("PC", "PE", "PG", "PI", "PS", "LPC", "LPE",
                                                                         "LPI", "dMePE", "CL"))~ "Glycerophospholipids",
                                           input_lipids$lipclass %in% c("TG", "DG", "MG") ~ "Neutral lipids",
                                           input_lipids$lipclass %in% c("SM", "So", "SoG1", "Cer", "CerG1","CerG2",
                                                                        "GM2","GM3", "CerG2GNAc1", "CerG3GNAc1", "CerG3NAc2") ~ "Sphingolipids",
                                           input_lipids$lipclass %in% c("ChE","Cholestoral") ~ "Sterols",
                                           TRUE ~ "Other lipids")
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
fc_points <- comparison %>% 
  rownames_to_column('lipid') %>% 
  filter(abs(logFC) > 5 & adj.P.Val>= 0.05 ) %>% 
  column_to_rownames('lipid')




ggplot(input_lipids, aes(logFC, -log10(adj.P.Val), color = sig)) +
  geom_point(alpha= 0.8) + geom_vline(xintercept = c(-1, 1), linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="black", size=0.2)+
  xlab(bquote("Fold change, " ~ log[2]~"("~textstyle(frac({x}, {y}))~scriptstyle(, AUC)~")")) +
  ylab(bquote(~-Log[10]~adjusted~italic(P))) +
  ggtitle("volcano plot") +
  theme_bw()+
  theme(line=element_blank(),
        axis.line = element_line(colour = "black", size = 1),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_color_manual(values = c("#67a9cf", "#ef8a62", "#2166ac","#b2182b", "#d1e5f0" ), drop=FALSE) +
  # geom_text_repel(data = subset(input_lipids, sig=="significant"), 
  #                 aes(label=lipclass, 
  #                    # color = c("#999999" ) 
  #                     ),
  #                 
  #                 size=1,
  #                #  force        = 0.5,
  #                  direction    = "both",
  #                  hjust        = 0,
  #                  segment.size = 0.2,
  #                 box.padding = unit(0.35, "lines"),
  #                 point.padding = unit(0.3, "lines"),
  #                 ) +
  # geom_label_repel(
  #   size = 1,
  #   #arrow = arrow(length = unit(0.03, "npc"), type = "closed", ends = "first"),
  #   direction = "both",
  #   label.size = 0.01,
  #   segment.size = 0.08
  #     ) +
  #geom_point(aes(col=lipclass, size=AveExpr)) +
  scale_size_continuous(range = c(0.25, 2.5), breaks=c(10, 20, 30))














EnhancedVolcano(comparison, lab = rownames(comparison), x = "logFC", y = "adj.P.Val", FCcutoff = 1,
                     transcriptLabSize = 2.0,
                     pCutoff = 0.05,
                     ylab = bquote(~-Log[10]~adjusted~italic(P)),
                     legendPosition = 'bottom',
                     legendLabSize = 10,
                     subtitle = NULL,
                     ylim= c(0,  max(-log10(comparison$adj.P.Val), na.rm=TRUE) + 2), 
               
                #shape = c(6, 4, 2, 11),
                #colCustom = keyvals.col,
                colAlpha = 1,
               
                #drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'grey50',
                gridlines.major = TRUE,
                gridlines.minor = FALSE,
                border = 'partial',
                borderWidth = 1.5,
                borderColour = 'black')

# cols <- brewer.pal(6, "Set2")
# keyvals.col <- rep(cols[1], nrow(input_lipids))
# 
# names(keyvals.col) <- rownames(input_lipids)
# 
# keyvals.col[rownames(input_lipids)[which(input_lipids$lipclass %in% "Glycerophospholipids")]] <- cols[1]
# names(keyvals.col)[which(input_lipids$lipclass %in% "Glycerophospholipids")] <- rownames(input_lipids)[which(input_lipids$lipclass %in% "Glycerophospholipids")]
# 
# keyvals.col[rownames(input_lipids)[which(input_lipids$lipclass %in% "Neutral lipids" )]] <- cols[2]
# names(keyvals.col)[which(input_lipids$lipclass %in% "Neutral lipids")] <- rownames(input_lipids)[which(input_lipids$lipclass %in% "Neutral lipids")]
# 
# keyvals.col[rownames(input_lipids)[which(input_lipids$lipclass %in% "Sphingolipids" )]] <- cols[3]
# names(keyvals.col)[which(input_lipids$lipclass %in% "Sphingolipids")] <- rownames(input_lipids)[which(input_lipids$lipclass %in% "Sphingolipids")]
# 
# keyvals.col[rownames(input_lipids)[which(input_lipids$lipclass %in% "Other lipids" )]] <- cols[5]
# names(keyvals.col)[which(input_lipids$lipclass %in% "Other lipids" )] <- rownames(input_lipids)[which(input_lipids$lipclass %in% "Other lipids" )]
# 
# keyvals.col[rownames(input_lipids)[which(input_lipids$lipclass %in% "n.s" )]] <- cols[6]
# names(keyvals.col)[which(input_lipids$lipclass %in% "n.s")] <- rownames(input_lipids)[which(input_lipids$lipclass %in% "n.s")]
# 
# keyvals.col[rownames(input_lipids)[which(input_lipids$lipclass %in% "Sterols" )]] <- cols[4]
# names(keyvals.col)[which(input_lipids$lipclass %in% "Sterols")] <- rownames(input_lipids)[which(input_lipids$lipclass %in% "Sterols")]
# 
# 
# filtered_data <- imputated_lipids
# #contmatrix
# # fit2 



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


path.class <- c('fatty acids', 'G3P', 'LCBs', 'Cer', 'SM', 'LPA', 'PA', 'DAG', 'TAG', 'PC',
                'LPC', 'PE', 'LPE', 'CDP-DAG', 'PI', 'LPI', 'PG', 'LPG', 'PS', 'LPS')
control_graph <- buildPathWay(path.class)
library("igraph")
library("visNetwork")
nodes <- data.frame(nodes = path.class, log2fc = rep(length(path.class), 1))
edges <- data.frame(from=c("fatty acids", "LCBs", "Cer", "fatty acids", "G3P", "LPA", "PA", 
                           "PA", "DAG", "DAG", "DAG","PC", "PE", "PS","CDP-DAG", "CDP-DAG",
                           "CDP-DAG","PI", "PG", "PS", "LPC","LPS", "LPG", "LPI", "PC"),
                    to=c("LCBs", "Cer", "SM", "LPA","LPA", "PA", "DAG", "CDP-DAG", "TAG", 
                         "PC", "PE", "LPC", "LPE", "PE", "PI", "PG", "PS", "LPI", "LPG", 
                         "LPS", "PC", "PS", "PG", "PI", "PE"))

g <- graph_from_data_frame(edges, directed = TRUE, vertices = nodes, remove.multiple=FALSE)

visIgraph(control_graph)

gg <- simplify(control_graph, remove.multiple = FALSE)


g <- plot(control_graph,
     remove.loops= FALSE, edge.curved=TRUE)
visIgraph(g)

initial_path <- createNetworkFromIgraph(control_graph,"Igraph")



message("Please type the control group and contrast groups for making pathway plot")
fc_data <- TotalNormalClass(foldchange_input, 1)
fc_data_long <- fc_data[[2]][[1]]
fc_class <- fc_data_long %>% spread(., key = experiment.group, value = value)

control.group <- fc_data_long$experiment.group[1]


exist.class <- intersect(fc_class$class, path.class)
detect.lipid <- fc_class %>% filter(class %in% exist.class)
empty.lipid <- path.class[!path.class %in% exist.class]
empty.class <- sapply(group_names, function(x) ifelse(x == control.group, return(rep(1, length(empty.lipid))), return(rep(0, length(empty.lipid)))))
empty.data <- data.frame(class=empty.lipid, empty.class)
path.data <- rbind(detect.lipid, empty.data)
write_csv(path.data, "data/path.csv")

contrast_group <- readline("Please input the contrast group name: ")
if(!contrast_group %in% group_names ){
  contrast_group <- readline("The group is not in the experiment. Please input the  name again: ")
}

size <- as.numeric(unlist(path.data[, contrast_group]))

df <- data.frame(nodes=path.data$class, log2fc= size, row.names = path.data$class, stringsAsFactors = FALSE)
loadTableData(df)
setNodeSizeMapping('log2fc', sizes=size)

df
# 
createNetworkFromIgraph()