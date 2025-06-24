#######################################################################################################################
#######################################################################################################################
################################################### June 12th, 2025 ###################################################
############################################# Author: Stefany Ayala-Montaño ###########################################
#######################################################################################################################
##################### This is a script that addresses contamination based on coverage, prevalence, ####################
##################### and 'suspected' contamination status based on microbiological reports. ##########################
########## It processes kraken-bracken biom format files provided by the user in '/$PATH/input_files/biom/' ###########
#######################################################################################################################
#######################################################################################################################

# Load necessary libraries for data manipulation, visualization, and data import/export

# Core tidyverse (automatically loads ggplot2, dplyr, tidyr, readr, purrr, etc.)
library(tidyverse)

# Additional tidyverse-friendly packages
library(magrittr)
library(ggthemes)
library(scales)
library(readxl)
library(writexl)

# File/path/config handling
library(here)   # for file paths
library(yaml)   # for config files

# Bioinformatics/data
library(phyloseq)
library(biomformat)
library(stringr) 


# Load config
params <- yaml::read_yaml(here("config.yaml"))

#call the 'functions' file
source(here::here("R", "functions.R"))

# Define base paths
data_dir <- here("input_files") 
output_dir <- here("ouput_files")

########## here to list all input files
input_path <- list(
  bracken_combined_select = here("input_files", "bracken_combine_select_2023_2024.tsv"), #bracken_combined_select.tsv
  negative = here("input_files", "Microbiologically_negative.xlsx"),
  microbiology = here("input_files", "Access_TAPIR_250325.xlsx"),
  coverage = here("input_files", "bracken_plus_coverage_v6.csv"),
  biom_folder = here("input_files", "biom"),
  tax_id = here("input_files", "taxid_from_cov_file.xlsx")
)

########## here to list all output files
output_path <- list(
  biom_go_folder = here("output_files", "biom_decontaminated"),
  raw_classification_system = here("output_files", "Samples_system_raw.xlsx"),
  samples_to_Go = here("output_files", "Samples_in_Go_after_Check.xlsx"),
  samples_to_Stop = here("output_files", "Samples_in_Stop_after_Check.xlsx")
)

###### Create folders 

if (!dir.exists(input_path$biom_folder)) {
  dir.create(input_path$biom_folder, recursive = TRUE)
}

if (!dir.exists(output_path$biom_go_folder)) {
  dir.create(output_path$biom_go_folder, recursive = TRUE)
}


################# Load input files

cov_file <- read.csv(file.path(input_path$coverage))
negative <- read_excel(file.path(input_path$negative))
bracken_combined <- read.delim(input_path$bracken_combined_select) 
microbiology <- read_excel(file.path(input_path$microbiology))
tax_id <- read_excel(file.path(input_path$tax_id))


# List all files in the "data/biom" folder
biom <- list.files(file.path(input_path$biom_folder), full.names = TRUE)

# how many files 
num_files <- length(biom)

# Print a message
message(paste0("📁 The user has provided ", num_files, " files.\n🗂️  File paths: ", paste(biom, collapse= ",\n ")))

###### weeks to process based on in the biom files provided in '/$PATH/input_files/biom/'

input_weeks <- stringr::str_extract(biom, "W\\d{2}") |> 
  stringr::str_remove("W") |> 
  as.integer()

# Print a message
message(paste0("📁 The weeks for downstream analyses are: ",  paste(input_weeks, collapse= "," )))

############################## bracken_parse_df function present in functions.R

bracken_select_2023_2024_modified_total_df <- bracken_parse_df(bracken_combined)

# Extract week number from Filename and append to Sample identifier

bracken_select_2023_2024_modified_total_df <- bracken_select_2023_2024_modified_total_df %>%
  mutate(week = ifelse(grepl("^[^-]+-[0-9]+[-_].*", Filename), 
                       as.numeric(gsub("^[^-]+-([0-9]+)[-_].*", "\\1", Filename)), 
                       NA))  %>% 
  mutate(Sample = paste0(Sample, '-', week)) %>%
  mutate(type = gsub("-[0-9]+", "", Sample)) %>%
  mutate(Type = case_when(type=="Swab" ~ "Swab",
                          type=="Water" ~ "Water",
                          TRUE ~ "Sample")) %>%
  mutate(Clase = case_when(str_detect(Filename, '[Ss]wab') ~ 'Swab',
                           str_detect(Filename, '[Ww]ater') ~ 'Water',
                           str_detect(Filename, 'A') ~ 'Anus',
                           str_detect(Filename, 'B|N') ~ 'Nose',
                           TRUE ~ 'Check')) %>%
  mutate(Read_category=case_when(Read_category == '**' ~ '> 10K',
                                 Read_category == '*' ~ '> 1K',
                                 TRUE ~ "< 1K"                              
  ))  %>% 
  mutate(name = as.character(name),   # Ensure name is character, not factor
         genus = sapply(strsplit(name, " "), `[`, 1)) %>%  # Extract genus
  filter(genus != "Others" & name != "unassigned" & new_est_reads > 0) %>%        # Filter low-quality
  mutate(sample_ID = sub("_.*", "", Filename))         # Simplify ID


# Assign contamination reasons based on read category and microbiological report
# > 10K means: higher than 10.000 reads; > 1K means: higher than 1.000 reads
# Site 'Neg_Ctrl' refers to the control (swab/water) from where we do not expect high number of reads


df_estrato_exclCriteria_file <- bracken_select_2023_2024_modified_total_df %>%
  mutate(microbio = case_when(
    sample_ID %in% negative ~ "negative",  # If sample ID is in 'negative' list → "negative"
    TRUE ~ "positive"                      # All other samples → "positive"
  )) %>%
mutate(cont_reason = case_when(
  Site == "Neg_Ctrl" & microbio == "positive" ~ "Cross-contamination",
  Site == "Neg_Ctrl" & microbio == "negative" & Read_category == "> 10K" ~ "Library-contamination-10K",
  Site == "Neg_Ctrl" & microbio == "negative" & Read_category == "> 1K" ~ "Library-contamination-1K",
  TRUE ~ ""  # No contamination reason
)) %>% 
  filter(week %in% input_weeks) ### subset of the weeks in process

# Separate control and biological sample types
controls <- df_estrato_exclCriteria_file %>% #bracken_exclCriteria
  filter(Clase %in% c("Water", "Swab"))  # Control samples

samples <- df_estrato_exclCriteria_file  %>% #bracken_exclCriteria
  filter(Clase %in% c("Anus", "Nose"))  # Biological samples

# Join controls and samples on 'name' and 'week' to assess contamination 
result <- samples %>%
  left_join(controls, by = c("name", "week"), suffix = c("_sample", "_control"), relationship = "many-to-many") %>%
  mutate(difference_fraction_reads = fraction_total_reads_control - fraction_total_reads_sample)  # Difference metric

##### subset of the coverage file for the weeks in process
select_cov_file <- cov_file %>% select(Filename, name, av_read_length, av_cov, week) %>% 
  filter(week %in% input_weeks)

# merge result with coverage data
threshold <- result %>%
  mutate(genus = sapply(strsplit(name, " "), `[`, 1)) %>%  # add genus from taxon name
  left_join(select_cov_file, by = c("Filename_sample" = "Filename", "name"), relationship = "many-to-many")

# Select relevant columns for downstream classification
df_threshold <- threshold %>% select(
  Filename_sample, name, genus, new_est_reads_sample, week.x, av_read_length, 
  Clase_sample, Site_sample, av_cov, fraction_total_reads_sample,
  fraction_total_reads_control, difference_fraction_reads, microbio_sample,
  microbio_control, cont_reason_control)

#######################################################################################################################
################################# classification of the system 'Stop-Check-Go' ########################################
#######################################################################################################################


df_system_definition <- df_threshold %>% mutate(coverage = case_when(av_cov < 20 ~ "low",
                                                             av_cov >= 20 ~ "ok",
                               ### low quality reads not assembled with FLYE will have no coverage
                                                             TRUE ~ "no coverage")) %>% 
  mutate(diff_fr = case_when(difference_fraction_reads < 0 | is.na(difference_fraction_reads) ~ "go",
                             difference_fraction_reads > 0 ~ "stop",
                             TRUE ~ "empty_dfr"))  %>%  mutate(system = case_when(
                               coverage=="low" & diff_fr=="stop" ~ "stop",
                               microbio_sample=="negative"  ~ "stop",
                               coverage=="low" & diff_fr=="go" & is.na(cont_reason_control) ~ "check", 
                               coverage=="low" & diff_fr=="go" & cont_reason_control=="Cross-contamination" ~ "stop",
                               coverage=="low" &  diff_fr=="go" & cont_reason_control=="Library-contamination-1K" ~ "check",
                               coverage=="low" &  diff_fr=="go" & cont_reason_control=="Library-contamination-10K" ~ "stop",
                               coverage=="ok" & diff_fr=="stop" ~ "stop",
                               coverage=="ok" & diff_fr=="go" & is.na(cont_reason_control) ~ "go",
                               coverage=="ok" & diff_fr=="go" & cont_reason_control=="Cross-contamination" ~ "check",
                               coverage=="ok" &  diff_fr=="go" & cont_reason_control=="Library-contamination-1K" ~ "go",
                               coverage=="ok" &  diff_fr=="go" & cont_reason_control=="Library-contamination-10K" ~ "check",
                              ###### reads without coverage have no system 'empty' 
                               TRUE ~ "empty_system")) %>% filter(name!="unassigned") ### exclude unassigned reads

write_xlsx(df_system_definition, output_path$raw_classification_system)

#######################################################################################################################
########################### In this section the samples in the category 'Check' are evaluated #########################
#######################################################################################################################

##################################### There are four cases based on the variables: #################################### 
###################### A. Coverage (low, go), B. difference of reads fraction -dfr- (stop, go), #######################
############ C. contamination reason (no reason 'NA', Cross-contamination, Library-contamination-1k/10k) ##############
########## We change samples from 'Check' to 'Go' if there is a positive-growth culture (Enterobacterales) ############
################ We change samples from 'Check' to 'Stop' if there is a negative-growth culture (Enterob.) ############
#######################################################################################################################

############## this file contains the reports from microbiology

df_microbiology <- microbiology %>% 
  mutate(genus=sub("_.*", "", name)) %>%
  select(sample_ID, genus, name)

############## system categories 

go <- df_system_definition %>% subset(system=="go")
check <- df_system_definition %>% subset(system=="check")
stop <- df_system_definition %>% subset(system=="stop")


###################################################### FIRST CASE #####################################################
######## coverage 'low', difference of reads fraction 'go' and no contamination reason

check_clow_frgo_nacont <- check %>% filter(is.na(cont_reason_control))
check_clow_frgo_nacont_enterob <- check_clow_frgo_nacont %>% mutate(sample_ID = sub("_.*", "", Filename_sample)) %>%
  filter(genus=="Escherichia" | genus=="Klebsiella" | genus=="Enterobacter" 
         | genus=="Citrobacter" | genus =="Acinetobacter")

check_clow_frgo_nacont_enterob2 <- left_join(check_clow_frgo_nacont_enterob, df_microbiology, by=c("sample_ID", "name"))

###### To 'Go' after 'Check'
check_clow_frgo_nacont_pass <- check_clow_frgo_nacont_enterob2 %>% filter(!is.na(genus.y))

## To 'Stop' after 'Check'

check_clow_frgo_nacont_stop <- check_clow_frgo_nacont_enterob2 %>% filter(is.na(genus.y))%>%mutate(system="stop after check")


################################################### SECOND CASE #######################################################
######## covergae 'low', dfr 'go', contamitation reason:'Library-contamination-1K'

check_clow_frgo_1k <- check %>% filter(cont_reason_control=="Library-contamination-1K") %>% mutate(sample_ID = sub("_.*", "", Filename_sample)) %>%
  filter(genus=="Escherichia" | genus=="Klebsiella" | genus=="Enterobacter" 
         | genus=="Citrobacter" | genus =="Acinetobacter")

check_clow_frgo_1k2 <- left_join(check_clow_frgo_1k, df_microbiology, by=c("sample_ID", "name"))

###### To 'Go' after 'Check' 
check_clow_frgo_1k_pass <- check_clow_frgo_1k2 %>% filter(!is.na(genus.y))

## To 'Stop' after 'Check'
check_clow_frgo_1k_stop  <- check_clow_frgo_1k2 %>% filter(is.na(genus.y))%>%mutate(system="stop after check")


################################################### THIRD CASE ########################################################
######## coverage 'ok', dfr 'go', contamination reason: 'Cross-contamination' >>> stop cross-contaminants <<<

check_cgo_frgo_cc <-  check %>% filter(cont_reason_control=="Cross-contamination") %>% 
  mutate(id_cc=case_when(week.x==2 & genus=="Staphylococcus" ~ "match", ## in this week the lab report at genus level
                         week.x==2 &  name=="Enterococcus faecalis" ~ "match",
                         week.x==4 &  name=="Staphylococcus aureus" ~ "match",
                         week.x==4 &  name=="Staphylococcus hominis" ~ "match",
                         week.x==12 &  name=="Escherichia coli" ~ "match",
                         week.x==13 &  name=="Staphylococcus epidermidis" ~ "match",
                         week.x==14 &  name=="Enterococcus faecalis" ~ "match",
                         week.x==15 &  name=="Klebsiella oxytoca" ~ "match",
                         week.x==16 &  name=="Klebsiella oxytoca" ~ "match",
                         week.x==16 &  name=="Staphylococcus aureus" ~ "match",
                         week.x==16 &  name=="Staphylococcus epidermidis" ~ "match",
                         week.x==20 &  name=="Staphylococcus epidermidis" ~ "match",
                         week.x==20 &  name=="Staphylococcus hominis" ~ "match",
                         week.x==29 &  name=="Enterococcus epidermidis" ~ "match",
                         week.x==48 &  name=="Staphylococcus hominis" ~ "match",
                         TRUE ~ "not match"))       

###### To 'Go' after 'Check' 
check_cgo_frgo_cc_pass <- check_cgo_frgo_cc %>% filter(id_cc=="not match") %>%
  filter(genus=="Escherichia" | genus=="Klebsiella" | genus=="Enterobacter" 
         | genus=="Citrobacter" | genus =="Acinetobacter")%>% mutate(system="go after check")

check_cgo_frgo_cc_pass2 <- check_cgo_frgo_cc_pass %>% mutate(sample_ID = sub("_.*", "", Filename_sample))
check_cgo_frgo_cc_pass2$id_cc <- NULL

## To 'Stop' after 'Check'

check_cgo_frgo_cc_stop <- check_cgo_frgo_cc %>% filter(id_cc=="match") %>%
  filter(genus=="Escherichia" | genus=="Klebsiella" | genus=="Enterobacter" 
         | genus=="Citrobacter" | genus =="Acinetobacter")%>%
  mutate(system="stop after check")%>%
  mutate(sample_ID = sub("_.*", "", Filename_sample))

check_cgo_frgo_cc_stop$id_cc <- NULL

#################################################### FOURTH CASE ######################################################

########## coverage 'ok', dfr 'go' and contamination reason: 'Library-contamination-10K'

check_clow_frgo_10k <- check %>% filter(cont_reason_control=="Library-contamination-10K") %>% mutate(sample_ID = sub("_.*", "", Filename_sample)) %>%
  filter(genus=="Escherichia" | genus=="Klebsiella" | genus=="Enterobacter" 
         | genus=="Citrobacter" | genus =="Acinetobacter")

check_clow_frgo_10k_enterob <- left_join(check_clow_frgo_10k, df_microbiology, by=c("sample_ID", "name"))

###### To 'Go' after 'Check' 
check_clow_frgo_10k_pass <- check_clow_frgo_10k_enterob %>% filter(!is.na(genus.y))

## To 'Stop' after 'Check'
check_clow_frgo_10k_stop <- check_clow_frgo_10k_enterob %>% filter(is.na(genus.y))%>% mutate(system="stop after check")


######################################## bind the four conditions (to 'Go') ###########################################

cond_one <- check_clow_frgo_nacont_pass
cond_one$genus.y <- NULL 
cond_one_df <- cond_one %>% dplyr::rename(genus= genus.x)

cond_two <- check_clow_frgo_1k_pass
cond_two$genus.y <- NULL 
cond_two_df <- cond_two %>% dplyr::rename(genus= genus.x)

cond_three <- check_clow_frgo_10k_pass
cond_three$genus.y <- NULL 
cond_three_df <- cond_three %>% dplyr::rename(genus= genus.x)

go_after_check <- rbind(cond_one_df,
                               cond_two_df,
                               check_cgo_frgo_cc_pass2,
                               cond_three_df) 

######### bind the system  to 'Go' after review to 'Check'
df_go <- go %>% mutate(sample_ID = sub("_.*", "", Filename_sample))
go_plus_conditions_in_check <- rbind(df_go, go_after_check)


### create the output file for future reference
write_xlsx(go_plus_conditions_in_check , output_path$samples_to_Go)

######################################## bind the four conditions (to 'Stop') ###########################################

cond_one_stop <- check_clow_frgo_nacont_stop
cond_one_stop$genus.y <- NULL 
cond_one_df_stop <- cond_one_stop %>% dplyr::rename(genus= genus.x)

cond_two_stop <- check_clow_frgo_1k_stop
cond_two_stop$genus.y <- NULL 
cond_two_df_stop <- cond_two_stop %>% dplyr::rename(genus= genus.x)

cond_three_stop <- check_clow_frgo_10k_stop
cond_three_stop$genus.y <- NULL 
cond_three_df_stop <- cond_three_stop %>% dplyr::rename(genus= genus.x)

stop_after_check <- rbind(cond_one_df_stop ,
                                 cond_two_df_stop ,
                                 check_cgo_frgo_cc_stop,
                                 cond_three_df_stop )

######### bind the system  to 'Stop' after review to 'Check'

df_stop <- stop %>% mutate(sample_ID = sub("_.*", "", Filename_sample))
stop_plus_conditions_in_check <- rbind(df_stop, stop_after_check)

### create the output file for future reference

write_xlsx(stop_plus_conditions_in_check, output_path$samples_to_Stop)


#########################################################################################################################
########################## Call the functions to process the bioligical data  (biom format files) #######################
#########################################################################################################################


# List all .biom files
biom_files <- biom

# Initialize list to store the biom objects
biom_list <- list()

# Loop through the biom files
for (file in biom_files) {
  # Extract week code from the filename 
  week <- tolower(stringr::str_extract(file, "W\\d{2}"))
  
  if (!is.na(week)) {
    # Construct the object name, e.g., "w06_biom"
    object_name <- paste0(week, "_biom")
    
    # Import the biom file
    biom_obj <- import_biom_data(file)
    
    # Assign it to the global environment 
    assign(object_name, biom_obj, envir = .GlobalEnv)
    
    # Save into the list with the same name
    biom_list[[object_name]] <- biom_obj
    
    message(paste("✔ Imported and saved to list as:", object_name))
  } else {
    warning(paste("❌ Could not process", basename(file)))
  }
}

#names(biom_list)
####### negative-growth IDs 

negative_list <- negative$sample_ID

########## data wrangling of biom files with microViz R package

for (week_name in names(biom_list)) { 
  
  
  # Assign biom object from list
  input_biom <- biom_list[[week_name]]
  
  # Desired output name
  output_name <- sub("_biom$", "_biom_m", week_name)
  
  # process itself
  processed <- process_biom_data(input_biom)
  
  # Assign it to the global environment 
  assign(output_name, processed, envir = .GlobalEnv)
  
  message(paste("✔ Processed and saved:", output_name))
}


# List all processed biom objects in the environment
list_biom_processed <- ls(pattern = "^w\\d{2}_biom_m")

weeks_to_process <- as.numeric(stringr::str_extract(list_biom_processed, "\\d{2}"))

###### Create the dataframe of contaminants

for (week in weeks_to_process) {
  week_str <- sprintf("%02d", week)  # Format week as two digits (e.g., 06, 11)
  result <- create_go_dataset(week)  # Generate dataset for the current week
  assign(paste0("for_go_system_w", week_str), result, envir = .GlobalEnv)  # Save as variable in global environment
}

# List all variables in the environment that match the naming pattern for the 'Go' datasets

databases <- ls(pattern = "^for_go_system_w")

# print message
message(paste0("📁 These are the databases with contaminants: ",  paste(databases, collapse= ", " )))


######## create the otu table for the biom files 

for (name in list_biom_processed) {
  if (exists(name)) {
    obj <- get(name)
    otu_obj <- otu_table(obj)
    
    # Desired output name for the file
    new_name <- sub("_biom_m$", "_otu_table", name)
    
    # Assign it to the global environment 
    assign(new_name, otu_obj, envir = .GlobalEnv)
    print(new_name)
  } else {
    warning(paste("Object", name, "not found"))
  }
}

#########################################################################################################################
############## Main function: function that assigns zero to the positions where it can be contamination #################
#### the files are not overwritten, so the user can compare both biom files 'w*_biom' or 'w*_biom_m' WITH 'w*_biom_go' ##
########## the difference between w*_biom' and 'w*_biom_m' is the number of variables (1 vs. 8, respectivelly) ##########
##################################### thanks to the 'process_biom_data' function ########################################

for (name in list_biom_processed) {
  
  # Extract the 2-digit week from the name 
  week_str <- str_extract(name, "w\\d{2}") |> str_remove("w")
  
  # Build object names
  otu_name <- paste0("w", week_str, "_otu_table")
  df_name  <- paste0("for_go_system_w", week_str)
  biom_name <- paste0("w", week_str, "_biom_m")
  output_name <- paste0("w", week_str, "_biom_go")
  
  
  # Check which objects exist
  missing <- c() ## empty list
  if (!exists(otu_name)) missing <- c(missing, otu_name)
  if (!exists(df_name)) missing <- c(missing, df_name)
  if (!exists(biom_name)) missing <- c(missing, biom_name)
  if (!exists(otu_name)) missing <- c(missing, otu_name)

  
  if (length(missing) == 0) {
    
    # All objects exist – proceed
    otu_matrix <- get(otu_name)
    df <- get(df_name)
    biom_obj <- get(biom_name)
    
    # Clean sample names
    ## if the pattern is not found, it does nothing
    sample_names(biom_obj) <- gsub(
      pattern = ".S.bracken.kraken.report",
      replacement = "",
      x = sample_names(biom_obj)
    )
    
    #Clean otu_table names
    sample_names(otu_matrix) <- gsub(
      pattern = ".S.bracken.kraken.report",
      replacement = "",
      x = sample_names(otu_matrix) 
    )
    
    ## function
    otu_matrix_modified <- assign_zero_to_positions(df, otu_matrix)
    otu_table_modified <- otu_table(otu_matrix_modified, taxa_are_rows = TRUE)
    
    biom_go <- phyloseq(
      otu_table_modified,
      sample_data(biom_obj),
      tax_table(biom_obj)
    )
    
    # Assign it to the global environment 
    assign(output_name, biom_go, envir = .GlobalEnv)
    
    output_file <- file.path(output_path$biom_go_folder, paste0(output_name, ".rds"))
    saveRDS(biom_go, output_file)
    
    message(paste("✔ Week", week_str, "processed and saved as", output_name))
    
  } else {
    #warning(paste("❌ Missing objects for week", week_str))
    warning(paste("❌ Missing objects for week", week_str, ":", paste(missing, collapse = ", ")))
  }
}

##########################################################################################################################
##########################################################################################################################
################################################## End of the script #####################################################
##########################################################################################################################
##########################################################################################################################