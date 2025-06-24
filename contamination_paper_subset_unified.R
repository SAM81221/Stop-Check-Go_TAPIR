
###################################### June 12th, 2025
################################### Author: Stefany Ayala-Montaño
######################################### This code has been commented with the use of ChatGPT (GPT-4.5-turbo model)
#########################################. prompt: This is a script that addresses contamination based on coverage, prevalence, and contamination status from microbiological reports. Comment out this script to make it easier to understand by a third party. Do not modify in any sense the code:  


# Load necessary libraries for data manipulation, visualization, and Excel

library(readr)
library(tidyr)
library(tidyverse)
library(readxl)
library(ggplot2)
library(writexl)
library(here) # use here() to build the file paths
library(yaml) # config. file with coverage threshold value and contamination definitions


# Load config
params <- yaml::read_yaml(here("config.yaml"))

#########################################################
#call the 'functions' file
source(here::here("R", "functions.R"))

# Define base paths
data_dir <- here("input_files") # or wherever your inputs live
output_dir <- here("ouput_files")

########## here to list all input files
input_path <- list(
  bracken_combined_select = here("input_files", "bracken_combine_select_2023_2024_stef_df_250124_woG230616_G230828.tsv"), #bracken_combined_select.tsv
  negative = here("input_files", "Microbiologically_negative.xlsx"),
#  negative_controls = here("input_files", "Negative_controls.xlsx"),
 # estrato_file = here("data", "estrato_file.csv"),
  coverage = here("input_files", "bracken_plus_coverage_v6.csv"),
  biom_folder = here("input_files", "biom"),
  tax_id = here("input_files", "taxid_from_cov_file_plus_manual108.xlsx")
)

########## here to list all output files
output_path <- list(
  biom_go_folder = here("output_files", "biom_decontaminated"),
  samples_to_Go = here("output_files", "Samples_in_Go_system.xlsx"),
  samples_to_Stop = here("output_files", "Samples_in_Stop_after_Check.xlsx")
)

# 2. Create folders 

if (!dir.exists(input_path$biom_folder)) {
  dir.create(input_path$biom_folder, recursive = TRUE)
}

if (!dir.exists(output_path$biom_go_folder)) {
  dir.create(output_path$biom_go_folder, recursive = TRUE)
}


################# Load files

# Use paths
cov_file <- read.csv(file.path(input_path$coverage))



#################

write_xlsx(df_threshold, file.path(output_dir, "Stop-Check-Go_system.xlsx"))
# Use them like this:
df_estrato_file <- read.csv(input_path$estrato_file)
write_xlsx(df_threshold, output_path$system_output)

write_xlsx(df_estrato_exclCriteria_file, output_path$df_excl)

df_estrato_file <- read.csv(input_path$estrato_file)

files <- list.files("Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/",
                    pattern = "combined_bracken_S_.*\\.txt$", full.names = TRUE)
bracken_data <- purrr::map_dfr(files, readr::read_tsv)
#####
#######################################
# weeks randomly selected (12): 6, 11, 20, 26, 27, 29, 34, 39, 41, 43, 50, 52

# Read in Bracken abundance reports 
# Files are named by sequencing batch and date codes

#bracken_G230119 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G230119.txt')
#bracken_G230126 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_GG230126.txt')
#bracken_G230202 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_GG230202.txt')
#bracken_G230208 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G230208.txt')
bracken_G230216 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G230216.txt')
#bracken_G230302 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G230302.txt')
#bracken_G230309 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G230309.txt')
#bracken_G230317 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G230317.txt')
bracken_G230324 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G230324.txt')
# bracken_G230332 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G230332.txt')
# bracken_G230406 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G230406.txt')
# bracken_G230414 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G230414.txt')
# bracken_G230428 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G230428.txt')
# bracken_G230505 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G230505.txt')
# bracken_G230511 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G230511.txt')
#bracken_G230519 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G230519.txt')
bracken_G230606 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G230606.txt')
#bracken_G230610 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G230610.txt')
# bracken_G230617 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G2306')
#bracken_G230626 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G230626.txt')
#bracken_G230632 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G230632.txt')
bracken_G230708 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G230708.txt')
bracken_G230713 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G230713.txt')
#bracken_G230826 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G230826.txt')
# bracken_G230828 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G230828.txt')
# bracken_G230902 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G230902.txt')
# bracken_G230908 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G230908.txt')
# bracken_G230917 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G230917.txt')
bracken_G240217 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G240217.txt')
# bracken_G240123 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G240123.txt')
# bracken_G240209 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G240209.txt')
bracken_G240222 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G240222.txt')
# 
# 
# bracken_G230420 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G230420.txt')
# bracken_G230616 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/G230616_combined_parsed_bracken.tsv') %>% 
#   select(1:8) %>% filter(!name %in% c("Others", "unassigned"))
# bracken_G230223 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/G230223_combined_parsed_bracken.tsv') %>% 
#   select(1:8) %>% filter(!name %in% c("Others", "unassigned"))
# bracken_G230607 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/G230606_week21_72hr_combined_parsed_bracken.tsv') %>% 
#   select(1:8) %>% filter(!name %in% c("Others", "unassigned"))
# bracken_G240424 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G240424.txt')
# bracken_G240424 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G240424.txt')
# bracken_G240515 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G240515.txt')
bracken_G240516 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G240516.txt')
# bracken_G240522 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G240522.txt')
bracken_G240528 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G240528.txt')
# bracken_G240612 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G240612.txt')
bracken_G240613 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G240613.txt')
# bracken_G240621 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G240621.txt')
# bracken_G240622 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G240622.txt')
# bracken_G240703 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G240703.txt')
# bracken_G240704 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G240704.txt')
# bracken_G240705 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G240705.txt')
# bracken_G240711 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G240711.txt')
# bracken_G240712 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G240712.txt')
 bracken_G240723 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G240723.txt')
# bracken_G240724 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G240724.txt')
 bracken_G240725 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G240725.txt')
# bracken_G240806 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G240806.txt')
# bracken_G240807 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G240807.txt')
# bracken_G240808 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G240808.txt')
# bracken_G240823 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G240823.txt')
# bracken_G240827 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G240827.txt')
# bracken_G240830 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G240830.txt')
# bracken_G240910 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G240910.txt')
# bracken_G240912 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G240912.txt')
# bracken_G240913 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G240913.txt')
# bracken_G240927 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G240927.txt')

 # Combine selected Bracken reports into one dataframe 
 
 bracken_combine_select_2023_2024_stef_df <- bind_rows(bracken_G230119, bracken_G230126, bracken_G230202, bracken_G230208,
                                                      bracken_G230216, bracken_G230223, bracken_G230302, bracken_G230309, 
                                                      bracken_G230317, bracken_G230324, bracken_G230332, bracken_G230406, 
                                                      bracken_G230414, bracken_G230420, bracken_G230428, bracken_G230505,   
                                                      bracken_G230511, bracken_G230519, bracken_G230606, bracken_G230607,
                                                      bracken_G230610, bracken_G230616, bracken_G230626, bracken_G230632,  
                                                      bracken_G230708, bracken_G230713, bracken_G230828, bracken_G240217, 
                                                      bracken_G240123, bracken_G240209, bracken_G230902, bracken_G230826, 
                                                      bracken_G240222, bracken_G240424, bracken_G230908, bracken_G230917,
                                                      bracken_G240515, bracken_G240516, bracken_G240522, bracken_G240528,
                                                      bracken_G240612, bracken_G240613, bracken_G240621, bracken_G240622,
                                                      bracken_G240703, bracken_G240704, bracken_G240705, bracken_G240711,
                                                      bracken_G240712, bracken_G240723, bracken_G240724, bracken_G240725,
                                                      bracken_G240806, bracken_G240807, bracken_G240808, bracken_G240823,
                                                      bracken_G240827, bracken_G240830, bracken_G240910, bracken_G240912,
                                                      bracken_G240913, bracken_G240927) %>% 
  filter(!str_detect(name, '       ')) # Remove entries with blank names
 
 # Save the combined Bracken data to a TSV file
 write.table(bracken_combine_select_2023_2024_stef_df, 
             "bracken_combine_select_2025.tsv", 
             sep = '\t', row.names = FALSE)
 
 # Load the file locally (used for further parsing/processing)

bracken_combined <- read.delim(input_path$bracken_combined_select) 


############################## bracken_parse_df function is in other R script called functions.R

bracken_select_2023_2024_modified_total_df <- bracken_parse_df(bracken_combined)#%>%

# Extract week number from Filename and append to Sample identifier

bracken_select_2023_2024_modified_total_df <- bracken_select_2023_2024_modified_total_df %>%
  mutate(week = ifelse(grepl("^[^-]+-[0-9]+[-_].*", Filename), 
                       as.numeric(gsub("^[^-]+-([0-9]+)[-_].*", "\\1", Filename)), 
                       NA))  %>% 
  mutate(Sample = paste0(Sample, '-', week)) 

##### subset

# Filter dataset to only include samples from randomly selected weeks
bracken_select_2023_2024_modified_ss_df <- bracken_select_2023_2024_modified_total_df %>% filter(week %in% subset_weeks)
#######


###subset 
# Get unique filenames from subset, in ascending order by week

sample_name_elements_ss <- bracken_select_2023_2024_modified_ss_df %>%
  # filter(grepl('Swab', Filename) | grepl('Water', Filename)) %>% 
  arrange((week))%>%   ######REMOVE DESC TO HAVE THE NORMAL ORDER  # Arrange first by prefix, then by numeric part
  distinct(Filename) %>% 
  pull() %>% 
  unique() 
#######

# Create a working dataframe with sample classifications and simplified types
df_ss <- bracken_select_2023_2024_modified_ss_df %>%filter(Filename %in% sample_name_elements_ss) %>%
  mutate(type = gsub("-[0-9]+", "", Sample)) %>%
  mutate(Type = case_when(type=="Swab" ~ "Swab",
                          type=="Water" ~ "Water",
                          TRUE ~ "Sample")) %>%
  mutate(Clase = case_when(str_detect(Filename, '[Ss]wab') ~ 'Swab',
                           str_detect(Filename, '[Ww]ater') ~ 'Water',
                           str_detect(Filename, 'A') ~ 'Anus',
                           str_detect(Filename, 'B|N') ~ 'Nose',
                           TRUE ~ 'Check'))

######### subset




##################### DOCUMENT WHY WE WANT TO REMOVE 'SWAB-26_G230708
######################################################################
# Final cleaning and transformation
df_estrato_file <- df_ss %>%
  mutate(Read_category=case_when(Read_category == '**' ~ '> 10K',
                                 Read_category == '*' ~ '> 1K',
                                 TRUE ~ "< 1K"                              
  ))  %>% 
mutate(name = as.character(name),   # Ensure name is character, not factor
       genus = sapply(strsplit(name, " "), `[`, 1)) %>%  # Extract genus
  filter(genus != "Others" & new_est_reads > 0) %>%        # Filter low-quality
  mutate(sample_ID = sub("_.*", "", Filename)) %>%         # Simplify ID
  filter(Filename != "Swab-26_G230708")                    # Remove a specific outlier

###### reports from microbiology
# weeks randomly selected (12): 6, 11, 20, 26, 27, 29, 34, 39, 41, 43, 50, 52
subset_weeks <-c(6, 11, 20, 26, 27, 29, 34, 39, 41, 43, 50, 52)

##### read microbiology negative samples list
#mnp <- read_excel("Q:/IUK-A-MIGE/PROJECTS/AYALAS/TAPIR/Contamination/Microbiologically_negative_patients.xlsx")
mnp <- read_excel("/Users/osas/Desktop/01_TAPIR/Contamination/v3/Microbiologically_negative_patients.xlsx")
subset_mnp <- mnp %>% filter(KW %in% subset_weeks)
#head(mnp)
list_mnp <- mnp$sample_ID
list_mnp_ss <- subset_mnp$sample_ID

##### read microbiology negative controls list
#################MBIO NEGATIVE CONTROLS (MNC)
#mnc <- read_excel("Q:/IUK-A-MIGE/PROJECTS/AYALAS/TAPIR/Contamination/Negative_controls.xlsx")
mnc <- read_excel("/Users/osas/Desktop/01_TAPIR/Contamination/v3/Negative_controls.xlsx")
mnc2 <- mnc %>% filter(sample_ID!="Swab-58" & sample_ID!="Water-58") ##### problem with week 6 in the year after the start, both are technically week 6
subset_mnc <- mnc2 %>% filter(KW %in% subset_weeks)
#head(mnc)
list_mnc_ss <- subset_mnc$sample_ID

################## create the list of negative samples/controls
negative_ss <- c(list_mnp_ss, list_mnc_ss)

pre_mnc <- subset_mnc %>% select(sample_ID, KW) %>% mutate(type="control")
pre_mnp <- subset_mnp %>% select(sample_ID, KW) %>% mutate(type="sample")
df_negative_ss <- rbind(pre_mnc, pre_mnp)

write_xlsx(df_negative_ss, "/Users/osas/Desktop/01_TAPIR/Contamination/paper/Github/negative_subset.xlsx")

############################ TO REVIEW 

# Assign microbiological status and contamination reasons based on conditions
df_estrato_exclCriteria_file <- df_estrato_file %>%
  mutate(microbio = case_when(
    sample_ID %in% negative ~ "negative",  # If sample ID is in 'negative' list → "negative"
    TRUE ~ "positive"                      # All other samples → "positive"
  )) %>%
mutate(cont_reason = case_when(
  Site == "Neg_Ctrl" & microbio == "positive" ~ "Cross-contamination",
  Site == "Neg_Ctrl" & microbio == "negative" & Read_category == "> 10K" ~ "Library-contamination-10K",
  Site == "Neg_Ctrl" & microbio == "negative" & Read_category == "> 1K" ~ "Library-contamination-1K",
  TRUE ~ ""  # No contamination reason
))

# Optional output of the preprocessed data
write.xlsx(df_estrato_exclCriteria_file, "/Users/osas/Desktop/01_TAPIR/Contamination/df_estrato_exclCriteria_ss.xlsx")

# Load coverage data and select key columns

cov_file <- read.csv(input_path$coverage)

#cov_file <- read.csv("/Users/osas/Desktop/01_TAPIR/Contamination/paper/RAISA/bracken_plus_coverage_v6.csv")
select_cov_file <- cov_file %>% select(Filename, name, av_read_length, av_cov, week)

# Filter for a subset of weeks of interest
cov_file_ss <- select_cov_file %>% filter(week %in% subset_weeks)

# Export filtered coverage file
write_xlsx(cov_file_ss, "/Users/osas/Desktop/01_TAPIR/Contamination/paper/Github/coverage_subset.xlsx")

# Load or reuse the cleaned bracken output, removing 'unassigned' taxa
bracken_exclCriteria <- df_estrato_exclCriteria_file %>% filter(name != "unassigned")

# Separate control and biological sample types
controls <- bracken_exclCriteria %>%
  filter(Clase %in% c("Water", "Swab"))  # Control samples

samples <- bracken_exclCriteria %>%
  filter(Clase %in% c("Anus", "Nose"))  # Biological samples

# Join controls and samples on 'name' and 'week' to assess contamination overlaps
result <- samples %>%
  left_join(controls, by = c("name", "week"), suffix = c("_sample", "_control"), relationship = "many-to-many") %>%
  mutate(difference_fraction_reads = fraction_total_reads_control - fraction_total_reads_sample)  # Difference metric

# Add genus and merge with coverage data
threshold <- result %>%
  mutate(genus = sapply(strsplit(name, " "), `[`, 1)) %>%  # Extract genus from taxon name
  left_join(select_cov_file, by = c("Filename_sample" = "Filename", "name"), relationship = "many-to-many")

# Select relevant columns for downstream classification
df_threshold <- threshold %>% select(
  Filename_sample, name, genus, new_est_reads_sample, week, av_read_length,
  seq_date_sample, Clase_sample, Site_sample, av_cov, fraction_total_reads_sample,
  fraction_total_reads_control, difference_fraction_reads, microbio_sample,
  microbio_control, cont_reason_control
)

# Use parameters in your code
df_threshold <- df_threshold %>%
  mutate(
    coverage = case_when(
      av_cov < ~ "low",
      av_cov >= ~ "ok",
      TRUE ~ "no coverage"
    )
  )


# Classify samples based on coverage and contamination difference
df_threshold <- df_threshold %>%
  mutate(coverage = case_when(
    av_cov < params$thresholds$low_coverage_threshold  ~ "low",       # Low sequencing coverage
    av_cov >= params$thresholds$low_coverage_threshold  ~ "ok",       # Acceptable coverage
    TRUE ~ "no coverage"
  )) %>%
  mutate(diff_fr = case_when(
    difference_fraction_reads < 0 | is.na(difference_fraction_reads) ~ "go",  # More reads in sample or NA
    difference_fraction_reads > 0 ~ "stop",                                   # More reads in control
    TRUE ~ "empty_dfr"
  )) %>%
  mutate(system = case_when(
    # Contamination flags based on combination of read difference, coverage, and known control issues
    coverage == "low" & diff_fr == "stop" ~ "stop",
    microbio_sample == "negative" ~ "stop",
    coverage == "low" & diff_fr == "go" & is.na(cont_reason_control) ~ "check",
    coverage == "low" & diff_fr == "go" & cont_reason_control == "Cross-contamination" ~ "stop",
    coverage == "low" & diff_fr == "go" & cont_reason_control == "Library-contamination-1K" ~ "check",
    coverage == "ok" & diff_fr == "stop" ~ "stop",
    coverage == "ok" & diff_fr == "go" & is.na(cont_reason_control) ~ "go",
    coverage == "ok" & diff_fr == "go" & cont_reason_control == params$contamination_levels$cross_contamination ~ "check",
    coverage == "ok" & diff_fr == "go" & cont_reason_control == params$contamination_levels$lib_contamination_1k ~ "go",
    coverage == "ok" & diff_fr == "go" & cont_reason_control == params$contamination_levels$lib_contamination_10k  ~ "check",
    coverage == "low" & diff_fr == "go" & cont_reason_control == params$contamination_levels$lib_contamination_10k  ~ "stop",
    TRUE ~ "unclassified_system" ## where the cases do not match with any condition in 'system'
  )) %>% 
  filter(name != "unassigned")  # Final filter to remove non-informative taxonomic labels


# Final output with system classification
write_xlsx(df_threshold, output_path$threshold)

##################################################

####### criteria for defining contamination 1. cross-contamination, 2. library contamination 

df_estrato_exclCriteria_file <- df_estrato_file %>%
  mutate(microbio = case_when(sample_ID %in% negative ~ "negative",
                              TRUE ~ "positive"
  )) %>%
  mutate(cont_reason = ifelse (Site =="Neg_Ctrl" & microbio=="positive", "Cross-contamination", 
                               #ifelse(Read_category=="> 1K" & Site == "Neg_Ctrl" & microbio=="positive", "Cross-contamination",
                               ifelse(Read_category=="> 10K" & Site == "Neg_Ctrl" & microbio=="negative", "Library-contamination-10K",
                                      ifelse(Read_category=="> 1K" & Site == "Neg_Ctrl" & microbio=="negative", "Library-contamination-1K", ""
                                      )
                               )
  )
  )

####optional output
write.xlsx(df_estrato_exclCriteria_file, "/Users/osas/Desktop/01_TAPIR/Contamination/df_estrato_exclCriteria_ss.xlsx")

###############################

cov_file <- read.csv("/Users/osas/Desktop/01_TAPIR/Contamination/paper/RAISA/bracken_plus_coverage_v6.csv")
#head(cov_file)
select_cov_file <- cov_file %>% select(Filename, name, av_read_length, av_cov, week)

cov_file_ss <- select_cov_file %>% filter(week %in% subset_weeks)

write_xlsx(cov_file_ss, "/Users/osas/Desktop/01_TAPIR/Contamination/paper/Github/coverage_subset.xlsx")


###### This file already contains the info related with the microbiological report (positive, negative)
#bracken_exclCriteria <- read_excel("/Users/osas/Desktop/01_TAPIR/Contamination/paper/combined_bracken_kraken_estrato_exclCriteria250204.xlsx")
bracken_exclCriteria <- df_estrato_exclCriteria_file %>% filter(name!="unassigned")

########## removed unassinged reads
#bracken_exclCriteria <- bracken_exclCriteria %>% filter(name!="unassigned")

controls <- bracken_exclCriteria %>%
  filter(Clase %in% c("Water", "Swab"))# %>%
#   mutate(status="Red flag")
# # & cont_reason == "" )

samples <- bracken_exclCriteria %>%
  filter(Clase %in% c("Anus", "Nose"))

# Hacemos un join entre los controls y las samples por 'name' y 'week'
result <- samples %>%
  left_join(controls, by = c("name", "week"), suffix = c("_sample", "_control"), relationship = "many-to-many") %>%
  mutate(difference_fraction_reads = fraction_total_reads_control - fraction_total_reads_sample)

threshold <- result %>%
  mutate(genus = sapply(strsplit(name, " "), `[`, 1)) %>%
  left_join(select_cov_file, by = c("Filename_sample" = "Filename", "name"), relationship = "many-to-many") #%>%filter(!is.na(av_cov))#& !is.na(Type_control)) %>%filter(av_cov < 20) %>% count(status) #!is.na(Type_control)

df_threshold <-threshold %>% select(Filename_sample, name, genus, new_est_reads_sample, week, av_read_length, seq_date_sample, Clase_sample, Site_sample, av_cov, fraction_total_reads_sample, fraction_total_reads_control, difference_fraction_reads,
                                    microbio_sample, microbio_control, cont_reason_control)

##### classification
df_threshold <- df_threshold %>% mutate(coverage = case_when(av_cov < 20 ~ "low",
                                                             av_cov >= 20 ~ "ok",
                                                             TRUE ~ "no coverage")) %>%
  mutate(diff_fr = case_when(difference_fraction_reads < 0 | is.na(difference_fraction_reads) ~ "go",
                             difference_fraction_reads > 0 ~ "stop",
                             TRUE ~ "empty_dfr"))  %>%  mutate(system = case_when(coverage=="low" & diff_fr=="stop" ~ "stop",
                                                                                  microbio_sample=="negative"  ~ "stop",
                                                                                  coverage=="low" & diff_fr=="go" & is.na(cont_reason_control) ~ "check", 
                                                                                  coverage=="low" & diff_fr=="go" & cont_reason_control=="Cross-contamination" ~ "stop",
                                                                                  coverage=="low" &  diff_fr=="go" & cont_reason_control=="Library-contamination-1K" ~ "check",### is it a low biomass sample ??
                                                                                  coverage=="ok" & diff_fr=="stop" ~ "stop",
                                                                                  coverage=="ok" & diff_fr=="go" & is.na(cont_reason_control) ~ "go",
                                                                                  # coverage=="ok" & diff_fr=="go" & cont_reason_control==" " ~ "go",
                                                                                  coverage=="ok" & diff_fr=="go" & cont_reason_control=="Cross-contamination" ~ "check",
                                                                                  # !is.na(cont_reason_control) ~ "stop", # coverage=="low" & 
                                                                                  # coverage=="low" & !is.na(cont_reason_control) ~ "stop",
                                                                                  coverage=="ok" &  diff_fr=="go" & cont_reason_control=="Library-contamination-1K" ~ "go",
                                                                                  coverage=="ok" &  diff_fr=="go" & cont_reason_control=="Library-contamination-10K" ~ "check",
                                                                                  coverage=="low" &  diff_fr=="go" & cont_reason_control=="Library-contamination-10K" ~ "stop",
                                                                                  # coverage=="ok" &  diff_fr=="go" & !is.na(cont_reason_control) ~ "check",
                                                                                  #coverage=="ok" &  diff_fr=="go" & !is.na(cont_reason_control) & genus=="Staphylococcus" ~ "stop", ### update because this only applies to w7, w48
                                                                                  #coverage=="ok" &  diff_fr=="go" & !is.na(cont_reason_control) & genus!="Staphylococcus" ~ "go",
                                                                                  TRUE ~ "empty_system")) %>% filter(name!="unassigned")

### optional output file
write_xlsx(df_threshold, "/Users/osas/Desktop/01_TAPIR/Contamination/paper/Stop-Check-Go_system.xlsx")