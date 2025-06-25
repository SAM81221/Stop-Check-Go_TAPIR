
###################################### June 12th, 2025
################################### Author: Ayorinde Afolayan
######################################### This code has been commented with the use of ChatGPT (GPT-4.5-turbo model)
#########################################. prompt: "Only comment out the following function to make it easier to understand (by a third person), do not modify the code"

# A. bracken_parse_df
bracken_parse_df <- function(bracken_combined_report){
  # Load necessary libraries
  # library(tidyverse)
  # library(magrittr)
  # library(ggthemes)
  # library(scales)
  # library(dplyr)
  # library(readxl)
  
  # Step 1: Process the combined Bracken report
  # - Calculate percentage of total reads
  # - Classify the sample site (Anus, Nose, Negative Control)
  # - Extract sample name
  # - Categorize estimated read count with symbols
  # - Replace low-percentage or irrelevant taxa with "Others"
  bracken_combined_report_df <- bracken_combined_report %>% 
    dplyr::mutate(percentage_total_reads = fraction_total_reads * 100) %>% 
    dplyr::mutate(Site = case_when(
      str_detect(Filename, '[Ss]wab|[Ww]ater') ~ 'Neg_Ctrl',
      str_detect(Filename, 'A') ~ 'Anus',
      str_detect(Filename, 'B|N') ~ 'Nose'
    )) %>% 
    dplyr::mutate(Sample = str_extract(Filename, "^[a-zA-Z0-9]+")) %>% 
    dplyr::mutate(Read_category = case_when(
      new_est_reads < 1000 ~ "",
      new_est_reads >= 1000 & new_est_reads <= 10000 ~ "*",
      new_est_reads > 10000 ~ "**"
    )) %>% 
    dplyr::mutate(name = if_else(percentage_total_reads >= 1, name, 'Others')) %>% 
    dplyr::mutate(name = if_else(!str_detect(name, 'virus|phage'), name, 'Others'))
  
  # Step 2: Capture and summarize the "Others" category
  capture_others_df <- bracken_combined_report_df %>% 
    dplyr::filter(name == 'Others') %>%
    group_by(Filename) %>%
    dplyr::mutate(
      sumpercentage = sum(percentage_total_reads),
      sumkrakenassignedreads = sum(kraken_assigned_reads),
      sumaddedreads = sum(added_reads),
      sumnewestreads = sum(new_est_reads),
      sumfractiontotalreads = sum(fraction_total_reads),
      sumpercentagetotalreads = sum(percentage_total_reads)
    ) %>%
    ungroup() %>% 
    dplyr::select(-c(
      percentage_total_reads, kraken_assigned_reads, added_reads, 
      new_est_reads, taxonomy_id, fraction_total_reads, 
      percentage_total_reads, Read_category
    )) %>% 
    distinct() %>% 
    dplyr::rename(
      kraken_assigned_reads = sumkrakenassignedreads,
      added_reads = sumaddedreads,
      new_est_reads = sumnewestreads,
      fraction_total_reads = sumfractiontotalreads,
      percentage_total_reads = sumpercentagetotalreads
    ) %>% 
    dplyr::mutate(taxonomy_id = 0) %>% 
    dplyr::mutate(Read_category = case_when(
      new_est_reads < 1000 ~ "",
      new_est_reads >= 1000 & new_est_reads <= 10000 ~ "*",
      new_est_reads > 10000 ~ "**"
    )) %>% 
    distinct() %>% 
    dplyr::select(Filename, name, taxonomy_id, taxonomy_lvl, kraken_assigned_reads, added_reads, new_est_reads, 
                  fraction_total_reads, percentage_total_reads, Site, Sample, Read_category)
  
  # Step 3: Merge "Others" with the main dataframe (excluding original 'Others')
  bracken_combined_report_semifinal_df <- bracken_combined_report_df %>%
    dplyr::filter(name != 'Others') %>% 
    dplyr::bind_rows(capture_others_df)
  
  # Step 4: Calculate and add 'unassigned' taxon data
  capture_unassigned_df <- bracken_combined_report_semifinal_df %>% 
    group_by(Filename) %>%
    dplyr::mutate(
      sumkrakenassignedreads = sum(kraken_assigned_reads),
      sumaddedreads = sum(added_reads),
      sumnewestreads = sum(new_est_reads),
      sumfractiontotalreads = sum(fraction_total_reads),
      sumpercentagetotalreads = sum(percentage_total_reads)
    ) %>% 
    dplyr::select(-c(
      name, percentage_total_reads, kraken_assigned_reads, added_reads, 
      new_est_reads, taxonomy_id, fraction_total_reads, 
      percentage_total_reads, Read_category
    )) %>% 
    distinct() %>% 
    dplyr::mutate(
      kraken_assigned_reads = round(sumkrakenassignedreads / sumpercentagetotalreads * 100 - sumkrakenassignedreads),
      added_reads = round(sumaddedreads / sumpercentagetotalreads * 100 - sumaddedreads),
      new_est_reads = round(sumnewestreads / sumpercentagetotalreads * 100 - sumnewestreads),
      fraction_total_reads = round(sumfractiontotalreads / sumpercentagetotalreads * 100 - sumfractiontotalreads),
      percentage_total_reads = 100 - sumpercentagetotalreads,
      taxonomy_id = 0,
      taxonomy_lvl = 'U',
      name = 'unassigned'
    ) %>% 
    distinct() %>% 
    dplyr::mutate(Read_category = case_when(
      new_est_reads < 1000 ~ "",
      new_est_reads >= 1000 & new_est_reads <= 10000 ~ "*",
      new_est_reads > 10000 ~ "**"
    )) %>% 
    distinct() %>% 
    ungroup() %>%
    dplyr::select(Filename, name, taxonomy_id, taxonomy_lvl, kraken_assigned_reads, added_reads, new_est_reads, 
                  fraction_total_reads, percentage_total_reads, Site, Sample, Read_category)
  
  # Step 5: Merge with the previous result to get the final dataframe
  bracken_combined_report_final_df <- bracken_combined_report_semifinal_df %>%
    dplyr::bind_rows(capture_unassigned_df)
  
  # Return the final parsed dataframe
  return(bracken_combined_report_final_df)
}

###################################### June 12th, 2025
################################### Author: Stefany Ayala-Montaño
######################################### This code has been commented with the use of ChatGPT (GPT-4.5-turbo model)
#########################################. prompt: "Only comment out the following function to make it easier to understand (by a third person), do not modify the code"

# Function to import biom files using a default taxonomy parsing function
# Define your import function

B. import_biom_data

import_biom_data <- function(data_input) {
  import_biom(data_input, parseFunction = parse_taxonomy_default)
}


# Function to process biom data and add metadata columns

C. process_biom_data
process_biom_data <- function(data) {
  # Predefined list of sample names considered as microbiologically negative 
 
  
  if (!is.character(negative_list)) {
    stop("❌ 'negative' must be character")
  }
  
  # Flatten in case it's a list of character vectors
  negative <- negative_list
  
  # Print a summary
  #cat("🧪 Number of negative samples provided:", length(negative), "\n")
  #cat("📋 Sample IDs:\n")
 # print(negative)
  
  data %>%
    microViz::ps_mutate(
      # Add 'Site' column based on pattern in sample ID
      Site = dplyr::case_when(
        stringr::str_detect(Id, '-A_') ~ 'Anus',
        stringr::str_detect(Id, 'B|N') ~ 'Nose',
        TRUE ~ 'Neg_Ctrl'
      ),
      # Extract patient identifier from sample ID
      Patient = dplyr::case_when(
        stringr::str_detect(Id, '^T\\d+') ~ stringr::str_extract(Id, '^T\\d+'),
        stringr::str_detect(Id, 'Wa') ~ 'Water',
        stringr::str_detect(Id, 'Sw') ~ 'Swab',
        TRUE ~ Id
      ),
      # Identify whether it's a true sample or a negative control
      sample_control = dplyr::case_when(
        stringr::str_detect(Id, '-A_') ~ 'True sample',
        stringr::str_detect(Id, 'B|N') ~ 'True sample',
        TRUE ~ 'Neg_Ctrl'
      ),
      # Extract filename from the sample ID
      Filename = str_extract(Id, "^[^_]+"),
      # Extract week number from the filename
      week = ifelse(grepl("^[^-]+-[0-9]+", Filename), 
                    as.numeric(gsub("^[^-]+-([0-9]+).*", "\\1", Filename)), 
                    NA),
      # Classify sample as positive or negative based on filename list
      microbio = dplyr::case_when(
        Filename %in% negative ~ "negative",
        TRUE ~ "positive"
      ),
      # Determine overall classification of sample/control status
      Sample_or_Control = dplyr::case_when(
        sample_control == "True sample" & microbio == "negative" ~ "Microbio neg. sample",
        sample_control == "Neg_Ctrl" & microbio == "negative" ~ "Control",
        sample_control == "Neg_Ctrl" & microbio == "positive" ~ "Cross contam.",
        TRUE ~ "True sample"
      )
    )
}


# List of weeks to be processed
#weeks_to_process <- c(6, 11, 20, 26, 27, 29, 34, 39, 41, 43, 50, 52)

#weeks_to_process <- as.numeric(stringr::str_extract(list_biom_processed, "\\d{2}"))

# D. create_go_dataset 

# Function that generates a dataset for a single week
create_go_dataset <- function(week_number) {
  # Filter the contamination data to only include entries from the current week
  system_contam <- stop_plus_conditions_in_check %>% 
    filter(week == week_number)
  
  # Join with taxonomy ID mapping based on the 'name' column
  for_go_dataset <- left_join(system_contam, tax_id, by = "name")
  
  # Identify entries where taxonomy_id was not found
  missing_tax_id <- for_go_dataset %>% filter(is.na(taxonomy_id))
  if (nrow(missing_tax_id) > 0) {
    message(paste("⚠️ Warning for week", week_number, ": Missing taxonomy_id"))
    print(missing_tax_id)  # Print missing entries for review
  }
  
  # Return a unique data frame with filenames and their corresponding taxonomy IDs
  data.frame(
    filename = for_go_dataset$Filename_sample,
    taxonomy_id = for_go_dataset$taxonomy_id
  ) %>% unique()
}


# E. assign_zero_to_positions

assign_zero_to_positions <- function(df, otu_matrix) {
  # Loop through each row of the data frame
  for (i in 1:nrow(df)) {
    # Find the position of the taxonomy_id in the row names of the otu_matrix
    pos_row <- which(rownames(otu_matrix) == df$taxonomy_id[i])
    # Find the position of the filename in the column names of the otu_matrix
    pos_col <- which(colnames(otu_matrix) == df$filename[i])
    
    # If both positions exist, assign zero to that position in the matrix
    if (length(pos_row) > 0 & length(pos_col) > 0) {
      otu_matrix[pos_row, pos_col] <- 0
    }
  }
  # Return the modified otu_matrix
  return(otu_matrix)
}
  

