
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
import_biom_data <- function(data_input) {
  import_biom(data_input, parseFunction = parse_taxonomy_default)
}


# Function to process biom data and add metadata columns
process_biom_data <- function(data) {
  # Predefined list of sample names considered as microbiologically negative 
 
  # negative <- c("T053-2-A", "T053-2-B", "T061-3-A", "T065-4-A", "T071-5-A", "T071-5-N", "T072-5-A", 
  #               "T072-6-B", "T074-6-A", "T078-6-A", "T084-8-A", "T084-8-B", "T088-9-A", "T089-9-A", 
  #               "T089-9-N", "T090-9-A", "T090-9-B", "T091-9-A", "T092-9-A", "T092-9-B", "T097-10-A", 
  #               "T099-10-A", "T103-11-A", "T112-14-A", "T112-14-B", "T116-14-A", "T116-14-B", "T117-14-A", 
  #               "T117-14-N", "T120-15-A", "T120-15-B", "T126-16-A", "T127-16-A", "T129-16-A", "T130-16-A", 
  #               "T130-16-B", "T131-17-A", "T131-17-B", "T132-17-A", "T137-17-A", "T138-17-A", "T138-17-B", 
  #               "T143-18-A", "T144-18-B", "T146-18-B", "T147-19-A", "T147-19-B", "T149-20-B", "T152-21-A", 
  #               "T152-21-N", "T155-21-A", "T158-22-A", "T158-22-B", "T159-22-A", "T159-22-B", "T160-22-A", 
  #               "T160-22-B", "T167-23-B", "T168-23-A", "T168-23-N", "T169-23-A", "T169-23-B", "T153-24-A", 
  #               "T153-24-B", "T155-24-A", "T170-24-A", "T171-24-A", "T174-24-A", "T153-25-A", "T187-27-B", 
  #               "T189-27-A", "T189-27-N", "T190-27-A", "T191-27-A", "T191-27-B", "T192-27-A", "T193-27-A", 
  #               "T193-27-B", "T196-27-A", "T196-27-B", "T198-27-A", "T198-27-B", "T203-29-A", "T205-29-A", 
  #               "T206-29-A", "T206-29-B", "T208-29-A", "T210-29-A", "T210-29-N", "T213-30-A", "T213-30-B", 
  #               "T215-30-A", "T216-30-B", "T213-31-B", "T216-31-A", "T217-31-A", "T218-31-A", "T218-31-B", 
  #               "T219-31-B", "T213-32-N", "T224-32-A", "T224-32-B", "T228-32-N", "T229-32-A", "T229-32-B", 
  #               "T233-34-N", "T235-34-N", "T238-35-A", "T238-35-B", "T244-36-A", "T244-36-B", "T244-37-A", 
  #               "T244-37-B", "T247-37-A", "T247-37-B", "T244-38-A", "T259-39-A", "T259-39-B", "T262-39-A", 
  #               "T262-39-N", "T266-40-A", "T268-40-A", "T269-40-A", "T270-41-A", "T270-41-B", "T271-41-A", 
  #               "T273-41-A", "T273-41-B", "T275-41-A", "T275-41-N", "T278-42-A", "T280-42-A", "T280-42-B", 
  #               "T286-43-A", "T294-45-A", "T294-45-N", "T295-45-A", "T295-45-B", "T299-46-A", "T299-46-B", 
  #               "T300-46-A", "T300-46-B", "T259-47-B", "T303-47-A", "T305-47-A", "T312-48-A", "T315-49-A", 
  #               "T315-49-B", "T316-49-A", "T316-49-B", "T319-50-A", "T320-50-A", "T320-50-B", "T321-50-A", 
  #               "T323-51-A", "T323-51-B", "T324-51-A", "T324-51-B", "T336-1-A", "T336-1-B", "T339-2-A", 
  #               "T339-2-B", "T340-2-A", "T340-2-B", "T342-3-B", "T346-3-A", "T346-3-B", "T347-3-A", "T347-3-B", 
  #               "T294-4-A", "T354-4-A", "T355-4-A", "T355-4-B", "Swab-3", "Water-3", "Swab-5", "Water-5", 
  #               "Swab-6", "Water-6", "Swab-7", "Water-7","Swab-8", "Water-8", "Swab-9", "Water-9", "Swab-10", 
  #               "Water-10", "Swab-11", "Water-11", "Swab-13", "Water-14", "Swab-17", "Water-17", "Swab-18", 
  #               "Water-18", "Swab-19", "Water-19", 
  #               "Swab-21", "Water-21", "Swab-22", "Water-22", "Swab-23", "Water-23", "Swab-24", "Water-24", 
  #               "Swab-25", "Water-25", "Swab-26", "Water-26", "Swab-27", "Water-27", "Swab-28", "Water-28", 
  #               "Swab-29", "Swab-30", "Water-30", "Swab-31", "Water-31", "Swab-32", "Water-32", "Swab-33", 
  #               "Water-33", "Swab-34", "Water-34", "Swab-35", "Water-35", "Swab-36", "Water-36", "Swab-37", 
  #               "Water-37", "Swab-38", "Water-38", "Swab-39", "Water-39", "Swab-40", "Water-40", "Swab-41", 
  #               "Water-41", "Swab-42", "Water-42", "Swab-43", "Swab-44", "Water-44", "Swab-45", 
  #               "Water-45", "Swab-46", "Water-46", "Swab-47", "Water-47", "Swab-49", "Water-49", "Swab-50", 
  #               "Water-50", "Swab-51", "Water-51", "Swab-52", "Water-52", "Swab-53", "Water-53", "Swab-54", 
  #               "Water-54", "Swab-55", "Water-55", "Swab-56", "Water-56", "Swab-57", "Water-57", "Swab-58", 
  #               "Water-58", "Swab-59", "Water-59", "Swab-60", "Water-60", "Swab-61", "Water-61", "Swab-62", 
  #               "Water-62")
  
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

# Loop through all weeks and process the biom data objects
# for (week in names(biom_list)) {
#   # Format week number with two digits (e.g., "06", "11")
#   #week_str <- sprintf("%02d", week)
#   
#   # Construct input and output variable names dynamically
#  # input_name <- paste0("w", week_str, "_biom")
#   output_name <- paste0("w", week_str, "_biom_m")
#   
#   # Check if the input object exists before processing
#   if (exists(input_name)) {
#     processed <- process_biom_data(get(input_name))  # Process the biom data
#     assign(output_name, processed)  # Save the result to a new object
#   } else {
#     print(paste("The object was not found:", input_name))  # Notify if the input doesn't exist
#   }
# }
# 


 

# Load the main input file containing sample contamination conditions
#red_plus_conditions_in_check <- read_excel("red_plus_conditions_in_check_250409.xlsx")


##### B. Create green dataset 

# List of weeks to be processed
#weeks_to_process <- c(6, 11, 20, 26, 27, 29, 34, 39, 41, 43, 50, 52)

#weeks_to_process <- as.numeric(stringr::str_extract(list_biom_processed, "\\d{2}"))

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


## D. assign_zero_to_positions

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
    ##
    
#     
#     # Modify OTU matrix
#     otu_matrix_modified <- assign_zero_to_positions(df, otu_matrix)
#     otu_table_modified <- otu_table(otu_matrix_modified, taxa_are_rows = TRUE)
#     
#     # Create new phyloseq object
#     biom_green <- phyloseq(
#       otu_table_modified,
#       sample_data(biom_obj),
#       tax_table(biom_obj)
#     )
#     
#     # Save to global environment and RDS
#     assign(output_name, biom_green, envir = .GlobalEnv)
#     output_file <- file.path(output_path$biom_go_folder, paste0(output_name, ".rds"))
#     saveRDS(biom_green, output_file)
#     
#     message(paste("✔ Week", week_str, "processed and saved as", output_name))
#     
#   } else {
#     warning(paste("❌ Missing objects for week", week_str, ":", paste(missing, collapse = ", ")))
#   }
# }

