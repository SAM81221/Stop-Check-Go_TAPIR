
---

# The TAPIR project

**License Notice**

This code is part of an unpublished manuscript and is currently under
embargo.
Please do not redistribute. The code will be made publicly available
upon article publication.

## Description

1.  bracken_S_sequencingDate.txt

The files are generated individually from Bracken Kraken command and
then combined into a single file in layman's terms in Nextflow.

*Basic example*

`process_combine_Files { input: file(files) from file_list`

`output: file("combined.txt") into combined_output`

`cat ${files.join(' ')} > combined.txt`

```  } ```

`workflow { file_list = Channel.fromPath("data/\*.txt")`
`combineFiles(file_list) }`

2.  Bracken_combine_select.tsv

The output files from every sequencing week are transferred locally and
used as input files in Rstudio. The files are combined into a single
file using rbind command, this single file is used as an input variable
in the function called bracken_parse_df found in functions.R

3.  Microbiologically negative patients and microbiologically negative
    controls

This file contains the reports from microbiology from culture agar
plates, where there was negative growth for bacteria in samples and
controls. For these negative sample types, we expected ‘too low’ in the
Qubit DNA concentration after DNA extraction.

This is an excel user-input file with two columns: sample_ID, and
sample_type (sample or control).

| sample_ID | sample_type |
|-----------|-------------|
| T251-3-A  | sample      |
| T251-3-B  | sample      |
| Water-8   | control     |
| Swab-48   | control     |

4.  Coverage file

The coverage file is generated from the coverage results of
assembly_info.txt from FLYE assembler and the Python script of
bracken-translate from Kraken. Assembly_info.txt provides the
information on coverage for each contig within the metagenome, and the
broken translate correlates each contig with tax ID and species name.

`kraken --db \$DBNAME sequences.fa \> sequences.kraken kraken-translate`
`--db \$DBNAME sequences.kraken \> sequences.labels`

<https://ccb.jhu.edu/software/kraken/MANUAL.html>

The coverage file should contain ‘Filename’ as the sample ID, ‘name’ as
the species name, ‘av_cov’ with the coverage values, and ‘week’.

| Filename | name             | av_cov | week |
|----------|------------------|--------|------|
| T251-3-A | Escherichia coli | 32     | 8    |

5.  Kraken-biom

The biom files are generated from the individual Kraken output files
with the command Kraken-biom

Installation:

`conda install bioconda::kraken-biom`

command:

`kraken-biom /$PATH/*_SeqDate_S_bracken_kraken.report.txt -o /$PATH/bracken_kraken_SeqDate.biom --fmt json`

6.  Tax_id

The purpose of this file is to query taxonomic IDs to species name,
therefore it contains two columns ‘name’, and ‘taxonomy_id’. We plan to
improve this file for general use through the generation of a
comprehensive file (bacteria) by using taxonKit.

<https://bioinf.shenwei.me/taxonkit/tutorial/>

7. Microbiological reports 

Contains the reports from positive growth species. In the context of this TAPIR project, the species grown on MacConkey agar plates were identified with MALDI-TOF. 

| sample_ID | name             | 
|-----------|------------------|
| T251-3-A  | Escherichia coli |


* File structure *
```
Stop-Check-Go_TAPIR/ 
├── main.R
├── R/
│ └── functions.R
├── input_files/
│ ├── taxid_from_cov_file.xlsx
│ ├── Microbiologically_negative.xlsx
│ ├── bracken_plus_coverage_v6.csv
│ ├── bracken_combine_select_2023_2024.tsv
│ ├── Access_TAPIR_250325.xlsx
│ └── biom/
├── output_files/
│ ├── Samples_system_raw.xlsx
│ ├── Samples_in_Stop_after_Check.xlsx
│ └── Samples_in_Go_system.xlsx
├── README.md
├── config.yaml
├── .Rhistory
└── .gitignore

```

** Note that the input files are:

2.  Bracken_combine_select.tsv
3.  Microbiologically_negative.xlsx
4.  bracken_plus_coverage_v6.csv
5.  bracken_kraken.biom
6.  taxid_from_cov_file.xlsx
7.  Access_TAPIR_250325.xlsx

Files 2, 3, 4, 5 and 7 are required for the script to work. File 6 is a fixed file. The location of these files should be`/$PATH/input_files/`. The script contains functions to process the biom format files from 5 by extracting the weeks of analysis and subsetting files 2, 3, 4 and 7. The output files produced are Excel metadata files associated with the
classification system ‘Samples_system_raw’, ‘Samples_in_Stop_after_Check’, ‘Samples_in_Go_after_Check’, and the phyloseq decontaminated objects saved as '.rds' in `/$PATH/output_files/biom_decontaminated/`

* Disclaimer: Use of AI Assistance *

I used assistance of ChatGPT (GPT-4o-mini model) for the code publication. The AI provided support based on the prompt:

"I have a project in R with Git and now I need to push it for publication. Guide me in every step."

While the AI’s suggestions helped with templates for publication and the version control processes, all final decisions, code implementations, and project content were written, reviewed and verified by the project author.