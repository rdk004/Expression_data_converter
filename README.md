# Expression_data_converter
This is a script written in R to convert count values (.tsv files) to TPM (Transcripts per million). This code has been structured by me but I don't claim it as IP since it has been sourced from Copilot, ChatGPT and other sources.
Note: This code is limited to human data, but it can be generalized for expression data of any species

Let me explain the flow and structure of the code

1) Section 1: Loading all the necessary libraries (Note: Please use install.packages('package_name') before loading them)
              If that doesn't work please use this command:
                  if (!requireNamespace("BiocManager", quietly = TRUE))
                      install.packages("BiocManager")
                      BiocManager::install('package_name')

2) Section 2: Importing the count values (.tsv) file from your system

3) Section 3: Antilog transformation ( 2 ^ (count value) + 1) and extract the Ensembl Gene IDs as a separate dataframe
Note: Web databases like Xenabrowser (TCGA), HPA, GEO etc. usually have log transformed count data. It needs to be antilog transformed before downstream manipulations

4) Section 4: Importing the GTF file and processing it for gene lengths
Note: 1) TPM conversion requires the gene lengths of all the genes in your dataset.
      2) You can download the GTF file from here: http://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz (you can find gtf #           files for other species as well (just modify the path to this: http://ftp.ensembl.org/pub/release-109/gtf
      3) Make note of the directory in which it gets stored post installation
The GTF file will most likely have gene lengths for almost all the human genes. The code allows you to compare your Ensembl ID dataframe with the GTF file t
to extract the required gene lengths

5) Section 5: TPM transformations
              Formula for conversion: TPM = x/sum(x) * 1e6 (where x = c/l; c is the count value (antilog transformed) and l is the length of the gene in #  #               kb; sum(x) is all x values summed for a sample)
              a) Let us apply the first transformation - normalize each element with corresponding gene length
              b) Let us apply the second   transformation - divide each every element by column_sum/10^6

6) Section 6: Export the final TPM dataframe as a (.csv) file

