# Masters_Thesis_Code
A smattering of scripts that were used for computational aspects of my Master's Thesis project. Below, the functions of the various scripts are briefly described.

For more context as to the reasons for these scripts or for curiosity about the larger Thesis project, please see EC_Thesis_Final.pdf.


Create_Dir.py - script used for organizing urinary E. coli sequences for annotation and phage prediction. Also contains methods for checking the redundency of a phage network and the parsing of taxonomic matches from BLAST+ output files from predicted bacteriphage sequences

diff_networks.py - script used to determine the individual predicted phage sequences that belong to particular culsters of the large phage network

Intact_phage_histogram.r - uses ggplot to produce a graphical representation of numerical data

edit_distance.py - a script that screens one directory of phage sequences against another to determine if there are idenitcal sequences between the two directories

integrase parse.py - parses PATRIC output files for "Integrase" genes, compares the various integrase sequences, and derives statistics about identified integrase genes

names.py - used to correct BLAST+ taxonomy outputs for more organized taxonomy grouping

netprep.r - script developed by Shapiro and Putonti to determine relatedness of bacteriophage sequences based on number of shared genes (doi: 10.1128/mBio.01870-17)

num_intact_phages.py - smattering of code for identification of predicted phage sequences from bacteria sequences screened with the PHASTER phage prediction tool. Also performed muliple reorganization, remaning, and reformatting functions

phage_network.py - produced a phage network via an edge list from an outputted pan-genome file produced by an r-script developed by Shapiro and Putonti (doi: 10.1128/mBio.01870-17) and the Anvi'o genomics tool

separate_phages.py - Organized workflow of producing a phage network from a directory of predicted phage sequences

search_sheet.py - short script used for parsing a large Excel doc for specific metadata entries

superinfection_blast.py - script for running BLAST+ on a remote machine and pulling the results to a local machine for analysis
