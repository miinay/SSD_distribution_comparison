# SSD_distribution_comparison
# What is this page?
This page contains the R script to compare several statistical distributions to estimate species sensitivity distributions (SSDs), using data retrieved from the database. (See the below reviewed paper for the details.)

# Objective
SSD is a promising approach to derive predicted no effect concentration (PNEC) in ecological risk assessment. In order to estimate robust SSD, the sample size (number of species) is required to be ≧ 5 to 10. However, performing chronic toxicity tests for a lot of species is challenging due to the extensive costs and labor, in contrast to acute toxicity tests. To address this challenge, it is useful to estimate SSD based on chronic data from that based on acute data.
Here we show the R code and an example dataset to analyze the relationshipe between acute and chronic SSDs.

# Files
1. Rcode.md
An example R code for analysis and visualization. This code requires an input dataset (e.g., "example.xlsx").

2. example.xlsx
This dataset "example.xlsx" includes 20,000 test records randomly selected from the "EnviroTox" database only for demonstration.
All the data used in the study was collected from the "EnviroTox" database (ver.1.2.0) (https://envirotoxdatabase.org/). Please contact us if you like to exactly reproduce our results.

3. SSD_chronic_figs
Figures generated from the example R code.

# Publication
