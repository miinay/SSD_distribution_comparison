# SSD_distribution_comparison
# What is this page?
This page contains the R script to compare four statistical distributions to estimate species sensitivity distributions (SSDs), using data retrieved from the database. 

# Objective
Species sensitivity distributions (SSDs) are crucial tools used for environmental risk assessment. Whereas several statistical distributions are available for SSD estimation, the fundamental question of which statistical distribution should be used has received limited systematic analysis. In this analysis, we perform SSD estimations using four statistical distributions (log-normal, log-logistics, Burr type III, and Weibull distributions) and subsequently compare the four distributions based on AICc and HC5 differences. 

# Files
1. Rcode.md  
An example R code for analysis and visualization. This code requires an input dataset (e.g., "example.xlsx").

2. example.xlsx  
This dataset "example.xlsx" includes 20,000 test records randomly selected from the "EnviroTox" database only for demonstration.
All the data used in the study was collected from the "EnviroTox" database (version 2.0.0) (https://envirotoxdatabase.org/).

3. Figures_example_data  
Figures generated from the example data using the R code.

# Publication
In preparation
