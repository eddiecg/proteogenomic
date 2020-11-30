# proteogenomic
Group of R functions to identify genes specific to a given cell type or condition. These functions compute specificity scores for each gene using all replicates of a given cell type/condition and use permutations to perform statistical testing.

To use these functions, install the proteogenomic package using the following commands:
  
    install.packages("devtools") #If not yet installed
    devtools::install_github("eddiecg/proteogenomic")

Proteogenomic contains two main functions:

    getSpecificities()

Computes specificity scores but does no statistical testing

    testSpecificities()

Computes specificity scores and infers statistical significance via permutations. Also, corrects for multiple hypothesis testing using the Benjamini-Hochberg method
