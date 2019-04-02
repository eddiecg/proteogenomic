# proteogenomic
Group of R functions to identify genes specific to a given cell type or condition. These functions compute specificity scores for each gene using all replicates of a given cell type/condition and use permutations to perform statistical testing.

There are two main functions:

getSpecificities() = Computes specificity scores but does no statistical testing
testSpecificities() = Computs specificity scores and infers statistical significance via permutations. Also, corrects for multiple hypothesis testing using the Benjamini-Hochberg method

To use these functions, clone this repository and install it as an R package.
