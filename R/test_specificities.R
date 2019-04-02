#' Specificity Score Calculation with Statistical Testing
#'
#' @author Eddie Cano-Gamez, \email{ecg@@sanger.ac.uk}
#' @usage testSpecificities(rna, protein, sample_groups)
#' @description  Given an RNA expression matrix (and optionally a matching protein expresison matrix), this function calculates a score reflecting how specific the expresison of each gene is for each sample category compared to the rest. Then, it uses permutations to tests if the gene is more specific than expected by chance.
#' @param rna RNA expression matrix. Must be a data.frame with columns named after sample category and rows named after genes. Expression values must be possitive and comparable across samples (ie. normalized for library size) before running this function.
#' @param protein Protein expression matrix (optional). Defaults to "none". Must be a data.frame with columns named after sample category and rows named after genes. Rows and sample categories must match those of the RNA matrix. Expression values must be possitive and comparable across samples.
#' @param sample.labels List of sample categories (eg. biological replicates, cell types, tissues, etc...). Must be a character vector. Its elements must match the column names of the RNA (and protein) matrices. Each category is listed only once.
#' @param weight.rna When considering both RNA and protein expression, weight assigned to the RNA data. Defaults to 0.5. Must be a number between 0 and 1.
#' @param weight.protein When considering both RNA and protein expression, weight assigned to the protein data. Defaults to 0.5. Must be a number between 0 and 1.
#' @param iter Number of permutations run when testing for statistical signifiance. Defaults to 1000.
#' @details
#' This function takes either one or two correlated expression matrices (eg. RNA and protein expression from the same set of samples) and calculates a specificity score for each sample category (eg. tissue, cell type or biological replicate).
#' Specificity score calculation is done using getSpecificities() (see documentation for this function).
#' After computing specificity scores, the function generates null distributions by randomly permuting the matrix sample names. This is done as many times as specified by the user. The observed specificity score is compared to the scores observed in the permuted data and a P value is calculated as the number of times the observed score is larger than the permuted score.
#' To account for multiple hypothesis testing, the function also implements P-value correction using the Benjamini-Hochberg method (see documentation for p.adjust()).
#' @export
#' @examples
#' ## USING ONE DATA SET ONLY (eg. RNA ONLY)
#'
#' # Simulating mock RNA data:
#' rna.example <- data.frame(matrix(rnorm(9000,mean=2000,sd=100),ncol=9,nrow=100))
#' sample_groups <- c("A","B","C")
#' gene_names <- paste("g",1:100,sep="")
#' colnames(rna.example) <- rep(sample_groups,each=3)
#' rownames(rna.example) <- gene_names
#'
#' # Simulating sets of highly expressed genes in each sample group only
#' rna.example[1:10,1:3] <- rna.example[1:10,1:3] + rnorm(1,mean=4000,sd=1000)
#' rna.example[20:30,4:6] <- rna.example[20:30,4:6] + rnorm(1,mean=4000,sd=1000)
#' rna.example[90:100,7:9] <- rna.example[90:100,7:9] + rnorm(1,mean=4000,sd=1000)
#'
#' # Running the function:
#' testSpecificities(rna.example, sample.labels = sample_groups, iter=1000)
#'
#' @examples
#' ## USING TWO MATCHING DATA SETS (eg. RNA AND PROTEIN)
#'
#' # Simulating matching mock Protein data:
#' prot.example <- data.frame(matrix(rnorm(9000,mean=7000,sd=100),ncol=9,nrow=100))
#' colnames(prot.example) <- rep(sample_groups,each=3)
#' rownames(prot.example) <- gene_names
#'
#' # Simulating sets of highly expressed proteins in each sample group only:
#' prot.example[1:10,1:3] <- prot.example[1:10,1:3] + rnorm(1,mean=1500,sd=1000)
#' prot.example[20:30,4:6] <- prot.example[20:30,4:6] + rnorm(1,mean=1500,sd=1000)
#' prot.example[90:100,7:9] <- prot.example[90:100,7:9] + rnorm(1,mean=1500,sd=1000)
#'
#' # Running the function:
#' testSpecificities(rna.example, prot.example, sample.labels = sample_groups, iter=1000)
testSpecificities <- function(rna.exp, prot.exp="none", sample.labels, weight.rna=0.5, weight.protein=0.5, iter=1000){
  if(sum(names(table(colnames(rna.exp))) %in% sample.labels) != length(sample.labels)){
    stop("RNA columns and sample labels do not match",call.=F)
  }
  ifelse(prot.exp=="none",{
    S <- getSpecificities(rna.exp, prot.exp, sample.labels, weight.rna, weight.protein)
    test.res <- matrix(0,nrow = dim(S)[1], ncol=dim(S)[2])
    foreach(icount(iter), .combine='c', .errorhandling='pass') %do% {
      rna.null <- rna.exp
      colnames(rna.null) <- sample(colnames(rna.null))
      S.null <- getSpecificities(rna.null, sample.labels=sample.labels)
      comparison <- S < S.null
      comparison <- comparison*1
      test.res <- test.res + comparison
    }
  },
  {
    if(nrow(prot.exp)!=nrow(rna.exp)){
      stop("RNA and protein matrices have a different number of rows",call.=F)
    }
    if(sum(rownames(prot.exp)!= rownames(rna.exp)) > 0){
      stop("Gene names do not match between RNA and protein",call.=F)
    }
    if(sum(names(table(colnames(prot.exp))) %in% sample.labels) != length(sample.labels)){
      stop("Protein columns and sample labels do not match",call.=F)
    }
    S <- getSpecificities(rna.exp, prot.exp, sample.labels, weight.rna, weight.protein)
    test.res <- matrix(0,nrow = dim(S)[1], ncol=dim(S)[2])
    foreach(icount(iter), .combine='c', .errorhandling='pass') %do% {
      rna.null <- rna.exp
      colnames(rna.null) <- sample(colnames(rna.null))
      protein.null <- prot.exp
      colnames(protein.null) <- colnames(rna.null)
      S.null <- getSpecificities(rna.null, protein.null, sample.labels=sample.labels, weight.rna, weight.protein)
      comparison <- S < S.null
      comparison <- comparison*1
      test.res <- test.res + comparison
    }
  })

  pvals <- as.data.frame((test.res+1)/iter)
  padj <-as.data.frame(apply(pvals, MARGIN=2, FUN=function(p){p.adjust(p, method="BH")}))
  res <- list(specificities=S,p.val=pvals, p.adj=padj)
  return(res)
}
