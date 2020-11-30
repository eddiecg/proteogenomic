#' Specificity Score Calculation
#'
#' @author Eddie Cano-Gamez, \email{ecg@@sanger.ac.uk}
#' @usage getSpecificities(rna, protein, sample_groups)
#' @description  Given an RNA expression matrix (and optionally a matching protein expresison matrix), this function calculates a score reflecting how specific the expresison of each gene is for each sample category compared to the rest.
#' @param rna RNA expression matrix. Must be a data.frame with columns named after sample category and rows named after genes. Expression values must be possitive and comparable across samples (ie. normalized for library size) before running this function.
#' @param protein Protein expression matrix (optional). Defaults to "none". Must be a data.frame with columns named after sample category and rows named after genes. Rows and sample categories must match those of the RNA matrix. Expression values must be possitive and comparable across samples.
#' @param sample.labels List of sample categories (eg. biological replicates, cell types, tissues, etc...). Must be a character vector. Its elements must match the column names of the RNA (and protein) matrices. Each category is listed only once.
#' @param weight.rna When considering both RNA and protein expression, weight assigned to the RNA data. Defaults to 0.5. Must be a number between 0 and 1.
#' @param weight.protein When considering both RNA and protein expression, weight assigned to the protein data. Defaults to 0.5. Must be a number between 0 and 1.
#' @details
#' This function takes either one or two correlated expression matrices (eg. RNA and protein expression from the same set of samples) and calculates a specificity score for each sample category (eg. tissue, cell type or biological replicate).
#' This specificity score goes from 0 to 1, where 1 means the gene is expressed only in that particular sample group/tissue/cell type but not in the others, while 0 means the gene is not expressed in that sample group.
#' To calcualte this score, the function averages all replicates gene-wise to obtain the mean expression value of each gene per sample group. Then, it performs Euclidean normalization to obtain a specificity score per gene.
#' When two matching data sets are included (eg. RNA and protein), normalization is done independently for each data set and then a combined socre is derived via a weighted sum of scores. These weights default to 0.5 but can be fine tuned by the user.
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
#' getSpecificities(rna.example, sample.labels = sample_groups)
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
#' getSpecificities(rna.example, prot.example, sample.labels = sample_groups)
getSpecificities <- function(rna, protein="none", sample.labels, weight.rna=0.5, weight.protein=0.5){
  if(sum(names(table(colnames(rna))) %in% sample.labels) != length(sample.labels)){
    stop("RNA columns and sample labels do not match",call.=F)
  }
  rna <- data.frame(sapply(sample.labels, FUN=function(x){
    exp = rna[,colnames(rna)==x]
    if(is.null(dim(exp))){
      exp = exp
    } else{
      exp = rowMeans(exp, na.rm=T)
    }
    return(exp)
  }))
  colnames(rna) <- sample.labels
  rna.specs <- data.frame(t(apply(rna, 1, FUN = function(r){
    s <- r/sqrt(sum(r^2, na.rm=T))
    return(s)
  })))
  ifelse(protein=="none", {
    S <- rna.specs
  },
  {
    if(nrow(protein)!=nrow(rna)){
      stop("RNA and protein matrices have a different number of rows",call.=F)
    }
    if(sum(rownames(protein)!= rownames(rna)) > 0){
      stop("Gene names do not match between RNA and protein",call.=F)
    }
    if(sum(names(table(colnames(protein))) %in% sample.labels) != length(sample.labels)){
      stop("Protein columns and sample labels do not match",call.=F)
    }
    protein <- data.frame(sapply(sample.labels, FUN=function(x){
      exp = protein[,colnames(protein)==x]
      if(is.null(dim(exp))){
         exp = exp
      } else{
        exp = rowMeans(exp, na.rm=T)
      }
    colnames(protein) <- sample.labels
    protein.specs <- data.frame(t(apply(protein, 1, FUN = function(r){
      s <- r/sqrt(sum(r^2, na.rm=T))
      return(s)})))
    S <- data.frame(sapply(sample.labels, FUN=function(x){
      weight.rna*rna.specs[,x] + weight.protein*protein.specs[,x]}))
  })
  rownames(S) <- rownames(rna.specs)
  colnames(S) <- sample.labels
  S <- na.omit(S)
  return(S)
}
