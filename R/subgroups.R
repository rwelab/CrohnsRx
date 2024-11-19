#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# # subgroup.R
# 
# Two tasks are performed:
# 1) The drug class order (from best performing:drug1 to worst performing:drug3)
#    for the three modeled drug classes - interleukin-12/23 (il12), integrin (intg), 
#    and tumor necrosis factor-alpha (tnfi) - are determined for each participant
#    (row). 
# 2) Two-sample t-tests are applied to drug pairs (drug1 vs drug2 and drug2 vs drug3)
#    to determine if a participants significantly prefers one or more drug classes
#    over another.
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

library(dplyr)

#------------------------------------------------------------------------------#

srs.subgroups <- function(data, model.list, placebo.label="Placebo", fit.label=".attrib", se.label=".se", df=NULL) {
  
  # columns
  col.rem   <- srs.begins_with(data, tolower(placebo.label))
  col.y     <- setdiff(srs.ends_with(data, fit.label), col.rem)
  col.se    <- setdiff(srs.ends_with(data, se.label),  col.rem)
  
  # vectors
  vec.lab   <- gsub(fit.label,"",col.y)  # get trt names
  if (is.null(df)) { vec.df <- srs.df(model.list) }
  if (!is.null(df)){ vec.df <- df }      # get model degrees of freedom
  
  # matrices
  mat.Yhat  <- unname(as.matrix(data[,col.y]))  # get fit columns
  mat.SE    <- unname(as.matrix(data[,col.se])) # get se columns
  N <- nrow(mat.Yhat)
  C <- ncol(mat.Yhat)
  mat.lab   <- matrix(rep(vec.lab, N), nrow=N, byrow=TRUE) # create label mat
  mat.df    <- matrix(rep(vec.df , N), nrow=N, byrow=TRUE) # create df mat
  
  # ranked matrices (Ex. Y2 > Y1 > Y3 == c(2, 1, 3))
  mat.idx   <- rank.index(mat.Yhat)
  rank.Yhat <- srs.rank(mat.Yhat, mat.idx)
  rank.SE   <- srs.rank(mat.SE  , mat.idx)
  rank.lab  <- srs.rank(mat.lab , mat.idx)
  rank.df   <- srs.rank(mat.df  , mat.idx)

  # two sample t-test (Y1 >= Y2; Y2 >= Y3, etc.)
  p.val <- man.t.test(rank.Yhat, rank.SE, rank.df)

  # subgroup formatting (ex. "Y2 = Y1 > Y3")
  p.val.ohe <- one.hot.str(p.val)
  mat.sbgrp <- mat.sandwich(rank.lab, p.val.ohe)
  sbgrp     <- row.str.concat(mat.sbgrp)
  colnames(rank.lab)  <- paste0("drug", 1:C)
  colnames(p.val)     <- paste0("p",    1:(C-1))
  colnames(p.val.ohe) <- paste0("p",    1:(C-1),".ohe")
  
  data <- cbind(data, rank.lab, p.val, p.val.ohe)
  data$Subgroup <- sbgrp
  return(data)
}

srs.summarise <- function(data) {
  data %>% 
    group_by(Subgroup) %>%
    summarize(count = n()) %>% 
    arrange(desc(count))
}

#------------------------------------------------------------------------------#
# Functions for ranking predicted treatment efficacy

rank.index <- function(mat) { 
  t(apply(-mat, 1, order)) 
}

srs.rank <- function(mat.data, mat.idx) {
  nr <- nrow(mat.data)
  nc <- ncol(mat.data)
  mat.rank <- matrix(NA, nrow=nr, ncol=nc)
  for (i in 1:nr) {
    mat.rank[i,] <- mat.data[i, mat.idx[i,]]
  }
  return(mat.rank)
}

#------------------------------------------------------------------------------#
# Functions for calculating treatment preference

srs.df <- function(model.list, placebo.label="Placebo") {
  trt.list <- setdiff(names(model.list), tolower(placebo.label))
  df <- c()
  for (trt in trt.list) {
    df <- c(df, nrow(model.list[[trt]]@frame) - (ncol(model.list[[trt]]@frame) - 1))
  }
  return(df)
}

man.t.test <- function(Y, SE, df) {
  
  if (nrow(Y) != nrow(SE)) {
    stop("Y and SE must have the same number of rows.")
  }
  
  nr <- nrow(Y); nc <- ncol(Y)
  p.vals <- matrix(NA, nrow = nr, ncol = nc - 1)
  
  for (i in 1:(nc - 1)) {
    N1  <-  Y[ ,i]; N2  <-  Y[ ,i+1]
    D1  <- SE[ ,i]; D2  <- SE[ ,i+1]
    DF1 <- df[ ,i]; DF2 <- df[ ,i+1]
    
    p.vals[,i] <- 2* pt(abs(N1 - N2) / sqrt(D1^2 + D2^2), df = DF1+DF2, lower.tail=F)
  }
  return(p.vals)
}

one.hot.str <- function(mat) {
  ifelse(mat < 0.05, ">", "=")
}

one.hot.num <- function(mat) {
  ifelse(mat < 0.05, 1, 0)
}

#------------------------------------------------------------------------------#
# Functions for formatting data

srs.ends_with <- function(data, pattern) {
  grep(paste0("\\",pattern,"$"), names(data), value = TRUE)
}

srs.begins_with <- function(data, pattern) {
  grep(paste0("^",pattern), names(data), value = TRUE)
}

mat.sandwich <- function(mat1, mat2) {
  
  nr1 <- nrow(mat1); nc1 <- ncol(mat1)
  nr2 <- nrow(mat2); nc2 <- ncol(mat2)
  
  if (nr1 != nr2) {
    stop("Input matrices must have the same number of rows.")
  }
  
  if (nc1 != nc2+1) {
    stop("Input matrices must have C and C-1 columns respectively.")
  }
  
  res <- matrix(NA, nrow = nr1, ncol = (nc1+nc2))
  for (i in 1:nc1) {
    res[,2*i-1] <- mat1[,i]
    if (i <= nc2) { res[,2*i] <- mat2[,i] }
  }
  return(res)
}

row.str.concat <- function(mat) {
  matrix(apply(mat, 1, paste, collapse = " "), nrow=nrow(mat), byrow=TRUE)
}

#------------------------------------------------------------------------------#

parse_ranking <- function(ranking) {
  # Split by ">" to handle greater relationships
  groups <- strsplit(ranking, ">")[[1]]
  parsed <- list()
  rank <- 0  # Highest rank starts at 0, increases with each ">" group
  
  for (group in groups) {
    line <- strsplit(trimws(group), "=")[[1]]  # Handle equal relationships
    for (drug in line) {
      parsed[[trimws(drug)]] <- rank  # Assign rank to each drug
    }
    rank <- rank + 1
  }
  
  return(parsed)
}

is_logically_equivalent <- function(ranking1, ranking2) {
  parsed1 <- parse_ranking(ranking1)
  parsed2 <- parse_ranking(ranking2)
  drugs <- names(parsed1)
  
  # Check pairwise consistency
  for (d1 in drugs) {
    for (d2 in drugs) {
      if (d1 != d2) {
        # Compare the relative ranking of (d1, d2) in both parsed results
        rel1 <- parsed1[[d1]] - parsed1[[d2]]
        rel2 <- parsed2[[d1]] - parsed2[[d2]]
        if (rel1 * rel2 < 0) { return(0) }
      }
    }
  }
  return(1)  # All pairwise comparisons are consistent
}

srs.isEqual <- function(target, current) {
  isEqual <- mapply(is_logically_equivalent, target, current)
  res <- data.frame(target, current, isEqual)
  return(res)
}

isEqual.pct <- function(target, current) {
  res <- srs.isEqual(target, current)
  mean(res$isEqual)
} 

noPref.pct <- function(data) {
  nrow(data %>% filter(p1.ohe=='=' & p2.ohe=='=')) / nrow(data)
}

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
