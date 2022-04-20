
## these functions are modified from P. Wilmart
## https://pwilmart.github.io/TMT_analysis_examples/IRS_validation.html

# function to do median intensity normalization
med_norm <- function(df, print_factors = TRUE) {
  # Normalizes each channel's sum to the average grand total
  # df: data frame of TMT data (one column for each channel)
  # print_factors: logical to control printing
  
  # compute norm factors and scale columns
  norm_facs <- mean(colMedians(as.matrix(df), na.rm=TRUE), na.rm=TRUE) / 
    colMedians(as.matrix(df),na.rm=TRUE)
  df_med  <- sweep(df, 2, norm_facs, FUN = "*")
  
  # print the normalization factors for QC check
  if (print_factors == TRUE) {
    cat("\nMedian Normalization Factors:\n ")
    cat(sprintf("%s - %0.3f\n", colnames(df), norm_facs))
  }
  
  df_med # return normalized data
}


## SL (sample loading) normalization
## Need to account for NA values in the following code
SL_Norm <- function(df, color = NULL, plot = FALSE) {
  # This makes each channel sum to the average grand total
  # df - data frame of TMT intensities
  # returns a new data frame with normalized values
  
  # compute scaling factors to make colsums match the average sum
  norm_facs <- mean(c(colSums(df, na.rm=TRUE)), na.rm=TRUE) / colSums(df, na.rm=TRUE)
  cat("SL Factors:\n", sprintf("%-5s -> %f\n", colnames(df), norm_facs))
  df_sl  <- sweep(df, 2, norm_facs, FUN = "*")
  
  # visualize results and return data frame
  if(plot == TRUE) {
    boxplot(log10(df_sl), col = color, notch = TRUE, main = "SL Normalized data")
  }
  df_sl
}


# median subtraction
med_subtraction <- function(df, print_factors = TRUE) {
  
  # compute global and column medians
  median_log2Int <- colMedians(as.matrix(df), na.rm = TRUE)
  median_global <- median(median_log2Int)
  median_diff <- median_global - median_log2Int 
  
  num_cols <- dim(df)[2]
  
  for (i in seq(1,num_cols)){
    df[,i] <- df[,i] + median_diff[i]
  }
  
  # print the normalization factors for QC check
  if (print_factors == TRUE) {
    cat("\nMedian Normalization Factors:\n ")
    cat(sprintf("%s : %0.3f\n", colnames(df), median_diff))}
  
  return(df)
}



