
# this function returns the input dataframe with additional columns for plotting
# data contains adj.value, log2FC, `Gene Symbol`, Protein
# Define FCvalue here for plotting Protein_label_select

plot_label <- function(data, FCvalue, FDRsig, FDRplot, direction){
  # FCvalue is not yet log-scaled
  
  #checking function arguments
  if(missing(data) | missing(FCvalue) | missing(FDRsig))
    stop('data, FCvalue, FDRsig is missing')
  
  # use FDRsig if FDRplot is not supplied
  if(missing(FDRplot))
    FDRplot <- FDRsig
  # use direction = both if not providing
  if(missing(direction))
    direction <- 'both'
  
  FCvalue <- log2(FCvalue)
  data_tmp <- data
  
  if (direction == 'both'){
    
    for (i in seq_along(data_tmp$Protein)){
      # adj.pvalue < FDR
      data_tmp$sig[i] <- ifelse(data_tmp$adj.pvalue[i] < FDRsig, 1, 0)
      
      # up or down regulation
      data_tmp$reg[i] <- ifelse(data_tmp$sig[i] == 1, ifelse(data_tmp$log2FC[i] > 0, 'up', 'down'), 'none')
      
      # get Gene.Symbol for significantly enriched protein.
      data_tmp$Protein_label[i] <- ifelse(data_tmp$sig[i] == 1, data_tmp$`Gene Symbol`[i], NA)
      
      # to only display certain foldchange (both + and -)
      # to only display above certain FDRplot
      data_tmp$Protein_label_select[i] <- 
        ifelse(data_tmp$adj.pvalue[i] < FDRplot & #adj.pvalue cutoff here
               abs(data_tmp$log2FC)[i] > FCvalue, #both +/- FC
               data_tmp$`Gene Symbol`[i], NA)
    }
    return(data_tmp)} 
  
  else if (direction == 'pos') {
    for (i in seq_along(data_tmp$Protein)){
      # adj.pvalue < FDR and log2FC > 0
      data_tmp$sig[i] <- ifelse(data_tmp$adj.pvalue[i] < FDRsig & data_tmp$log2FC[i] > 0, 1, 0)
      
      # up or down regulation
      data_tmp$reg[i] <- ifelse(data_tmp$sig[i] == 1, ifelse(data_tmp$log2FC[i] > 0, 'up', 'down'), 'none')
      
      # get Gene.Symbol for significantly enriched protein.
      data_tmp$Protein_label[i] <- ifelse(data_tmp$sig[i] == 1, data_tmp$`Gene Symbol`[i], NA)
      
      # to only display certain foldchange (both + and -)
      # to only display above certain FDRplot
      data_tmp$Protein_label_select[i] <- 
        ifelse(data_tmp$adj.pvalue[i] < FDRplot & #adj.pvalue cutoff here
                 abs(data_tmp$log2FC)[i] > FCvalue, #both +/- FC
               data_tmp$`Gene Symbol`[i], NA)
    }    
    return(data_tmp)
  }
  
  else if (direction == 'neg') {
    for (i in seq_along(data_tmp$Protein)){
      # adj.pvalue < FDR and log2FC > 0
      data_tmp$sig[i] <- ifelse(data_tmp$adj.pvalue[i] < FDRsig & data_tmp$log2FC[i] < 0, 1, 0)
      
      # up or down regulation
      data_tmp$reg[i] <- ifelse(data_tmp$sig[i] == 1, ifelse(data_tmp$log2FC[i] > 0, 'up', 'down'), 'none')
      
      # get Gene.Symbol for significantly enriched protein.
      data_tmp$Protein_label[i] <- ifelse(data_tmp$sig[i] == 1, data_tmp$`Gene Symbol`[i], NA)
      
      # to only display certain foldchange (both + and -)
      # to only display above certain FDRplot
      data_tmp$Protein_label_select[i] <- 
        ifelse(data_tmp$adj.pvalue[i] < FDRplot & #adj.pvalue cutoff here
                 abs(data_tmp$log2FC)[i] > FCvalue, #both +/- FC
               data_tmp$`Gene Symbol`[i], NA)
    }    
    return(data_tmp)
  }
}

rank_label <- function(data, GOterm, FDRsig, FDRcutoff, skip){
  # GOterm for annotating cellular components
  # FDRcutoff can be 'pos' or 'neg' referring to plotting log2FC > 0 or < 0
  # skip.... annotate every skip protein (e.g., 4)
  
  data_tmp <- data
  
  data_tmp$GOterm <- ifelse(grepl(GOterm,
                               data_tmp$`Cellular Component`, ignore.case = TRUE),1, 0)
  
  # rank based on log2FC
  data_tmp <- data_tmp %>% arrange(log2FC) %>%
                           #arrange(desc(log2FC))
    dplyr::filter(adj.pvalue < FDRsig)
  
  if (FDRcutoff == 'pos') {
    data_tmp <- data_tmp %>% dplyr::filter(log2FC > 0)
  } else if (FDRcutoff == 'neg') {
    data_tmp <- data_tmp %>% dplyr::filter(log2FC < 0)    
  } else {return('FDRcutoff must be pos or neg')}

  data_tmp$Rank <- seq_along(data_tmp$Protein)
  
  for (i in seq_along(data_tmp$Protein)){
    # simplify by taking the first entry of Gene.Symbol
    data_tmp$Protein_label[i] <- 
      str_extract(data_tmp$Protein_label[i],"^([^;])+")
  }

  for (i in seq_along(data_tmp$Protein)){
    # simplify by taking the first entry of Gene.Symbol
    data_tmp$Protein_label_select2[i] <- 
      ifelse(
        (data_tmp$GOterm[i] == 0 & data_tmp$Rank[i] > 50), 
        data_tmp$Protein_label[i],
        ifelse( (data_tmp$GOterm[i] == 1 & data_tmp$Rank[i] %% skip == 0), #plot every skip
                data_tmp$Protein_label[i], '')
      )
  }
  
  #the first protein show up
  if (FDRcutoff == 'pos') {
    data_tmp$Protein_label_select2[length(seq_along(data_tmp$Protein))] <- 
      data_tmp$Protein_label[length(seq_along(data_tmp$Protein))]
  } else if (FDRcutoff == 'neg') {
    data_tmp$Protein_label_select2[1] <- 
      data_tmp$Protein_label[1]    
  } else {return('FDRcutoff must be pos or neg')}
  
  return(data_tmp)
}
