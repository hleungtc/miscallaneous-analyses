library(plyr)
library(dplyr)
library(tidyr)
library(readxl)

# TABULATES ALL SCFA DATA
# Specify location of .xls chromatograms using path argument

parse_scfa <- function(subset_pattern=NULL, path = "./"){
  file_list <- list.files(pattern="*.xls", path)
  sequence_list <- list()
  for(i in 1:(length(file_list))){
    sequence_list[[file_list[i]]] <- data.frame(read_xls(path=paste(path, file_list[[i]], sep=""), skip = 39,sheet = "Integration"))
  }
  samples <- list()
  for(j in 1:length(file_list)){
    samples[[file_list[j]]] <- sequence_list[[file_list[j]]][-c(1,2),c("Peak.Name","Amount")]
    samples[[file_list[j]]]$Amount <- as.numeric(as.character(samples[[file_list[j]]]$Amount))
    colnames(samples[[file_list[j]]])[2] <- file_list[j]}
  for(l in 1:length(file_list)){
    samples[[file_list[l]]] <- samples[[file_list[l]]][which(!is.na(samples[[file_list[l]]]$Peak.Name)),]
  }
  if(is.null(subset_pattern)){subset_samples <- samples}else{
    subset_samples <- samples[which(grepl(pattern = subset_pattern, x = names(samples)))]}
  table <- subset_samples[[1]]
  for(k in 2:length(subset_samples)){
    table <- join(table, subset_samples[[k]], "Peak.Name")
  }
  
  # REMOVES UNKNOWN PEAKS AND INTERNAL STANDARD
  table <- table[which(!table$Peak.Name %in% ""),]
  colnames(table) <- gsub(pattern = ".xls", replacement = "", x = colnames(table))
  table <- table[which(!table$Peak.Name %in% c("Component 2", "2-Ethylbutyric Acid")),]
  
  rownames(table) <- table$Peak.Name
  table$Peak.Name <- NULL
  parsed_table <- data.frame(t(table))
  parsed_table$Sample <- row.names(parsed_table)
  parsed_table[is.na(parsed_table)] <- 0
  return(parsed_table)
}

# NORMALIZES THE REPORTED mM CONCENTRATIONS (GENOTEK)
normalize_scfa <- function(table, full_tube, tube_subtract, pipet_vol = 0.2, tube_w_buffer = 13.74, tube_wt = 11.80, buffer_wt = 1.94, buffer_vol = 2.9){
  feces_wt <- full_tube - tube_w_buffer
  # PIPET_WT = WT OF SAMPLE USED FOR SCFA ANALYSIS
  pipet_wt <- c()
  for(i in 1:length(feces_wt)){
    # NORMAL SAMPLE, WEIGHT OF TUBE > WEIGHT OF STOCK TUBE + 0.25 g
    if(full_tube[i] > tube_w_buffer + 0.25){
      pipet_wt[i] <- (full_tube[i] - tube_subtract[i])*(feces_wt[i]/(buffer_wt+feces_wt[i]))/1000
    # COMPRIMISED SAMPLE (WITHOUT BUFFER), WEIGHT OF TUBE < WEIGHT OF STOCK TUBE + 0.25 g
      }else{
        pipet_wt[i] <- (full_tube[i] - tube_subtract[i])/1000}
  }
  
  normalized_table <- table * (buffer_vol/1000)/pipet_wt
  return(normalized_table)
}