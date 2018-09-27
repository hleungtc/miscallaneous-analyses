require(readxl)
require(ggplot2)
require(plyr)
require(dplyr)
require(reshape2)

lm_eqn <- function(df){
  m <- lm(CT ~ log10(Quantity), df);
  a = format(coef(m)[1], digits = 2)
  b = format(coef(m)[2], digits = 2)
  r2 = format(summary(m)$r.squared, digits = 3)
  eqn <- paste("CT = ", a, "+", "(", b, ")","*log10[Copy No.]", ", R2 = ", r2, sep="")
  return(eqn)
}

convert_copy <- function(value, curve){
  copy <- 10^((value - coef(lm(CT ~ log10(Quantity), curve))[1])/coef(lm(CT ~ log10(Quantity), curve))[2])
  return(copy)
}

import_qpcr <- function(path = "./", sheet = "Results", dilution = 1){
  RAWFILE <- data.frame(read_xls(path = path, sheet = sheet, skip = 7, col_names = T))
  RAWFILE <- RAWFILE[,c(1:5,7,10)]
  colnames(RAWFILE) <- c("Well", "Sample", "Target", "Task", "Reporter", "CT", "Quantity")
  RAWFILE[,"CT"][RAWFILE[,"CT"] %in% "Undetermined"] <- NA
  RAWFILE[,"CT"] <- as.numeric(RAWFILE[,"CT"])
  QPCR_DATA <- list()
  QPCR_DATA[["Raw"]] <- RAWFILE
  QPCR_DATA[["Standard"]] <- RAWFILE[which(RAWFILE$Task %in% "STANDARD"),]
  
  QPCR_DATA[["Samples"]] <- RAWFILE[which(RAWFILE$Task %in% "UNKNOWN"),] %>%
    group_by(Sample) %>%
    summarise(CT.mean=mean(CT))
  QPCR_DATA[["Samples"]] <- data.frame(QPCR_DATA[["Samples"]])
  QPCR_DATA[["Samples"]]$Quantity <- convert_copy(QPCR_DATA[["Samples"]]$CT.mean, QPCR_DATA[["Standard"]])
  QPCR_DATA[["Samples"]]$Quantity.undiluted <- QPCR_DATA[["Samples"]]$Quantity/dilution
  
  QPCR_DATA[["Curve"]] <- ggplot(data = QPCR_DATA[["Standard"]], aes(x = Quantity, y = CT)) + geom_point() +
    scale_x_log10() + geom_smooth(method="lm") +
    ggtitle(label=lm_eqn(QPCR_DATA[["Standard"]])) + theme_bw() 
  return(QPCR_DATA)
}

normalize_qpcr <- function(QPCR_DATA, weight, 
                           elution_vol = 100, 
                           lysate_proportion = 0.436, 
                           irt_proportion = 0.474){
  QPCR_DATA[["Samples"]]$Normalized.quantity <- QPCR_DATA[["Samples"]]$Quantity.undiluted * 
    lysate_proportion * irt_proportion/(elution_vol * weight)
  return(QPCR_DATA[["Samples"]])
}