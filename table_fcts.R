#### functions toCatTab and from2CatString for reporting on categorical variables ####

toCatTab <- function(x, ncat, bycat)
{

  total <- sum(table(x))
  toShow <- ""
  for (k in (1:ncat)) 
  {
    cnts.x <- sum(table(x)[((k-1)*bycat + 1):(k*bycat)]) 
    perc.x <- round(sum(table(x)[((k-1)*bycat + 1):(k*bycat)])*100/total, digits = 0) 
    if (bycat > 1 ) 
    {
      toShow <- paste0(toShow, names(table(x))[(k-1)*bycat + 1], "-", 
                       names(table(x))[k*bycat], 
                       ": ", cnts.x, " (", perc.x, "%); ")
    } else 
    {
      toShow <- paste0(toShow, names(table(x))[(k-1)*bycat + 1], 
                       ": ", cnts.x, " (", perc.x, "%);")
    }
  }
  return(toShow)
}

toCatLines <- function(x, ncat, bycat, vName)
{
  
  total <- sum(table(x))
  toShow <- c()
  if (grepl("charlson", vName) ) bycat <- 3
  if (bycat == 1)
  {
    for (k in (1:ncat)) 
    {
      cnts.x <- sum(table(x)[((k-1)*bycat + 1):(k*bycat)]) 
      perc.x <- round(sum(table(x)[((k-1)*bycat + 1):(k*bycat)])*100/total, 
                      digits = 0) 
      cnts.x[is.na(cnts.x)] <- perc.x[is.na(perc.x)] <- 0
      
      toShow <- rbind(toShow, 
                      c(paste0(vName, ":", names(table(x))[(k-1)*bycat + 1]), 
                          paste0(cnts.x, " (", perc.x, "%)") ))
      
    }
    
  } else
    { for (k in (1:floor(ncat/bycat)))
    {
      cnts.x <- sum(table(x)[((k-1)*bycat + 1):(k*bycat)]) 
      perc.x <- round(sum(table(x)[((k-1)*bycat + 1):(k*bycat)])*100/total, 
                      digits = 0) 
      cnts.x[is.na(cnts.x)] <- perc.x[is.na(perc.x)] <- 0
      toShow <- rbind(toShow,
                      c(paste0(vName, ":",
                               names(table(x))[(k-1)*bycat + 1], "-", 
                               names(table(x))[k*bycat]), 
                        paste0(cnts.x, " (", perc.x, "%) ")))
      }} 
      
  toShow <- as.data.frame(toShow)
  colnames(toShow) <- c("varName", "median")
  
  return(toShow)
}

from2CatString <- function(x, ncat, bycat, principalLevel)
  # here ncat = 2 and bycat = 1
{
  UseTab <- table(x)
  total <- sum(UseTab)
  toShow <- ""
  posPrinciplaLevel <- which(names(UseTab) == principalLevel)
  
  numToShow <-  UseTab[posPrinciplaLevel]
  percentToShow <- round(100*numToShow/total, digits = 0)
  
  if (length(percentToShow)==0) {percentToShow <- 0; numToShow <- 0}
  
  if (tolower(principalLevel) %in% c("female", "male", "woman", "man", "women", "men")){
    toShow <- paste0(#╧ principalLevel, ": ",
                     numToShow, " (",percentToShow,  "% )")
    
  } else toShow <- paste0(numToShow, " (",  percentToShow, "% )")
  return(toShow) 
}

descrTab <- function(df, quVec, CIcov = 0.95, 
                     allVars = NULL, catVars, binVars_yes, 
                     binVars_no, binVars_sex)
{
  library(psych)
  library(gtools)
  
  if (is.null(allVars)) allVars <- colnames(df)
  
  Tb <- psych::describe(df[ ,allVars], IQR = TRUE, quant = quVec)
  Tb <- as.data.frame(Tb)
  Tb$varName <- Tb$varname <- gsub("\\*", "", rownames(Tb))

  
  Tb$mean <- round(Tb$mean, digits = 2)
  Tb$sd <- round(Tb$sd, digits = 2)
  Tb$median <- round(Tb$median, digits = 2)
  Tb$IQR <- round(Tb$IQR, digits = 2)
  Tb$min <- round(Tb$min, digits = 2)
  Tb$max <- round(Tb$max, digits = 2)
  
  quNames <- paste0("Q", quVec)
  Tb[ ,quNames] <- round(Tb[ ,quNames], digits = 2)
  
  normQ <- qnorm(1 - (1-CIcov)/2)
  Tb$mean <- paste0(Tb$mean, " (", 
                    round(Tb$mean - normQ*Tb$sd, digits = 2), ", ", 
                    round(Tb$mean + normQ*Tb$sd, digits = 2), ")")
  Tb$median <- paste0(Tb$median, " [",Tb[ ,quNames[1]],", ",
                      Tb[ ,quNames[2]], "]")
  
  for (v in c(catVars) ){
    Tb$mean[which(Tb$varname == v)] <-  
      Tb$median[which(Tb$varname == v)] <- 
      Tb$min[which(Tb$varname == v)] <- 
      Tb$max[which(Tb$varname == v)] <- ""
      Tb <- smartbind(Tb,
                      toCatLines(df[ ,v], ncat = length(table(df[ ,v])), 
                                 bycat = 1, vName = v))
    
  }
  
  for (v in c(binVars_yes) ){
    Tb$mean[which(Tb$varname == v)] <-  
      Tb$median[which(Tb$varname == v)] <- from2CatString(df[ ,v], ncat = length(table(df[ ,v])), 
                                                    bycat = 1, principalLevel = "yes")
      Tb$min[which(Tb$varname == v)] <- 
      Tb$max[which(Tb$varname == v)] <- ""
  }
  for (v in c(binVars_no) ){
    Tb$mean[which(Tb$varname == v)] <-  
      Tb$median[which(Tb$varname == v)] <- from2CatString(df[ ,v], ncat = length(table(df[ ,v])), 
                                                          bycat = 1, principalLevel = "no")
    Tb$min[which(Tb$varname == v)] <- 
      Tb$max[which(Tb$varname == v)] <- ""
  }
  for (v in c(binVars_sex) ){
    Tb$mean[which(Tb$varname == v)] <-  
      Tb$median[which(Tb$varname == v)] <- from2CatString(df[ ,v], ncat = length(table(df[ ,v])), 
                                                          bycat = 1, principalLevel = "male")
    Tb$min[which(Tb$varname == v)] <- 
      Tb$max[which(Tb$varname == v)] <- ""
  }
  
  
  rownames(Tb) <- Tb$varname <- Tb$varName
  Tb$varName <- gsub(".*:","", Tb$varName)
  Tb$n[which(is.na(Tb$n))] <- ""
  return(Tb)
}





