### DESCREAPTION:
# expSummary() receives GSE experiment id from user and return experiment_summary dataframe object.
# This script export:
# 1. The minimal required fields data for a microarray experiment_summary as a csv.
# 2. RNA sequencing datasets as a csv in case the given GSE contains the dataset.


library(rentrez)
library(GEOquery)
library(magrittr) 
library(dplyr)    
library(XML)
library(jsonlite)
library(purrr)
library(data.table)

## Gets GSE experiment id from user via console
var = readline("Please enter a 'GSE' experiment id: ")
if (! startsWith(var,"GSE")) {
  var = readline("This is not a valied input. Please enter a 'GSE' experiment id: ")
}


expSummary <- function(var) {
  ##  microarray experiment_summary
  gds_search <- entrez_search(db="gds", term=paste(var,"[GEO Accession]",sep=""), retmax = 1)
  s<-entrez_summary(db="gds", id=gds_search$ids)
  
  df_exp<- data.frame(experiment_id=s[["accession"]],
                   platform_id=paste("GPL",s[["gpl"]],sep = ""), 
                   suppfile=s["suppfile"], 
                   ftplink=s["ftplink"])
  affy = "FALSE"
  if (grepl("CEL", s[["suppfile"]], fixed=TRUE)) {
    affy = "TRUE"
  }
  df_exp$affy<-affy
  
  rnaseq = "FALSE"
  if (length(s[["gdstype"]]) > 0) {
    if (grepl("Expression profiling by high throughput sequencing", s[["gdstype"]], fixed=TRUE)) {
      rnaseq = "TRUE"
    }
  }
  df_exp$rnaseq<-rnaseq
  
  
  
  ## rna sequencing (deep sequencing) datasets if exist in the given GSE
  sra_search <- entrez_search(db="sra", term=paste(s[["extrelations"]][["targetobject"]],"[sra]",sep=""),retmax=gds_search$count)
  
  if (length(sra_search$ids) > 0) {
    ef <-entrez_fetch(db = "sra", id=sra_search$ids, rettype = "runinfo",retmode = "xml",parsed = TRUE)
    x <- xmlToList(ef)
    
    dt_list <- map(x, as.data.table)
    dt <- rbindlist(dt_list, fill = TRUE, idcol = T)
    dt<-select(dt,-.id)
    write.csv(dt,"metadata.csv", row.names = TRUE)
    
  }
  return(df_exp)
}

# experiment_summary object
experiment_summary<-expSummary(var)

write.csv(experiment_summary,"experiment_summary.csv", row.names = TRUE)


