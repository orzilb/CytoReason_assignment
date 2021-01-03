### DESCREAPTION:
# module m includes 2 functions:
# 1. m$expSummary() 
#    input: GSE experiment id
#    return The minimal required fields data for a microarray experiment_summary as a dataframe.
# 2. m$metadata() 
#    input: GSE experiment id
#    return RNA sequencing datasets as a dataframe in case the given GSE contains the dataset, else return NULL.



m <- module ({
  
  import(rentrez)
  import(GEOquery)
  import(magrittr) 
  import(dplyr)    
  import(XML)
  import(jsonlite)
  import(purrr)
  import(data.table)
  
  
  ##  microarray experiment_summary
  
  expSummary <- function(gse) {
    gds_search <- entrez_search(db="gds", term=paste(gse,"[GEO Accession]",sep=""), retmax = 1)
    s<-entrez_summary(db="gds", id=gds_search$ids)
    
    summarry_table<- data.frame(experiment_id=s[["accession"]],
                        platform_id=paste("GPL",s[["gpl"]],sep = ""), 
                        suppfile=s["suppfile"], 
                        ftplink=s["ftplink"])
    affy = "FALSE"
    if (grepl("CEL", s[["suppfile"]], fixed=TRUE)) {
      affy = "TRUE"
    }
    summarry_table$affy<-affy
    
    rnaseq = "FALSE"
    if (length(s[["gdstype"]]) > 0) {
      if (grepl("Expression profiling by high throughput sequencing", s[["gdstype"]], fixed=TRUE)) {
        rnaseq = "TRUE"
      }
    }
    summarry_table$rnaseq<-rnaseq
    
    
    return(summarry_table)
      
    }
    
  ## rna sequencing (deep sequencing) datasets if exist in the given GSE
  
  metadata <- function(gse) {
    gds_search <- entrez_search(db="gds", term=paste(gse,"[GEO Accession]",sep=""), retmax = 1)
    s<-entrez_summary(db="gds", id=gds_search$ids)
    sra_search <- entrez_search(db="sra", term=paste(s[["extrelations"]][["targetobject"]],"[sra]",sep=""),retmax=gds_search$count)
    
    # if Metadata exist
    if (length(sra_search$ids) > 0) {
      ef <-entrez_fetch(db = "sra", id=sra_search$ids, rettype = "runinfo",retmode = "xml",parsed = TRUE)
      
      x <- xmlToList(ef)
      dt_list <- map(x, as.data.table)
      sra_meta <- rbindlist(dt_list, fill = TRUE, idcol = T)
      sra_meta<-select(sra_meta,-.id)
      return(sra_meta)
    }
    # if Metadata not exist
    return(NULL)
  }
})
