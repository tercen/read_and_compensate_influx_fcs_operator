library(BiocManager)
library(tercen)
library(dplyr)
library(flowCore)
library(flowWorkspace)
library(fuzzyjoin)
library(stringr)

# http://localhost:5402/admin/w/9bc1fd64ee4d8642eb4c61d22c1d7e8b/ds/c40928c8-b24e-4c51-a9eb-6f94e49f7d75
options("tercen.workflowId" = "9bc1fd64ee4d8642eb4c61d22c1d7e8b")
options("tercen.stepId" =     "c40928c8-b24e-4c51-a9eb-6f94e49f7d75")
##

sort_to_data = function(filename, 
                        comp=FALSE, comp_df=NULL,
                        transform="none") {
  INCLUDE <- c("Well", "FSC", "SSC", "*", "TIME", "Tray", "Well")
  
  # Read FCS file using flowCore::read.FCS
  flowfile <- read.FCS(filename,
                       transformation=FALSE) 
  
  if (keyword(flowfile, "INDEXSORTPOSITIONS") %in% keyword(flowfile) == TRUE){
    
    #Creates a DataFrame of the sorted events
    flowdata <- as.data.frame(exprs(flowfile))
    indexed_flowdata <- flowdata[flowdata$'Sort Result Bits' >0,]
    
    # Filter out correct rows
    indexed_flowdata = indexed_flowdata %>% 
      select(contains(INCLUDE)) %>% 
      rename_all(~str_replace_all(., "\\*",""))
    
    # Perform transformation if needed
    if (transform == "biexponential") {
      trans_f = flowWorkspace::flowjo_biexp()
      trans_flow_data = indexed_flowdata %>% select(-contains(c('Well', 
                                                              'TIME', 
                                                              'Tray')))
      
      for (c in colnames(trans_flow_data)) {
        indexed_flowdata[, c] = trans_f(indexed_flowdata[, c])
      }
    }
    
    # Re-create flowFrame object
    indexed_fcs = flowFrame(exprs = as.matrix(indexed_flowdata))
    
    # Perform compensation
    if (comp) {
      if (is.null(comp_df)) {
      indexed_fcs = compensate(indexed_fcs, spillover(flowfile)$SPILL)
      } else {
        colnames(comp_df) = colnames(spillover(flowfile)$SPILL)
        indexed_fcs = compensate(indexed_fcs, comp_df)
      }
    }
    
    # Final DF
    indexed = as.data.frame(exprs(indexed_fcs))
    
    #Creates a dataframe of the well positions
    wellorder<-strsplit(keyword(flowfile,"INDEXSORTPOSITIONS")[[1]], ",")
    wellID <- wellorder[[1]][seq(1, length(wellorder[[1]]), 3)]
    x<-wellorder[[1]][seq(2, length(wellorder[[1]]), 3)]
    y<-wellorder[[1]][seq(3, length(wellorder[[1]]), 3)]
    
    
    well_df<-data.frame(wellID,x,y)
    well_df[, 2] <- as.numeric(as.character(well_df[, 2]))
    well_df[, 3] <- as.numeric(as.character(well_df[, 3]))
    colnames(well_df) <- c("Well","Tray X", "Tray Y")
    
    #Software does not store the index positions correctly so we use dplyr and fuzzylogic packages
    indexdata<-difference_inner_join(well_df,indexed,
                                     by=c("Tray X", "Tray Y"),
                                     max_dist = 3)
    
    colnames(indexdata)<-c('Well','Tray X (kw)','Tray Y (kw)',
                           colnames(indexed)) 
    #Reading the fcs file again to get the column names slows the script down by 1/3.  
    #To avoid this pre populate this with the col names.  
    #Doing it this way makes the script more resistant to changes in the .fcs files
    
    # Add document ID and CI
    indexdata %>%
      mutate_if(is.logical, as.character) %>%
      mutate_if(is.integer, as.double) %>%
      mutate(.ci = rep_len(0, nrow(.))) %>%
      mutate(filename = rep_len(basename(filename), nrow(.)))
    
  } else {
    stop("FCS file is not a BD InFlux file.")
  }
  
}

ctx = tercenCtx()

if (!any(ctx$cnames == "documentId")) stop("Column factor documentId is required") 

# Setup operator properties
compensation <- TRUE
if(!is.null(ctx$op.value("compensation"))) type <- ctx$op.value("compensation")

transformation <- "biexponential"
if(!is.null(ctx$op.value("transformation"))) comparison <- ctx$op.value("transformation")

#1. extract files
df <- ctx$cselect()

docId = df$documentId[1]
doc = ctx$client$fileService$get(docId)
filename = tempfile()
writeBin(ctx$client$fileService$download(docId), filename)
on.exit(unlink(filename))

# unzip if archive
if(length(grep(".zip", doc$name)) > 0) {
  tmpdir <- tempfile()
  unzip(filename, exdir = tmpdir)
  f.names <- list.files(tmpdir, full.names = TRUE, 
                        pattern="\\.fcs$", ignore.case=TRUE)
  csv.names <- list.files(tmpdir, full.names = TRUE, 
                          pattern="\\.csv$", ignore.case=TRUE)
  
  if (length(csv.names) == 0) { 
    comp.df <- NULL
  } else {
    comp.df <- read.csv(csv.names[1], check.names=FALSE)[-1]
  }
  
} else {
  f.names <- filename
  comp.df <- NULL
}

# check FCS
if(any(!isFCSfile(f.names))) stop("Not all imported files are FCS files.")

assign("actual", 0, envir = .GlobalEnv)
task = ctx$task

#2. convert them to FCS files
f.names %>%
  lapply(function(filename){
    # pass CSV compensation matrix or NULL
    data = sort_to_data(filename, 
                        comp=compensation, comp_df=comp.df,
                        transform=transformation)
    if (!is.null(task)) {
      # task is null when run from RStudio
      actual = get("actual",  envir = .GlobalEnv) + 1
      assign("actual", actual, envir = .GlobalEnv)
      evt = TaskProgressEvent$new()
      evt$taskId = task$id
      evt$total = length(f.names)
      evt$actual = actual
      evt$message = paste0('processing FCS file ' , filename)
      ctx$client$eventService$sendChannel(task$channelId, evt)
    } else {
      cat('processing FCS file ' , filename)
    }
    data
  }) %>%
  bind_rows() %>%
  ctx$addNamespace() %>%
  ctx$save()