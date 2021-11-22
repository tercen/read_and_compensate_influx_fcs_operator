library(BiocManager)
library(tercen)
library(dplyr)
library(flowCore)
library(flowWorkspace)
library(fuzzyjoin)
library(stringr)

sort_to_data = function(path, display_name = "", 
                        compensation = FALSE, comp_df = NULL,
                        transformation = "none") {

  INCLUDE <- c("Well", "FSC", "SSC", "TIME", "Tray", "Well")
  
  # Read FCS file using flowCore::read.FCS
  flowfile <- read.FCS(path,
                       transformation=TRUE) 
  
  if (keyword(flowfile, "INDEXSORTPOSITIONS") %in% keyword(flowfile) == TRUE){
    
    #Creates a DataFrame of the sorted events
    flowdata <- as.data.frame(exprs(flowfile))
    indexed_flowdata <- flowdata[flowdata$'Sort Result Bits' >0,]
    
    # Filter cols
    filter_cols = colnames(indexed_flowdata %>% 
                             select(contains(c("*"))) %>% 
                             rename_all(~str_replace_all(., "\\*","")))

    # Filter out correct rows
    indexed_flowdata = indexed_flowdata %>% 
      select(-contains(c("*"))) %>%
      select(contains(c(INCLUDE, filter_cols)))
      
      #rename_all(~str_replace_all(., "\\*",""))
  
    # Perform transformation if needed
    if (transformation == "biexponential") {
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
    if (compensation) {
      if (is.null(comp_df)) {
        indexed_fcs = compensate(indexed_fcs, spillover(flowfile)$SPILL)
      } else {
        #comp_df = t(inv(as.matrix(comp_df)))
        colnames(comp_df) = colnames(spillover(flowfile)$SPILL)
        #print(comp_df)
        indexed_fcs = compensate(indexed_fcs, compensation(as.matrix(comp_df)))
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
      mutate(filename = rep_len(basename(display_name), nrow(.)))
    
  } else {
    stop("FCS file is not a BD InFlux file.")
  }
  
}

ctx = tercenCtx()

if (!any(ctx$cnames == "documentId")) stop("Column factor documentId is required") 

# Setup operator properties
compensation_param <- "none"
if(!is.null(ctx$op.value("compensation"))) compensation_param <- ctx$op.value("compensation")

transformation_param <- "biexponential"
if(!is.null(ctx$op.value("transformation"))) transformation_param <- ctx$op.value("transformation")

#1. extract files
df <- ctx$cselect()

docId = df$documentId[1]
doc = ctx$client$fileService$get(docId)
display_name = doc$name
filename = tempfile()
writeBin(ctx$client$fileService$download(docId), filename)
on.exit(unlink(filename))

# unzip if archive
if(length(grep(".zip", doc$name)) > 0) {
  tmpdir <- tempfile()
  unzip(filename, exdir = tmpdir)
  f.names <- list.files(tmpdir, full.names = TRUE, 
                        pattern="\\.fcs$", ignore.case=TRUE)
  comp.names <- list.files(tmpdir, full.names = TRUE, 
                           pattern="\\.comp$", ignore.case=TRUE)
  
  fcs_files = c()
  comp_files = c()
  
  for (f in f.names) {
    f_name = gsub(pattern = "\\.fcs$", "", basename(f), ignore.case=TRUE)
    comp_file = "none"
    
    for (c in comp.names) {
      c_name = gsub(pattern = "\\.comp$", "", basename(c), ignore.case=TRUE)
      if (f_name == c_name) comp_file = c_name
    }
    
    fcs_files = c(fcs_files, f)
    comp_files = c(comp_files, c)
    
  }
  
  data <- data.frame(fcs=fcs_files, comp=comp_files)
  
  display_name <- fcs_files

} else {
  data <- data.frame(fcs=c(filename), comp=c(""))
}

# check FCS
if(any(!isFCSfile(data$fcs))) stop("Not all imported files are FCS files.")

assign("actual", 0, envir = .GlobalEnv)
task = ctx$task

#2. convert them to FCS files
data %>%
  apply(1, function(row){
    fcs = row[1]
    comp = row[2]
    
    if (compensation_param != "none") {
      comp.df <- read.csv(comp, check.names=FALSE)[-1]
      # pass CSV compensation matrix or NULL
      data = sort_to_data(path = fcs, display_name = fcs,
                          compensation = TRUE, comp_df = comp.df,
                          transformation = transformation_param)
    } else {
      data = sort_to_data(path = fcs, display_name = fcs ,
                          compensation = FALSE, transform = transformation_param)
    }
    
    if (!is.null(task)) {
      # task is null when run from RStudio
      actual = get("actual",  envir = .GlobalEnv) + 1
      assign("actual", actual, envir = .GlobalEnv)
      evt = TaskProgressEvent$new()
      evt$taskId = task$id
      evt$total = length(data$fcs)
      evt$actual = actual
      evt$message = paste0('processing FCS file ' , fcs)
      ctx$client$eventService$sendChannel(task$channelId, evt)
    } else {
      cat('processing FCS file ' , fcs)
    }
    data
  }) %>%
  bind_rows() %>%
  ctx$addNamespace() %>%
  ctx$save()