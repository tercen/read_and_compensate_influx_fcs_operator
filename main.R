library(BiocManager)
library(tercen)
library(dplyr)
library(flowCore)
library(fuzzyjoin)

fcs_sort_to_data = function(filename) {
  # Read FCS file using flowCore::read.FCS
  flowfile <- read.FCS(filename, alter.names = TRUE)
  flowfile <- compensate(flowfile, spillover(flowfile)$SPILL)
  
  if (keyword(flowfile, "INDEXSORTPOSITIONS") %in% keyword(flowfile) == TRUE) {
    #Creates a dataframe of the sorted events
    ofintrest <- exprs(flowfile)
    dataframe <- as.data.frame(ofintrest)
    indexed <- dataframe[dataframe$'Sort.Result.Bits' > 0, ]
    
    #Creates a dataframe of the well positions
    wellorder <-
      strsplit(keyword(flowfile, "INDEXSORTPOSITIONS")[[1]], ",")
    wellID <- wellorder[[1]][seq(1, length(wellorder[[1]]), 3)]
    x <- wellorder[[1]][seq(2, length(wellorder[[1]]), 3)]
    y <- wellorder[[1]][seq(3, length(wellorder[[1]]), 3)]
    
    
    well_df <- data.frame(wellID, x, y)
    well_df[, 2] <- as.numeric(as.character(well_df[, 2]))
    well_df[, 3] <- as.numeric(as.character(well_df[, 3]))
    colnames(well_df) <- c("Well", "Tray.X", "Tray.Y")
    
    #Sortware does not store the index positions correctly so we use thr dplyr and fuzzylogic packages to correct for this
    indexdata <-
      difference_inner_join(well_df,
                            indexed,
                            by = c("Tray.X", "Tray.Y"),
                            max_dist = 3)
    colnames(indexdata) <- c('Well', 'Tray X (kw)', 'Tray Y (kw)',
                             colnames(read.FCS(filename))) #Reading the fcs file again to get the column names slows the script down by 1/3.  To avoid this pre populate this with the col names.  Doing it this way makes the script more resistant to changes in the .fcs files
    
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

if (!any(ctx$cnames == "documentId"))
  stop("Column factor documentId is required")

#1. extract files
df <- ctx$cselect()

docId = df$documentId[1]
doc = ctx$client$fileService$get(docId)
filename = tempfile()
writeBin(ctx$client$fileService$download(docId), filename)
on.exit(unlink(filename))

# unzip if archive
if (length(grep(".zip", doc$name)) > 0) {
  tmpdir <- tempfile()
  unzip(filename, exdir = tmpdir)
  f.names <- list.files(tmpdir, full.names = TRUE)
} else {
  f.names <- filename
}

# check FCS
if (any(!isFCSfile(f.names)))
  stop("Not all imported files are FCS files.")

assign("actual", 0, envir = .GlobalEnv)
task = ctx$task


#2. convert them to FCS files
f.names %>%
  lapply(function(filename) {
    data = fcs_sort_to_data(filename)
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