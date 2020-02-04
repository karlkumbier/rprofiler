#' Generate KS profile
#'
#' Core function for generating KS profile for all wells on a plate.
#'
#' @param plate (character) identifier for plate to generate KS profiles for.
#'   Must match a plateID in xmeta.
#' @param xmeta (data.table) well-level metadata as returned by loadMeta.
#' @param plate.dir (character) path to plate directory, containing the plate
#'   directory with image features for the selected plate.
#' @param control.cl (character) cell lines to use as KS reference population.
#' @param control.cpd (character) compound to use as KS reference
#'   population.
#' @param control.usg (character) cell line usage to set as KS reference
#'   population.
#' @param n.core (int) number of cores to parallelize over.
#' @param type (character) type of features to use for KS profile generation
#'   (one of: cbfeature, operetta)
#'
#' @return None. Writes KS profiles to plate directory.
#'
#' @export
#'
#' @importFrom data.table data.table rbindlist
#' @importFrom stringr str_c str_remove_all
#' @importFrom dplyr filter mutate
#' @importFrom parallel mclapply
#' @importFrom R.matlab readMat
generateProfile <- function(plate, xmeta, plate.dir,
                            control.cl=NULL,
                            control.cpd=NULL,
                            control.usg=NULL,
                            n.core=1, 
                            type='cbfeature') {
  
  # Check for valid input
  if (is.null(c(control.cl, control.cpd, control.usg)))
    stop('Specifiy a reference population -- cell line(s), compound(s), or usage')
  if (!type %in% c('cbfeature', 'operetta'))
    stop('type must be one of cbfeature or operetta')

  # Filter xmeta to current plate
  xmeta <- filter(xmeta, PlateID == plate)

  # Load in data for controls to generate KS reference distribution 
  xcontrol <- loadControl(xmeta, plate.dir, type, 
                          control.cl, control.cpd, control.usg)
  
  # Set path to write each cell/compound combination
  write.path <- str_c(plate.dir, '/', plate, '/ks_profiles/')
  dir.create(write.path, recursive=TRUE, showWarnings=FALSE)
  well.ids <- unique(xmeta$WellID)

  # Load in each well of plate
  out <- mclapply(well.ids, function(w) {
    cat(str_c('PROCCESSING WELL ', w, '...\n'))
   
    tryCatch({    
      xtreat <- loadTreatment(xmeta, plate.dir, w, type)
   
      # Generate KS statistics, dropping na values 
      out.raw <- lapply(xtreat, wellKS, xcontrol=xcontrol)
      out.raw <- cleanListNA(out.raw)

      # Format data for return with well and plate ids
      out <- rbindlist(out.raw)
      colnames(out) <- c(colnames(xcontrol), 'Ncells')
      id <- names(xtreat)
      out <- data.frame(ID=id, PlateID=plate, out)
      
      # Write to file for cell line/ compound combination
      cat('DONE!!!\n')
      return(out)
    
    }, error=function(e) {
      warning(str_c('Error processing well: ', w, '. Skipping well')) 
      return(NULL)
    })
  }, mc.cores=n.core)
  
  save(file='~/test.Rdata', out)
  out <- rbindlist(cleanListNull(out))
  write.csv(file=str_c(write.path, 'profiles.csv'), out, 
            row.names=FALSE, quote=FALSE)
}

#' Load control data
#' @export
loadControl <- function(xmeta, plate.dir,
                        type='cbfeature', 
                        control.cl=NULL,
                        control.cpd=NULL,
                        control.usg=NULL) {
  # Load in image features for specified control wells 

  # Filter metadata to control cell lines
  if (!is.null(control.cl)) {
    xmeta <- filter(xmeta, CellLine %in% control.cl)
    if (!is.null(c(control.cpd, control.usg)))
      warning('Multiple control types specified, using cell line')

  } else if (!is.null(control.cpd)) {
    xmeta <- filter(xmeta, Compound %in% control.cpd)
    if (!is.null(c(control.cpd, control.usg)))
      warning('Multiple control types specified, using cell line')
  
  } else if (!is.null(control.usg)) {
    xmeta <- filter(xmeta, Usage %in% control.usg) 

  } else {
    stop('Specifiy control cell line, compound, or usage')
  }

  # Load in selected wells and format for return
  wells <- loadWells(xmeta, plate.dir, type)
  wells <- unlist(wells, recursive=FALSE)
  wells <- do.call(rbind, wells)
  return(wells)
}

#' Load treatment data
#' @export
loadTreatment <- function(xmeta, plate.dir, well.id, type='cbfeature') {
  # Load in raw image features for specified well on given plate

  # Check that metadata has been filtered to specific plate
  plate <- unique(xmeta$plate)
  if (length(plate) > 1) stop()

  # Filter metadata to selcted cell line
  xmeta <- filter(xmeta, WellID %in% well.id)
    
  # Load in selected wellsand format for return
  wells <- loadWells(xmeta, plate.dir, type)
  wells <- unlist(wells, recursive=FALSE) 
  return(wells)
}


loadWells <- function(xmeta, plate.dir, type='cbfeature') {
  # Wrapper function to load data from multiple wells based on filtered metadata

  # Check that metadata has been filtered to specific plate
  plate <- as.character(unique(xmeta$PlateID))
  if (length(plate) > 1) stop('Multiple plates')

  # Get plate directory based on filtered metadata
  plates <- list.dirs(plate.dir, recursive=FALSE)
  plates <- plates[str_detect(plates, plate)]

  if (type == 'cbfeature') {
    feature.dir <- str_c(plates, '/features/cbfeatures/')
    ext <- '.mat'
    well.id <- str_c(xmeta$RowID, xmeta$ColID, sep='-')
    well.files <- str_c('cbfeatures-', well.id, ext)
    wells <- lapply(well.files, loadWellMat, feature.dir=feature.dir)
    wels <- cleanListNull(wells)
  } else {
    #TODO: check operetta formatting
    ext <-  '.csv'
  }

  return(wells)
} 


loadWellMat <- function(feature.dir, well.file) {
  # Function to load data from a selected well 
  
  tryCatch({  
    # Load in feature matrix for an individual well
    cat(str_c('Loading well ', well.file, '...\n'))
    well.id <- str_remove_all(well.file, '(cbfeatures-|\\.mat)')
    xmat <- readMat(str_c(feature.dir, '/', well.file))
    x <- xmat$features['data',,][[1]]

    # Get chanel and category metadata for each feature
    x.featnames <- unlist(xmat$features['name',,])
    id.feat <- 3:length(x.featnames)

    # Add channel id for unique feature name
    chan.id <- lapply(xmat$features['channel.id',,][[1]], unlist)
    chan.id <- sapply(chan.id, str_c, collapse='')
    x.featnames[id.feat] <- str_c('X', x.featnames[id.feat], chan.id[id.feat], sep='..')
    colnames(x) <- x.featnames
    
    # Subset to feature data
    x <- list(x[,id.feat])
    names(x) <- well.id
    return(x)
  }, error=function(e) {
    warning(str_c('Error loading well: ', well.file, '. Well has beend dropped'))
    return(NULL)
  })

}


wellKS <- function(xwell, xcontrol, id.feat=NULL) {
  # Generate KS statistics for indicated features.

  if (is.null(id.feat)) id.feat <- 1:ncol(xcontrol)
  suppressWarnings(
    out <- sapply(id.feat, function(i) signedKS(xwell[,i], xcontrol[,i]))
  )

  n <- nrow(xwell)
  out <- data.frame(matrix(out, nrow=1)) %>% mutate(NCells=n)
  return(out)
}


signedKS <- function(x, y, min.n=5) {
  # Cacluates a signed KS statistic between two empirical distirbutions x and y.
  # Drop NA values in x, y
  x <- c(na.omit(x))
  y <- c(na.omit(y))

  if (length(x) < min.n | length(y) < min.n) {
    warning('Fewer than min.n cells in well, skipping...')  
    return(NA)
  }

  ks1 <- ks.test(x, y, alternative='greater')$statistic
  ks2 <- ks.test(x, y, alternative='less')$statistic
  ks.sign <- ifelse(ks1 > ks2, 1, -1)
  ks <- unname(max(ks1, ks2) * ks.sign)
  return(ks)
}

################################################################################
# General utility functions for loading data from indicated files
################################################################################
cleanListNA <- function(x) {
  # Drop NA values from list
  id.drop <- sapply(x, is.na)
  x <- x[!id.drop]
  return(x)
}

cleanListNA <- function(x) {
  # Drop NA values from list
  id.drop <- sapply(x, is.na)
  x <- x[!id.drop]
  return(x)
}

cleanListNull <- function(x) {
  # Drop null values from list
  id.drop <- sapply(x, is.null)
  x <- x[!id.drop]
  return(x)
}
