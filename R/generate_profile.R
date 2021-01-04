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
#' @importFrom dplyr filter mutate select
#' @importFrom parallel mclapply mcmapply
#' @importFrom R.matlab readMat
generateProfile <- function(plate, xmeta, plate.dir, control.variable, controls,
                            type='cbfeature',
                            n.bs=0,
                            n.core=1) {
  
  # Check for valid input
  if (!type %in% c('cbfeature', 'operetta'))
    stop('type must be one of cbfeature or operetta')

  # Filter xmeta to current plate
  xmeta <- filter(xmeta, PlateID == plate)

  print('Loading Reference')
  # Load in data for controls to generate KS reference distribution 
  xcontrol <- loadControl(xmeta, plate.dir, control.variable, 
                          controls, type, n.core)
  print(str_c('N REFERENCE: ', nrow(xcontrol)))

  # Set path to write each cell/compound combination
  write.path <- str_c(plate.dir, '/', plate, '/ks_profiles/')
  dir.create(write.path, recursive=TRUE, showWarnings=FALSE)
  well.ids <- unique(xmeta$WellID)

  # Load in each well of plate
  out <- mclapply(well.ids, function(w) {
   
    tryCatch({    
      # Load in raw, cell-level feature data 
      xtreat <- loadTreatment(xmeta, plate.dir, w, type)
     
      # Generate KS statistics, dropping na values 
      out.raw <- lapply(xtreat, wellKS, xcontrol=xcontrol)
      
      # Format data for return with well and plate ids
      out <- rbindlist(out.raw)
      
      if (n.bs > 0) {
          xtreat.bs <- loadTreatment(xmeta, plate.dir, w, type, n.bs)
          out.bs <- lapply(xtreat.bs, wellKS, xcontrol=xcontrol)
          out.bs <- rbindlist(out.bs)
          
          out <- mutate(out, Bootstrap=FALSE)
          out.bs <- mutate(out.bs, Bootstrap=TRUE)
          out <- rbind(out, out.bs)
          
          id <- c(names(xtreat), names(xtreat.bs))
          colnames(out) <- c(colnames(xcontrol), 'Ncells', 'Bootstrap')
      } else {
          id <- names(xtreat)
          colnames(out) <- c(colnames(xcontrol), 'Ncells')
      }

      out <- data.frame(ID=id, PlateID=plate, out)
      
      # Write to file for cell line/ compound combination
      return(out)
    
    }, error=function(e) {
      warning(str_c('Error processing well: ', w, '. Skipping well'))
      print(w)
      print(e) 
      return(NULL)
    })
  }, mc.cores=n.core)
  
  out <- rbindlist(cleanListNull(out))
  write.csv(file=str_c(write.path, type, '_profiles.csv'), out, 
            row.names=FALSE, quote=FALSE)
}


loadControl <- function(xmeta, plate.dir, control.variable, controls,
                        type='cbfeature',
                        n.core=1) { 
  # Load in image features for specified control wells 

  # Filter metadata to control cell lines
  if (!all(unlist(control.variable) %in% colnames(xmeta))) {
    stop('Control variable not in metadata')
  }

  for (i in 1:length(control.variable)) {
    id.keep <- xmeta[[control.variable[[i]]]] %in% controls[[i]]
    xmeta <- xmeta[id.keep,]
  }

  # Load in selected wells and format for return
  print(str_c('Number metadata controls: ', nrow(xmeta)))
  wells <- loadWells(xmeta, plate.dir, type, n.core=n.core)
  wells <- unlist(wells, recursive=FALSE)
  wells <- rbindlist(wells)
  return(wells)
}

loadTreatment <- function(xmeta, plate.dir, well.id, 
                          type='cbfeature', n.bs=0) {
  # Load in raw image features for specified well on given plate

  # Check that metadata has been filtered to specific plate
  plate <- unique(xmeta$plate)
  if (length(plate) > 1) stop()

  # Filter metadata to selcted cell line
  xmeta <- filter(xmeta, WellID %in% well.id)
    
  # Load in selected wellsand format for return
  wells <- loadWells(xmeta, plate.dir, type, n.bs=n.bs)
  wells <- unlist(wells, recursive=FALSE) 
  return(wells)
}


loadWells <- function(xmeta, plate.dir, type='cbfeature', n.bs=0, n.core=1) {
  # Wrapper function to load data from multiple wells based on filtered metadata

  # Check that metadata has been filtered to specific plate
  plate <- as.character(unique(xmeta$PlateID))
  if (length(plate) > 1) stop('Multiple plates')

  # Get plate directory based on filtered metadata
  plate.path <- list.dirs(plate.dir, recursive=FALSE)
  plate.path <- plate.path[str_detect(plate.path, plate)]

  if (type == 'cbfeature') {
    feature.dir <- str_c(plate.path, '/features/cbfeatures/')
    ext <- '.mat'
    well.id <- str_c(xmeta$RowID, xmeta$ColID, sep='-')
    well.files <- str_c('cbfeatures-', well.id, ext)
    
    wells <- mclapply(well.files, loadWellMat, 
                      feature.dir=feature.dir,
                      n.bs=n.bs, 
                      mc.cores=n.core)
    wells <- cleanListNull(wells)
  } else {
    feature.dir <- str_c(plates, '/features/opfeatures/')
    
    input.file <- list.files(feature.dir)
    input.file <- input.file[str_detect(input.file, 'Objects_Population')]
    if (length(input.file) > 1) stop('Multiple feature files detected')
    if (length(input.file) == 0) stop('No feature files detected')
    
    feature.id <- str_remove_all(input.file, 'Objects_Population - ') %>%
      str_remove_all('\\..*')
    
    xop <- fread(str_c(feature.dir, input.file), skip='Row')
    wells <- mcmapply(function(row, col) {
                      loadWellOp(row, col, x=xop, feature.id=feature.id) 
                      }, xmeta$RowID, xmeta$ColID, 
                      mc.cores=n.core, SIMPLIFY=FALSE)
    wells <- cleanListNull(wells)
  }

  return(wells)
} 


loadWellMat <- function(feature.dir, well.file, n.bs=0) {
  # Function to load data from a selected well 
  
  tryCatch({  
    # Load in feature matrix for an individual well
    well.id <- str_remove_all(well.file, '(cbfeatures-|\\.mat)')
    xmat <- readMat(str_c(feature.dir, '/', well.file))
    x <- xmat$features['data',,][[1]]

    # Get chanel and category metadata for each feature
    x.featnames <- unlist(xmat$features['name',,])
    id.feat <- 3:length(x.featnames)

    # Add channel id for unique feature name
    chan.id <- lapply(xmat$features['channel.id',,][[1]], unlist)
    chan.id <- sapply(chan.id, str_c, collapse='')
    x.featnames[id.feat] <- str_c('X', x.featnames[id.feat], 
                                  chan.id[id.feat], sep='..')
    
    # Subset to feature data
    x <- matrix(x[,id.feat], ncol=length(id.feat)) 
    colnames(x) <- x.featnames[id.feat]
    x <- data.table(x)

    # Bootstrap sample cells from each well
    if (n.bs > 0) {
        bsid <- bsSample(nrow(x), n.bs)
        x <- lapply(bsid, function(i) x[i,])
        names(x) <- rep(well.id, n.bs)
    } else {
        x <- list(x)
        names(x) <- well.id
    }
    
    return(x)
  }, error=function(e) {
    warning(str_c('Error loading well: ', well.file, 
                  '. Well has beend dropped'))
    return(NULL)
  })

}

loadWellOp <- function(x, row, col, feature.id) {
  tryCatch({
    xop <- filter(x, Row == row) %>%
      filter(Column == col) %>%
      select(contains(feature.id))

    x <- list(xop)
    names(x) <- str_c(row, '-', col)
    return(x)
  }, error=function(e) {
    warning(str_c('Error loading well: ', row, '-', col, 
                  '. Well has beend dropped'))
    return(NULL)
  })
}

#' Evaluate KS given well level and control matrices
#'
#' @export
#'
wellKS <- function(xwell, xcontrol, id.feat=NULL) {
  # Generate KS statistics for indicated features.

  xwell <- as.data.frame(xwell)
  xcontrol <- as.data.frame(xcontrol)
  if (is.null(id.feat)) id.feat <- 1:ncol(xcontrol)
  out <- sapply(id.feat, function(i) signedKS(xwell[,i], xcontrol[,i]))

  n <- nrow(xwell)
  out <- data.frame(matrix(out, nrow=1)) %>% mutate(NCells=n)
  return(out)
}


signedKS <- function(x, y, min.n=5) {
  # Cacluates a signed KS statistic between two empirical distirbutions x and y.
  # Drop NA values in x, y
  x <- c(na.omit(x))
  y <- c(na.omit(y))

  ks <- ksTest(x, y)
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

bsSample <- function(n, n.bs) {
  # Take k bootstrap samples of length n
  genSample <- function() sample(1:n, n, replace=TRUE)
  id <- replicate(n.bs, genSample(), simplify=FALSE) 
  return(id)
}

