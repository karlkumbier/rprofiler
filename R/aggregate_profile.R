#' Aggregate Profiles
#' 
#' Function to load in and aggregate profiles generated for multiple plates.
#'
#' @param plates (character) identifiers of the plates whose ks profiles will be
#'   loaded and aggregated.
#' @param plate.dir (character) path to plate directory, containing the plate
#'   directory with image features for the selected plate.
#' @param write.path (character) path to write aggregated ks profiles to.
#'
#' @return None. Will save a csv file containing KS profiles for all selected
#'   plates.
#'
#' @export
#'
#' @importFrom stringr str_c str_detect
#' @importFrom data.table fread rbindlist
aggregate_profiles <- function(plate.ids, plate.dir, write.dir, type) {
  
  # Check for valid profile files in each plate directory
  plates <- list.files(plate.dir)
  plates.select <- lapply(plate.ids, function(z) plates[str_detect(plates, z)])
  id.drop <- sapply(plates.select, length) != 1
  plates.select <- plates.select[!id.drop]

  profile.paths <- str_c(plate.dir, '/', plates.select, '/', 
                         'ks_profiles/', type, '_profiles.csv')
  meta.paths <- str_c(plate.dir, '/', plates.select, '/metadata.csv')
  
  id.exist.p <- sapply(profile.paths, file.exists)
  id.exist.m <- sapply(meta.paths, file.exists)
  id.exist <- id.exist.p & id.exist.m
  if (sum(!id.exist) > 0) {
    missing <- profile.paths[!id.exist]
    for (m in missing) {
      warning(str_c('Missing KS profiles for: ', m, '. Plate dropped'))
    }
  }

  if (sum(id.exist) == 0) {
    warning('All KS profiles missing')
  } else {
    
    # Load in profiles for all available plates
    profile.paths <- profile.paths[id.exist]
    profiles <- lapply(profile.paths, fread)
    profiles <- rbindlist(profiles, fill=TRUE)

    meta <- lapply(meta.paths[id.exist], fread)
    meta <- rbindlist(meta, fill=TRUE)
    
    meta.path <- str_c(write.dir, '/metadata.csv')
    profile.path <- str_c(write.dir, '/', type, '_profiles.csv')
    dir.create(write.dir, recursive=TRUE, showWarnings=FALSE)
    write.csv(file=profile.path, profiles, row.names=FALSE, quote=FALSE)
    write.csv(file=meta.path, meta, row.names=FALSE, quote=FALSE)
  }
}

#' Aggregate feature descriptors
#' 
#' Function to load in and aggregate profiles generated for multiple plates.
#'
#' @param plates (character) identifiers of the plates whose ks profiles will be
#'   loaded and aggregated.
#' @param plate.dir (character) path to plate directory, containing the plate
#'   directory with image features for the selected plate.
#' @param write.path (character) path to write aggregated ks profiles to.
#'
#' @return None. Will save a csv file containing KS profiles for all selected
#'   plates.
#'
#' @export
#'
#' @importFrom stringr str_c str_detect
#' @importFrom data.table fread rbindlist
aggregate_features <- function(plate.ids, plate.dir, write.dir) {

    # Check for valid feature files in each plate directory
    plates <- list.files(plate.dir)
    plates.select <- lapply(plate.ids, function(z) plates[str_detect(plates, z)])
    id.drop <- sapply(plates.select, length) != 1
    plates.select <- plates.select[!id.drop]

    feature.paths <- str_c(plate.dir, '/', plates.select, '/feature_descriptors.csv')
    id.exist <- sapply(feature.paths, file.exists)

    if (all(id.exist)) {
        out <- rbindlist(lapply(feature.paths, fread))
        feature.path <- str_c(write.dir, '/feature_descriptors.csv')
        write.csv(file=feature.path, out, row.names=FALSE, quote=FALSE)

    } else {
        stop(str_c('Missing feature files for plate: ', feature.paths[!id.exist]))
    }

    return(out)
}
