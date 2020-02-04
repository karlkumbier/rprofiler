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
#' @importFrom stringr str_c
#' @importFrom data.table fread
aggregate_profiles <- function(plates, plate.dir, write.path) {
  
  # Check for valid profile files in each plate directory
  profile.paths <- str_c(plate.dir, '/', plates, '/', 
                         'ks_profiles/profiles.csv')
  
  id.exist <- sapply(profile.paths, file.exists)
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
    
    write.dir <- str_remove_all(write.path, '([^\\/]+$)')
    dir.create(write.dir, recursive=TRUE, showWarnings=FALSE)
    write.csv(file=write.path, profiles, row.names=FALSE, quote=FALSE)
  }
}
