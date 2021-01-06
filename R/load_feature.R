#' Load feature descriptors
#'
#' Core function for generating a table of feature descriptors for a select
#' plate.
#'
#' @param plate (character) path to plate directory
#'
#' @return data.frame indicating feature name, channel, and category
#'
#' @export
#'
#' @importFrom stringr str_c str_subset
#' @importFrom R.matlab readMat
#' @importFrom dplyr mutate
loadFeature <- function(plate) {
    # Generate data frame of feature descriptors
    feature.dir <- str_c(plate, 'features/cbfeatures/')
    files <- list.files(feature.dir) %>% str_subset('.*mat$')

    # Read feature descriptors and format for table
    xfeat <- readMat(str_c(feature.dir, files[1]))$features
    features <- clean_unlist(xfeat[2,,])[-(1:2)]
    channel <- lapply(clean_unlist(xfeat[3,,], FALSE), unlist)
    channel <- sapply(channel, str_c, collapse='')[-(1:2)]
    cat <- clean_unlist(xfeat[5,,])[clean_unlist(xfeat[4,,])]

    out <- data.frame(Feature=features) %>%
        mutate(Channel=channel) %>%
        mutate(Category=cat)

    return(out)
}

clean_unlist <- function(x, recursive=TRUE) {
    return(unname(unlist(x, recursive)))
}

