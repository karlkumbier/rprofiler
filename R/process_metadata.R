#' Load metadata
#'
#' Load in metadata, stored as excel spreadsheet with standard formatting,  
#' for a select plate. TODO: Note on standard formatting.
#'
#' @param meta.file (character): path to file containing plate metadata.
#'
#' @return a data table with one entry per well, indicating metadata associated
#'   with each well.
#'
#' @export
#'
#' @importFrom stringr str_c str_split str_remove_all
#' @importFrom dplyr mutate filter
#' @importFrom readxl read_excel
#' @importFrom data.table data.table
loadMeta <- function(meta.file) {

  # Check whether file exists at indicated path
  if (!file.exists(meta.file)) {
    warning(str_c(meta.file, 'not found'))
    return(NULL)
  }

  # Read in info sheet to determine meta data pages
  # TODO: remove hard coding for excel file info
  workbook.info <- read_excel(meta.file, col_names=FALSE, skip=8, n_max=2)
  sheets <- unname(unlist(workbook.info[1, -1]))
  markers <- str_split(unlist(workbook.info[2, 2]), ', ')[[1]]

  # Read sheets of plate metadata
  workbook <- lapply(sheets, read_excel, path=meta.file)

  # Extract metadata from each sheet
  # TODO: add checks for empty rows
  n.row <- sum(!is.na(workbook[[1]][,1]))
  n.col <- sum(!is.na(workbook[[1]][1,])) - 1
  xmeta <- sapply(workbook, wellMeta)
  col.id <- rep(1:n.col, each=n.row)
  row.id <- rep(1:n.row, times=n.col)

  # Filter metadata to used wells and set well/plate.fileentifiers
  colnames(xmeta) <- sheets
  xmeta <- data.table(xmeta) %>%
    mutate(RowID=row.id, ColID=col.id, WellID=str_c(row.id, '-', col.id)) %>%
    mutate(UsedWells=str_remove_all(UsedWells, '(^\ |\ $)')) %>%
    filter(UsedWells == 1) %>%
    mutate(Markers=str_c(markers, collapse='_')) %>%
    mutate_if(is.character, str_replace_all, pattern=',', replacement='\\.')
  return(xmeta)
}

wellMeta <- function(x) {
  # Extract well-specific metadata from excel sheet
  id.drop <- is.na(x[,1])
  x <- as.matrix(x)[1:nrow(x), 2:ncol(x)]
  x <- x[!id.drop,]
  rownames(x) <-  NULL
  colnames(x) <-  NULL
  return(c(x))
}

