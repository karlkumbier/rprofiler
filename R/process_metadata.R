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
  workbook.info <- read_excel(meta.file, col_names=FALSE, skip=8, n_max=2)
  sheets <- unname(unlist(workbook.info[1, -1]))

  if (!'UsedWells' %in% sheets) stop('Specify UsedWells')
  markers <- str_split(unlist(workbook.info[2, 2]), ', *')[[1]]

  # Read sheets of plate metadata
  workbook <- lapply(sheets, function(s) {
                       tryCatch({
                         read_excel(s, path=meta.file, trim_ws=TRUE, range='B1:Y17')
                       }, error=function(e) {
                         print(str_c('Missing sheet: ', s))
                         return(NULL)
                       })
  })

  # Drop missing sheets
  id.drop <- sapply(workbook, is.null)
  workbook <- workbook[!id.drop]
  sheets <- sheets[!id.drop]

  # Drop empty sheets
  id.drop <- sapply(workbook, function(z) mean(is.na(z)) == 1)
  workbook <- workbook[!id.drop]
  sheets <- sheets[!id.drop]

  # Extract metadata from each sheet
  id.used <- which(sheets == 'UsedWells')
  n.row <- nrow(workbook[[id.used]]) #sum(rowMeans(is.na(workbook[[id.used]])) != 1)
  n.col <- ncol(workbook[[id.used]]) #sum(colMeans(is.na(workbook[[id.used]])) != 1)
  xmeta <- sapply(workbook, wellMeta)

  col.id <- rep(1:n.col, each=n.row)
  row.id <- rep(1:n.row, times=n.col)

  # Filter metadata to used wells and set well/plate.fileentifiers
  colnames(xmeta) <- sheets
  xmeta <- data.table(xmeta) %>%
    mutate(RowID=row.id, ColID=col.id, WellID=str_c(row.id, '-', col.id)) %>%
    mutate(UsedWells=str_remove_all(UsedWells, '(^\ |\ $)')) %>%
    filter(UsedWells == 1) %>%
    mutate_if(is.character, str_replace_all, pattern=',', replacement='\\.')
  
  if (!all(is.na(markers))) {
    xmeta <- mutate(xmeta, Markers=str_c(markers, collapse='_')) %>%
      mutate(Markers=str_remove_all(Markers, '\\((.+?)\\)')) %>%
      mutate(Markers=str_replace_all(Markers, '/', '-'))
  }

  return(xmeta)
}

wellMeta <- function(x) {
  # Extract well-specific metadata from excel sheet
  #id.drop <- is.na(x[,1])
  #x <- as.matrix(x)[1:nrow(x), 2:ncol(x)]
  #x <- x[!id.drop,]
  rownames(x) <-  NULL
  colnames(x) <-  NULL
  return(c(as.matrix(x)))
}


