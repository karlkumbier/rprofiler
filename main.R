library(rprofiler)
library(stringr)
library(dplyr)
library(data.table)
options(stringsAsFactors=FALSE)
args <- R.utils::commandArgs(asValues=TRUE)
save(file='args.Rdata', args)

# Set paths for project, platesmetadata, code, and output
meta.dir <- args$META_DIR
plate.dir <- args$PLATE_DIR

plate.ids <- args$PLATE_ID
if (str_detect(plate.ids, 'csv')) {
  plate.ids <- as.character(fread(plate.ids, header=FALSE)$V1)
} else {
  plate.ids <- str_split(plate.ids, ',')[[1]]
}

write.dir <- args$WRITE_DIR

# Check that plate IDs exist in plate dir and are unique
plates <- list.files(plate.dir)
if (length(plates) == 0) {
  stop('No plates found in plate directory')
}

counts <- colSums(as.matrix(sapply(plate.ids, str_detect, string=plates)))
if (any(counts > 1) | any(counts == 0)) {
  warning('Missing or non-unique plate ids. Subsetting to unique plates')
  plate.ids <- plate.ids[counts == 1]
}

# Format controls for profile generation
if (is.null(args$CONTROLS))
  stop('Specify KS controls')
if (is.null(args$CONTROL_VARIABLE))
  stop('Specify a variable to use for KS controls')

controls <- str_split(args$CONTROLS, '\\|')[[1]]
controls <- str_split(controls, ',')

control.variable <- str_split(args$CONTROL_VARIABLE, '\\|')[[1]]
control.variable <- str_split(control.variable, ',')

n.bs <- as.numeric(ifelse(is.null(args$NBS), 0, args$NBS))
n.core <- as.numeric(args$N_CORE)
print(n.core)
type <-  args$TYPE
if (is.null(n.core)) n.core <- 1

# Check for valid unput arguments
if (is.null(plate.ids))
  stop('Specify plate')

if (is.null(meta.dir))
  stop('Specify metadata directory')

if (is.null(plate.dir))
  stop('Specify plate directory')

if (is.null(type))
  stop('Specify type of features (i.e. cbfeature or operetta)')


###############################################################################
# Iterate over all plates, generating metadata files and profiles
###############################################################################
meta.files <- str_c(meta.dir, '/', plate.ids, '.xlsx')
for (i in 1:length(plate.ids)) {
  tryCatch({
  plate <- plates[str_detect(plates, plate.ids[i])]

  print(str_c('Processing plate: ', plate.ids[i]))
  # Load in metadata for each plate
  xmeta <- loadMeta(meta.files[i]) %>% 
    mutate(PlateID=plate) %>%
    mutate(ID=str_c(PlateID, '_', WellID))

  # Clean marker set names for output csv
  plate <- plates[str_detect(plates, plate.ids[i])]
  meta.output <- str_c(plate.dir, '/', plate, '/metadata.csv')
  write.csv(file=meta.output, xmeta, quote=FALSE, row.names=FALSE)
 
  # Generate KS profiles for plate
  generateProfile(plate, 
                  xmeta=xmeta, 
                  plate.dir=plate.dir, 
                  controls=controls,
                  control.variable=control.variable,
                  n.bs=n.bs,
                  n.core=n.core,
                  type=type)
  }, error=function(e) {
    print(str_c('Error processing plate ', plate.ids[i]))
    print(e)
  })
}

# Aggregate profiles if processing multiple plates
if (!is.null(write.dir)) {
  aggregate_profiles(plate.ids, plate.dir, write.dir, type)
}

print('PROFILES DONE!!!')
