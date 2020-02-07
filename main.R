library(rprofiler)
library(stringr)
library(dplyr)
options(stringsAsFactors=FALSE)
args <- R.utils::commandArgs(asValues=TRUE)

# Set paths for project, platesmetadata, code, and output
meta.files <- str_split(args$META_FILE, ',')[[1]]
plate.dir <- args$PLATE_DIR
plate.ids <- str_split(args$PLATE_ID, ',')[[1]]
write.path <- args$WRITE_PATH

if (is.null(args$CONTROLS))
  stop('Specify KS controls')
if (is.null(args$CONTROL_VARIABLE))
  stop('Specify a variable to use for KS controls')

controls <- str_split(args$CONTROLS, ',')[[1]]
control.variable <- str_split(args$CONTROL_VARIABLE, ',')[[1]]

n.core <- as.numeric(args$N_CORE)
type <-  args$TYPE
if (is.null(n.core)) n.core <- 1

# Check for valid unput arguments
if (is.null(plate.ids))
  stop('Specify plate')

if (is.null(meta.files))
  stop('Specify metadata file')

if (length(plate.ids) != length(meta.files))
  stop('Specify metadata file for each plate')

if (is.null(plate.dir))
  stop('Specify plate directory')

if (is.null(type))
  stop('Specify type of features (i.e. cbfeature or operetta)')

###############################################################################
# Iterate over all plates, generating metadata files and profiles
###############################################################################
for (i in 1:length(plate.ids)) {
  
  # Load in metadata for each plate
  xmeta <- loadMeta(meta.files[i]) %>% 
    mutate(PlateID=plate.ids[i]) %>%
    mutate(ID=str_c(PlateID, '_', WellID))

  # Clean marker set names for output csv
  meta.output <- str_c(plate.dir, '/', plate.ids[i], '/metadata.csv')
  write.csv(file=meta.output, xmeta, quote=FALSE, row.names=FALSE)

  # Generate KS profiles for plate
  generateProfile(plate.ids[i], 
                  xmeta=xmeta, 
                  plate.dir=plate.dir, 
                  controls=controls,
                  control.variable=control.variable,
                  n.core=n.core,
                  type=type)
}

# Aggregate profiles if processing multiple plates
if (!is.null(write.path)) aggregate_profiles(plate.ids, plate.dir, write.path)
print('PROFILES DONE!!!')
