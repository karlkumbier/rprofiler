library(rprofiler)
library(tidyverse)
options(stringsAsFactors=FALSE)
args <- R.utils::commandArgs(asValues=TRUE)

# Set paths for project, platesmetadata, code, and output
meta.file <- args$META_FILE
plate.dir <- args$PLATE_DIR
plate.id <- args$PLATE_ID
controls.cl <- str_split(args$CONTROL_CL, '\\|')[[1]]
controls.cpd <- str_split(args$CONTROL_CPD, '\\|')[[1]]
controls.usg <- str_split(args$CONTROL_USG, '\\|')[[1]]
meta.output <- str_c(plate.dir, '/', plate.id, '/metadata.csv')
n.core <- as.numeric(args$N_CORE)
type < <-  args$TYPE

# Load in plates for selected experiments
if (!file.exists(plate.id.file)) stop('plate id file not found')
plate.ids <- as.character(unlist(fread(plate.id.file)))

# Load in metadata for each plate
xmeta <- loadMeta(meta.file) %>% 
  mutate(PlateID=plate.id, ID=str_c(plate.id, '_', WellID))

# Clean marker set names for output csv
xmeta$Markers <- str_remove_all(xmeta$Markers, '\\((.+?)\\)')
xmeta$Markers <- str_replace_all(xmeta$Markers, '/', '-')
write.csv(file=meta.output, xmeta, quote=FALSE, row.names=FALSE)

# Generate KS profiles for plate
generateProfile(plate.id, 
                xmeta=xmeta, 
                plate.dir=plate.dir, 
                control.cl=control.cl,
                control.cpd=control.cpd,
                control.usg=control.usg,
                n.core=n.core,
                type=type)
