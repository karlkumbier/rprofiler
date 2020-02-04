library(rprofiler)
library(tidyverse)
options(stringsAsFactors=FALSE)

# Set paths for project, platesmetadata, code, and output
meta.file <- '~/github/Rpackages/rprofiler/data/metadata/2020019002.xlsx'
plate.dir <- '~/github/Rpackages/rprofiler/data/plates/'
plate <- 'MAR_20200115_DensityDilutions_plate_2020019002'
control.cl <- c('ALS15', 'ALS24', 'ALS45')
control.cpd <- NULL
control.usg <- NULL
meta.output <- str_c(plate.dir, '/', plate, '/metadata.csv')
n.core <- as.numeric(2)
type <- 'cbfeature'

# Load in metadata for each plate
xmeta <- loadMeta(meta.file) %>% 
  mutate(PlateID=plate, ID=str_c(plate, '_', WellID))

# Clean marker set names for output csv
xmeta$Markers <- str_remove_all(xmeta$Markers, '\\((.+?)\\)')
xmeta$Markers <- str_replace_all(xmeta$Markers, '/', '-')
write.csv(file=meta.output, xmeta, quote=FALSE, row.names=FALSE)

# Generate KS profiles for plate
generateProfile(plate=plate, 
                xmeta=xmeta, 
                plate.dir=plate.dir, 
                control.cl=control.cl,
                control.cpd=control.cpd,
                control.usg=control.usg,
                n.core=n.core,
                type=type)
