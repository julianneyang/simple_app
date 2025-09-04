renv::init(bioconductor="3.20")
renv::install("rsconnect")
renv::snapshot()
rsconnect::writeManifest(
  appFiles= files,
  #appFileManifest = "src/app.R",
  appMode= "shiny",
  envManagementR = TRUE,
  verbose = TRUE
)


files <- list.files(
  path = ".", 
  recursive = TRUE,   # look into subfolders
  full.names = TRUE   # keep relative paths
)

