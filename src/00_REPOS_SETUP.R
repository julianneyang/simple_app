renv::init(bioconductor="3.20")
renv::install("rsconnect")
renv::snapshot()
rsconnect::writeManifest(
  #appFileManifest = "src/app.R",
  appMode= "shiny",
  envManagementR = TRUE,
  verbose = TRUE
)
