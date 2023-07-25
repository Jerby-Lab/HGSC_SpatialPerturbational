# files<-list.files(paste0("/local/path/to/this/repository/Code"),
files <- list.files(paste0("~/Projects/HGSC_SpatialPerturbational/Code"),
                  include.dirs = F,
                  pattern = ".R",
                  full.names = T)
lapply(files, source)
start <- Sys.time()
HGSC_main()
print(Sys.time() - start)
