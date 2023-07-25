files<-list.files(paste0("/path/to/your/local/clone/of/this/repo/Code"),
                  include.dirs = F,
                  pattern = ".R",
                  full.names = T)
lapply(files, source)
start <- Sys.time()
HGSC_main()
print(Sys.time() - start)
