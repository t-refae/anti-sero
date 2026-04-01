library(targets)

t0 <- Sys.time()
tar_make()
t1 <- Sys.time()

runtime <- difftime(t1, t0)

writeLines(
  paste("Total runtime:", round(runtime, digits=3), units(runtime)),
  "outputs/runtime.txt"
)