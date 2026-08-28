library(targets)

t0 <- Sys.time()
tar_make()
# tar_make(names = !tidyselect::matches("_loo|loo_|_omega"))
t1 <- Sys.time()

runtime <- difftime(t1, t0)

writeLines(
  paste("Total runtime:", round(runtime, digits=3), units(runtime)),
  "outputs/runtime.txt"
)

store_report <- function(path = "_targets/objects", top = 20) {
  f <- list.files(path, full.names = TRUE)
  i <- file.info(f)
  d <- data.frame(name = basename(f), mb = i$size / 1024^2)
  d <- d[order(-d$mb), ]
  cat(sprintf("total: %.2f GB across %d objects\n", sum(d$mb) / 1024, nrow(d)))
  cat(sprintf("scratch: %.2f GB\n",
              sum(file.info(list.files("_targets/scratch", full.names = TRUE,
                                       recursive = TRUE))$size, na.rm = TRUE) / 1024^3))
  # group by target family
  d$family <- sub("_(CVA6|EV71|EV68).*$", "", d$name)
  fam <- aggregate(mb ~ family, d, sum)
  print(head(fam[order(-fam$mb), ], 15))
  head(d, top)
}
store_report()
