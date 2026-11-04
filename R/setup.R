# setup.R — source every module in order. Scripts start with: source("R/setup.R").
local({
  here <- tryCatch({
    a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(a)) dirname(sub("^--file=", "", a[1])) else "R"
  }, error = function(e) "R")
  # resolve to the R/ directory whether called from root or scripts/
  cand <- c("R", file.path("..", "R"), here, file.path(here, "..", "R"))
  rdir <- cand[which(file.exists(file.path(cand, "config.R")))[1]]
  if (is.na(rdir)) stop("cannot locate R/ modules")
  for (f in c("config.R","params.R","model.R","arms.R","metrics.R"))
    source(file.path(rdir, f))
})
