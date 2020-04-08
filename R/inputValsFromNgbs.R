inputValsFromNgbs <- function(ras, cells, aggregateFun = "mean",
                              saveRas = TRUE, overwrite = TRUE) {
  rasFileName <- filename(ras)

  aggregateFun <- get(aggregateFun)
  ngbs <- as.data.table(adjacent(ras, cells, directions = 8, sorted = TRUE))

  ngbsVals <- data.table(to = ngbs$to, toVals = as.numeric(ras[ngbs$to]))
  ngbs <- ngbsVals[ngbs, on = .(to)]
  ngbs[, fromAvgVal := aggregateFun(toVals, na.rm = TRUE), by = from]

  inputVals <- unique(ngbs[, .(from, fromAvgVal)])
  inputVals[is.na(fromAvgVal), fromAvgVal := 0]

  ras[inputVals$from] <- inputVals$fromAvgVal

  ## overwrite raster and reload (otherwise we lose the filename for later)
  if (saveRas) {
    writeRaster(ras, filename = rasFileName, overwrite = overwrite)
    ras <- raster(rasFileName)
  }
  ras
}
