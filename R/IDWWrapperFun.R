#' Inverse distance-weighted interpolation wrapper function
#'
#' @param field a variable (i.e. column name in \code{data}) to be
#'   interpolated across \code{rasTemplate}
#' @param data a data.table with variables to be interpolated
#' @param rasTemplate a \code{RasterLayer} or \code{SpatRaster} on which interpolation occurs
#' @param ... passed to \code{gstat::gstat}
#'
#' @importFrom terra interpolate mask
#' @importFrom gstat gstat

IDWWrapperFun <- function(field, data, rasTemplate, ...) {
  form <- as.formula(paste(field, "~ 1"))
  interpModel <- gstat(formula = as.formula(form), data = data, ...)
  # interpModel <- gstat(formula = form, locations = weatherDataShort)   ## for IDW interpolation
  if (is(rasTemplate, "RasterLayer")) {
    weatherDataIntrpRas <- interpolate(object = rasTemplate, model = interpModel)  ## interpolate on RTM
  }
  if (is(rasTemplate, "SpatRaster")) {
    weatherDataIntrpRas <- interpolate(object = rasTemplate, model = interpModel,
                                       debug.level = 0, index = 1)  ## interpolate on RTM
  }
  weatherDataIntrpRas <- mask(weatherDataIntrpRas, rasTemplate)
  weatherDataIntrpRas
}

#' Convert raster to points, adding pix IDs
#'
#' @param ras a \code{RasterLayer} or \code{SpatRaster} to be
#'   converted to points
#' @param rasterToMatch a \code{RasterLayer} or \code{SpatRaster} that
#'   will supply pixelIDs (so needs to be identical to \code{ras}, except in values)
#' @param crs passed to \code{st_transform}, the final projection
#'
#' @importFrom terra as.points
#' @importFrom sf st_as_sf st_transform
.rasterToPoints <- function(ras, rasterToMatch, crs) {
  x <- st_as_sf(as.points(ras))
  x$pixelIndex <- which(!is.na(rasterToMatch[]))
  x <- st_transform(x, crs = crs)
  x
}
