# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects and functions, use sim$xxx.
defineModule(sim, list(
  name = "Biomass_fireProperties",
  description = paste("Complement to fire spread that calculates fire (behaviour) properties in function of vegetation (fuels),",
                      "climate and topography conditions, using the Canadian Forest Fire Behaviour Prediction System"),
  keywords = c("fire behaviour", "fuels", "fire-vegetation feedbacks", "fire-climate feedbacks", "FBP system", "topography"),
  authors = person("Ceres", "Barros", email = "cbarros@mail.ubc.ca", role = c("aut", "cre")),
  childModules = character(0),
  version = list(Biomass_fireProperties = numeric_version("0.2.0"),
                 Biomass_core = numeric_version("1.3.2"),
                 LandR = "0.0.3.9000", SpaDES.core = "0.2.7",
                 raster = "3.1-5"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "Biomass_fireProperties.Rmd"),
  reqdPkgs = list("R.utils", "raster", "data.table", "dplyr", "gdalUtilities",
                  "sp", "sf", "cffdrs", "amc", "fasterize", "gstat", "crayon",
                  "PredictiveEcology/LandR@development",
                  "PredictiveEcology/SpaDES.core@development",
                  "PredictiveEcology/SpaDES.tools@development",
                  "PredictiveEcology/reproducible@development"),
  parameters = rbind(
    defineParameter("fireWeatherMethod", "character", "sample", NA, NA,
                    desc = paste("How fire weather is summarized by location. When 'weatherData' contains",
                                 "several fire weather observations (i.e. fire-days' weather values) per point,",
                                 "they will be summarized by either calculating average fire-day weather per",
                                 "location ('average'), or by sampling a fire-day ('sample'). If sampling,",
                                 "at each fire year a new sample will be drawn. Defaults to 'average' for speed.")),
    defineParameter("vegFeedback", "logical", TRUE, NA, NA, desc = "Should vegetation feedbacks unto fire be simulated? Defaults to TRUE"),
    defineParameter(name = "fireInitialTime", class = "numeric", default = 2L,
                    desc = "The event time that the first fire disturbance event occurs"),
    defineParameter(name = "fireTimestep", class = "numeric", default = 2L,
                    desc = "The number of time units between successive fire events in a fire module"),
    defineParameter(".plotMaps", "logical", FALSE, NA, NA, "This describes whether maps should be plotted or not"),
    defineParameter(".plotInterval", "numeric", 1, NA, NA, "This describes the simulation time interval between plot events"),
    defineParameter(".studyAreaName", "character", NA, NA, NA,
                    "Human-readable name for the study area used. If NA, a hash of studyArea will be used."),
    defineParameter(".useCache", "logical", "init", NA, NA,
                    desc = "use caching for the spinup simulation?")
  ),
  inputObjects = bind_rows(
    expectsInput(objectName = "DEMRas", objectClass = "RasterLayer",
                 desc = paste("Digital elevation model (DEM) raster used to make 'elevationRas', 'slopeRas' and 'aspectRas'.",
                              "Ideally the same used by the weather generator to make 'weatherData'."),
                 sourceURL = "https://drive.google.com/file/d/1Fosf9xfD4UmljwZCxH7MHqsO9EtK4nvp/view?usp=sharing"),
    expectsInput(objectName = "fuelTypesMaps", objectClass = "list",
                 desc = paste("List of RasterLayers of fuel types ('finalFuelType'), coniferDominance ('coniferDom'), and optionally",
                              "the degree of curing ('curing') per pixel. 'finalFuelType' needs to be a ratified raster (raster with levels)")),
    expectsInput(objectName = "FWIinit", objectClass = "data.frame",
                 desc = paste("Initalisation parameter values for FWI calculations. Defaults to default values in cffdrs::fwi.
                 This table should be updated every year)")),
    expectsInput(objectName = "pixelGroupMap", objectClass = "RasterLayer",
                 desc = "updated community map at each succession time step"),
    expectsInput(objectName = "pixelNonForestFuels", objectClass = "data.table",
                 desc = paste("Table of non forest fuel attributes (pixel ID, land cover, fuel type",
                              "name and code, and degree of curing) for each pixel with non-forest fuels.
                               Defaults to NULL, if not provided by another module")),
    expectsInput("rasterToMatch", "RasterLayer",
                 desc = "a raster of the studyArea in the same resolution and projection as rawBiomassMap",
                 sourceURL = NA),
    expectsInput(objectName = "rasterToMatchFBPPoints", objectClass = "sf",
                 desc = paste("the spatial points version of rasterToMatch reprojected",
                              "to FBP-compatible lat/long projection:",
                              "'+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'"),
                 source = NA),
    expectsInput("studyArea", "SpatialPolygonsDataFrame",
                 desc = paste("Polygon to use as the study area.",
                              "Defaults to  an area in Southwestern Alberta, Canada."),
                 sourceURL = NA),
    expectsInput(objectName = "studyAreaFBP", objectClass = "SpatialPolygonsDataFrame",
                 desc = paste("same as studyArea, but on FBP-compatible lat/long projection:",
                              "'+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'"),
                 sourceURL = NA),
    expectsInput(objectName = "topoData", objectClass = "data.table",
                 desc = paste("Table with slope and aspect values extracted from 'DEMRas', by pixelIndex",
                              "and on FBP compatible projection")),
    expectsInput(objectName = "weatherData", objectClass = "sf",
                 desc = paste("Weather point data, to be used to identify fire days.",
                              "Needs to have the following columns: 'month' (optional), 'day' (optional), 'temperature',",
                              "'relativeHumidity', 'windSpeed' (optional), 'precipitation'. 'temperature' refers to average air temperature",
                              "in Celsius, 'windSpeed' should be average wind speed (m/s) at 10m height, and 'precipitation' is",
                              "total precipitation in mm. These data can be daily, monthly or yearly averages.",
                              "Defaults to rasters of annual climate data downloaded from ClimateNA for 2011 using CanESM2_RCP45_r11i1p1_2011MSY"),
                 sourceURL = "https://drive.google.com/open?id=12iNnl3P7VjisVKC0vatSrXyhYtl6w-D1"),
    expectsInput(objectName = "weatherDataCRS", objectClass = "character",
                 desc = paste("The original projection of the climate data table. Must be supplied if",
                              "weatherData is supplied by the user or a module."))
  ),
  outputObjects = bind_rows(
    # createsOutput(objectName = "FBPinputs", objectClass = "RasterLayer",
    #               desc = "Fire behaviour prediction system inputs table"),
    # createsOutput(objectName = "FBPoutputs", objectClass = "list",
    #               desc = "Fire weather outputs table"),
    # createsOutput(objectName = "FWIinputs", objectClass = "RasterLayer",
    #               desc = "Fire weather inputs table"),
    # createsOutput(objectName = "FWIoutputs", objectClass = "list",
    #               desc = "Fire weather outputs table"),
    # createsOutput(objectName = "aspectRas", objectClass = "RasterLayer",
    #               desc = "Raster of aspect values extracted from 'DEMRas'"),
    # createsOutput(objectName = "slopeRas", objectClass = "RasterLayer",
    #               desc = "Raster of slope values extracted from 'DEMRas'")
    createsOutput(objectName = "fireCFBRas", objectClass = "RasterLayer",
                  desc = "Raster of crown fraction burnt"),
    createsOutput(objectName = "fireIntRas", objectClass = "RasterLayer",
                  desc = "Raster of equilibrium head fire intensity [kW/m]"),
    createsOutput(objectName = "fireROSRas", objectClass = "RasterLayer",
                  desc = "Raster of equilibrium rate of spread [m/min]"),
    createsOutput(objectName = "fireRSORas", objectClass = "RasterLayer",
                  desc = "Critical spread rate for crowning [m/min]"),
    createsOutput(objectName = "fireTFCRas", objectClass = "RasterLayer",
                  desc = "Raster of total fuel consumed [kg/m^2]"),
    createsOutput(objectName = "fireYear", objectClass = "numeric", desc = "Next fire year"),
    createsOutput(objectName = "weatherDataShort", objectClass = "sf",
                  desc = paste("Weather point data summarized according to 'fireWeatherMethod', by pixelIndex",
                               "and on FBP compatible projection")),
    createsOutput(objectName = "topoData", objectClass = "data.table",
                  desc = paste("Table with slope and aspect values extracted from 'DEMRas', by pixelIndex",
                               "and on FBP compatible projection"))
  )
))

doEvent.Biomass_fireProperties = function(sim, eventTime, eventType, debug = FALSE) {
  switch(
    eventType,
    init = {
      ## Initialise module
      sim <- firePropertiesInit(sim)

      ## schedule future event(s)
      sim <- scheduleEvent(sim, P(sim)$fireInitialTime, "Biomass_fireProperties",
                           "calcFireProperties", eventPriority = 2.25) ## always calculate fire parameters before the first fire time
    },
    calcFireProperties = {
      ## in the first year of fire always calculate parameters
      if (time(sim) == P(sim)$fireInitialTime) {
        ## calculate fire parameters
        sim <- calcFBPProperties(sim)
      }

      ## in subsequent years evaluate if parameters are to be calculated again (veg feedbacks = TRUE)
      if (time(sim) == sim$fireYear) {
        if (P(sim)$vegFeedback) {
          ## calculate fire parameters
          sim <- calcFBPProperties(sim)

          ## schedule future event(s)
          ## only calculate parameters in fire years.
          sim <- scheduleEvent(sim, time(sim) + P(sim)$fireTimestep, "Biomass_fireProperties",
                               "calcFireProperties", eventPriority = 2.25)
        }
      }

    },
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}

### module initialization
firePropertiesInit <- function(sim) {
  message(blue("Processing climate and topo. data for fire weather and fuel calculation"))
  cacheTags <- c(currentModule(sim), "firePropertiesInit")
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)

  ## checks
  if (start(sim) == P(sim)$fireInitialTime)
    warning(red("start(sim) and P(sim)$fireInitialTime are the same.\nThis may create bad scheduling with init events"))

  if (!P(sim)$fireWeatherMethod %in% c("average", "sample"))
    stop("'fireWeatherMethod'can only be of type 'average' or 'sample'.")

  ## define first fire year
  sim$fireYear <- as.integer(P(sim)$fireInitialTime)

  ## MAKE TOPO DATA ------------------------------------------
  ## extract slope and aspect from DEM raster -
  ## use gdalUtils::gdaldem instead of raster::terrain which was not giving consistent no. of NAs across machines
  ## don't cache, because Cache won't be able to tell if the input raster  changed if the name hasn't
  slopeRas <- gdaldem(mode = "slope",
                      input_dem = filename(sim$DEMRas),
                      output_map = normalizePath(file.path(inputPath(sim), "slopeRas.tif")),
                      compute_edges = TRUE,
                      p = TRUE)
  if (file.exists(slopeRas)) {
    slopeRas <- raster(slopeRas)
  } else {
    stop("Could not calculate/save slope raster from DEM")
  }

  aspectRas <- gdaldem(mode = "aspect",
                       input_dem = filename(sim$DEMRas),
                       output_map = file.path(inputPath(sim), "aspectRas.tif"),
                       compute_edges = TRUE,
                       trigonometric = TRUE)
  if (file.exists(aspectRas)) {
    aspectRas <- raster(aspectRas)
  } else {
    stop("Could not calculate/save aspect raster from DEM")
  }

  ## if they differ in number of NAs, input data from neighbours
  if (sum(!is.na(slopeRas[])) < sum(!is.na(sim$DEMRas[]))) {
    pixMismatch <- which(is.na(slopeRas[]) & !is.na(sim$DEMRas[]))
    slopeRas <- inputValsFromNgbs(ras = slopeRas, cells = pixMismatch)
  }

  if (sum(!is.na(aspectRas[])) < sum(!is.na(sim$DEMRas[]))) {
    pixMismatch <- which(is.na(aspectRas[]) & !is.na(sim$DEMRas[]))
    aspectRas <- inputValsFromNgbs(ras = aspectRas, cells = pixMismatch)
  }

  ## make points and reproject
  slopePoints <- st_as_sf(rasterToPoints(slopeRas, spatial = TRUE))
  slopePoints <- st_transform(slopePoints, crs = crs(sim$studyAreaFBP))
  aspectPoints <- st_as_sf(rasterToPoints(aspectRas, spatial = TRUE))
  aspectPoints <- st_transform(aspectPoints, crs = crs(sim$studyAreaFBP))

  names(slopePoints) <- c("slope", "geometry")
  names(aspectPoints) <- c("aspect", "geometry")

  ## join to get IDs
  slopePoints <- st_join(slopePoints, sim$rasterToMatchFBPPoints,
                         join = st_nearest_feature)
  aspectPoints <- st_join(aspectPoints, sim$rasterToMatchFBPPoints,
                          join = st_nearest_feature)
  slopePoints$rasterToMatch <- NULL
  aspectPoints$rasterToMatch <- NULL

  ## make a table
  sim$topoData <- data.table(pixelIndex = slopePoints$pixelIndex,
                             longitude = st_coordinates(slopePoints)[,"X"],
                             latitude = st_coordinates(slopePoints)[,"Y"],
                             slope = slopePoints$slope,
                             aspect = aspectPoints$aspect)

  ## check that the number of points matches
  if (getOption("LandR.assertions")) {
    if (length(unique(c(crs(slopePoints), crs(aspectPoints), crs(sim$rasterToMatchFBPPoints)))) > 1)
      stop("Reprojecting topography data to FBP-compatible lat/long projection failed.\n
             Please inspect Biomass_fireProperties::firePropertiesInit")

    if (length(unique(c(nrow(slopePoints), nrow(aspectPoints), nrow(sim$rasterToMatchFBPPoints)))) > 1)
      stop("Topography data layers and rasterToMatch differ in number of points in FBP-compatible lat/long projection.\n
             Please inspect Biomass_fireProperties::firePropertiesInit")
  }

  ## MAKE FIRE WEATHER --------------------------------------
  ## check if there are month and day columns, if not add them with defaults
  addCols <- c("month", "day")
  addCols <- addCols[!c("month", "day") %in% names(sim$weatherData)]
  if (length(addCols) == 2) {
    sim$weatherData$month <- 7  ## default to July
    sim$weatherData$day <- 1
  } else {
    if (length(addCols) == 1) {
      if (addCols == "month")
        sim$weatherData$month <- 7 else
          sim$weatherData$day <- 1
    }
  }

  ## check if there's only one weather value per point
  coords <- data.table(st_coordinates(sim$weatherData))
  coords <- unique(coords)
  if (!nrow(coords) == nrow(sim$weatherData)) {
    ## METHOD 1 - AVERAGE ACROSS FIRE DAYS -----
    if (P(sim)$fireWeatherMethod == "average") {
      message(blue("Calculating average weather data per point in 'weatherData'"))
      weatherDataShort <- data.table(st_drop_geometry(sim$weatherData))

      ## add point IDS
      weatherDataShort[, `:=`(X = st_coordinates(sim$weatherData)[,"X"],
                              Y = st_coordinates(sim$weatherData)[,"Y"])]
      coords <- coords[, ID := 1:nrow(coords)]
      weatherDataShort <- coords[weatherDataShort, on = .(X, Y)]

      weatherDataShort <- weatherDataShort[, list(longitude = X,
                                                  latitude = Y,
                                                  month = round(mean(month)),
                                                  day = round(mean(day)),
                                                  temperature = mean(temperature),
                                                  precipitation = mean(precipitation),
                                                  relativeHumidity = mean(relativeHumidity),
                                                  windSpeed = mean(windSpeed)),
                                           by = ID] %>%
        unique(.)
    }

    if (P(sim)$fireWeatherMethod == "sample") {
      ## METHOD 2 - SAMPLE FIRE DAYS RANDOMLY -----
      ## (needs to be done once each fire year)
      message(blue("Sampling weather data (one fire day) per point in 'weatherData'"))
      weatherDataShort <- data.table(st_drop_geometry(sim$weatherData))

      ## add point IDS
      weatherDataShort[, `:=`(X = st_coordinates(sim$weatherData)[,"X"],
                              Y = st_coordinates(sim$weatherData)[,"Y"])]
      coords <- coords[, ID := 1:nrow(coords)]
      weatherDataShort <- coords[weatherDataShort, on = .(X, Y)]

      ## sample a fire day per point
      weatherDataShort[, sampDay := sample(1:length(.I), 1), by = ID]
      weatherDataShort <- weatherDataShort[, list(longitude = X[sampDay],
                                                  latitude = Y[sampDay],
                                                  month = month[sampDay],
                                                  day = day[sampDay],
                                                  temperature = temperature[sampDay],
                                                  precipitation = precipitation[sampDay],
                                                  relativeHumidity = relativeHumidity[sampDay],
                                                  windSpeed = windSpeed[sampDay]),
                                           by = ID] %>%
        unique(.)
    }

    ## first convert to sf
    weatherDataShort <- st_as_sf(weatherDataShort, coords = c("longitude", "latitude"),
                                 crs = sim$weatherDataCRS, agr = "constant")
    ## export to sim - don't replace weatherData - if sampling we need the full data intact
    sim$weatherDataShort <- weatherDataShort
    rm(weatherDataShort)
    .gc()
  } else {
    message(blue("Only one weather value found per spatial point. No sampling or averaging of weather data will be done"))
    sim$weatherDataShort <- sim$weatherData
  }

  ## INTERPOLATE WEATHER TO RASTERTOMATCH AND RE-MAKE SHORT WEATHER TABLE
  ## to avoid data losses all re-projections need to be made in vector format:
  ## 1. re-project points to match original RTM crs
  ## 2. interpolate using points
  ## 3. make a table using pixIDs
  ## 4. add coordinates
  message(blue("Interpolating weather data for extraction of weather values for each simulated pixel"))
  weatherDataShort <- sim$weatherDataShort
  weatherDataShort <- st_transform(weatherDataShort, crs = st_crs(sim$rasterToMatch))
  weatherDataShort <- as_Spatial(weatherDataShort) ## convert to sp for gstat

  ## nearest neighbour interpolation using 8 neighbours
  ## although interpolating dates is very weird, so is climate data from different dates.
  fields <- c("precipitation", "temperature", "windSpeed", "relativeHumidity", "month", "day")
  weatherDataIntrpList <- lapply(fields, FUN = function(x) {
    form <- as.formula(paste(x, "~ 1"))
    interpModel <- gstat(formula = form, data = weatherDataShort, set = list(idp = 0),
                         nmax = 8)
    # interpModel <- gstat(formula = form, locations = weatherDataShort)   ## for IDW interpolation
    weatherDataIntrpRas <- interpolate(object = sim$rasterToMatch, model = interpModel)  ## interpolate on RTM
    mask(weatherDataIntrpRas, sim$rasterToMatch)
  })

  names(weatherDataIntrpList) <- fields

  weatherDataIntrpPointsList <- lapply(weatherDataIntrpList, FUN = function(x) {
    x <- st_as_sf(rasterToPoints(x, spatial = TRUE))
    x$pixelIndex <- which(!is.na(getValues(sim$rasterToMatch)))
    x <- st_transform(x, crs = crs(sim$studyAreaFBP))
    x
  })

  names(weatherDataIntrpPointsList) <- fields

  weatherDataShort2 <- data.table(pixelIndex = st_drop_geometry(sim$rasterToMatchFBPPoints)$pixelIndex,
                                  longitude = st_coordinates(sim$rasterToMatchFBPPoints)[, "X"],
                                  latitude = st_coordinates(sim$rasterToMatchFBPPoints)[, "Y"])
  weatherDataShort3 <- data.table(pixelIndex = weatherDataIntrpPointsList$precipitation$pixelIndex,
                                  month = as.integer(weatherDataIntrpPointsList$month$var1.pred),
                                  day = as.integer(weatherDataIntrpPointsList$day$var1.pred),
                                  precipitation = weatherDataIntrpPointsList$precipitation$var1.pred,
                                  temperature = weatherDataIntrpPointsList$temperature$var1.pred,
                                  windSpeed = weatherDataIntrpPointsList$windSpeed$var1.pred,
                                  relativeHumidity = weatherDataIntrpPointsList$relativeHumidity$var1.pred)
  ## join and make sure no pixels were lost
  weatherDataShort3 <- weatherDataShort3[weatherDataShort2, on = .(pixelIndex)]

  if (getOption("LandR.assertions")) {
    if (nrow(weatherDataShort2) != nrow(weatherDataShort3))
      stop("Some pixels have no weather data")
  }

  ## export to sim
  sim$weatherDataShort <- weatherDataShort3

  return(invisible(sim))
}

## Derive fire parameters from FBP system - rasters need to be in lat/long
calcFBPProperties <- function(sim) {
  cacheTags <- c(currentModule(sim), "FBPPercParams")

  if (getOption("LandR.assertions")) {
    if (length(unique(c(ncell(sim$fuelTypesMaps$finalFuelType), ncell(sim$fuelTypesMaps$coniferDom),
                        ncell(sim$fuelTypesMaps$curing), ncell(sim$rasterToMatch)))) > 1)
      stop("One or more 'fuelTypesMaps' rasters have different number of cells from 'rasterToMatch'")
  }

  ## FUEL TYPES ------------------------------
  ## reproject to fuelTypesMaps to FBP-compatible this results in the loss of many pixels
  ## so use spatial points for projection. pixelIDS need to be made into a column, since
  ## there may be fewer points than in RTMFBPPoints
  ## each object is then converted to a DT before joining
  fuelTypePoints <- st_as_sf(rasterToPoints(sim$fuelTypesMaps$finalFuelType, spatial = TRUE))
  fuelTypePoints$pixelIndex <- which(!is.na(getValues(sim$fuelTypesMaps$finalFuelType)))
  fuelTypePoints <- st_transform(fuelTypePoints, crs = crs(sim$studyAreaFBP))

  coniferDomPoints <- st_as_sf(rasterToPoints(sim$fuelTypesMaps$coniferDom, spatial = TRUE))
  coniferDomPoints$pixelIndex <- which(!is.na(getValues(sim$fuelTypesMaps$coniferDom)))
  coniferDomPoints <- st_transform(coniferDomPoints, crs = crs(sim$studyAreaFBP))

  curingPoints <- st_as_sf(rasterToPoints(sim$fuelTypesMaps$curing, spatial = TRUE))
  curingPoints$pixelIndex <- which(!is.na(getValues(sim$fuelTypesMaps$curing)))
  curingPoints <- st_transform(curingPoints, crs = crs(sim$studyAreaFBP))

  if (getOption("LandR.assertions")) {
    if (length(unique(c(crs(fuelTypePoints), crs(coniferDomPoints),
                        crs(curingPoints), crs(sim$rasterToMatchFBPPoints)))) > 1)
      stop("Reprojecting climate data to FBP-compatible lat/long projection failed.\n
             Please inspect Biomass_fireProperties::firePropertiesInit")
  }

  fuelTypeDT <- as.data.table(st_drop_geometry(fuelTypePoints))
  coniferDomDT <- as.data.table(st_drop_geometry(coniferDomPoints))
  curingDT <- as.data.table(st_drop_geometry(curingPoints))
  RTMpixelIndex <- data.table(pixelIndex = sim$rasterToMatchFBPPoints$pixelIndex)

  ## make table of final fuel types by merging
  FTs <- Reduce(function(x,y) merge.data.table(x, y, by = "pixelIndex", all = TRUE),
                list(fuelTypeDT, coniferDomDT, curingDT, RTMpixelIndex))
  setnames(FTs, old = names(FTs), new = c("pixelIndex", "FuelType", "coniferDom", "curing"))

  ## remove pixels with no fuels
  pixNoFuels <- FTs[is.na(FuelType) & is.na(coniferDom) & is.na(curing), pixelIndex]
  FTs <- FTs[!pixelIndex %in% pixNoFuels]

  ## add FBP forest fuel type names
  FTequivTable <- as.data.table(raster::levels(sim$fuelTypesMaps$finalFuelType)[[1]])
  setnames(FTequivTable, "ID", "FuelType")
  FTs <- unique(FTequivTable[, .(FuelTypeFBP, FuelType)])[FTs, on = "FuelType"]

  ## add FBP non-forest fuel type names if running fires in non-forested pixels
  if (!is.null(sim$pixelNonForestFuels)) {
    ## unfortunately some pixels are lost here.
    FTs <- unique(sim$nonForestFuelsTable[,.(FuelTypeFBP, FuelType)])[FTs, on = "FuelType"]
    ## a duplicated column was created, containing the missing  fuel type names
    FTs[is.na(i.FuelTypeFBP) & !is.na(FuelTypeFBP), i.FuelTypeFBP := FuelTypeFBP]

    ## remove duplicated column and rename
    FTs[, FuelTypeFBP := NULL]
    setnames(FTs, old = "i.FuelTypeFBP", new = "FuelTypeFBP")
  }

  FTs <- FTs[!is.na(FuelType)]

  ## replace all D2s for D1s - cffdrs::fbp does not include D2
  ## need to make a new FuelType if D1 doesn't exist already
  if (any(FTequivTable$FuelTypeFBP == "D2")) {
    message(blue("Replacing all D2 fuel types by D1, as cffdrs::fbp does not include D2"))
    D1No <- unique(FTequivTable[FuelTypeFBP == "D1", FuelType])
    if (length(D1No)) {
      FTs[FuelTypeFBP == "D2", `:=` (FuelTypeFBP = "D1", FuelType = D1No)]
    } else {
      D1No <- max(FTequivTable$FuelType) + 1
      FTs[FuelTypeFBP == "D2", `:=` (FuelTypeFBP = "D1", FuelType = D1No)]

      rclmat <- matrix(c(FTequivTable[FuelTypeFBP == "D2", FuelType], D1No),
                       nrow = 1, ncol = 2, byrow = TRUE)
      sim$fuelTypesMaps$finalFuelType <- reclassify(sim$fuelTypesMaps$finalFuelType, rclmat)
      sim$fuelTypesMaps$finalFuelType <- ratify(sim$fuelTypesMaps$finalFuelType)
      levs <- as.data.table(raster::levels(sim$fuelTypesMaps$finalFuelType)[[1]])
      levs2 <- unique(FTs[, .(FuelType, FuelTypeFBP)])
      setnames(levs2, "FuelType", "ID")
      levs <- levs2[levs, on = "ID"]
      levels(sim$fuelTypesMaps$finalFuelType) <- as.data.frame(levs)
      setColors(sim$fuelTypesMaps$finalFuelType, n = nrow(levs)) <- brewer.pal(n = nrow(levs), "Accent")
    }
  }

  ## check for duplicates (there shouldn't be any)
  if (getOption("LandR.assertions"))
    if (any(duplicated(FTs))) {
      stop("Duplicated pixels found in fuel types table.",
           " Please debug Biomass_fireProperties calcFBPProperties() event function")
    }

  ## check all pixel indices are part of RTM
  if (getOption("LandR.assertions"))
    if (!all(FTs$pixelIndex %in% which(!is.na(getValues(sim$rasterToMatch)))) &
        !all(FTs$pixelIndex %in% sim$rasterToMatchFBPPoints$pixelIndex)) {
      stop("Wrong pixel indices attributed to fuel types.\n",
           "Please debug Biomass_fireProperties::calcFireProperties")
    }

  ## FWI ------------------------------
  ## make/update table of FWI inputs
  FWIinputs <- data.frame(id = sim$weatherDataShort$pixelIndex,
                          lat = sim$weatherDataShort$latitude,
                          long = sim$weatherDataShort$longitude,
                          day = sim$weatherDataShort$day,
                          mon = sim$weatherDataShort$month,
                          temp = sim$weatherDataShort$temperature,
                          rh = sim$weatherDataShort$relativeHumidity,
                          ws = sim$weatherDataShort$windSpeed,
                          prec = sim$weatherDataShort$precipitation)

  if (getOption("LandR.assertions"))
    if (!all(FTs$pixelIndex %in% FWIinputs$id)) {
      warning("Some pixels with fuels have no climate data.\n",
              "Please debug Biomass_fireProperties::calcFireProperties")
    }

  ## calculate FW indices
  FWIoutputs <- suppressWarnings({
    cffdrs::fwi(input = FWIinputs,
                init = na.omit(sim$FWIinit),
                batch = FALSE, lat.adjust = TRUE)
  })
  FWIoutputs <- data.table(FWIoutputs)

  ## FBP -----------------------------
  ## make inputs dataframe for FBI
  ## add fuel types and conifer dominance to FWIOutputs
  ## note that because climate/topo data is "larger" there are pixels that have no fuels - these are removed.

  ## check pixelIndex matches (there shouldn't be any)

  if (getOption("LandR.assertions"))
    if (!all(FTs$pixelIndex %in% FWIoutputs$ID)) {
      warning("Some pixels with fuels have no fire weather indices.\n",
              "Please debug Biomass_fireProperties::calcFireProperties")
    }
  setnames(FWIoutputs, "ID", "pixelIndex")
  FWIoutputs <- FTs[FWIoutputs, on = "pixelIndex", nomatch = 0]

  ## add slope, aspect and day/month (convert date to julian day)
  ## again, only keep pixels that have fuels
  FWIoutputs <- sim$topoData[, .(pixelIndex, slope, aspect)][FWIoutputs, on = "pixelIndex", nomatch = 0]
  FWIinputs <- data.table(FWIinputs)
  setnames(FWIinputs, "id", "pixelIndex")
  FWIoutputs <- FWIinputs[, .(pixelIndex, day, mon)][FWIoutputs, on = "pixelIndex", nomatch = 0]
  FWIoutputs[, `:=`(mon = as.character(mon),
                    day = as.character(day))]
  FWIoutputs[nchar(mon) == 1, mon := paste0("0", mon)]
  FWIoutputs[nchar(day) == 1, day := paste0("0", day)]
  FWIoutputs[, julDay := julian(as.Date(paste(1900, mon, day, sep = "-")), origin = as.Date("1900-01-01"))]

  FBPinputs <- data.frame(id = FWIoutputs$pixelIndex,
                          FuelType = FWIoutputs$FuelTypeFBP,
                          LAT = FWIoutputs$LAT,
                          LONG = FWIoutputs$LONG,
                          FFMC = FWIoutputs$FFMC,
                          BUI = FWIoutputs$BUI,
                          WS = FWIoutputs$WS,
                          GS = FWIoutputs$slope,
                          Dj = FWIoutputs$julDay,
                          Aspect = FWIoutputs$aspect,
                          PC = FWIoutputs$coniferDom,
                          cc = FWIoutputs$curing)

  if (all(is.na(FBPinputs$cc))) {
    FBPinputs$cc <- NULL
  }

  FBPoutputs <- suppressWarnings({
    cffdrs::fbp(input = FBPinputs, output = "All")
  })
  FBPoutputs <- data.table(FBPoutputs)

  ## RASTERIZE FBP OUTPUTS BY PIXEL INDEX
  ## Crown fraction burnt
  sim$fireCFBRas <- suppressWarnings(setValues(sim$rasterToMatch, rep(NA, ncell(sim$rasterToMatch))))
  sim$fireCFBRas[FBPoutputs$ID] <- FBPoutputs$CFB

  ## Head fire intensity
  sim$fireIntRas <- suppressWarnings(setValues(sim$rasterToMatch, rep(NA, ncell(sim$rasterToMatch))))
  sim$fireIntRas[FBPoutputs$ID] <- FBPoutputs$HFI

  ## Rate of spread
  sim$fireROSRas <- suppressWarnings(setValues(sim$rasterToMatch, rep(NA, ncell(sim$rasterToMatch))))
  sim$fireROSRas[FBPoutputs$ID] <- FBPoutputs$ROS

  ## Critical spread rate for crowning
  sim$fireRSORas <- suppressWarnings(setValues(sim$rasterToMatch, rep(NA, ncell(sim$rasterToMatch))))
  sim$fireRSORas[FBPoutputs$ID] <- FBPoutputs$RSO

  ## Total fuel consumption
  sim$fireTFCRas <- suppressWarnings(setValues(sim$rasterToMatch, rep(NA, ncell(sim$rasterToMatch))))
  sim$fireTFCRas[FBPoutputs$ID] <- FBPoutputs$TFC

  return(invisible(sim))
}

## OTHER INPUTS AND FUNCTIONS --------------------------------
.inputObjects = function(sim) {
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  cacheTags = c(currentModule(sim), "function:.inputObjects")

  ## project to Lat/Long (decimal degrees) for compatibility with FBP system
  latLong <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

  if (!suppliedElsewhere("studyAreaFBP", sim)) {
    if (!suppliedElsewhere("studyArea", sim)) {
      stop("Please provide a 'studyArea' polygon")
      # message("'studyArea' was not provided by user. Using a polygon (6250000 m^2) in southwestern Alberta, Canada")
      # sim$studyArea <- randomStudyArea(seed = 1234, size = (250^2)*100)  # Jan 2021 we agreed to force user to provide a SA/SAL
    }
    sim$studyAreaFBP <- sim$studyArea
  }

  ## if necessary reproject to lat/long - for compatibility with FBP
  if (!compareCRS(latLong, crs(sim$rasterToMatch))) {
    sim$studyAreaFBP <- spTransform(sim$studyAreaFBP, latLong) #faster without Cache
  }

  if (is.na(P(sim)$.studyAreaName)) {
    params(sim)[[currentModule(sim)]][[".studyAreaName"]] <- reproducible::studyAreaName(sim$studyAreaLarge)
    message("The .studyAreaName is not supplied; derived name from sim$studyAreaLarge: ",
            params(sim)[[currentModule(sim)]][[".studyAreaName"]])
  }

  ## RASTER TO MATCH
  needRTM <- FALSE
  if (is.null(sim$rasterToMatch)) {
    if (!suppliedElsewhere("rasterToMatch", sim)) {      ## if one is not provided, re do both (safer?)
      needRTM <- TRUE
      message("There is no rasterToMatch supplied; will attempt to use rawBiomassMap")
    } else {
      stop("rasterToMatch is going to be supplied, but ", currentModule(sim), " requires it ",
           "as part of its .inputObjects. Please make it accessible to ", currentModule(sim),
           " in the .inputObjects by passing it in as an object in simInit(objects = list(rasterToMatch = aRaster)",
           " or in a module that gets loaded prior to ", currentModule(sim))
    }
  }

  if (needRTM) {
    if (!suppliedElsewhere("rawBiomassMap", sim) ||
        !compareRaster(sim$rawBiomassMap, sim$studyArea, stopiffalse = FALSE)) {
      rawBiomassMapURL <- paste0("http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
                                 "canada-forests-attributes_attributs-forests-canada/",
                                 "2001-attributes_attributs-2001/",
                                 "NFI_MODIS250m_2001_kNN_Structure_Biomass_TotalLiveAboveGround_v1.tif")
      rawBiomassMapFilename <- "NFI_MODIS250m_2001_kNN_Structure_Biomass_TotalLiveAboveGround_v1.tif"
      rawBiomassMap <- Cache(prepInputs,
                             targetFile = rawBiomassMapFilename,
                             url = rawBiomassMapURL,
                             destinationPath = dPath,
                             studyArea = sim$studyArea,
                             rasterToMatch = NULL,
                             maskWithRTM = FALSE,
                             useSAcrs = FALSE,     ## never use SA CRS
                             method = "bilinear",
                             datatype = "INT2U",
                             filename2 = NULL,
                             userTags = c(cacheTags, "rawBiomassMap"),
                             omitArgs = c("destinationPath", "targetFile", "userTags", "stable"))
    } else {
      rawBiomassMap <- Cache(postProcess,
                             x = sim$rawBiomassMap,
                             studyArea = sim$studyArea,
                             useSAcrs = FALSE,
                             maskWithRTM = FALSE,   ## mask with SA
                             method = "bilinear",
                             datatype = "INT2U",
                             filename2 = NULL,
                             overwrite = TRUE,
                             userTags = cacheTags,
                             omitArgs = c("destinationPath", "targetFile", "userTags", "stable"))
    }

    ## if we need rasterToMatchLarge, that means a) we don't have it, but b) we will have rawBiomassMap
    if (is.null(sim$rasterToMatch))
      warning(paste0("rasterToMatchLarge is missing and will be created \n",
                     "from rawBiomassMap and studyAreaLarge.\n
              If this is wrong, provide raster"))

    sim$rasterToMatch <- rawBiomassMap
    RTMvals <- getValues(sim$rasterToMatch)
    sim$rasterToMatch[!is.na(RTMvals)] <- 1

    sim$rasterToMatch <- Cache(writeOutputs,
      sim$rasterToMatch,
      filename2 = .suffix(file.path(dPath, "rasterToMatch.tif"), paste0("_", P(sim)$.studyAreaName)),
      datatype = "INT2U",
      overwrite = TRUE,
      userTags = c(cacheTags, "rasterToMatch"),
      omitArgs = c("userTags"))
  }

  if (!compareCRS(sim$studyArea, sim$rasterToMatch)) {
    warning(paste0("studyArea and rasterToMatch projections differ.\n",
                   "studyArea will be projected to match rasterToMatchLarge"))
    sim$studyArea <- spTransform(sim$studyArea, crs(sim$rasterToMatch))
    sim$studyArea <- fixErrors(sim$studyArea)
  }

  ## convert all inputs to Lat/Long (decimal degrees)
  ## for compatibility with FBP system
  ## convert to points before projecting to avoid data loss.

  ## convert rasterToMatch to spatial points and project
  ## to prevent data loss
  ## note: don't mask to area until the end.
  sim$rasterToMatchFBPPoints <- st_as_sf(rasterToPoints(sim$rasterToMatch, spatial = TRUE))
  sim$rasterToMatchFBPPoints$pixelIndex <- which(!is.na(getValues(sim$rasterToMatch)))
  sim$rasterToMatchFBPPoints <- st_transform(sim$rasterToMatchFBPPoints, crs = crs(sim$studyAreaFBP))

  ## DEFAULT WEATHER DATA
  ## these defaults are only necessary if weatherData is not supplied by another module
  ## climate defaults to Climate NA Data, year 2011, RCP4.5
  ## note that some Climate NA data were multiplied by 10
  if (!suppliedElsewhere("weatherData", sim)) {
    message(blue("Getting default 'temperatureRas' to make default 'weatherData'.",
                 "If this is not correct, make sure Biomass_fireProperties can detect 'weatherData' is supplied"))
    ## get default temperature values, summer average
    temperatureRas <- Cache(prepInputs, targetFile = "Tave_sm.asc",
                            url = extractURL("temperatureRas", sim),
                            archive = "CanESM2_RCP45_r11i1p1_2011MSY.zip",
                            alsoExtract = NA,
                            destinationPath = dPath,
                            fun = "raster",
                            filename2 = NULL,
                            userTags = cacheTags)
    temperatureRas <- temperatureRas/10  ## back transform temp values

    ## add the original CRS if it's not defined
    if (is.na(crs(temperatureRas)))
      crs(temperatureRas) <- "+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"

    temperatureRas <- Cache(postProcess,
                            x = temperatureRas,
                            rasterToMatch = sim$rasterToMatch,
                            maskWithRTM = TRUE,
                            method = "bilinear",
                            filename2 = NULL,
                            userTags = c(cacheTags, "temperatureRas"))

    message(blue("Getting default 'precipitationRas' to make default 'weatherData'.",
                 "If this is not correct, make sure Biomass_fireProperties can detect 'weatherData' is supplied"))
    ## get default precipitation values, summer cummulative precipitation
    precipitationRas <- Cache(prepInputs, targetFile = "PPT_sm.asc",
                              url = extractURL("precipitationRas", sim),
                              archive = "CanESM2_RCP45_r11i1p1_2011MSY.zip",
                              alsoExtract = NA,
                              destinationPath = dPath,
                              fun = "raster",
                              filename2 = NULL,
                              userTags = cacheTags)

    ## add the original CRS if it's not defined
    if (is.na(crs(precipitationRas)))
      crs(precipitationRas) <- "+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"

    precipitationRas <- Cache(postProcess,
                              x = precipitationRas,
                              rasterToMatch = sim$rasterToMatch,
                              maskWithRTM = TRUE,
                              method = "bilinear",
                              filename2 = NULL,
                              userTags = c(cacheTags, "precipitationRas"))

    message(blue("Getting default 'relativeHumRas' to make default 'weatherData'.",
                 "If this is not correct, make sure Biomass_fireProperties can detect 'weatherData' is supplied"))
    ## get default precipitation values, summer average
    relativeHumRas <- Cache(prepInputs, targetFile = "RH_sm.asc",
                            url = extractURL("relativeHumRas", sim),
                            archive = "CanESM2_RCP45_r11i1p1_2011MSY.zip",
                            alsoExtract = NA,
                            destinationPath = dPath,
                            fun = "raster",
                            filename2 = NULL,
                            userTags = cacheTags)

    ## add the original CRS if it's not defined
    if (is.na(crs(relativeHumRas)))
      crs(relativeHumRas) <- "+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"

    relativeHumRas <- Cache(postProcess,
                            x = relativeHumRas,
                            rasterToMatch = sim$rasterToMatch,
                            maskWithRTM = TRUE,
                            method = "bilinear",
                            filename2 = NULL,
                            userTags = c(cacheTags, "relativeHumRas"))

    ## PROJECT CLIMATE/TOPO RASTERS AS POINTS
    temperaturePoints <- st_as_sf(rasterToPoints(temperatureRas, spatial = TRUE))
    temperaturePoints <- st_transform(temperaturePoints, crs = crs(sim$studyAreaFBP))

    precipitationPoints <- st_as_sf(rasterToPoints(precipitationRas, spatial = TRUE))
    precipitationPoints <- st_transform(precipitationPoints, crs = crs(sim$studyAreaFBP))

    relativeHumPoints <- st_as_sf(rasterToPoints(relativeHumRas, spatial = TRUE))
    relativeHumPoints <- st_transform(relativeHumPoints, crs = crs(sim$studyAreaFBP))

    if (getOption("LandR.assertions")) {
      if (length(unique(c(crs(temperaturePoints), crs(precipitationPoints),
                          crs(relativeHumPoints), crs(sim$rasterToMatchFBPPoints)))) > 1)
        stop("Reprojecting climate data to FBP-compatible lat/long projection failed.\n
             Please inspect Biomass_fireProperties::firePropertiesInit")

      if (length(unique(c(nrow(temperaturePoints), nrow(precipitationPoints),
                          nrow(relativeHumPoints), nrow(sim$rasterToMatchFBPPoints)))) > 1)
        stop("Climate data layers and rasterToMatch differ in number of points in FBP-compatible lat/long projection.\n
             Please inspect Biomass_fireProperties::firePropertiesInit")
    }

    message(blue(currentModule(sim), " is making 'weatherData' from default temperature, precipitation and relative humidity raster layers"))
    sim$weatherDataCRS <- crs(sim$rasterToMatchFBPPoints)
    weatherData <- data.table(temperature = temperaturePoints[, 1, drop = TRUE],
                              precipitation = precipitationPoints[, 1, drop = TRUE],
                              relativeHumidity = relativeHumPoints[, 1, drop = TRUE],
                              month = 7,  ## consider July.
                              day = 1,
                              windSpeed = 0, ## consider no wind
                              latitude = st_coordinates(sim$rasterToMatchFBPPoints)[,2],
                              longitude = st_coordinates(sim$rasterToMatchFBPPoints)[,1])
    weatherData <- na.omit(weatherData)
    weatherData <- st_as_sf(weatherData, coords = c("longitude", "latitude"),
                            crs = sim$weatherDataCRS, agr = "constant")

    ## this is no longer necessary as ClimateNA has relative humidity data
    ## relative humidity
    ## using dew point between -3 and 20%, quarterly seasonal for Jun 2013
    ## https://calgary.weatherstats.ca/metrics/dew_point.html
    # weatherData[, relHum := RH(t = weatherData$temp, Td = runif (nrow(weatherData), -3, 20), isK = FALSE)]

    ## export to sim
    sim$weatherData <- weatherData
  } else {
    if (!suppliedElsewhere("weatherDataCRS", sim))
      stop(red("'weatherData' appears to be supplied to Biomass_fireProperties,",
               "but not weatherDataCRS. Please make sure 'weatherDataCRS' is also provided."))
  }

  ## DEM RASTER
  if (!suppliedElsewhere("DEMRas", sim)) {
    sim$DEMRas <- Cache(prepInputs, targetFile = "DEM1kmRes.tif",
                        alsoExtract = "DEM1kmRes.prj",
                        archive = "DEMrasters.zip",
                        destinationPath = dPath,
                        url = extractURL("DEMRas", sim),
                        rasterToMatch = sim$rasterToMatch,
                        maskWithRTM = TRUE,
                        method = "bilinear",
                        filename2 = .suffix("DEMRas.tif", paste0("_", P(sim)$.studyAreaName)),
                        overwrite = TRUE,
                        useSAcrs = FALSE,
                        userTags = c(cacheTags, "DEMRas"),
                        omitArgs = "userTags")

    ## the DEM used in BioSIM loses pixels (has more NA's) after reproj.
    ## this is a quickfix - best solution is to get another, larger DEM:
    ## input average ngb values
    if (sum(!is.na(sim$DEMRas[])) < sum(!is.na(sim$rasterToMatch[]))) {
      pixMismatch <- which(is.na(sim$DEMRas[]) & !is.na(sim$rasterToMatch[]))
      sim$DEMRas <- inputValsFromNgbs(ras = sim$DEMRas, cells = pixMismatch)
    }
  }

  ## FWI INITIALISATION DATAFRAME
  ## TODO: FWIinit should be updated every year from previous year's/days/months results
  if (!suppliedElsewhere("FWIinit", sim)) {
    sim$FWIinit <- data.frame(ffmc = 85, dmc = 6, dc = 15)
  }

  if (!suppliedElsewhere("pixelNonForestFuels", sim)) {
    sim$pixelNonForestFuels <- NULL
  }

  return(invisible(sim))
}
