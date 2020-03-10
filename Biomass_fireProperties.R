# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects and functions, use sim$xxx.
defineModule(sim, list(
  name = "Biomass_fireProperties",
  description = "Complement to fire spread that calculates fire (behaviour) properties in fucntion of vegetation, climate and topography conditions
  usign the Canadian Forest Fire Behaviour Prediction System",
  keywords = c("fire behaviour", "fuels", "fire-vegetation feedbacks", "fire-climate feedbacks", "FBP system", "topography"),
  authors = person("Ceres", "Barros", email = "cbarros@mail.ubc.ca", role = c("aut", "cre")),
  childModules = character(0),
  version = list(Biomass_fireProperties = numeric_version("0.1.0"),
                 Biomass_core = numeric_version("1.3.2"),
                 LandR = "0.0.3.9000", SpaDES.core = "0.2.7"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "Biomass_fireProperties.Rmd"),
  reqdPkgs = list("R.utils", "raster", "data.table", "dplyr",
                  "sp", "sf", "cffdrs", "amc",
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
    defineParameter(".useCache", "logical", "init", NA, NA,
                    desc = "use caching for the spinup simulation?")
  ),
  inputObjects = bind_rows(
    expectsInput(objectName = "ForestFuelTypes", objectClass = "data.table",
                 desc = "Table of Fuel Type parameters, with  base fuel type, species (in LANDIS code), their - or + contribution ('negSwitch'),
                 min and max age for each species"),
    expectsInput(objectName = "fuelTypesMaps", objectClass = "list",
                 desc = "List of RasterLayers of fuel types and coniferDominance per pixel."),
    expectsInput(objectName = "FWIinit", objectClass = "data.frame",
                 desc = "Initalisation parameter values for FWI calculations. Defaults to default values in cffdrs::fwi.
                 This table should be updated every year"),
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
                 sourceURL = ""),
    expectsInput(objectName = "studyAreaFBP", objectClass = "SpatialPolygonsDataFrame",
                 desc = paste("same as studyArea, but on FBP-compatible lat/long projection:",
                              "'+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'"),
                 sourceURL = ""),
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
                               "and on FBP compatible projection")),
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
  cacheTags <- c(currentModule(sim), "firePropertiesInit")
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)

  ## checks
  if (start(sim) == P(sim)$fireInitialTime)
    warning(red("start(sim) and P(sim)$fireInitialTime are the same.\nThis may create bad scheduling with init events"))

  ## define first fire year
  sim$fireYear <- as.integer(P(sim)$fireInitialTime)

  ## MAKE TOPO DATA ------------------------------------------
  ## extract slope and aspect from DEM raster - assume NAs are 0s to remove border effect
  DEMRas <- sim$DEMRas
  DEMRas[is.na(DEMRas)] <- 0
  tempBrick <- terrain(DEMRas, opt = c("slope", "aspect"))
  tempBrick <- mask(tempBrick, sim$DEMRas)
  slopeRas <- tempBrick$slope
  aspectRas <- tempBrick$aspect
  rm(tempBrick, DEMRas);
  .gc()

  ## make points and reproject
  slopePoints <- st_as_sf(rasterToPoints(slopeRas, spatial = TRUE))
  slopePoints <- st_transform(slopePoints, crs = crs(sim$studyAreaFBP))
  aspectPoints <- st_as_sf(rasterToPoints(aspectRas, spatial = TRUE))
  aspectPoints <- st_transform(aspectPoints, crs = crs(sim$studyAreaFBP))

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
  ## check first if there's only one weather value per point
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


  ## check if IDs match rasterToMatchFBPPoints (they may not, if data is supplied from user/another module)
  ## if they dont, the data needs to be matched
  if (!all(sim$topoClimData$ID %in% sim$rasterToMatchFBPPoints$pixelIndex)) {
    ## TODO: CHECK PROJECTION
    browser()

    message(blue("'topoClimData' does not conform to 'rasterToMatchFBPPoints'.",
                 "Matching the two geographically"))
    ## add a column of point ID
    sim$topoClimData[, pointID := as.factor(paste(longitude, latitude, sep = "_"))]
    sim$topoClimData[, pointID := as.numeric(pointID)]

    ## make a points sf object
    topoClimDataPoints <- unique(sim$topoClimData[, .(pointID, longitude, latitude)])
    coordinates(topoClimDataPoints) <- c("longitude", "latitude")
    topoClimDataPoints <- st_as_sf(topoClimDataPoints)
    st_crs(topoClimDataPoints) <- sim$topoClimCRS

    if (st_crs(topoClimDataPoints) != st_crs(sim$rasterToMatchFBPPoints)) {
      message(blue("Projecting 'topoClimData points to rasterToMatchFBPPoints projection"))
      topoClimDataPoints <- st_transform(topoClimDataPoints, crs = st_crs(rasterToMatchFBPPoints))
    }

  }


  return(invisible(sim))
}

## Derive fire parameters from FBP system - rasters need to be in lat/long
calcFBPProperties <- function(sim) {
  cacheTags <- c(currentModule(sim), "FBPPercParams")
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

  ## make table of final fuel types by joining
  FTs <- Reduce(function(x,y) merge.data.table(x, y, by = "pixelIndex", all = TRUE),
                list(fuelTypeDT, coniferDomDT, curingDT, RTMpixelIndex))
  setnames(FTs, old = names(FTs), new = c("pixelIndex", "FuelType", "coniferDom", "curing"))

  ## remove pixels with no fuels
  pixNoFuels <- FTs[is.na(FuelType) & is.na(coniferDom) & is.na(curing), pixelIndex]
  FTs <- FTs[!pixelIndex %in% pixNoFuels]

  ## add FBP forest fuel type names
  FTs <- unique(sim$ForestFuelTypes[, .(FuelTypeFBP, FuelType)])[FTs, on = "FuelType"]

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
  D1No <- unique(sim$ForestFuelTypes[FuelTypeFBP == "D2", FuelType]) - 1
  FTs[FuelTypeFBP == "D2", `:=` (FuelTypeFBP = "D1",
                                 FuelType = D1No)]


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
  FWIinputs <- data.frame(id = sim$topoClimData$ID,
                          lat = sim$topoClimData$lat,
                          long = sim$topoClimData$long,
                          mon = sim$topoClimData$mon,  ## consider July.
                          temp = sim$topoClimData$temp,
                          rh = sim$topoClimData$relHum,
                          ws = sim$topoClimData$ws,
                          prec = sim$topoClimData$precip)

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
  setnames(FTs, "pixelIndex", "ID")
  FWIoutputs <- FTs[FWIoutputs, on = "ID", nomatch = 0]

  ## add slope and aspect
  ## again, only keep pixels that have fuels
  FWIoutputs <- sim$topoClimData[, .(ID, slope, aspect)][FWIoutputs, on = "ID", nomatch = 0]

  FBPinputs <- data.frame(id = FWIoutputs$ID,
                          FuelType = FWIoutputs$FuelTypeFBP,
                          LAT = FWIoutputs$LAT,
                          LONG = FWIoutputs$LONG,
                          FFMC = FWIoutputs$FFMC,
                          BUI = FWIoutputs$BUI,
                          WS = FWIoutputs$WS,
                          GS = FWIoutputs$slope,
                          Dj = rep(180, nrow(FWIoutputs)),
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
      if (getOption("LandR.verbose", TRUE) > 0)
        message("'studyArea' was not provided by user. Using a polygon (6250000 m^2) in southwestern Alberta, Canada")
      sim$studyArea <- randomStudyArea(seed = 1234, size = (250^2)*100)
    }
    sim$studyAreaFBP <- sim$studyArea
  }

  ## if necessary reproject to lat/long - for compatibility with FBP
  if (!identical(latLong, crs(sim$studyAreaFBP))) {
    sim$studyAreaFBP <- spTransform(sim$studyAreaFBP, latLong) #faster without Cache
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
    if (!suppliedElsewhere("rawBiomassMap", sim)) {
      rawBiomassMapURL <- paste0("http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
                                 "canada-forests-attributes_attributs-forests-canada/",
                                 "2001-attributes_attributs-2001/",
                                 "NFI_MODIS250m_2001_kNN_Structure_Biomass_TotalLiveAboveGround_v1.tif")

      rawBiomassMap <- Cache(prepInputs,
                             targetFile = rawBiomassMapFilename,
                             url = rawBiomassMapURL,
                             destinationPath = dPath,
                             studyArea = sim$studyArea,
                             rasterToMatch = if (!needRTM) sim$rasterToMatch else NULL,
                             # maskWithRTM = TRUE,    ## if RTM not supplied no masking happens (is this intended?)
                             maskWithRTM = if (!needRTM) TRUE else FALSE,
                             ## TODO: if RTM is not needed use SA CRS? -> this is not correct
                             # useSAcrs = if (!needRTM) TRUE else FALSE,
                             useSAcrs = FALSE,     ## never use SA CRS
                             method = "bilinear",
                             datatype = "INT2U",
                             filename2 = NULL,
                             userTags = c(cacheTags, "rawBiomassMap"),
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
    sim$rasterToMatch <- Cache(writeOutputs, sim$rasterToMatch,
                               filename2 = file.path(cachePath(sim), "rasters", "rasterToMatch.tif"),
                               datatype = "INT2U", overwrite = TRUE,
                               userTags = c(cacheTags, "rasterToMatch"),
                               omitArgs = c("userTags"))
  }

  if (!identical(crs(sim$studyArea), crs(sim$rasterToMatch))) {
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

  ## DEFAULT TOPO, TEMPERATURE AND PRECIPITATION
  ## these defaults are only necessary if topoClimData is not supplied by another module
  ## climate defaults to Climate NA Data, year 2011, RCP4.5
  ## note that some Climate NA data were multiplied by 10
  if (!suppliedElsewhere("topoClimData", sim)) {
    message(blue("Getting default 'temperatureRas' to make default 'topoClimData'.",
                 "If this is not correct, make sure Biomass_fireProperties can detect 'topoClimData' is supplied"))
    ## get default temperature values, summer average
    temperatureRas <- Cache(prepInputs, targetFile = "Tave_sm.asc",
                            url = extractURL("temperatureRas", sim),
                            archive = "CanESM2_RCP45_r11i1p1_2011MSY.zip",
                            alsoExtract = NA,
                            destinationPath = dPath,
                            fun = "raster",
                            filename2 = FALSE,
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
                            filename2 = FALSE,
                            userTags = c(cacheTags, "temperatureRas"))

    message(blue("Getting default 'precipitationRas' to make default 'topoClimData'.",
                 "If this is not correct, make sure Biomass_fireProperties can detect 'topoClimData' is supplied"))
    ## get default precipitation values, summer cummulative precipitation
    precipitationRas <- Cache(prepInputs, targetFile = "PPT_sm.asc",
                              url = extractURL("precipitationRas", sim),
                              archive = "CanESM2_RCP45_r11i1p1_2011MSY.zip",
                              alsoExtract = NA,
                              destinationPath = dPath,
                              fun = "raster",
                              filename2 = FALSE,
                              userTags = cacheTags)

    ## add the original CRS if it's not defined
    if (is.na(crs(precipitationRas)))
      crs(precipitationRas) <- "+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"

    precipitationRas <- Cache(postProcess,
                              x = precipitationRas,
                              rasterToMatch = sim$rasterToMatch,
                              maskWithRTM = TRUE,
                              method = "bilinear",
                              filename2 = FALSE,
                              userTags = c(cacheTags, "precipitationRas"))

    message(blue("Getting default 'relativeHumRas' to make default 'topoClimData'.",
                 "If this is not correct, make sure Biomass_fireProperties can detect 'topoClimData' is supplied"))
    ## get default precipitation values, summer average
    relativeHumRas <- Cache(prepInputs, targetFile = "RH_sm.asc",
                            url = extractURL("relativeHumRas", sim),
                            archive = "CanESM2_RCP45_r11i1p1_2011MSY.zip",
                            alsoExtract = NA,
                            destinationPath = dPath,
                            fun = "raster",
                            filename2 = FALSE,
                            userTags = cacheTags)

    ## add the original CRS if it's not defined
    if (is.na(crs(relativeHumRas)))
      crs(relativeHumRas) <- "+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"

    relativeHumRas <- Cache(postProcess,
                            x = relativeHumRas,
                            rasterToMatch = sim$rasterToMatch,
                            maskWithRTM = TRUE,
                            method = "bilinear",
                            filename2 = FALSE,
                            userTags = c(cacheTags, "relativeHumRas"))

    ## TODO defaults of slope/aspect should cover whole of Canada
    ## get default slope values
    slopeRas <- Cache(prepInputs, targetFile = "dataset/SLOPE.tif",
                      archive = "DEM_Foothills_study_area.zip",
                      alsoExtract = NA,
                      destinationPath = getPaths()$inputPath,
                      rasterToMatch = sim$rasterToMatch,
                      maskWithRTM = TRUE,
                      method = "bilinear",
                      datatype = "FLT4S",
                      filename2 = FALSE,
                      userTags = c(cacheTags, "slopeRas"))

    ## get default aspect values
    aspectRas <- Cache(prepInputs, targetFile = "dataset/ASPECT.tif",
                       archive = "DEM_Foothills_study_area.zip",
                       alsoExtract = NA,
                       destinationPath = getPaths()$inputPath,
                       rasterToMatch = sim$rasterToMatch,
                       maskWithRTM = TRUE,
                       method = "bilinear",
                       datatype = "FLT4S",
                       filename2 = FALSE,
                       userTags = c(cacheTags, "aspectRas"))

    ## PROJECT CLIMATE/TOPO RASTERS AS POINTS
    message(blue("Processing climate data for fire weather and fuel calculation"))
    slopePoints <- st_as_sf(rasterToPoints(slopeRas, spatial = TRUE))
    slopePoints <- st_transform(slopePoints, crs = crs(sim$studyAreaFBP))

    aspectPoints <- st_as_sf(rasterToPoints(aspectRas, spatial = TRUE))
    aspectPoints <- st_transform(aspectPoints, crs = crs(sim$studyAreaFBP))

    temperaturePoints <- st_as_sf(rasterToPoints(temperatureRas, spatial = TRUE))
    temperaturePoints <- st_transform(temperaturePoints, crs = crs(sim$studyAreaFBP))

    precipitationPoints <- st_as_sf(rasterToPoints(precipitationRas, spatial = TRUE))
    precipitationPoints <- st_transform(precipitationPoints, crs = crs(sim$studyAreaFBP))

    relativeHumPoints <- st_as_sf(rasterToPoints(relativeHumRas, spatial = TRUE))
    relativeHumPoints <- st_transform(relativeHumPoints, crs = crs(sim$studyAreaFBP))

    if (getOption("LandR.assertions")) {
      if (length(unique(c(crs(slopePoints), crs(aspectPoints), crs(sim$rasterToMatchFBPPoints)))) > 1)
        stop("Reprojecting topography data to FBP-compatible lat/long projection failed.\n
             Please inspect Biomass_fireProperties::firePropertiesInit")

      if (length(unique(c(nrow(slopePoints), nrow(aspectPoints), nrow(sim$rasterToMatchFBPPoints)))) > 1)
        stop("Topography data layers and rasterToMatch differ in number of points in FBP-compatible lat/long projection.\n
             Please inspect Biomass_fireProperties::firePropertiesInit")

      if (length(unique(c(crs(temperaturePoints), crs(precipitationPoints),
                          crs(relativeHumPoints), crs(sim$rasterToMatchFBPPoints)))) > 1)
        stop("Reprojecting climate data to FBP-compatible lat/long projection failed.\n
             Please inspect Biomass_fireProperties::firePropertiesInit")

      if (length(unique(c(nrow(temperaturePoints), nrow(precipitationPoints),
                          nrow(relativeHumPoints), nrow(sim$rasterToMatchFBPPoints)))) > 1)
        stop("Climate data layers and rasterToMatch differ in number of points in FBP-compatible lat/long projection.\n
             Please inspect Biomass_fireProperties::firePropertiesInit")
    }

    message(blue(currentModule(sim), " is making 'topoClimData' from default temperature, precipitation and relative humidity raster layers"))
    sim$topoClimDataCRS <- as.character(st_crs(sim$rasterToMatchFBPPoints))
    topoClimData <- data.table(temperature = temperaturePoints[, 1, drop = TRUE],
                               precipitation = precipitationPoints[, 1, drop = TRUE],
                               relativeHumidity = relativeHumPoints[, 1, drop = TRUE],
                               slope = slopePoints[, 1, drop = TRUE],
                               aspect = aspectPoints[, 1, drop = TRUE],
                               month = 7,  ## consider July.
                               day = 1,
                               windSpeed = 0, ## consider no wind
                               latitude = st_coordinates(sim$rasterToMatchFBPPoints)[,2],
                               longitude = st_coordinates(sim$rasterToMatchFBPPoints)[,1])
    topoClimData <- na.omit(topoClimData)
    topoClimData <- st_as_sf(topoClimData, coords = c("longitude", "latitude"),
                             crs = sim$topoClimDataCRS, agr = "constant")


    ## this is no longer necessary as ClimateNA has relative humidity data
    ## relative humidity
    ## using dew point between -3 and 20%, quarterly seasonal for Jun 2013
    ## https://calgary.weatherstats.ca/metrics/dew_point.html
    # topoClimData[, relHum := RH(t = topoClimData$temp, Td = runif (nrow(topoClimData), -3, 20), isK = FALSE)]

    ## export to sim
    sim$topoClimData <- topoClimData

  } else {
    if (!suppliedElsewhere("topoClimDataCRS", sim))
      stop(red("'topoClimData' appears to be supplied to Biomass_fireProperties,",
               "but not topoClimDataCRS. Please make sure 'topoClimDataCRS' is also provided."))
  }


  ## FWI INITIALISATION DATAFRAME
  ## TODO: FWIinit should be updated every year from previous year's/days/months results
  if (!suppliedElsewhere("FWIinit", sim)) {
    sim$FWIinit <- data.frame(ffmc = 85, dmc = 6, dc = 15)
  }

  if(!suppliedElsewhere("pixelNonForestFuels", sim)) {
    sim$pixelNonForestFuels <- NULL
  }

  return(invisible(sim))
}
