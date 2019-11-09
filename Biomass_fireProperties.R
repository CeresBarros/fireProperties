# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects and functions, use sim$xxx.
defineModule(sim, list(
  name = "Biomass_fireProperties",
  description = "Complement to fire spread that calculates fire (behaviour) properties in fucntion of vegetation, climate and topography conditions
  usign the Canadian Forest Fire Behaviour Prediction System",
  keywords = c("fire behaviour", "fuels", "fire-vegetation feedbacks", "fire-climate feedbacks", "FBP system", "topography"),
  authors = person("Ceres", "Barros", email = "cbarros@mail.ubc.ca", role = c("aut", "cre")),
  childModules = character(0),
  version = numeric_version("0.1.0"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "Biomass_fireProperties.Rmd"),
  reqdPkgs = list("R.utils", "raster", "data.table", "dplyr",
                  "sf", "cffdrs",
                  "PredictiveEcology/SpaDES.core@development",
                  "PredictiveEcology/SpaDES.tools@development",
                  "PredictiveEcology/reproducible@development"),
  parameters = rbind(
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
    expectsInput(objectName = "aspectRas", objectClass = "RasterLayer",
                 desc = "Raster of aspect values - needs to be previously downloaded at this point"),
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
    expectsInput(objectName = "precipitationRas", objectClass = "RasterLayer",
                 desc = paste0("Raster of summer average precipitation values.",
                               "Defaults data downloaded from Climate NA for 2011 using: CanESM2_RCP45_r11i1p1_2011MSY"),
                 sourceURL = "https://drive.google.com/open?id=12iNnl3P7VjisVKC0vatSrXyhYtl6w-D1"),
    expectsInput(objectName = "rasterToMatchFBP", objectClass = "RasterLayer",
                  desc = "a rasterToMatch reprojected to FBP-compatible projection"),
    expectsInput(objectName = "relativeHumRas", objectClass = "RasterLayer",
                 desc = paste0("Raster of summer average relative humidity values.",
                               "Defaults data downloaded from Climate NA for 2011 using: CanESM2_RCP45_r11i1p1_2011MSY"),
                 sourceURL = "https://drive.google.com/open?id=12iNnl3P7VjisVKC0vatSrXyhYtl6w-D1"),
    expectsInput(objectName = "simulatedBiomassMap", objectClass = "RasterLayer",
                 desc = "Biomass map at each succession time step. Default is Canada national biomass map",
                 sourceURL = "http://tree.pfc.forestry.ca/kNN-StructureBiomass.tar"),
    expectsInput(objectName = "slopeRas", objectClass = "RasterLayer",
                 desc = "Raster of slope values - needs to be previously downloaded at this point"),
    expectsInput("studyArea", "SpatialPolygonsDataFrame",
                 desc = paste("Polygon to use as the study area.",
                              "Defaults to  an area in Southwestern Alberta, Canada."),
                 sourceURL = ""),
    expectsInput(objectName = "studyAreaFBP", objectClass = "SpatialPolygonsDataFrame",
                 desc = "same as studyArea,  but on FBP-compatible projection", sourceURL = ""),
    expectsInput(objectName = "temperatureRas", objectClass = "RasterLayer",
                 desc = paste0("Raster of summer average temperature values.",
                               "Defaults data downloaded from Climate NA for 2011 using: CanESM2_RCP45_r11i1p1_2011MSY"),
                 sourceURL = "https://drive.google.com/open?id=12iNnl3P7VjisVKC0vatSrXyhYtl6w-D1"),
    expectsInput(objectName = "topoClimData", objectClass = "data.table",
                 desc = "Climate data table with temperature, precipitation and relative humidity for each pixelGroup")
  ),
  outputObjects = bind_rows(
    # createsOutput(objectName = "FBPinputs", objectClass = "RasterLayer",
    #               desc = "Fire behaviour prediction system inputs table"),
    # createsOutput(objectName = "FBPoutputs", objectClass = "list",
    #               desc = "Fire weather outputs table"),
    createsOutput(objectName = "aspectRas", objectClass = "RasterLayer",
                  desc = "Raster of aspect values - reprojected/cropped"),
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
    # createsOutput(objectName = "FWIinputs", objectClass = "RasterLayer",
    #               desc = "Fire weather inputs table"),
    # createsOutput(objectName = "FWIoutputs", objectClass = "list",
    #               desc = "Fire weather outputs table"),
    createsOutput(objectName = "precipitationRas", objectClass = "RasterLayer",
                  desc = "Raster of precipitation values - reprojected/cropped"),
    createsOutput(objectName = "rasterToMatchFBP", objectClass = "RasterLayer",
                  desc = "a rasterToMatch reprojected to FBP-compatible projection"),
    createsOutput(objectName = "slopeRas", objectClass = "RasterLayer",
                  desc = "Raster of slope values - reprojected/cropped"),
    createsOutput(objectName = "temperatureRas", objectClass = "RasterLayer",
                  desc = "Raster of temperature values - reprojected/cropped"),
    createsOutput(objectName = "topoClimData", objectClass = "data.table",
                  desc = "Climate data table with temperature, precipitation and relative humidity for each pixelGroup")
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

  ## project all inputs to Lat/Long (decimal degrees)
  ## for compatibility with FBP system

  ## buffer rasterToMatch to prevent data loss, using twice the pixel size and distance
  ## then reproject to FBP compatible projection and mask buffer out (later)
  ## note: don't mask to area until the end.
  rasterToMatchFBP <- buffer(sim$rasterToMatch,
                             width = res(sim$rasterToMatch)[1]*2)
  rasterToMatchFBP <- projectRaster(rasterToMatchFBP, method = "ngb",
                                    crs = crs(sim$studyAreaFBP))

  ## PROJECT CLIMATE/TOPO RASTERS
  message(blue("Processing climate data for fire weather and fuel calculation"))
  sim$temperatureRas <- Cache(postProcess,
                              x = sim$temperatureRas,
                              rasterToMatch = rasterToMatchFBP,
                              maskWithRTM = TRUE,
                              method = "bilinear",
                              filename2 = NULL,
                              userTags = c(cacheTags, "topoClimRas"),
                              omitArgs = c("userTags"))
  sim$precipitationRas <- Cache(postProcess,
                                x = sim$precipitationRas,
                                rasterToMatch = rasterToMatchFBP,
                                maskWithRTM = TRUE,
                                method = "bilinear",
                                filename2 = NULL,
                                userTags = c(cacheTags, "topoClimRas"),
                                omitArgs = c("userTags"))
  sim$relativeHumRas <- Cache(postProcess,
                              x = sim$relativeHumRas,
                              rasterToMatch = rasterToMatchFBP,
                              maskWithRTM = TRUE,
                              method = "bilinear",
                              filename2 = NULL,
                              userTags = c(cacheTags, "topoClimRas"),
                              omitArgs = c("userTags"))
  sim$slopeRas <- Cache(postProcess,
                        x = sim$slopeRas,
                        rasterToMatch = rasterToMatchFBP,
                        maskWithRTM = TRUE,
                        method = "bilinear",
                        filename2 = NULL,
                        userTags = c(cacheTags, "topoClimRas"),
                        omitArgs = c("userTags"))
  sim$aspectRas <- Cache(postProcess,
                         x = sim$aspectRas,
                         rasterToMatch = rasterToMatchFBP,
                         maskWithRTM = TRUE,
                         method = "bilinear",
                         filename2 = NULL,
                         userTags = c(cacheTags, "topoClimRas"),
                         omitArgs = c("userTags"))

  ## TOPOCLIMDATA TABLE ----------------------
  ## TODO: change to draw from fire weather distributions
  topoClimData <- data.table(ID = seq_len(ncell(rasterToMatchFBP)),
                             pixelGroup = getValues(rasterToMatchFBP),
                             temp = getValues(sim$temperatureRas), precip = getValues(sim$precipitationRas),
                             relHum = getValues(sim$relativeHumRas),
                             slope = getValues(sim$slopeRas), aspect = getValues(sim$aspectRas),
                             lat = coordinates(rasterToMatchFBP)[,2],
                             long = coordinates(rasterToMatchFBP)[,1])

  ## this is no longer necessary as ClimateNA has relative humidity data
  ## relative humidity
  ## using dew point between -3 and 20%, quarterly seasonal for Jun 2013
  ## https://calgary.weatherstats.ca/metrics/dew_point.html
  # topoClimData[, relHum := RH(t = topoClimData$temp, Td = runif (nrow(topoClimData), -3, 20), isK = FALSE)]

  ## export to sim
  sim$rasterToMatchFBP <- rasterToMatchFBP
  sim$topoClimData <- topoClimData

  return(invisible(sim))
}

## Derive fire parameters from FBP system - rasters need to be in lat/long
calcFBPProperties <- function(sim) {
  cacheTags <- c(currentModule(sim), "FBPPercParams")

  ## FUEL TYPES ------------------------------
  ## reproject to fuelTypesMaps to FBP-compatible crs
  fuelTypeRas <- Cache(postProcess,
                       x = sim$fuelTypesMaps$finalFuelType,
                       rasterToMatch = sim$rasterToMatchFBP,
                       maskWithRTM = TRUE,
                       method = "ngb",
                       filename2 = NULL,
                       userTags = c(cacheTags, "fuelTypeRas"),
                       omitArgs = c("userTags"))

  coniferDomRas <- Cache(postProcess,
                         x = sim$fuelTypesMaps$coniferDom,
                         rasterToMatch = sim$rasterToMatchFBP,
                         maskWithRTM = TRUE,
                         method = "bilinear",
                         filename2 = NULL,
                         userTags = c(cacheTags, "coniferDomRas"),
                         omitArgs = c("userTags"))

  if (!is.null(sim$pixelNonForestFuels)) {
    curingRas <- Cache(postProcess,
                       x = sim$fuelTypesMaps$curing,
                       rasterToMatch = sim$rasterToMatchFBP,
                       maskWithRTM = TRUE,
                       method = "bilinear",
                       filename2 = NULL,
                       userTags = c(cacheTags, "curingRas"),
                       omitArgs = c("userTags"))
  } else {
    ## post-process fails with a NA raster
    curingRas <- cropInputs(sim$fuelTypesMaps$curing, rasterToMatch = sim$rasterToMatchFBP)
    curingRas <- projectInputs(curingRas, rasterToMatch = sim$rasterToMatchFBP)
    curingRas <- maskInputs(curingRas, rasterToMatch = sim$rasterToMatchFB0, maskWithRTM = TRUE)

    if (getOption("LandR.assertions"))
      if (!compareRaster(curingRas, sim$rasterToMatchFBP))
        stop("Can't reproject NA curing raster to rasterToMatchFBP",
             " Please debug Biomass_fireProperties calcFBPProperties() event function")
  }


  ## make table of final fuel types
  FTs <- data.table(ID = seq_len(ncell(sim$rasterToMatchFBP)),
                    FuelType = getValues(fuelTypeRas),
                    coniferDom = getValues(coniferDomRas),
                    curing = getValues(curingRas))

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

  ## check for duplicates (there shouldn't be any)
  if (getOption("LandR.assertions"))
    if (any(duplicated(FTs))) {
      stop("Duplicated pixels found in fuel types table.",
           " Please debug Biomass_fireProperties calcFBPProperties() event function")
  }

  ## FWI ------------------------------
  ## make/update table of FWI inputs
  FWIinputs <- data.frame(id = sim$topoClimData$ID,
                          lat = sim$topoClimData$lat,
                          long = sim$topoClimData$long,
                          mon = 7,  ## consider July.
                          temp = sim$topoClimData$temp,
                          rh = sim$topoClimData$relHum,
                          ws = 0,
                          prec = sim$topoClimData$precip)

  ## calculate FW indices
  FWIoutputs <- suppressWarnings({
    cffdrs::fwi(input = na.omit(FWIinputs),
                init = na.omit(sim$FWIinit),
                batch = FALSE, lat.adjust = TRUE)
  })
  FWIoutputs <- data.table(FWIoutputs)

  ## FBP -----------------------------
  ## make inputs dataframe for FBI
  ## add fuel types and conifer dominance to FWIOutputs
  ## note that because climate/topo data is "larger" there are pixels that have no fuels - these are removed.
  FWIoutputs <- FTs[FWIoutputs, on = "ID", nomatch = 0,
                    allow.cartesian = TRUE]

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
    cffdrs::fbp(input = na.omit(FBPinputs), output = "All")
  })
  FBPoutputs <- data.table(FBPoutputs)

  ## FBP OUTPUTS TO SPATIALPOINTS
  FBPOutputsPts <- FBPoutputs[, .(ID, CFB, ROS, RSO, HFI, TFC)]
  FBPOutputsPts <- FBPOutputsPts[FWIoutputs[, .(ID, LAT, LONG)],
                                 on = "ID", nomatch = 0][, ID := NULL]

  FBPOutputsSf <- st_as_sf(FBPOutputsPts, coords = c("LONG", "LAT"),
                           crs =  as.character(crs(sim$studyAreaFBP)),
                           agr = "constant")

  ## REPROJECT TO ORIGINAL CRS without data loss and convert to raster
  ## note that after rasterizing it is safer to mask, in case some of
  ## the buffer artefacts come through
  FBPOutputsSf <- st_transform(FBPOutputsSf,
                               crs = as.character(crs(sim$rasterToMatch)))
  FBPOutputsPoly <- as_Spatial(FBPOutputsSf)

  ## Crown fraction burnt
  sim$fireCFBRas <- rasterize(FBPOutputsPoly, sim$rasterToMatch,
                              field = "CFB", fun = function(x, na.rm = TRUE) max(x))
  sim$fireCFBRas <- mask(sim$fireCFBRas, sim$rasterToMatch)

  ## Head fire intensity
  sim$fireIntRas <- rasterize(FBPOutputsPoly, sim$rasterToMatch,
                              field = "HFI", fun = function(x, na.rm = TRUE) max(x))
  sim$fireIntRas <- mask(sim$fireIntRas, sim$rasterToMatch)

  ## Rate of spread
  sim$fireROSRas <- rasterize(FBPOutputsPoly, sim$rasterToMatch,
                              field = "ROS", fun = function(x, na.rm = TRUE) max(x))
  sim$fireROSRas <- mask(sim$fireROSRas, sim$rasterToMatch)

  ## Critical spread rate for crowning
  sim$fireRSORas <- rasterize(FBPOutputsPoly, sim$rasterToMatch,
                              field = "RSO", fun = function(x, na.rm = TRUE) max(x))
  sim$fireRSORas <- mask(sim$fireRSORas, sim$rasterToMatch)

  ## Total fuel consumption
  sim$fireTFCRas <- rasterize(FBPOutputsPoly, sim$rasterToMatch,
                              field = "TFC", fun = function(x, na.rm = TRUE) max(x))
  sim$fireTFCRas <- mask(sim$fireTFCRas, sim$rasterToMatch)

  ## export to sim
  # sim$FWIinputs <- FWIinputs
  # sim$FWIoutputs <- FWIoutputs
  # sim$FBPinputs <- FBPinputs
  # sim$FBPoutputs <- FBPoutputs

  return(invisible(sim))
}

## OTHER INPUTS AND FUNCTIONS --------------------------------
.inputObjects = function(sim) {
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  cacheTags = c(currentModule(sim), "function:.inputObjects")

  ## project to Lat/Long (decimal degrees) for compatibility with FBP system
  ## TODO: this results in data loss - but LandR doesn't deal well with lat/long
  ## need to find long term solution
  latLong <- "+proj=longlat +datum=WGS84"

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

  ## DEFAULT TOPO, TEMPERATURE AND PRECIPITATION
  ## these defaults are only necessary if the rasters are not supplied by another module
  ## climate defaults to Climate NA Data, year 2011, RCP4.5
  ## note that some Climate NA data were multiplied by 10
  if (!suppliedElsewhere("temperatureRas", sim)) {
    ## get default temperature values, summer average
    sim$temperatureRas <- Cache(prepInputs, targetFile = "Tave_sm.asc",
                                url = extractURL("temperatureRas", sim),
                                archive = "CanESM2_RCP45_r11i1p1_2011MSY.zip",
                                alsoExtract = NA,
                                destinationPath = dPath,
                                fun = "raster",
                                filename2 = FALSE,
                                userTags = cacheTags)
    sim$temperatureRas <- sim$temperatureRas/10  ## back transform temp values

    ## add the original CRS if it's not defined
    if (is.na(crs(sim$temperatureRas)))
      crs(sim$temperatureRas) <- "+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
  }

  if (!suppliedElsewhere("precipitationRas", sim)) {
    ## get default precipitation values, summer cummulative precipitation
    sim$precipitationRas <- Cache(prepInputs, targetFile = "PPT_sm.asc",
                                  url = extractURL("precipitationRas", sim),
                                  archive = "CanESM2_RCP45_r11i1p1_2011MSY.zip",
                                  alsoExtract = NA,
                                  destinationPath = dPath,
                                  fun = "raster",
                                  filename2 = FALSE,
                                  userTags = cacheTags)

    ## add the original CRS if it's not defined
    if (is.na(crs(sim$precipitationRas)))
      crs(sim$precipitationRas) <- "+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
  }

  if (!suppliedElsewhere("relativeHumRas", sim)) {
    ## get default precipitation values, summer average
    sim$relativeHumRas <- Cache(prepInputs, targetFile = "RH_sm.asc",
                                url = extractURL("relativeHumRas", sim),
                                archive = "CanESM2_RCP45_r11i1p1_2011MSY.zip",
                                alsoExtract = NA,
                                destinationPath = dPath,
                                fun = "raster",
                                filename2 = FALSE,
                                userTags = cacheTags)
    if (is.na(crs(sim$relativeHumRas)))
      crs(sim$relativeHumRas) <- "+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
  }

  if (!suppliedElsewhere("slopeRas", sim)) {
    ## TODO defaults of slope/aspect should cover whole of Canada
    ## get default slope values
    sim$slopeRas <- Cache(prepInputs, targetFile = "dataset/SLOPE.tif",
                          archive = "DEM_Foothills_study_area.zip",
                          alsoExtract = NA,
                          destinationPath = getPaths()$inputPath,
                          datatype = "FLT4S",
                          filename2 = FALSE,
                          userTags = cacheTags)
  }

  if (!suppliedElsewhere("aspectRas", sim)) {
    ## get default aspect values
    sim$aspectRas <- Cache(prepInputs, targetFile = "dataset/ASPECT.tif",
                           archive = "DEM_Foothills_study_area.zip",
                           alsoExtract = NA,
                           destinationPath = getPaths()$inputPath,
                           datatype = "FLT4S",
                           filename2 = FALSE,
                           userTags = cacheTags)
  }

  ## FWI INITIALISATION DATAFRAME
  ## TODO:FWIinit should be updated every year from previous year's/days/months results
  if (!suppliedElsewhere("FWIinit", sim)) {
    sim$FWIinit = data.frame(ffmc = 85,
                             dmc = 6,
                             dc = 15)
  }

  if(!suppliedElsewhere("pixelNonForestFuels", sim)) {
    sim$pixelNonForestFuels <- NULL
  }

  return(invisible(sim))
}
