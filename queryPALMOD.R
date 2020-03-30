# function to query the PALMOD 130k marine palaeoclimate data synthesis
queryPALMOD <- function(sites, Parameter, ParameterDetail = NULL, SensorSpecies = NULL, lat_min = -90, lat_max = 90, lon_min = -180, lon_max = 180, age_min = 0, age_max = Inf, min_parameter = 1, min_tiepoint = 0, createLists = TRUE){
  query <- function(x){
    out <- tryCatch({
      # check input
      if(!Parameter %in% allParameter){
        stop('parameter not included in database. Check spelling?')
      }
      if(lat_min >= lat_max | lon_min >= lon_max){
        stop('region not correctly defined.')
      }
      if(age_min < 0 | age_min > age_max){
        stop('check age range.')
      }
      if(!is.null(ParameterDetail)){
        if(!ParameterDetail %in% allParameterDetail){
          stop('parameter detail not in list.')
        }
      }
      if(!is.null(SensorSpecies)){
        if(!SensorSpecies %in% allSensorSpecies){
          stop('sensor species not in list.')
        }
      }
      # read data
      dat <- readRDS(x)
      # in regional range
      require(sp)
      inRegion <- point.in.polygon(dat$site$SiteLon, dat$site$SiteLat, c(lon_min, lon_max, lon_max, lon_min), c(lat_max, lat_max, lat_min, lat_min)) > 0
      # 
      if(inRegion & Parameter %in% dat$meta$Parameter & is.null(ParameterDetail) & is.null(SensorSpecies)){ # no ParameterDetail, no SensorSpecies
        # how many data points within age range
        pIndx <- which(dat$meta$Parameter %in% Parameter)
        pData <- lapply(pIndx, function(x){
          tem <- dat$data[, c(1, x+1)]
          names(tem)[1] <- 'depth_m'
          tem
        })
        subData <- lapply(pData, function(x) inner_join(dat$AgeModel[,1:2], x, by = 'depth_m'))
        subData <- lapply(subData, function(x) x[complete.cases(x),])
        nParameter <- max(sapply(subData, function(x) sum(x$meanAge_kaBP >= age_min & x$meanAge_kaBP <= age_max)))
        # how many tiepoints within age range
        agesUsed <- dat$BACON$AgeData$UseFlag
        agesUsed[grep('BACON_ON', dat$BACON$AgeData$Comment)] <- 1
        allTiepoints <- round(dat$BACON$AgeData$CalendarAge_kaBP[agesUsed == 1], 3)
        nTiepoints <- sum(allTiepoints >= age_min & allTiepoints <= age_max)
        nParameter >= min_parameter & nTiepoints >= min_tiepoint
      } else if(inRegion & Parameter %in% dat$meta$Parameter & !is.null(ParameterDetail) & is.null(SensorSpecies)){ # no SensorSpecies
        # how many data points within age range
        pIndx <- which(dat$meta$Parameter %in% Parameter & dat$meta$Material %in% ParameterDetail)
        if(length(pIndx)> 0){
          pData <- lapply(pIndx, function(x){
            tem <- dat$data[, c(1, x+1)]
            names(tem)[1] <- 'depth_m'
            tem
          })
          subData <- lapply(pData, function(x) inner_join(dat$AgeModel[,1:2], x, by = 'depth_m'))
          subData <- lapply(subData, function(x) x[complete.cases(x),])
          nParameter <- max(sapply(subData, function(x) sum(x$meanAge_kaBP >= age_min & x$meanAge_kaBP <= age_max)))
          # how many tiepoints within age range
          agesUsed <- dat$BACON$AgeData$UseFlag
          agesUsed[grep('BACON_ON', dat$BACON$AgeData$Comment)] <- 1
          allTiepoints <- round(dat$BACON$AgeData$CalendarAge_kaBP[agesUsed == 1], 3)
          nTiepoints <- sum(allTiepoints >= age_min & allTiepoints <= age_max)
          nParameter >= min_parameter & nTiepoints >= min_tiepoint
        } else {
          FALSE
        }
      } else if(inRegion & Parameter %in% dat$meta$Parameter & is.null(ParameterDetail) & !is.null(SensorSpecies)){ # no ParameterDetail
        # how many data points within age range
        pIndx <- which(dat$meta$Parameter %in% Parameter & dat$meta$Species %in% SensorSpecies)
        if(length(pIndx)> 0){
          pData <- lapply(pIndx, function(x){
            tem <- dat$data[, c(1, x+1)]
            names(tem)[1] <- 'depth_m'
            tem
          })
          subData <- lapply(pData, function(x) inner_join(dat$AgeModel[,1:2], x, by = 'depth_m'))
          subData <- lapply(subData, function(x) x[complete.cases(x),])
          nParameter <- max(sapply(subData, function(x) sum(x$meanAge_kaBP >= age_min & x$meanAge_kaBP <= age_max)))
          # how many tiepoints within age range
          agesUsed <- dat$BACON$AgeData$UseFlag
          agesUsed[grep('BACON_ON', dat$BACON$AgeData$Comment)] <- 1
          allTiepoints <- round(dat$BACON$AgeData$CalendarAge_kaBP[agesUsed == 1], 3)
          nTiepoints <- sum(allTiepoints >= age_min & allTiepoints <= age_max)
          nParameter >= min_parameter & nTiepoints >= min_tiepoint
        } else {
          FALSE
        } 
      } else if(inRegion & Parameter %in% dat$meta$Parameter & !is.null(ParameterDetail) & !is.null(SensorSpecies)){ # all specified
        # how many data points within age range
        pIndx <- which(dat$meta$Parameter %in% Parameter & dat$meta$Material %in% ParameterDetail & dat$meta$Species %in% SensorSpecies)
        if(length(pIndx)> 0){
          pData <- lapply(pIndx, function(x){
            tem <- dat$data[, c(1, x+1)]
            names(tem)[1] <- 'depth_m'
            tem
          })
          subData <- lapply(pData, function(x) inner_join(dat$AgeModel[,1:2], x, by = 'depth_m'))
          subData <- lapply(subData, function(x) x[complete.cases(x),])
          nParameter <- max(sapply(subData, function(x) sum(x$meanAge_kaBP >= age_min & x$meanAge_kaBP <= age_max)))
          # how many tiepoints within age range
          agesUsed <- dat$BACON$AgeData$UseFlag
          agesUsed[grep('BACON_ON', dat$BACON$AgeData$Comment)] <- 1
          allTiepoints <- round(dat$BACON$AgeData$CalendarAge_kaBP[agesUsed == 1], 3)
          nTiepoints <- sum(allTiepoints >= age_min & allTiepoints <= age_max)
          nParameter >= min_parameter & nTiepoints >= min_tiepoint
        } else {
          FALSE
        }
      } else {
        FALSE
      }
    },
    error = function(cond){
      message(cond)
      return(NULL)
    })
    return(out)
  }
  # create these first before filtering
  if(createLists){
    allParameter <- sort(unique(unlist(sapply(sites, function(x) readRDS(x)$meta$Parameter))))
    allParameterDetail <- sort(unique(unlist(sapply(sites, function(x) readRDS(x)$meta$Material))))
    allSensorSpecies <- sort(unique(unlist(sapply(sites, function(x) readRDS(x)$meta$Species))))
  }
  # apply filter
  which(sapply(sites, query))
}
