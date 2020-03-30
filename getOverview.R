# function to create an overview of the data in  the PALMOD 130k marine palaeoclimate data synthesis similar to Table 5 in
# Jonkers, L., Cartapanis, O., Langner, M., McKay, N., Mulitza, S., Strack, A., and Kucera, M.: Integrating palaeoclimate time series with rich metadata for uncertainty modelling: strategy and documentation of the PALMOD 130k marine palaeoclimate data synthesis, Earth Syst. Sci. Data Discuss., https://doi.org/10.5194/essd-2019-223, in review, 2020.
getOverview <- function(x){
  require(tidyverse)
  require(reshape2)
  extract <- function(site){
    out <- tryCatch({
      dat <- readRDS(site)
      # get indices of proxies and extract metadata
      proxies <- c('benthic.d18O', 'benthic.d13C', 'planktonic.d18O', 'planktonic.d13C', 'surface.temp', 'deep.temp', 'CaCO3', 'TOC', 'BSi')
      pIndx <- which(dat$meta$Parameter %in% proxies)
      subMeta <- dat$meta[pIndx, ]
      overview <- cbind.data.frame(Site = dat$site$SiteName,
                                   Longitude = dat$site$SiteLon,
                                   Latitude = dat$site$SiteLat,
                                   Elevation = dat$site$SiteDepth_m,
                                   Parameter = subMeta$Parameter,
                                   ParameterDetail = subMeta$Material,
                                   SensorSpecies = subMeta$Species,
                                   stringsAsFactors = FALSE)
      # get minimum, maximum age and resolution
      pData <- dat$data[, c(1, pIndx+1)]
      names(pData)[1] <- 'depth_m'
      subData <- inner_join(dat$AgeModel[,1:2], pData, by = 'depth_m')
      # ageInfo
      ageInfo <- if(ncol(subData)>3){
        sapply(subData[, 3:ncol(subData)], function(x){
          c(minAge_kyr = min(subData$meanAge_kaBP[!is.na(x)]),
          maxAge_kyr = max(subData$meanAge_kaBP[!is.na(x)]),
          medianResolution_kyr = median(diff(subData$meanAge_kaBP[!is.na(x)])))
        })} else {
          c(minAge_kyr = min(subData$meanAge_kaBP[!is.na(subData[,3])]),
          maxAge_kyr = max(subData$meanAge_kaBP[!is.na(subData[,3])]),
          medianResolution_kyr = median(diff(subData$meanAge_kaBP[!is.na(subData[,3])])))
        }
      
      agesUsed <- dat$BACON$AgeData$UseFlag
      agesUsed[grep('BACON_ON', dat$BACON$AgeData$Comment)] <- 1
      allTiepoints <- round(dat$BACON$AgeData$Depth_m[agesUsed == 1], 3)
      
      depthRanges <- data.frame(n = round(sapply(subData, function(x) range(subData$depth_m[!is.na(x)]))[, -c(1, 2)], 3))
      nTiepoints <- sapply(depthRanges, function(x) sum(allTiepoints >= x[1] & allTiepoints <= x[2]))
      
      data.frame(overview, t(ageInfo), nTiepoints, stringsAsFactors = FALSE)
    },
    error = function(cond){
      message(paste(site, 'gives error'))
      message(cond)
      return(NULL)
    },
    finally = message(paste('Processed', site))
    )
    return(out)
  }
  bind_rows(sapply(x, extract, simplify = FALSE))
}
