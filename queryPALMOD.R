# function to query the PALMOD 130k marine palaeoclimate data synthesis
# extracting the age and/or sst ensembles requires custom code
# the option to filter by the minimum number of tie points could 
#  be replaced with an age uncertainty criterion

queryPALMOD <- function(sites,
                        Parameter, 
                        ParameterDetail = NULL,
                        SensorSpecies = NULL, 
                        lat_min = -90, 
                        lat_max = 90, 
                        lon_min = -180,
                        lon_max = 180,
                        ocean_basin = NULL,
                        age_min = 0,
                        age_max = Inf, 
                        min_duration = 0, # in kyr
                        min_n_obs = 1, # evaluated within min duration window
                        max_sample_interval = Inf, # in kyr evaluated within min duration window
                        min_n_tiepoint = 0,
                        extract_data = FALSE
                        ){
  require(tidyverse)
  query <- function(x){
    out <- tryCatch({
      # check input
      parameters <- c(
        "age", "Age_kaBP_2.5percentile", "Age_kaBP_97.5percentile",
        "benthic.d13C", "benthic.d13C.error", "benthic.d18O",
        "benthic.d18O.error", "benthic.MgCa", "BIT",
        "BSi", "C37.concentration", "CaCO3", "DBD", 
        "deep.temp", "depth", "depth_m", "depth_merged", "depth.bottom",
        "depth.composite", "depth.top", "IRD", "label",
        "LDI", "meanAge_kaBP", "medianAge_kaBP", "note",                      
        "notes", "planktonic.d13C", "planktonic.d13C.error",
        "planktonic.d18O", "planktonic.d18O.error",
        "planktonic.foram.abundance", "planktonic.MgCa",
        "planktonic.MgCa.error", "seaice.concentration","surface.temp",
        "surface.temp.error","surface.temp.hi", "surface.temp.lo",
        "TEX86", "TOC", "TOC.error", "UK37"
      )
      
      if(!Parameter %in% parameters){
        stop(paste("Parameter must be one of:", paste(parameters, collapse = ", ")))
      }
      
      if(lat_min >= lat_max | lon_min >= lon_max){
        stop('region not correctly defined. Check latitude and longitude')
      }
      
      oceans <- c("Southern Ocean",
                  "South Atlantic Ocean",
                  "South Pacific Ocean",
                  "North Pacific Ocean",
                  "South China and Easter Archipelagic Seas",
                  "Indian Ocean",
                  "Mediterranean Region",
                  "Baltic Sea",
                  "North Atlantic Ocean",
                  "Arctic Ocean")
      
      if(!is.null(ocean_basin)){
        if(!all(ocean_basin %in% oceans)
           ){
          stop(paste("ocean_basin must be one of:", paste(oceans, collapse = ", ")))
        }
      }
      
      if(age_min < -0.1 | age_min > age_max){
        stop('check age range.')
      }
      
      details <- c("aquatic palynomorphs assemblages, benthic foraminifera calcite, bulk sediment, CaCO3 and TOC free 2-63 fraction, coarse fraction, diatom assemblages, dinocyst assemblages, fine fraction, LDI, MgCa, planktonic foraminifera, planktonic foraminifera assemblages, planktonic foraminifera calcite, pteropod aragonite, radiolaria assemblages, TEX86, TOC, UK37, unknown")
      
      if(!is.null(ParameterDetail)){
        if(!ParameterDetail %in% details){
          stop('parameter detail not in list.')
        }
      }
      
      species <- c("A. angulosa, B. aculeata, B. marginata, C. acicula, C. kullenbergi, C. lobatulus, C. mabahethi, C. mckannai, C. mundulus, C. pachyderma, C. pseudoungeriana, C. robertsonianus, C. teretis, C. wuellerstorfi, Cibicides sp, Cibicides spp, Cibicidoides sp, Cibicidoides spp, D. anfracta, E. bartletti, E. batialis, E. clavatum, E. excavatum, E. excavatum forma clavata, E. exigua, E. groenlandicum, G. aequilateralis, G. affinis, G. bradyi, G. bulloides, G. calida, G. cassaformis, G. conglobatus, G. crassaformis, G. digitata, G. falconensis, G. glutinata, G. hexagonus, G. hirsuta, G. humilis, g. inflata, G. inflata, G. menardii, G. menardii (including G. flexuosa and tumida), G. pacifica, G. rubenscens, G. ruber, G. ruber pink, G. ruber white, G. ruber white s.l., G. ruber white s.s., G. rubescens, G. sacculifer no sac, G. sacculifer sac, G. scitula, G. siphonifera, G. tenella, G. truncatulinoides, G. truncatulinoides d, G. truncatulinoides dextral, G. truncatulinoides s, G. truncatulinoides sinistral, G. tumida, Globigerinita sp., Globobulimina spp, H. balthica, H. elegans, H. orbiculare, L. inflata, M. affinis, M. barleanum, M. barleeanus, M. zaandami, Melonis sp, mixed, Mixed, N. dutertrei, N. incompta, N. incompta and N. dutertrei integrade, N. labradorica, N. pachyderma, N. pachyderma + N. dutertrei, N. pachyderma sinistral, N. umbonifera, O. umbonatus, O. universa, P. ariminensis, P. obliquiloculata, PDintergrade, Pyrgo sp, T_quinqueloba, T. iota, T. quinqueloba, T. quinquiloba, T. sacculifer, U. akita, U. auberiana, U. bifurcata, U. costata, U. excellens, U. hispida, U. hollicki, U. mediterranea, U. peregrina, U. senticosa, unknown, Uvigerina sp, Uvigerina spp")
      
      if(!is.null(SensorSpecies)){
        if(!SensorSpecies %in% species){
          stop('sensor species not available')
        }
      }
      if(min_duration > age_max - age_min){
        stop('min duration longer than time window')
      }
      # read data
      dat <- readRDS(x)
      
      # in regional range
      inBasin <- if(is.null(ocean_basin)){
        TRUE
      } else {
        dat$site$OceanBasin %in% ocean_basin
      }
      
      require(sp)
      inPolygon <- point.in.polygon(dat$site$SiteLon, dat$site$SiteLat, c(lon_min, lon_max, lon_max, lon_min), c(lat_max, lat_max, lat_min, lat_min)) > 0
      
      inRegion <- all(inBasin, inPolygon)
      
      # enough tiepoints?
      enoughTiepoint <- dat$chron %>%
        filter(
          ChronAgeRejected == FALSE,
          ChronCalendarAge_kaBP >= age_min,
          ChronCalendarAge_kaBP <= age_max
          ) %>%
        reframe(enoughTiepoint = n()>= min_n_tiepoint) %>%
        unlist(use.names = FALSE)
  
      aIndx <- grep("medianAge_kaBP", dat$meta$Parameter)
      dIndx <- grep("depth_merged", dat$meta$Parameter)
      
      test <- if(inRegion & Parameter %in% dat$meta$Parameter & enoughTiepoint){
        pIndx <- if(is.null(ParameterDetail) & is.null(SensorSpecies)){
          which(dat$meta$Parameter %in% Parameter)
        } else if(!is.null(ParameterDetail) & is.null(SensorSpecies)){
          which(dat$meta$Parameter %in% Parameter & dat$meta$Material %in% ParameterDetail)
        } else if (is.null(ParameterDetail) & !is.null(SensorSpecies)){
          which(dat$meta$Parameter %in% Parameter & dat$meta$Species %in% SensorSpecies)
        } else if (!is.null(ParameterDetail) & !is.null(SensorSpecies)){
          which(dat$meta$Parameter %in% Parameter & dat$meta$Material %in% ParameterDetail & dat$meta$Species %in% SensorSpecies)
        }
        
        subData <- map(pIndx, function(x){
          dat$data[, c(aIndx, dIndx, x)] %>%
            drop_na() %>%
            filter(medianAge_kaBP >= age_min & medianAge_kaBP <= age_max)
        })
        
        subIndx <- map_lgl(subData, ~nrow(.) > 1)
        
        pIndx <- pIndx[subIndx]
        
        if(any(subIndx)){
          subData[subIndx] %>%
            map_lgl(function(x){
              max(x$medianAge_kaBP) - min(x$medianAge_kaBP) >= min_duration &
                max(lag(x$medianAge_kaBP), na.rm = TRUE) < max_sample_interval &
                nrow(x) >= min_n_obs
            }) %>%
            any()
        } else {
          FALSE
        }
      } else {
        FALSE
      }
      
      if(extract_data & test){
        subData <- map2(subData[subIndx], dat$meta$Parameter[pIndx], function(x, y){
          names(x)[3] <- y
          x
        })
        tibble(
          dat$site,
          dat$meta[pIndx,],
          data = subData
        )
      } else if (extract_data & !test){
        NULL
      } else {
        test
      }

    },
    error = function(cond){
      message(cond)
      return(NULL)
    })
    return(out)
  }
  # apply query function
  if(extract_data){
    map_df(sites, query) %>%
      select(where(~sum(!is.na(.)) > 0))
  } else {
    map_lgl(sites, query)
  }
}
