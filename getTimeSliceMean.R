# extract mean values for time slice
# use query first to determine which sites meet the criteria, this means that only sites that have at least one data point
# on the average age model within the time slice are included.
# I could consider changing this and search within the ensembles for time series that have at least n ages within the time slice
# or avoid estimates based on a few age ensembles only (as age distribution may be non-normal or other time series may have fewer data points
# and only a few ensembles within the time slice)
# have not done that yet to allow for more flexibility with querying


getTimeSliceMean <- function(x, files, Parameter, age_min, age_max, defaultErr = 0.1){
  require(tidyverse)

  # read data
  dat <- readRDS(files[x])
  
  # parameter index
  Pindx <- which(dat$meta$Parameter %in% Parameter)
  
  # data vs depth
  Pdat <- dat$data[ , c(1,Pindx)]
  
  # get uncertainty from measurements
  # use analytical reproducibility if available, if not
  # analytical uncertainty, and if that not, the assumed default error
  
  measU <- tibble(ParameterOriginal = dat$meta$ParameterOriginal[Pindx], 
                  measU = map_dbl(Pindx, function(y){
                    if('ParameterReproducibility' %in% names(dat$meta)){
                      if(!is.na(dat$meta$ParameterReproducibility[y])){
                        as.numeric(dat$meta$ParameterReproducibility[y])
                      }
                      else {
                        defaultErr
                      }
                    } 
                    else if ('ParameterAnalyticalError' %in% names(dat$meta)){
                      if(!is.na(dat$meta$ParameterAnalyticalError[y])){
                        as.numeric(dat$meta$ParameterAnalyticalError[y])
                      }
                      else {
                        defaultErr
                      }
                    }
                    else {
                      defaultErr
                    }
                  })
  ) %>%
    mutate(measU = case_when(is.na(measU) ~ defaultErr,
                             TRUE ~ measU))
  
  
  # time slice mean including uncertainty from age
  TSavg <- tryCatch({
    tibble(dat$AgeEnsembles) %>%
      rename(depthMerged_m = depth_m) %>%
      inner_join(Pdat, ., by = 'depthMerged_m') %>%
      pivot_longer(-(1:(length(Pindx)+1)), names_to = 'Ens', values_to = 'age') %>%
      filter(age >= age_min & age < age_max) %>%
      pivot_longer(dat$meta$ParameterOriginal[Pindx], names_to = 'ParameterOriginal', values_to = 'Value') %>%
      drop_na() %>%
      # group_by(ParameterOriginal) %>%
      # mutate(nEms = n_distinct(Ens)) %>%
      # filter(nEms >= minEns) %>%
      left_join(., measU, by = 'ParameterOriginal') %>%
      group_by(ParameterOriginal) %>%
      mutate(Noise = rnorm(n(), mean = 0, sd = measU),
             ValuePlusNoise = Value + Noise) %>%
      group_by(ParameterOriginal, Ens) %>%
      summarise(EnsAvgValue = mean(ValuePlusNoise), n = n(), .groups = 'drop') %>%
      group_by(ParameterOriginal) %>%
      summarise(TSValAvg = mean(EnsAvgValue), TSValSD = sd(EnsAvgValue), AvgN = mean(n)) %>%
      mutate(extractNote = 'OK')
  },
  error = function(cond){
    message(paste("File error:", x[1]))
    message(cond)
    return(tibble(ParameterOriginal = NA_character_,
                  TSValAvg = NA_real_,
                  TSValSD = NA_real_,
                  AvgN = NA_real_,
                  extractNote = 'oops, something went wrong'))
  }) %>%
    mutate(TSStdErr = TSValSD/sqrt(AvgN),
           age_min = age_min,
           age_max = age_max)

  # output
  tibble(dat$site) %>%
    mutate(SiteNotes = as.character(SiteNotes),
           version = dat$PALMODversion,
           extractDate = Sys.Date()) %>%
    bind_cols(.,
              left_join(TSavg,
                        tibble(dat$meta[Pindx, ] %>%
                                 select(-RetrievalNumber)),
                        by = 'ParameterOriginal'))
}
