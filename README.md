## PALMODutils
example R code to interact with the PALMOD 130k marine palaeoclimate data synthesis 
The contents of the data product are described in detail in https://doi.org/10.5194/essd-2019-223. And the data are available for download at https://doi.pangaea.de/10.1594/PANGAEA.908831.

The function getOverview creates an overview of the contents of the data base, similar to Table 5 in https://doi.org/10.5194/essd-2019-223.

The function queryPALMOD allows querying of the data by:
* geographical range
* parameter
* parameter detail
* sensor species
* age range
* minimum number of data points
* minimum number of age control points

The function searches with calendar age. A complete parameter list can be found in Table 4 of https://doi.org/10.5194/essd-2019-223. Specifying parameter detail makes most sense when filtering near surface temperature estimates (using for instance: MgCa, UK37 or planktonic foraminifera assemblage).

Both getOverview and queryPALMOD  should be applied over a vector that contains the location of the individual files, which can be created with:

`PALMOD130k_V1_0_1 <- list.files('PALMOD_RDS', full.names = TRUE)`

### create the overview
```
source('getOverview.R')
overview <- getOverview(PALMOD130k_V1_0_1)
```

### query the data
Note that this function return the indices of sites that meet the criteria, which can then be used the extract the exact data you need.
The query function creates a list of parameters, parameter details and species. To speed things up one can create these once and set createLists to FALSE.
```
allParameter <- sort(unique(unlist(sapply(PALMOD130k_V1_0_1, function(x) readRDS(x)$meta$Parameter))))
allParameterDetail <- sort(unique(unlist(sapply(PALMOD130k_V1_0_1, function(x) readRDS(x)$meta$Material))))
allSensorSpecies <- sort(unique(unlist(sapply(PALMOD130k_V1_0_1, function(x) readRDS(x)$meta$Species))))

source('queryPALMOD.R')
meetsCrit <- queryPALMOD(sites = PALMOD130k_V1_0_1,
                          Parameter = 'benthic.d18O',
                          ParameterDetail = NULL,
                          SensorSpecies = NULL,
                          lat_min = 0,
                          lat_max = 90,
                          lon_min = -180,
                          lon_max = 180,
                          age_min = 0,
                          age_max = 130,
                          min_parameter = 100,
                          min_tiepoint = 10,
                          createLists = FALSE)
```
