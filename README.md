## PALMODutils
example R code to interact with the PALMOD 130k marine palaeoclimate data synthesis version 2 (to be released soon). 
The contents of version 1 of the data product are described in detail in https://doi.org/10.5194/essd-2019-223. Version 1 data are available for download at https://doi.pangaea.de/10.1594/PANGAEA.908831. The links to the data and the data descriptor for version 2 will be mentioned here in due time. This code cannot be used to anymore to query previous versions of the data product (use the initial commit instead).

The function `queryPALMOD()` allows querying of the data by:
* geographical range
* ocean basin
* parameter
* parameter detail
* sensor species
* age range
* minimum lenght (within the age range)
* minimum number of data points
* maximum time between data points
* minimum number of age control points

The function searches within calendar age expressed as ka BP or kyr. 

A complete parameter list can be found in Table 4 of https://doi.org/10.5194/essd-2019-223. Specifying parameter detail makes most sense when filtering near surface temperature estimates (using for instance: MgCa, UK37 or planktonic foraminifera assemblage).

The function either returns `TRUE/FALSE` indicating whether the query criteria are met when `extract_data = FALSE` (default) or a tibble with the (meta)data when `extract_data = TRUE`

`queryPALMOD()`  should be applied over a vector that contains the location of the individual files, which can be created with:

`PALMOD130k_V2_0_0 <- list.files("path", full.names = TRUE)`

### example to query the data
```
source('queryPALMOD.R')
meetsCrit <- queryPALMOD(sites = PALMOD130k_V2_0_0,
                         Parameter = "planktonic.d18O",
                         ParameterDetail = NULL,
                         SensorSpecies = NULL,
                         lat_min = -90,
                         lat_max = 90,
                         lon_min = -180,
                         lon_max = 180,
                         ocean_basin = NULL,
                         age_min = 0,
                         age_max = 2,
                         min_n_obs = 50,
                         min_duration = 0, 
                         max_sample_interval = Inf, 
                         min_n_tiepoint = 3,
                         extract_data = TRUE
                         )
```
