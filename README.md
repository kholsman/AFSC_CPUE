<!-- default, tango, pygments, kate, monochrome, espresso, zenburn, haddock, and textmate. -->

#### [**AFSC Survey CPUE data: github.com/kholsman/AFSC_CPUE**](https://github.com/kholsman/AFSC_CPUE "AFSC_CPUE code Repo")

Repo maintained by:  
Kirstin Holsman  
Alaska Fisheries Science Center  
NOAA Fisheries, Seattle WA  
**[kirstin.holsman@noaa.gov](kirstin.holsman@noaa.gov)**  
*Last updated: Jan 05, 2023*

# Overview

The below scripts return a list object cpue_data saved as a compressed
Rdata file with the naming ‘reg.srvy#.spp.cpue_data.Rdata’ such as
“ebs.srvy98.plk.cpue_data.Rdata”. Each cpue_data list contains 8
data.frames:

``` r
load("data/out/2023_01_03/cpue/ebs/ebs.srvy98.plk.cpue_data.Rdata")

names(cpue_data)
```

    ## [1] "totalB_N_SEBS_NBS"      "totalB_N_SEBS"          "mnCPUE_strata_yr"      
    ## [4] "total_bin_B_N_SEBS"     "total_bin_B_N_SEBS_NBS" "mnCPUE_strata_bin_yr"  
    ## [7] "CPUE_station_bin_yr"    "CPUE_station_yr"

The data.frames are

1.  **totalB_N\_SEBS_NBS**: Total biomass (kg) or abundance (# of fish)
    for the species in each year for NEBS + SEBS survey results  

2.  **totalB_N\_SEBS**: Total biomass (kg) or abundance (# of fish) for
    the species in each year for only SEBS survey results

3.  **mnCPUE_strata_yr** : Average survey CPUE (kg per Km2) or abundance
    (# per Km2) for the species in each strata and year

4.  **total_bin_B\_N_SEBS**: Total biomass (kg) or abundance (# of fish)
    for each bin (10 mm) for the species in each year for NEBS + SEBS
    survey results

5.  **total_bin_B\_N_SEBS_NBS** : Total biomass (kg) or abundance (# of
    fish) for each bin (10 mm) for the species in each year for only the
    SEBS survey results

6.  **mnCPUE_strata_bin_yr**: Average survey CPUE (kg per Km2) or
    abundance (# per Km2) for each size bin for the species in each
    strata and year

7.  **CPUE_station_bin_yr**: Station specific survey CPUE (kg per Km2)
    or abundance (# per Km2) for each size bin for the species in
    eachyear

8.  **CPUE_station_yr**: Station specific survey CPUE (kg per Km2) or
    abundance (# per Km2) for the species in each year

These are calculated from the RACEBASE data tables for survey results
where total CPUE was recorded for the species *s* (location_catch) at
each haul, expanded to include stations *i* where CPUE=0 (location) and
expanded to each size bin *l* using the proportional subset of frequency
of fish of given length (mm), binned into 10 mm bins (*l*) and predicted
weight (*Ŵ*) for each size bin *l* at each station *i*:

$$B\_{s,y} =  \\bar{CPUE\_{s,k,y}} \\dot{}A\_{k}$$
where *A*<sub>*k*</sub> is the area of the strata *k* in
*K**m*<sup>2</sup> and $\\bar{CPUE\_{s,k,y}}$ is the strata specific
average CPUE (kg per *K**m*<sup>2</sup> or number per
*K**m*<sup>2</sup>) of all stations *i* in strata *k*:
$$\\bar{CPUE\_{s,k,y}} = \\frac{1}{n_k}\\dot{}\\sum\_{n_k}{CPUE\_{s,k,y,i}}$$
where $ CPUE\_{s,k,y,i} $ is the station specific CPUE (saves as the
object `cpue_data$CPUE_station_yr`).

To obtain population level estimates of the biomass or abundance of fish
by size bin *l*, we used a length weight regression to esimate the
weight of each size fish *j* measured (*Ŵ*) to calculate the proportion
by weight or frequency at each station where
*Ŵ* = *α*<sub>*s*</sub> + *L*<sub>*j*</sub><sup>*β*<sub>*s*</sub></sup>

and
$$p^w\_{l,i} = \\frac{N\_{l,i}\\dot{}\\hat{\\bar{W\_{l,i}}}}{\\sum\_{}{N\_{l,i}\\dot{}\\hat{\\bar{W\_{l,i}}}}}$$
and
$$p^N\_{l,i} = \\frac{N\_{l,i}}{\\sum\_{}{N\_{l,i}}}$$
This was then multiplied by the CPUE at each station
(*C**P**U**E*<sub>*s*, *k*, *y*, *i*</sub>) to obtain a station estimate
of CPUE by size bin *l*
$$CPUE\_{s,k,y,l,i} = p^N\_{l,i}\\dot{}CPUE\_{s,k,y,i}$$

Finally, the average strata CPUE ($\\bar{CPUE\_{s,k,y,l}}$) and whole of
EBS biomass by size bin (*B*<sub>*s*, *y*, *l*</sub>) was calculated as:
$$\\bar{CPUE\_{s,k,y,l}} = \\frac{1}{n_k}\\dot{}\\sum\_{n_k}{CPUE\_{s,k,y,l,i}}$$
and
$$B\_{s,y,l}= \\bar{CPUE\_{s,k,y,l}}\\dot{}A\_{k}$$
\# Code  
.

## Step 0: Set up the R workspace

The first step is to set up the switches for what files to update and
create in the file `R/setup.R`. The code below then loads these settings
as well as base data, functions, and packages.

``` r
  # get everything set up:
  #----------------------------------------
    # rm(list=ls())
# this uses the password saved in R/password.R
    suppressMessages(source("R/make.R"))
```

## Step 1: Update SQL queries (level1)

This step must be run on a computer that has access to RACEBASE. The
code below will generate the base files for steps 2 and 3 below,and will
save them in the folder data/in/2023_01_03 under subfolders for each
region in `srvys$reg` and each species in `splist` (see `R/setup.R` to
change these settings).

**IMPORTANT:**

-   **This step must be run in 32 bit R (RODBC doesn’t run in 64 bit R)
    and must be connected to the RACEBASE SQL database**

-   **To change R studio from the default 64 bit to 32 bit go to
    Tools>Global options and select the 32 bit version of R.**

-   **The code will connect to the SQL database using your password and
    username. Remember to update the `username_path` in the first line
    of the `R/setup.R` file and corresponding `username` and `password`
    under `username_password.R`. A template is available under `R/`.**

<!-- ![Header of `setup.R` where `username_path` can be adjusted. This file also is where species, regions, and bins are specified.](figs/setup.jpg){width=80%} -->

``` r
  # update the SQL queries
  #---------------------------------------------  
  source(file.path(code.path,"R/sub_scripts/runRACE_qrys.R"))
```

## Step 2: Update the LWA regressions

The default code for RACEBASE uses set LW relationships, however we
prefer to update the LW regressions using glms. Depending on how many
observations exist the LW relationships can be region specific or use
data across all regions.The default below is all regions combined. This
code generates two outputs in data/out/2023_01_03,
LWGlms_srvy_all_regs.Rdata and `LW_SmryTable.Rdata`. It also updates the
`species_lkup$LW_a` and species_lkup$LW_b\` parms used in Step 3.

``` r
  # update the LW regressions 
  #---------------------------------------------  

  if(update_LWdata==1){     
     source(file.path(code.path,"R/sub_scripts/updateLW.R"))
     # reload with updated data:
     source(file.path(code.path,"R/load_data.R"))
  }
  species_lkup
```

## Step 3: Get CPUE data from the surveys

This code is the core script for generating the CPUE_NUMKM2 and
CPUE_BIOMKM2 values by size bin, region, and species.

``` r
  nreg <- length(srvys$reg)
  nspp <- length(species_lkup$sp)
  
  for (r in 1:nreg){
    for(s in 1:nspp){
      flnm <- paste0(srvys[r,]$reg,".srvy",
                     srvys[r,]$num,".",
                     species_lkup[s,]$sp)
      cat("now getting data for: ",flnm,"\n")
      cpue_data <- suppressMessages(
        get_CPUE_DATA(
        datapath   = data.path,
        out_dir    = file.path(data.out),
        flnm       = flnm,
        species    = species_lkup[s,]$SPECIES_CODE,
        survey     = srvys[r,]$num,
        includeNBS = TRUE,
        saveit     = T,
        bins       = sp_bins[[ species_lkup[s,]$sp ]]))
     
      # # check the data :
      if(1==10){
      tt <- cpue_data%>%
            group_by(YEAR,REGION,STATIONID,SN)%>%
            filter(BIN ==400)%>%
            summarize(cnt =length(STATIONID))
       max(tt$cnt)  #Should be 1
      #this looks to be a duplicate sampling...
      #mis-entry or code error ?
       cpue_data%>%filter(YEAR==1988,STATIONID=="J-13")
      }
       rm(cpue_data)

    }
  }
```

*The cpue files are now saved in the directory data/out/2023_01_03/..*

## Step 4: calc propB for use in biomass weighting age or diet data

TBA

<img src="figs/out_dir.jpg" style="width:90.0%" /> ## Step 5: Get annual
biomass by age and CEATTLE bins:

## Appendix 1: `R/setup.R`primary setup script

<img src="figs/setup_large.jpg" style="width:90.0%" />
