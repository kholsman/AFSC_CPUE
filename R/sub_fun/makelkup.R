#'
#'
#' 
makelkup<-function(conIN,outfl = "data/in/lookup_files"){

  cat(" -- Make lookup tables\n")

  # Get lookup tables
  #---------------------------------------------------
  cat(" --   NODC\n")
  NODC     <- sqlQuery(conIN,"SELECT * FROM foodlab.nodc")
  cat(" --   stations\n")
  stations <- sqlQuery(conIN,"SELECT * FROM foodlab.STATION_LOOKUP")
  cat(" --   Haul\n")
  HAUL      <- sqlQuery(conIN,"SELECT * FROM foodlab.HAUL")
  cat(" --   HAUL_Race\n")
  HAUL_Race <- sqlQuery(conIN,"SELECT * FROM racebase.HAUL")
  HAUL_Race$YEAR <- as.numeric(substr(HAUL_Race$START_TIME,1,4))
  #HAUL_RACE     <- sqlQuery(conIN,"SELECT * FROM RACEBASE.HAUL", max = 2)
  # cat(" --HAUL_RB\n")
  # HAUL_RB  <- sqlQuery(conIN,"SELECT * FROM racebase.HAUL")
  cat(" --   SPECIES\n")
  SPECIES  <- sqlQuery(conIN,"SELECT * FROM RACEBASE.SPECIES")
  cat(" --   surveys\n")
  surveys  <- sqlQuery(conIN,"SELECT * FROM RACE_DATA.V_CRUISES")
  cat(" --   STRATA_AREA\n")
  STRATA_AREA  <- (sqlQuery(con,"SELECT * FROM racebase.STRATUM"))

  cat(" --   saving files in ", outfl,"\n")
  
  save(HAUL,     file=file.path(outfl,"HAUL.Rdata"))
  save(HAUL_Race, file=file.path(outfl,"HAUL_Race.Rdata"))
  save(NODC,     file=file.path(outfl,"NODC.Rdata"))
  save(SPECIES,  file=file.path(outfl,"SPECIES.Rdata"))
  save(stations, file=file.path(outfl,"stations.Rdata"))
  save(surveys,  file=file.path(outfl,"surveys.Rdata"))
  save(STRATA_AREA,file=file.path(outfl,"STRATA_AREA.Rdata"))
  


}