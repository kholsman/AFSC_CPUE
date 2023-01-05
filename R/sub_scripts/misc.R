

# fill in the zeros for missing bins



# get the biomass by size bin



#count number of fish in each length bin:
# data1      <- length%>%
#   group_by(REGION,YEAR,STRATUM,STATIONID,LAT,LON,VESSEL,CRUISE,HAUL,HAULJOIN,
#            BIN,BIN_mm,CPUE_NUMKM2)%>%
#   summarize(FREQ = sum(FREQUENCY))%>%ungroup()
# 
#data1.w      <- length%>%
data1 <- length%>%
  group_by(REGION,YEAR,STRATUM,STATIONID,LAT,LON,VESSEL,CRUISE,HAUL,HAULJOIN,
           BIN,BIN_mm,CPUE_NUMKM2,CPUE_KGKM2)%>%
  summarize(FREQ  = sum(FREQUENCY),
            SUB_W = sum(W_hat))%>%ungroup()
# 
# #create sum of all fish measured at the station:
# #data2.w      <- length%>%
# data2      <- length%>%
#   group_by(REGION,YEAR,STRATUM,STATIONID,LAT,LON,VESSEL,CRUISE,HAUL,HAULJOIN,
#            BIN,BIN_mm,CPUE_NUMKM2,CPUE_KGKM2)%>%
#   summarize(SUM   = sum(FREQUENCY),
#             SUM_W = sum(W_hat))%>%ungroup()
# # data2      <- length%>%
# #   group_by(VESSEL,CRUISE,HAUL,HAULJOIN,YEAR,LAT,LON,BIN,BIN_mm,CPUE_NUMKM2,CPUE_KGKM2)%>%
# #   summarize(SUM = sum(FREQUENCY))%>%ungroup()
# 
# data3       <- data1%>%full_join(data2)
# data3.w     <- data1.w%>%full_join(data2.w)
# 
# sum is station total, freq is the bin total - prop is proportion by length bin
#data3$PROPN <- data3$FREQ/data3$SUM 
data1$PROPN <- data1$FREQ/data1$SUM 
# sum is station total, freq is the bin total - prop is proportion by length bin
data1$PROPW <- data1$SUB_W/data1$SUM_W  


data3$NUM_KM2    <- data3$PROPN*data3$CPUE_NUMKM2 # BIOM is the catch by length bin , is the number of fish per haul divided by the area swept so number/km^2
data3$BIOM_KGKM2 <- data3$PROPW*data3$CPUE_KGKM2 # BIOM is the catch by length bin , is the number of fish per haul divided by the area swept so number/km^2

data3_0    <- data3%>%
  group_by(VESSEL,CRUISE,HAUL,HAULJOIN,YEAR,LAT,LON,BIN)%>%
  summarize(
    BIOM_KGKM2 = sum(BIOM_KGKM2),
    NUM_KM2 = sum(NUM_KM2))%>%ungroup()

dataL     <- location%>%
  left_join(data3,by = c("VESSEL","CRUISE","HAUL","HAULJOIN","YEAR","LAT","LON"))
dataL     <- dataL[order(dataL$YEAR,
                         dataL$BIN,
                         dataL$HAULJOIN,
                         dataL$CRUISE,
                         dataL$VESSEL,
                         dataL$HAUL),]

nn        <- which(is.na(dataL$BIN))

dataL$BIOM_KGKM2[is.na(dataL$BIOM_KGKM2)] <-0   # empty catches
dataL$NUM_KM2[is.na(dataL$NUM_KM2)]       <-0
dataL$SPECIES_CODE<-species
dataL$SN      <- SN
dataL$CN      <- CN
dataL$LW_a    <- LW_a
dataL$LW_b    <- LW_b
dataL$REGION  <- region
dataL$reg     <- reg
dataL$srvy_nm <- survey_name

cpue_data <- dataL