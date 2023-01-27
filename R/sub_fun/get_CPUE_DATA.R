# species=spp$num[s]; survey = 98;bins = bins1; LW_a=LW_a_in[s]; LW_b=LW_b_in[s];species=10120

get_CPUE_DATA     <- function(
  datapath  = data.path,
  out_dir   = file.path(data.out,"../"),
  cpue_dir  = "cpue",
  saveit    = T,
  flnm      = flnm,
  species   = 10120,
  survey    = 98,
  includeNBS = TRUE,
  bins = c(0,100,200,300,400,500,600,700,800,900,1500)){
  
  reg    <- srvys[srvys$num ==survey,"reg"]
  rnum   <- srvys[srvys$num ==survey,"num"]
  
  sub <- surveys%>%filter(SURVEY_DEFINITION_ID==survey)
  region           <- as.character(sub$REGION[1])
  survey_name      <- as.character(sub$SURVEY_NAME[1])
  rm(sub)
  
  sub <- species_lkup%>%filter(SPECIES_CODE==species)
  sp    <- sub$sp
  spnum <- sub$SPECIES_CODE
  SN    <- sub$SPECIES_NAME
  CN    <- sub$COMMON_NAME
  LW_a  <- sub$LW_a
  LW_b  <- sub$LW_b
  rm(sub)
  
  sub_dir  <- file.path(datapath,reg)
  sp_dir   <- file.path(sub_dir,sp)
  
  load(file.path(sp_dir,"location.Rdata"))
  load(file.path(sp_dir,"location_catch.Rdata"))
  load(file.path(sp_dir,"length.Rdata"))
  
  length         <- length%>%filter(!is.na(STRATUM))
  location       <- location%>%filter(!is.na(STRATUM))
  location_catch <- location_catch%>%filter(!is.na(STRATUM),!is.na(CPUE_NUMKM2))
  
  if(location_catch$SPECIES_CODE[1]!=spnum)
    stop(cat("species codes do not match!location_cat species code = ",location_catch$SPECIES_CODE[1],
               ", while species num =", spnum ))
  if(length$SPECIES_CODE[1]!=spnum)
    stop(cat("species codes do not match!length species code = ",length$SPECIES_CODE[1],
             ", while species num =", spnum ))
  
  ## if survey is EBS, exclude far north regions.
  if(survey==98){
    if(!includeNBS){
      length         <- length%>%filter(!is.na(STRATUM),STRATUM<63)
      location       <- location%>%filter(!is.na(STRATUM),STRATUM<63)
      location_catch <- location_catch%>%filter(!is.na(STRATUM),STRATUM<63)
    }
  }
  
  ## for SLope survey exclude 2000 from all plots
  if(survey==78){
    # LWA_table     <- subset(LWA_table,lenLWA_tablegth$YEAR!=2000)
    length         <- subset(length,length$YEAR!=2000)
    location       <- subset(location,location$YEAR!=2000)
    location_catch <- subset(location_catch,location_catch$YEAR!=2000)
  }
  
  ##Exclude null locations and transform to all positive longitudes
  location       <- location%>%filter(!is.na(location$LON))
  location_catch <- location_catch%>%filter(!is.na(location_catch$LON))
  length         <- length%>%filter(!is.na(length$LON))

  # inner_join includes all row in x and y, 
  # full_join includes all rows in x or y
  length          <- length%>%
    inner_join(location_catch,
               by = c("VESSEL","CRUISE","HAUL","STRATUM",
                      "HAULJOIN","YEAR","LON","LAT",
                      "SURVEY_DEFINITION_ID", "SPECIES_CODE"))%>%
    ungroup()
  length          <- length%>%
    left_join(location,by = c("VESSEL","CRUISE",
                              "HAUL","STRATUM","TEMP","DEPTH",
                              "HAULJOIN","YEAR","LON","LAT"))%>%
    mutate(REGION=region)%>%ungroup()

  # ggplot(sub)+geom_point(aes(x=LENGTH,y=WEIGHT))
  # sub$W_hat <- exp(predict(LW.glm$WALLEYE_POLLOCK,
  #                          newdata=data.frame(L=sub$LENGTH/10)))
  # ggplot(sub)+
  #   geom_point(aes(x=LENGTH,y=WEIGHT))+
  #   geom_point(aes(x=LENGTH,y=W_hat),color="red")
    
  # Step 1. get W hat using LW_a and LW_b
  # -----------------------------------------
  if(dim(length)[1]>0){
    length$W_hat   <- LW_a*((length$LENGTH/10)^LW_b)
    length$W_hat   <- exp(log(LW_a)+log(length$LENGTH/10)*LW_b)
    length$qrydate <- as.character(qrydate)
  }
  cat(paste("LW_a=",round(LW_a*(10^3),4)," x10-3 and LW_b=",round(LW_b,4),"\n"))
  
  cat("\n assigning bins....\n"  )
 
  length$BIN    <- NA
  length$BIN_mm <- ""
  
  tmp    <- getBIN2(x=bins,BIN_IN = bins,divid=1,type=3)
  binALL <- data.frame(BIN = unlist(tmp[1,]),BIN_mm = unlist(tmp[2,]))
  rm(tmp)
 
  length <-length%>%mutate(BIN = getBIN2(x=LENGTH, BIN_IN=bins,divid=1,type=1),
                  BIN_mm =  getBIN2(x=LENGTH, BIN_IN=bins,divid=1,type=2))
  # rm(tmp)
  #     
  # pb <- txtProgressBar(min = 0, max = length(bins), style = 3)
  # 
  # if(dim(length)[1]>0){
  #   i<-1
  #   length$BIN[ length$LENGTH<=bins[i] & length$LENGTH >=0]  <- mean(c(0,bins[i]))
  #   length$BIN_mm[ length$LENGTH<=bins[i] & length$LENGTH >=0]  <- paste0("[0,",bins[i],")")
  #   BIN=mean(c(0,bins[i]))
  #   BIN_mm =paste0("[0,",bins[i],")") 
  #   
  #   for ( i in 2:length(bins)){
  #     # update progress bar
  #     setTxtProgressBar(pb, i)
  #     
  #     length$BIN[length$LENGTH<= bins[i] & length$LENGTH > bins[(i-1)] ] <-  mean(bins[(i-1):i])
  #     length$BIN_mm[length$LENGTH<= bins[i] & length$LENGTH > bins[(i-1)] ] <-   paste0("[",bins[(i-1)],",",bins[i],")")
  #  
  #     BIN    <-  c(BIN,mean(bins[(i-1):i]))
  #     BIN_mm <-  c(BIN_mm,paste0("[",bins[(i-1)],",",bins[i],")"))
  # 
  #   }
  #   #plus group
  #   length$BIN[length$LENGTH> bins[i] ] <- mean(c(bins[i], bins[i]+mean(bins[-1]-bins[1:(length(bins)-1)])))
  #   length$BIN_mm[length$LENGTH>bins[i] ] <-   paste0(bins[(i)],"+")
  #    BIN    <-  c(BIN,mean(c(bins[i], bins[i]+mean(bins[-1]-bins[1:(length(bins)-1)]))))
  #    BIN_mm <-  c(BIN_mm, paste0(bins[(i)],"+"))
  #    binALL <- data.frame(BIN,BIN_mm)
  #   
  #   close(pb)
  # }
  # cat("complete\n")
  
  # Step 2. Get annual biomass and abundance estimates
  # -----------------------------------------
  
  # get Strata AREA estimates from 2019 - they were updated in 2022 but are larger??
  sub_SA <- STRATA_AREA%>%filter(YEAR>2009,YEAR<2022)%>%
    group_by(REGION,STRATUM)%>%
    summarize(AREA = mean(AREA, na.rm=T))%>%ungroup()
  sub_SA <- STRATA_AREA%>%filter(YEAR==2022)%>%
    group_by(REGION,STRATUM)%>%
    summarize(AREA = mean(AREA, na.rm=T))%>%ungroup()
  
  # create full matrix of stations and populate 0s for where haul took place 
  # but CPUE == 0 using location (all stations sampled) and location_catch (with CPUE>0)
  station <- location%>%
    group_by(YEAR, STATIONID,HAULJOIN, VESSEL, 
             CRUISE, HAUL,MONTH, DAY,LAT, LON)%>%
    summarize(STRATUM = mean(STRATUM,na.rm=T),
              TEMP    = mean(TEMP,na.rm=T),
              SST     = mean(SST,na.rm=T),
              DEPTH   = mean(DEPTH,na.rm=T))

  station_catch <- location_catch%>%
    group_by(YEAR,HAULJOIN, VESSEL, CRUISE, HAUL,LAT, LON)%>%
    summarize(STRATUM = mean(STRATUM,na.rm=T),
              CPUE_NUMKM2  = mean(CPUE_NUMKM2,na.rm=T),
              CPUE_KGKM2    = mean(CPUE_KGKM2,na.rm=T),
              SURVEY_DEFINITION_ID = mean(SURVEY_DEFINITION_ID, na.rm=T))
  
  # get station specific data and join with strata AREA estimates from 2019
  CPUE_station_yr <- 
    station%>% 
    full_join(station_catch) %>%
    distinct(HAULJOIN,.keep_all = TRUE)%>%
    mutate(REGION=region)%>%ungroup()
  
  nn <- which(is.na(CPUE_station_yr$CPUE_KGKM2))
  if(any(CPUE_station_yr$HAULJOIN[nn]%in%station_catch$HAULJOIN))
    stop("error with duplicated HAULJOINS in output")
  CPUE_station_yr$CPUE_NUMKM2[nn]<- CPUE_station_yr$CPUE_KGKM2[nn] <-0
  
  CPUE_station_yr <- CPUE_station_yr%>%left_join(sub_SA)%>%
    mutate(SURVEY_DEFINITION_ID =survey,
           SPECIES_CODE = spnum)
  
  se<-function(x,na.rm=T){
    if(na.rm==T){
      if(any(is.na(x)))
        x<-as.numeric(na.omit(x))
    }
    return(sd(x)/sqrt(length(x)))
  }
  # get mean CPUE by strata
  mnCPUE_strata_yr <- CPUE_station_yr%>%
    group_by(REGION,YEAR,STRATUM,SPECIES_CODE)%>%
    summarize(
      AREA        = mean(AREA, na.rm=T),
      nobs_yk  = length(CPUE_KGKM2),
      mnCPUE_KGKM2_yk  = mean(CPUE_KGKM2,na.rm=T),
      mnCPUE_NUMKM2_yk = mean(CPUE_NUMKM2,na.rm=T),
      sdCPUE_KGKM2_yk  = sd(CPUE_KGKM2,na.rm=T),
      sdCPUE_NUMKM2_yk = sd(CPUE_NUMKM2,na.rm=T),
      seCPUE_KGKM2_yk  = se(CPUE_KGKM2,na.rm=T),
      seCPUE_NUMKM2_yk = se(CPUE_NUMKM2,na.rm=T))%>%
    mutate(
      # seCPUE_KGKM2_yk  = sdCPUE_KGKM2_yk/sqrt(nCPUE_KGKM2_yk),
           # seCPUE_NUMKM2_yk = sdCPUE_NUMKM2_yk/sqrt(nCPUE_NUMKM2_yk),
           B_KG_yk  = mnCPUE_KGKM2_yk*AREA,
           N_yk     = mnCPUE_NUMKM2_yk*AREA,
           sdB_KG_yk  = sdCPUE_KGKM2_yk*AREA,
           sdN_yk     = sdCPUE_NUMKM2_yk*AREA,
           seB_KG_yk  = seCPUE_KGKM2_yk*AREA,
           seN_yk     = seCPUE_NUMKM2_yk*AREA)%>%
    left_join(species_lkup%>%
                rename(CN=COMMON_NAME,SN = SPECIES_NAME)%>%
                select(SPECIES_CODE,CN,SN,sp,num,NAME))%>%
    ungroup()
  
  #Total Biomass and abundance for NBS+SEBS annually
  totalB_N_SEBS_NBS <- mnCPUE_strata_yr%>%
    group_by(REGION, YEAR,SPECIES_CODE,CN,SN,sp,num,NAME)%>%
    summarize(nstrata = length(mnCPUE_KGKM2_yk),
              nobs      = sum(nobs_yk),
              B_KG_y    = sum(B_KG_yk,na.rm=T),
              N_y       = sum(N_yk,na.rm=T),
              sdB_KG_y  = sum(sdB_KG_yk,na.rm=T),
              sdN_y     = sum(sdN_yk,na.rm=T),
              seB_KG_y  = sum(seB_KG_yk,na.rm=T),
              seN_y     = sum(seN_yk,na.rm=T))%>%ungroup()
  
  #Total Biomass and abundance for SEBS only annually
  totalB_N_SEBS <- mnCPUE_strata_yr%>%filter(STRATUM<63)%>%
    group_by(REGION,YEAR,SPECIES_CODE,CN,SN,sp,num,NAME)%>%
    summarize(nstrata = length(mnCPUE_KGKM2_yk),
              nobs      = sum(nobs_yk),
              B_KG_y    = sum(B_KG_yk,na.rm=T),
              N_y       = sum(N_yk,na.rm=T),
              sdB_KG_y  = sum(sdB_KG_yk,na.rm=T),
              sdN_y     = sum(sdN_yk,na.rm=T),
              seB_KG_y  = sum(seB_KG_yk,na.rm=T),
              seN_y     = sum(seN_yk,na.rm=T))%>%ungroup()
  # source("data/in/lookup_files/2022_assessmentB.R")
  # # load obs data
  # tmp<- obs%>%filter(species==species_lkup[species_lkup$SPECIES_CODE==spnum,]$sp)
  # ggplot()+
  #   geom_line(data=totalB_N_SEBS_NBS,aes(x=YEAR,y=totB),size=.8)+
  #   geom_line(data=totalB_N_SEBS,aes(x=YEAR,y=totB),color="red",size=.8,linetype="dashed")+
  #   geom_errorbar(data=totalB_N_SEBS_NBS,
  #                 aes(x=YEAR+.3,ymin=(totB-totB_se),ymax=(totB+totB_se)),color="black")+
  #   geom_point(data=tmp,
  #              aes(x=YEAR,y=totB*1000),color="blue",size=2)+
  #   geom_errorbar(data=tmp,
  #              aes(x=YEAR,ymin=(totB-totBse)*1000,ymax=(totB+totBse)*1000),color="blue")

 # Step 3. Get proportion of total biomass by bin by station i 
 # ----------------------------------------- 

  station$ID <- 1:dim(station)[1]
  bins_dat   <- expand.grid(ID=station$ID,BIN=binALL$BIN)
  bins_dat   <- bins_dat%>%left_join(station)%>%left_join(binALL)
  bins_dat   <- bins_dat%>%
    mutate(LAB = paste0(BIN,"_",YEAR,"_",STRATUM,"_",STATIONID))
  
  # collapse freq and mean weight into bins by station and year
  station_bin_yr <- length%>%
    group_by(qrydate,YEAR,STATIONID,HAULJOIN, 
             VESSEL, CRUISE, HAUL,
             MONTH, DAY,LAT, LON, STRATUM,
             BIN,BIN_mm )%>%
    summarize(STRATUM=mean(STRATUM,na.rm=T),
              TEMP = mean(TEMP,na.rm=T),
              SST  = mean(SST,na.rm=T),
              DEPTH = mean(DEPTH,na.rm=T),
              SPECIES_CODE = mean(SPECIES_CODE,na.rm=T),
              SURVEY_DEFINITION_ID = mean(SURVEY_DEFINITION_ID,na.rm=T),
              CPUE_NUMKM2  = mean(CPUE_NUMKM2,na.rm=T),
              CPUE_KGKM2   = mean(CPUE_KGKM2,na.rm=T),
              LENGTH       = mean(LENGTH,na.rm=T),
              FREQUENCY    = sum(FREQUENCY,na.rm=T),
              W_hat        = mean(W_hat,na.rm=T))%>%
    mutate(LAB = paste0(BIN,"_",YEAR,"_",STRATUM,"_",STATIONID))%>%ungroup()
  
  
  # create empty matrix of unmeasured bins to add to each station
    bin_ADD <- bins_dat%>%left_join(CPUE_station_yr)
    
    bin_ADD$FREQUENCY    <- 
      bin_ADD$W_hat      <-  
      bin_ADD$LENGTH     <- 0
    bin_ADD$REGION       <- region
    bin_ADD$SPECIES_CODE <- spnum
    bin_ADD$qrydate      <- qrydate
  
    CPUE_station_bin_yr  <- 
      bin_ADD%>% 
      full_join(station_bin_yr)%>%
      distinct%>%
      mutate(REGION=region)%>%ungroup()
    
  
    data.frame(CPUE_station_bin_yr%>%filter(LAB == bin_ADD$LAB[1]))
    dim(CPUE_station_bin_yr)
    
    nn <- which(is.na(CPUE_station_bin_yr$CPUE_KGKM2))
    
    if(length(nn)>0)
      stop("error with CPUE_station_bin_yr CPUE == NA")
   
  # collapse freq and mean weight into station and year
  data1 <- CPUE_station_bin_yr%>%
    group_by(REGION,YEAR,STRATUM,STATIONID,LAT,LON,VESSEL,CRUISE,HAUL,HAULJOIN,
             CPUE_NUMKM2,CPUE_KGKM2)%>%
    summarize(sumFREQ  = sum(FREQUENCY,na.rm=T),
              sumW_hat = sum(W_hat*FREQUENCY,na.rm=T))%>%ungroup()
  
  # get propW and propN by station year and stratum
  CPUE_station_bin_yr <- CPUE_station_bin_yr%>%
    left_join(data1)%>%
    mutate(LAB   = paste0(BIN,"_",YEAR,"_",STRATUM,"_",STATIONID),
           propN = FREQUENCY/sumFREQ,
           propW = (W_hat*FREQUENCY)/sumW_hat)%>%ungroup()
  
  # check the data:
  datacheck <- CPUE_station_bin_yr%>%
    group_by(REGION,YEAR,STRATUM,STATIONID,LAT,LON,VESSEL,CRUISE,HAUL,HAULJOIN,
             CPUE_NUMKM2,CPUE_KGKM2)%>%
    summarize(sum_propN  = sum(propN,na.rm=T),
              sum_propW = sum(propW,na.rm=T))%>%ungroup()
  
  if(any(unique(datacheck$sum_propW)<.999&unique(datacheck$sum_propW)>0))
    stop("problem with datacheck")
  
  CPUE_station_bin_yr <- CPUE_station_bin_yr%>%
    mutate(bin_CPUE_NUMKM2  = CPUE_NUMKM2*propN,
           bin_CPUE_KGKM2   = CPUE_KGKM2*propW)%>%
    rename(haul_CPUE_NUMKM2 = CPUE_NUMKM2,
           haul_CPUE_KGKM2  = CPUE_KGKM2)
  
 
  # Step 3 get ebs wide biomass by bin:
  # get mean CPUE by strata and bin
  mnCPUE_strata_bin_yr <- CPUE_station_bin_yr%>%
    group_by(REGION,YEAR,STRATUM,SPECIES_CODE,BIN,BIN_mm)%>%
    summarize(
      AREA              = mean(AREA, na.rm=T),
      nobs_ykl          = length(bin_CPUE_KGKM2),
      mnCPUE_KGKM2_ykl  = mean(bin_CPUE_KGKM2,na.rm=T),
      mnCPUE_NUMKM2_ykl = mean(bin_CPUE_NUMKM2,na.rm=T),
      sdCPUE_KGKM2_ykl  = sd(bin_CPUE_KGKM2,na.rm=T),
      sdCPUE_NUMKM2_ykl = sd(bin_CPUE_NUMKM2,na.rm=T),
      seCPUE_KGKM2_ykl  = se(bin_CPUE_KGKM2,na.rm=T),
      seCPUE_NUMKM2_ykl = se(bin_CPUE_NUMKM2,na.rm=T))%>%
    mutate(B_KG_ykl     = mnCPUE_KGKM2_ykl*AREA,
           N_ykl        = mnCPUE_NUMKM2_ykl*AREA,
           sdB_KG_ykl   = sdCPUE_KGKM2_ykl*AREA,
           sdN_ykl      = sdCPUE_NUMKM2_ykl*AREA,
           seB_KG_ykl   = seCPUE_KGKM2_ykl*AREA,
           seN_ykl      = seCPUE_NUMKM2_ykl*AREA)%>%left_join(sub_SA%>%rename(AREA2=AREA))%>%
    left_join(species_lkup%>%
                rename(CN=COMMON_NAME,SN = SPECIES_NAME)%>%
                select(SPECIES_CODE,CN,SN,sp,num,NAME))%>%ungroup()
  
  tmp<-mnCPUE_strata_bin_yr%>%ungroup()%>%
    group_by(REGION,YEAR,SPECIES_CODE,CN,SN,sp,num,NAME)%>%
    summarize(nobs_y   = sum(nobs_ykl),
              sumAREA_y   = sum(AREA,na.rm=T),
              sumB_KG_yl    = sum(B_KG_ykl,na.rm=T),
              sumN_yl       = sum(N_ykl,na.rm=T))%>%
    left_join(totalB_N_SEBS_NBS)%>%
    ungroup()
  
  mnCPUE_strata_bin_yr_SEBS_NBS<-mnCPUE_strata_bin_yr%>%
    left_join(tmp)%>%
    mutate(B_KG_ykl    = (B_KG_y/sumB_KG_yl)*B_KG_ykl,
           N_ykl       = (N_y/sumN_yl)*N_ykl,
           sdB_KG_yl   = (B_KG_y/sumB_KG_yl)*sdB_KG_ykl,
           sdN_yl      = (N_y/sumN_yl)*sdN_ykl,
           seB_KG_yl   = (B_KG_y/sumB_KG_yl)*seB_KG_ykl,
           seN_yl      = (N_y/sumN_yl)*seN_ykl)%>%
    ungroup()
  
  
  #Total Biomass and abundance for NBS+SEBS annually
  total_bin_B_N_SEBS_NBS <- mnCPUE_strata_bin_yr_SEBS_NBS%>%
    group_by(REGION,YEAR,SPECIES_CODE,CN,SN,sp,num,NAME,BIN,BIN_mm)%>%
    summarize(nobs_yl   = sum(nobs_ykl),
              sumAREA_yl = sum(AREA,na.rm=T),
              B_KG_yl    = sum(B_KG_ykl,na.rm=T),
              N_yl       = sum(N_ykl,na.rm=T),
              sdB_KG_yl  = sum(sdB_KG_ykl,na.rm=T),
              sdN_yl     = sum(sdN_ykl,na.rm=T),
              seB_KG_yl  = sum(seB_KG_ykl,na.rm=T),
              seN_yl     = sum(seN_ykl,na.rm=T))%>%
    ungroup()
  
  rm(tmp)
  
  tmp<-mnCPUE_strata_bin_yr%>%ungroup()%>%filter(STRATUM<63)%>%
    group_by(REGION,YEAR,SPECIES_CODE,CN,SN,sp,num,NAME)%>%
    summarize(nobs_y   = sum(nobs_ykl),
              sumAREA_y   = sum(AREA,na.rm=T),
              sumB_KG_yl    = sum(B_KG_ykl,na.rm=T),
              sumN_yl       = sum(N_ykl,na.rm=T))%>%
    left_join(totalB_N_SEBS)%>%
    ungroup()
  
  mnCPUE_strata_bin_yr_SEBS<-mnCPUE_strata_bin_yr%>%
    left_join(tmp)%>%
    mutate(B_KG_ykl    = (B_KG_y/sumB_KG_yl)*B_KG_ykl,
           N_ykl       = (N_y/sumN_yl)*N_ykl,
           sdB_KG_yl   = (B_KG_y/sumB_KG_yl)*sdB_KG_ykl,
           sdN_yl      = (N_y/sumN_yl)*sdN_ykl,
           seB_KG_yl   = (B_KG_y/sumB_KG_yl)*seB_KG_ykl,
           seN_yl      = (N_y/sumN_yl)*seN_ykl)%>%
    ungroup()
  
  
  #Total Biomass and abundance for SEBS only annually
  total_bin_B_N_SEBS <- mnCPUE_strata_bin_yr_SEBS%>%filter(STRATUM<63)%>%
    group_by(REGION,YEAR,SPECIES_CODE,CN,SN,sp,num,NAME,BIN,BIN_mm)%>%
    summarize(nobs_yl   = sum(nobs_ykl),
              sumAREA_yl = sum(AREA,na.rm=T),
              B_KG_yl    = sum(B_KG_ykl,na.rm=T),
              N_yl       = sum(N_ykl,na.rm=T),
              sdB_KG_yl  = sum(sdB_KG_ykl,na.rm=T),
              sdN_yl     = sum(sdN_ykl,na.rm=T),
              seB_KG_yl  = sum(seB_KG_ykl,na.rm=T),
              seN_yl     = sum(seN_ykl,na.rm=T))%>%
    ungroup()
  rm(mnCPUE_strata_bin_yr)
  
  # k = strata; l = bin; y = year
  # use propB_ykl If you want annual sums  use  and summarize as : val*propB_ykl & group_by(YEAR,REGION,SN,CN)%>%
  # use propB_yl  If you want annual and strata sums 
  #     (across bins) : val*propB_yl & group_by(YEAR,REGION,STRATUM,SN,CN)%>%
  # use propB_yk  If you want annual and bin sums 
  #     (across across strata) : val*propB_yk &  group_by(YEAR,REGION,BIN_cm_mid,SN,CN)%>%
  
  propByStrataBin_SEBS <- mnCPUE_strata_bin_yr_SEBS%>%
    select("REGION","YEAR","STRATUM",SPECIES_CODE,CN,SN,sp,num,NAME,
           "BIN","BIN_mm", "B_KG_ykl"  ,"N_ykl",
           "nobs_ykl","sdB_KG_ykl","sdN_ykl","seB_KG_ykl","seN_ykl") %>%
    left_join(totalB_N_SEBS%>%
                select("REGION","YEAR","SPECIES_CODE", "B_KG_y"  ,"N_y",
                       "nobs","sdB_KG_y","sdN_y","seB_KG_y","seN_y"))
  
  tmpP<-propByStrataBin_SEBS
  # "prop B in a given strata and length bin": use this when collapsing dbin and strata specific values up to the regional level:
  tmpP$propB_ykl     <- 0; cc<- which(tmpP$B_KG_y>0)
  tmpP$propB_ykl[cc] <- tmpP$B_KG_ykl[cc]/tmpP$B_KG_y[cc]
  tmpP$propN_ykl     <- 0;cc<- which(tmpP$N_y>0)
  tmpP$propN_ykl[cc] <- tmpP$N_ykl[cc]/tmpP$N_y[cc] 
  propByStrataBin_SEBS<-tmpP%>%ungroup()
  rm(tmpP)
  
  propByStrataBin_SEBS_NBS <- mnCPUE_strata_bin_yr_SEBS_NBS%>%
    select("REGION","YEAR","STRATUM",SPECIES_CODE,CN,SN,sp,num,NAME,
           "BIN","BIN_mm", "B_KG_ykl"  ,"N_ykl",
           "nobs_ykl","sdB_KG_ykl","sdN_ykl","seB_KG_ykl","seN_ykl") %>%
    left_join(totalB_N_SEBS_NBS%>%
                select("REGION","YEAR","SPECIES_CODE", "B_KG_y"  ,"N_y",
                       "nobs","sdB_KG_y","sdN_y","seB_KG_y","seN_y"))
  
  tmpP<-propByStrataBin_SEBS_NBS
  # "prop B in a given strata and length bin": use this when collapsing dbin and strata specific values up to the regional level:
  tmpP$propB_ykl     <- 0; cc<- which(tmpP$B_KG_y>0)
  tmpP$propB_ykl[cc] <- tmpP$B_KG_ykl[cc]/tmpP$B_KG_y[cc]
  tmpP$propN_ykl     <- 0;cc<- which(tmpP$N_y>0)
  tmpP$propN_ykl[cc] <- tmpP$N_ykl[cc]/tmpP$N_y[cc] 
  propByStrataBin_SEBS_NBS<-tmpP%>%ungroup()
  rm(tmpP)
    
  propByStrata <- mnCPUE_strata_yr%>%
    select("REGION","YEAR","STRATUM",SPECIES_CODE,CN,SN,sp,num,NAME,
           "B_KG_yk"  ,"N_yk",
           "nobs_yk",
           "sdB_KG_yk","sdN_yk",
           "seB_KG_yk","seN_yk") %>%
    left_join(totalB_N_SEBS_NBS%>%
                select("REGION","YEAR","SPECIES_CODE",
                       "B_KG_y"  ,"N_y",
                       "nobs",
                       "sdB_KG_y","sdN_y",
                       "seB_KG_y","seN_y"))
  tmpP<-propByStrata
  # "prop B in a given strata": use this when collapsing strata specific values up to the regional level:
  tmpP$propB_yk     <- 0; cc<- which(tmpP$B_KG_y>0)
  tmpP$propB_yk[cc] <- tmpP$B_KG_yk[cc]/tmpP$B_KG_y[cc]
  tmpP$propN_yk     <- 0;cc<- which(tmpP$N_y>0)
  tmpP$propN_yk[cc] <- tmpP$N_yk[cc]/tmpP$N_y[cc] 
  propByStrata<-tmpP%>%ungroup()
  rm(tmpP)
  
  propByBin<- total_bin_B_N_SEBS_NBS%>%
    select("REGION","YEAR",SPECIES_CODE,CN,SN,sp,num,NAME,"BIN","BIN_mm",
           "nobs_yl",
           "B_KG_yl"  ,"N_yl",
           "nobs_yl",
           "sdB_KG_yl","sdN_yl",
           "seB_KG_yl","seN_yl") %>%
    left_join(totalB_N_SEBS_NBS%>%
                select("REGION","YEAR","SPECIES_CODE",
                       "B_KG_y"  ,"N_y",
                       "nobs",
                       "sdB_KG_y","sdN_y",
                       "seB_KG_y","seN_y"))%>%ungroup()
  tmpP<-propByBin
  # "prop B in a given bin": use this when collapsing bin specific values up to the regional level:
  tmpP$propB_yl     <- 0; cc<- which(tmpP$B_KG_y>0)
  tmpP$propB_yl[cc] <- tmpP$B_KG_yl[cc]/tmpP$B_KG_y[cc]
  tmpP$propN_yl     <- 0;cc<- which(tmpP$N_y>0)
  tmpP$propN_yl[cc] <- tmpP$N_yl[cc]/tmpP$N_y[cc] 
  propByBin<-tmpP%>%ungroup()
  rm(tmpP)
  

  
  #double check the results
  cnt_ByStrataBin <- propByStrataBin_SEBS_NBS%>%
    select("REGION","YEAR","STRATUM",BIN,BIN_mm,
           SPECIES_CODE,CN,SN,sp,num,
           "propB_ykl","propN_ykl")%>%
    group_by(YEAR,REGION,SN,CN)%>%
    summarise(sum_propB_ykl=sum(propB_ykl,na.rm=T),
              sum_propN_ykl=sum(propN_ykl,na.rm=T))
  
  cnt_ByStrataBin_SEBS <- propByStrataBin_SEBS%>%
    select("REGION","YEAR","STRATUM",BIN,BIN_mm,
           SPECIES_CODE,CN,SN,sp,num,
           "propB_ykl","propN_ykl")%>%
    group_by(YEAR,REGION,SN,CN)%>%
    summarise(sum_propB_ykl=sum(propB_ykl,na.rm=T),
              sum_propN_ykl=sum(propN_ykl,na.rm=T))
  
  cnt_ByStrata <- propByStrata%>%
    select("REGION","YEAR","STRATUM",
           SPECIES_CODE,CN,SN,sp,num,
           "propB_yk","propN_yk")%>%
    group_by(YEAR,REGION,SN,CN)%>%
    summarise(sum_propB_yk=sum(propB_yk,na.rm=T),
              sum_propN_yk=sum(propN_yk,na.rm=T))
  
  
  cnt_ByBin <- propByBin%>%
    select("REGION","YEAR",BIN,BIN_mm,num,
           SPECIES_CODE,CN,SN,sp,
           "propB_yl","propN_yl")%>%
    group_by(YEAR,REGION,SN,CN)%>%
    summarise(sum_propB_yl=sum(propB_yl ,na.rm=T),
              sum_propN_yl=sum(propN_yl ,na.rm=T))

  
  
  
  rm(list=c("location","location_catch","length"))
  cat("saving results...")
  if(!dir.exists(out_dir)) dir.create(out_dir)
  cpue_dir <- file.path(out_dir,"cpue")
  if(!dir.exists(cpue_dir)) dir.create(cpue_dir)
  cpue_dir <- file.path(cpue_dir,reg)
  if(!dir.exists(cpue_dir)) dir.create(cpue_dir)
  # 
  # save(totalB_N_SEBS_NBS,file = file.path(cpue_dir,paste0(flnm,".totalB_N_SEBS_NBS.Rdata")))
  # save(totalB_N_SEBS,file = file.path(cpue_dir,paste0(flnm,".totalB_N_SEBS.Rdata")))
  # save(mnCPUE_strata_yr,file = file.path(cpue_dir,paste0(flnm,".mnCPUE_strata_yr.Rdata")))
  # 
  # save(total_bin_B_N_SEBS,file = file.path(cpue_dir,paste0(flnm,".total_bin_B_N_SEBS.Rdata")))
  # save(total_bin_B_N_SEBS_NBS,file = file.path(cpue_dir,paste0(flnm,".total_bin_B_N_SEBS_NBS.Rdata")))
  # save(mnCPUE_strata_bin_yr,file = file.path(cpue_dir,paste0(flnm,".mnCPUE_strata_bin_yr.Rdata")))
  # 
  # save(CPUE_station_bin_yr,file = file.path(cpue_dir,paste0(flnm,".CPUE_station_bin_yr.Rdata")))
  # save(CPUE_station_yr,file = file.path(cpue_dir,paste0(flnm,".CPUE_station_yr.Rdata")))
 
  cpue_data <-list(totalB_N_SEBS_NBS = totalB_N_SEBS_NBS,
                   totalB_N_SEBS = totalB_N_SEBS,
                   mnCPUE_strata_yr = mnCPUE_strata_yr,
                   total_bin_B_N_SEBS = total_bin_B_N_SEBS,
                   total_bin_B_N_SEBS_NBS = total_bin_B_N_SEBS_NBS,
                   mnCPUE_strata_bin_yr_SEBS = mnCPUE_strata_bin_yr_SEBS%>%
                     select(-"sumAREA_y",-"sumB_KG_yl",-"sumN_yl",-AREA2),
                   mnCPUE_strata_bin_yr_SEBS_NBS = mnCPUE_strata_bin_yr_SEBS_NBS%>%
                     select(-"sumAREA_y",-"sumB_KG_yl",-"sumN_yl",-AREA2),
                   CPUE_station_bin_yr = CPUE_station_bin_yr,
                   CPUE_station_yr = CPUE_station_yr,
                   propByBin = propByBin%>%
                     select("REGION","YEAR",BIN,BIN_mm,num,
                            SPECIES_CODE,CN,SN,sp,"propB_yl","propN_yl"),
                   propByStrata = propByStrata%>%
                     select("REGION","YEAR","STRATUM",
                            SPECIES_CODE,CN,SN,sp,num,
                            "propB_yk","propN_yk"),
                   propByStrataBin_SEBS = propByStrataBin_SEBS%>%
                     select("REGION","YEAR","STRATUM",BIN,BIN_mm,
                            SPECIES_CODE,CN,SN,sp,num,
                            "propB_ykl","propN_ykl"),
                   propByStrataBin_SEBS_NBS = propByStrataBin_SEBS_NBS%>%
                     select("REGION","YEAR","STRATUM",BIN,BIN_mm,
                            SPECIES_CODE,CN,SN,sp,num,
                            "propB_ykl","propN_ykl")
                   )
  save(cpue_data,file = file.path(cpue_dir,paste0(flnm,".cpue_data.Rdata")))
  
  cat("files saved \n")
  return(cpue_data)
  
}