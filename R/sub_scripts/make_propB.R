# make_propB:
# ------------------------------------
# k = strata; l = bin; y = year
# use propB_ykl If you want annual sums  use  and summarize as : val*propB_ykl & group_by(YEAR,REGION,SN,CN)%>%
# use propB_yl  If you want annual and strata sums 
#     (across bins) : val*propB_yl & group_by(YEAR,REGION,STRATUM,SN,CN)%>%
# use propB_yk  If you want annual and bin sums 
#     (across across strata) : val*propB_yk &  group_by(YEAR,REGION,BIN_cm_mid,SN,CN)%>%

propB_tmp <- makePropB2(dat=datUSE,binlimtIN = binlimt,BigbinIN=BigBIN)

# GET PROP B and fill in missing years with annual averages:
propB     <- propB_tmp$propB%>%
  select("YEAR","REGION","BIN_cm_mid","STRATUM","CN","SN",
         "propB_ykl","propN_ykl",
         "propB_yk","propN_yk",
         "propB_yl","propN_yl")

# fill in missing years with average across years:
YEARS   <- sort(unique(unique(c_USE$YEAR),unique(propB$YEAR)))
STRATA  <- sort(unique(unique(c_USE$STRATUM),unique(propB$STRATUM)))
BINS    <- sort(unique(unique(c_USE$BIN_cm_mid),unique(propB$BIN_cm_mid)))
SPECIES <- sort(unique(unique(c_USE$CN),unique(propB$CN)))
missing <- list()
missing[["YRS"]]    <- setdiff(YEARS,unique(propB$YEAR))
missing[["STRATA"]] <- setdiff(STRATA,unique(propB$STRATUM))
missing[["BINS"]]   <- setdiff(BINS,unique(propB$BIN_cm_mid))

# mean propB by year(y) strata(K) bin(l) 
mn_propB_ykl  <- propB%>%
  group_by(REGION,BIN_cm_mid,STRATUM,CN,SN)%>%
  summarise_at(c("propB_ykl","propN_ykl"),mean,na.rm=T)
# now recenter using sumB
mn_propB_ykl$propB_ykl <-  mn_propB_ykl$propB_ykl/sum(mn_propB_ykl$propB_ykl)
mn_propB_ykl$propN_ykl <-  mn_propB_ykl$propN_ykl/sum(mn_propB_ykl$propN_ykl)

# mean propB by year(y) and strata (k) (uses to get bin mean values)
mn_propB_yk  <- propB%>%
  group_by(REGION,STRATUM,CN,SN)%>%
  summarise_at(c("propB_yk","propN_yk"),mean,na.rm=T)
mn_propB_yk$propB_yk <-  mn_propB_yk$propB_yk/sum(mn_propB_yk$propB_yk)
mn_propB_yk$propN_yk <-  mn_propB_yk$propN_yk/sum(mn_propB_yk$propN_yk)

# mean propB by year(y) and length (l) (used to get strata mean values)
mn_propB_yl  <- propB%>%
  group_by(REGION,BIN_cm_mid,CN,SN)%>%
  summarise_at(c("propB_yl","propN_yl"),mean,na.rm=T)
mn_propB_yl$propB_yl <-  mn_propB_yl$propB_yl/sum(mn_propB_yl$propB_yl)
mn_propB_yl$propN_yl <-  mn_propB_yl$propN_yl/sum(mn_propB_yl$propN_yl)

n_strata     <- length(unique(propB$STRATUM))
n_bin        <- length(unique(propB$BIN_cm_mid))
n_stratabin  <- 1

# if values are missing from the diet database 
# (i.e., bin, strata year combo not sampled, add NA for propB)

# only for missing years: fill in missing stratum with averages across years:
if(length(missing[["YRS"]])>0){
  
  #create prop B for missing year:
  tmp_grid <- expand.grid(YEAR         = missing[["YRS"]],
                          STRATUM      = unique(propB$STRATUM),
                          BIN_cm_mid   = unique(propB$BIN_cm_mid))
  tmp      <- merge(tmp_grid,mn_propB_ykl,by = c("BIN_cm_mid","STRATUM"),all.x=T)
  tmp      <- merge(tmp,mn_propB_yk,      by = c("STRATUM","REGION","CN","SN"),all.x=T)
  tmp      <- merge(tmp,mn_propB_yl,      by = c("BIN_cm_mid","REGION","CN","SN"),all.x=T)
  
  sum_yl   <- sum(tmp["propB_yl"], na.rm=T)
  sum_yk   <- sum(tmp["propB_yk"], na.rm=T)
  sum_ykl  <- sum(tmp["propB_ykl"],na.rm=T)
  propB    <- merge(propB, tmp, 
                    by=c("YEAR",
                         "REGION",
                         "STRATUM",
                         "BIN_cm_mid",
                         "CN","SN",
                         "propB_ykl",
                         "propN_ykl",
                         "propB_yk",
                         "propN_yk",
                         "propB_yl",
                         "propN_yl"),all.x=T,all.y=T)
}

