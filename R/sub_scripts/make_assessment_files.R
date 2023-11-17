#'
#'
#'
#'make_assessment_files.R
#'
#'



thisYr <- format(Sys.time(), "%Y")
today  <- format(Sys.time(), "%b %d, %Y")
source("R/make.R") 
ss     <- 1
bins   <- sp_bins[[ species_lkup[ss,]$sp ]]
reg <- c(BS = "ebs")
load(paste0("data/out/",qrydate,"/cpue/",reg,"/",reg,".srvy98.",species_lkup[ss,]$sp,".cpue_data.Rdata"))
load(paste0("data/in/",qrydate,"/LWA_srvy_noObs.Rdata"))
head(LWA_srvy)
rm(s)

get_time<- function(x, nm = "year",formatIN = "%Y-%m-%d %H:%M:%S"){
  a <- strptime(x, format = formatIN)
  if(is.null(nm)){
    out <- a
  }else{
    out <- eval(parse(text = paste0("a$",nm)))
  }
  return(out)
}


LWA <- LWA_srvy%>%filter(
  SPECIES_CODE == species_lkup[ss,]$SPECIES_CODE,
  REGION  == names(reg))%>%
  mutate(
  TIME = get_time(x= START_TIME,nm=NULL),
  YEAR = get_time(x= START_TIME,nm="year")+1900,
  MONTH = get_time(x= START_TIME,nm="mon")+1,
  DAY = get_time(x= START_TIME,nm="mday"),
  BIN = getBIN2(x=LENGTH, BIN_IN=bins,divid=1,type=1),
  BIN_mm =  getBIN2(x=LENGTH, BIN_IN=bins,divid=1,type=2))

# now get pdf of Age Length by strata and year


age_comp_yr_strata_bin <- LWA%>%
  group_by(REGION,STRATUM,SPECIES_CODE,YEAR,AGE,BIN, BIN_mm)%>%
  summarize(nLen  = length.na(LENGTH))

age_comp_yr_strata <- LWA%>%
  group_by(REGION,STRATUM,SPECIES_CODE,YEAR,BIN, BIN_mm)%>%
  summarize(tot  = length.na(LENGTH))

age_comp_yr_strata_bin <- age_comp_yr_strata_bin%>%
  left_join(age_comp_yr_strata)%>%
  mutate(age_prop = nLen/tot)

test<-age_comp_yr_strata_bin%>%
  group_by(REGION,STRATUM,SPECIES_CODE,YEAR,BIN, BIN_mm)%>%
  summarize(sum  = sum(age_prop,na.rm=T))
unique(test$sum)
# now combine with propW by strata 

age_comp_yr_strata_bin <- age_comp_yr_strata_bin%>%
  left_join(cpue_data$propByStrataBin)%>%
  filter(!is.na(SN ))



# now roll up into year basin age_comp
age_comp_yr_bin <- age_comp_yr_strata_bin%>%
  group_by(REGION, SPECIES_CODE,  YEAR,AGE,BIN,BIN_mm,CN,SN,sp)%>%
  summarize(age_prop = sum(age_prop *propB_ykl ),
            sumpropB = sum(propB_ykl ))%>%
  mutate(age_len_comp = age_prop/sumpropB)%>%
  arrange(REGION, SN, YEAR, AGE, BIN)


age_comp_bin <- age_comp_yr_bin%>%
  group_by(REGION, SPECIES_CODE,YEAR,BIN,BIN_mm,CN,SN,sp)%>%
  summarize(tot = sum(age_len_comp, na.rm=T))


age_comp_yr_bin<- age_comp_yr_bin%>%left_join(age_comp_bin)%>%
  group_by(REGION, SPECIES_CODE,YEAR,AGE,BIN,BIN_mm,CN,SN,sp)%>%
  mutate(prop_age = age_len_comp/tot)


mn_age_comp_bin <- age_comp_yr_bin%>%
  group_by(REGION, SPECIES_CODE,AGE,BIN,BIN_mm,CN,SN,sp)%>%
  summarize(prop_age = mean(prop_age, na.rm=T))

mn_comp_bin <- mn_age_comp_bin%>%
  group_by(REGION, SPECIES_CODE,BIN,BIN_mm,CN,SN,sp)%>%
  summarize(tot = sum(prop_age, na.rm=T))

mn_age_comp_bin<- mn_age_comp_bin%>%left_join(mn_comp_bin)%>%
  group_by(REGION, SPECIES_CODE,AGE,BIN,BIN_mm,CN,SN,sp)%>%
  mutate(prop_age = prop_age/tot)


ggplot(data = age_comp_yr_bin)+
  geom_point(aes(x=BIN,y= AGE+prop_age,color=factor(AGE)))+
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3),aes(x=BIN,y= AGE+prop_age,fill =factor(AGE), color=factor(AGE)))+
  theme_minimal()

ggplot(data = age_comp_yr_bin)+
  geom_point(aes(x=BIN,y= prop_age,color=factor(AGE)))+
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3),aes(x=BIN,y= prop_age,fill =factor(AGE), color=factor(AGE)))+
  theme_minimal()


test<-age_comp_yr_strata_bin%>%
  group_by(REGION,STRATUM,SPECIES_CODE,YEAR,AGE)%>%
  summarize(sum  = sum(age_prop,na.rm=T))



tmp    <- getBIN2(x=bins,BIN_IN = bins,divid=1,type=3)
binALL <- data.frame(BIN = unlist(tmp[1,]),BIN_mm = unlist(tmp[2,]))
rm(tmp)




age_comp <- LWA%>%group_by(REGION,STRATUM,SPECIES_CODE,YEAR,AGE)%>%
  summarize(mnLen = mean(LENGTH, na.rm=T),
            sdLen = sd(LENGTH, na.rm=T),
            nLen  = length.na(LENGTH))


dat <- cpue_data$CPUE_station_bin_yr

dat%>%
  filter(YEAR==1982,STATIONID=="A-05")%>%
  group_by(YEAR,STATIONID)%>%
  summarise(sum(propW))

dat%>%
  filter(YEAR==1982,STATIONID=="A-05")%>%
  group_by(YEAR,STATIONID,haul_CPUE_KGKM2 )%>%
  summarise(sum(propW),sum(bin_CPUE_KGKM2))

ggplot(dat%>%filter(YEAR==1982,STATIONID=="A-05"))+  
  geom_line(aes(x=BIN, y = propW, color = "A-05"))+
  geom_point(aes(x=BIN, y = propW, color = "A-05"))+ theme_minimal()


ggplot(dat%>%filter(YEAR==1982,STATIONID=="A-05"))+  
  geom_line(aes(x=BIN, y = bin_CPUE_KGKM2, color = "A-05"))+
  geom_point(aes(x=BIN, y = bin_CPUE_KGKM2, color = "A-05"))+ theme_minimal()

