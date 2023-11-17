#'
#'
#'
#'check_CPUE_dat.R
#'
#'


rm(list=ls())
qrydate <-  "2023_03_02"
library(dplyr)
library(ggplot2)
thisYr <- format(Sys.time(), "%Y")
today  <- format(Sys.time(), "%b %d, %Y")
#source("R/make.R") 
#load(paste0("data/out/",qrydate,"/cpue/ebs/ebs.srvy98.plk.cpue_data.Rdata"))

load("C:/Users/kirstin.holsman/Downloads/ebs.srvy98.plk.cpue_data.Rdata")

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



#Maurice;s code;


survey_file <- cpue_data
#load(survey_file)


load(paste0("data/in/",qrydate,"/ebs/plk/location_catch.Rdata"))
raw <- location_catch

load(paste0("data/in/",qrydate,"/ebs/plk/length.Rdata"))
raw_len <- length

## Filter to hauls where CPUE is non-zero but there are NA bin cpues
na_data <- cpue_data$CPUE_station_bin_yr |> arrange(YEAR, HAULJOIN) |> filter(haul_CPUE_KGKM2 > 0 & is.na(bin_CPUE_KGKM2))

na_data |> filter(BIN==25)


unique(na_data$propW) # In all cases, when bin CPUE is NA, propW is NA too

unique(table(na_data$HAULJOIN)) ## 151 --> all bins are NA for these hauls


unique(na_data$HAULJOIN)

stat <- "K-14"
stat <- "A-05"
stat <- "K-06"

subd <- dat%>%
  filter(YEAR==1982,STATIONID==stat)

na_sub <- unique(na_data$HAULJOIN)


head(subd)
raw%>%filter(YEAR==1982,HAULJOIN%in%unique(na_data$HAULJOIN))
raw_len%>%filter(YEAR==1982,HAULJOIN%in%unique(na_data$HAULJOIN))

dat%>%
  filter(YEAR==1982,STATIONID==stat)%>%
  group_by(YEAR,STATIONID)%>%
  summarise(sum(propW))

dat%>%
  filter(YEAR==1982,STATIONID==stat)%>%
  group_by(YEAR,STATIONID,haul_CPUE_KGKM2 )%>%
  summarise(sum(propW),sum(bin_CPUE_KGKM2))

ggplot(dat%>%filter(YEAR==1982,STATIONID==stat))+  
  geom_line(aes(x=BIN, y = propW, color = stat))+
  geom_point(aes(x=BIN, y = propW, color = stat))+ theme_minimal()


ggplot(dat%>%filter(YEAR==1982,STATIONID==stat))+  
  geom_line(aes(x=BIN, y = bin_CPUE_KGKM2, color = stat))+
  geom_point(aes(x=BIN, y = bin_CPUE_KGKM2, color = stat))+ theme_minimal()
