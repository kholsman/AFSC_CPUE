#load_data.R
# 

# Load base files:
#----------------------
if(!file.exists("data/in/lookup_files/SPECIES.Rdata")|update_lkups){
  source(file.path(username_path,"username_password.R"))
  if(R.Version()$arch=="i386"){
    con <-  odbcConnect("AFSC",username,password)
  }else{
    con <-  odbcConnect("AFSC",username,password,believeNRows=FALSE)
  }
  makelkup(conIN = con, outfl = "data/in/lookup_files")
  close(con)
}

tmpfl <- "data/in/lookup_files"
load(file.path(tmpfl,"HAUL.Rdata"))
load(file.path(tmpfl,"NODC.Rdata"))
load(file.path(tmpfl,"SPECIES.Rdata"))
load(file.path(tmpfl,"stations.Rdata"))
load(file.path(tmpfl,"surveys.Rdata"))
load(file.path(tmpfl,"STRATA_AREA.Rdata"))

#load(file.path("data/in/lookup_files","STATION_LKUP_noObs.Rdata"))
#load(file.path("data/in/lookup_files","LWA_srvy.Rdata"))

sub_surveys     <- surveys%>%filter(SURVEY_DEFINITION_ID%in%srvys$num)
# region          <-  as.character(surveys[surveys$SURVEY_DEFINITION_ID==survey,][1,]$REGION[1])
# survey_name     <-  as.character(surveys[surveys$SURVEY_DEFINITION_ID==survey,][1,]$SURVEY_NAME[1])

if(update_lkups)
  source("R/sub_scripts/make_species_lkup.R")
load(file.path("data/in/lookup_files","species_lkup.Rdata"))

if(update_LWdata){
  source("R/sub_scripts/updateLW.R")
  source("R/sub_scripts/make_species_lkup.R")
}
load(file.path("data/in/lookup_files",LWname))



#splist_n <- species_lkup$SPECIES_CODE



