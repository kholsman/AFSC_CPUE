##################################################
# Rcode to query RACE database
# Kirstin.holsman@noaa.gov
# updated 2023
##################################################


# ************************
# ************************

require(RODBC)
source(file.path(username_path,"username_password.R"))

#require(data.frame)
qrydate  <-  format(Sys.time(), "%Y_%m_%d")

pswdTEST  <-  FALSE

if(R.Version()$arch=="i386"){
  con <-  odbcConnect("AFSC",username,password)
}else{
  con <-  odbcConnect("AFSC",username,password,believeNRows=FALSE)
}

testT  <-  odbcGetInfo(con)

if(any(ls()=="testT")) 
  pswdTEST<-TRUE

if(pswdTEST){ 
  cat("Successfully tested connection to RACEBASE\n")
  }else{
  stop("could not connect to RACEBASE\n")
}

# if(pswdTEST)
#   source(file.path(code.path,"R/sub_scripts/sub_runRACE_qrys.R"))
# 


cat("______________________________________________\n")
cat("______________________________________________\n")
cat("running SQL Queries - takes 10-15 mins to run\n")
cat("______________________________________________\n")
cat("______________________________________________\n")


####################################################################
## connect to database
####################################################################

con <- odbcConnect("AFSC",uid=username,pwd=password)   # update this in setup!
odbcGetInfo(con)

cat("--> connected to AFSC Database <--\n")


  tables   <- sqlTables(con,schema="FOODLAB")$TABLE_NAME
  # tables2   <- sqlTables(con,schema=RACE_schema)$TABLE_NAME
 
  cat("   -- Create data tables\n")
 
  # Create data tables (not saved to Oracle for now)
  #---------------------------------------------------

  for(r in 1:dim(srvys)[1]){
    sub_dir <- file.path(data.path,srvys$reg[r])
    if(!dir.exists(sub_dir)) dir.create(sub_dir)
    
    for(spn in 1:length(splist)){
      splk <- species_lkup[species_lkup$spnm==splist[spn],]
      cat(" -- ",srvys$reg[r],": ",splk$sp[1],"  \n")
      sp_dir <- file.path(sub_dir,splk$sp[1])
      if(!dir.exists(sp_dir)) dir.create(sp_dir)
      
      out_nm         <- paste0(".",splk$sp[1],".",srvys$reg[r])
      location       <- data.frame(sqlQuery(con,location_sql(surveyIN = srvys$num[r] )))
      location_catch <- data.frame(sqlQuery(con,location_catch_sql(surveyIN = srvys$num[r],
                                                               speciesIN = splk$SPECIES_CODE)))
      length         <- data.frame(sqlQuery(con,length_sql(surveyIN = srvys$num[r],
                                                           speciesIN = splk$SPECIES_CODE )))
      save(location,       file=file.path(sp_dir,"location.Rdata"))
      save(location_catch, file=file.path(sp_dir,"location_catch.Rdata"))
      save(length,         file=file.path(sp_dir,"length.Rdata"))
      rm(list=c("location","location_catch","length"))
     
    }
  }

  if(1==10){
    SRVY_BIOM  <-  sqlQuery(con,qry_SRVY_BIOM)
    save(SRVY_BIOM,file=file.path(data.path,"SRVY_BIOM.Rdata"))
    saveAndSub("SRVY_BIOM",saveDB=TRUE,path1=file.path(data.path))
  }
  
  # LengthComp<-(sqlQuery(con,"SELECT * FROM racebase.LENGTH"))	

  ####################################################################
  ## Now copy files
  ####################################################################
  cat(" -- Copy n save \n")
  cat("--> queries complete; copying files and wrapping up <--\n")
  cat("______________________________________________\n")
  cat("______________________________________________\n")
  
  unlink("data/in/Newest", recursive = TRUE)
  # 
  dir.create("data/in/Newest")
  # dir.create(file.path(data.path,"Data_qrys","Newest","AFSC_DB_QRYS"))
  # 
  file.copy(file.path(data.path), 
            "data/in/Newest", 
            overwrite = TRUE, recursive=TRUE)
close(con)


cat("--> all done <--\n")
cat("______________________________________________\n")
cat("______________________________________________\n")