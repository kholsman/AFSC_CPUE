# updateLW.R

#_____________________________________
#  update LWA regressions:
#_____________________________________

  source(file.path(username_path,"username_password.R"))

  if(R.Version()$arch=="i386"){
    con <-  odbcConnect("AFSC",username,password)
  }else{
    con <-  odbcConnect("AFSC",username,password,believeNRows=FALSE)
  }
# 
#  tt2 <- sqlTables(channel = con, schema = "FOODLAB")
#  tt  <- sqlFetch(channel = con, schema = "RACEBASE",sqtable = "HAUL", max = 2, rows_at_time = 1)
#  tt  <- sqlTables(channel = con, schema = "RACEBASE")
#  tt  <- sqlFetch(channel = con, schema = "foodlab",sqtable = "HAUL", max = 2, rows_at_time = 1)
#  

  cat("\n-- Get LWA_ALL\n")
  LWA_ALL     <- sqlQuery(con,qry_LWA_all)
  LW_ALL      <- sqlQuery(con,qry_LW_all)
  LWA_srvy    <- LWA_ALL%>%filter(CRUISE_TYPE =="Race_Groundfish")
  LW_srvy     <- LW_ALL%>%filter(CRUISE_TYPE =="Race_Groundfish")
  save(LWA_ALL,  file=file.path(data.path,"LWA_ALL.Rdata"))
  save(LWA_srvy, file=file.path(data.path,"LWA_srvy_noObs.Rdata"))
  save(LW_ALL,  file=file.path(data.path,"LW_ALL.Rdata"))
  save(LW_srvy, file=file.path(data.path,"LW_srvy_noObs.Rdata"))
  
  close(con)
  
  for(reg in c("BS","GOA","AI")){
    cat("-- create LWA_glms: ",reg,"\n")
    LWdata <-  createLWA_glms(
      data.path1    = data.path,
      #data.pathout1 = data.out,
      data.pathout1 = "data/in/lookup_files",
      reg_list      = reg,
      spcode_IN     = species_lkup$SPECIES_CODE,
      LWDATA        = "LW_srvy_noObs.Rdata",
      LWDATA_nm     = "LW_srvy",
      flname        = paste0("LWGlms_srvy_",reg,".Rdata"),
      saveit        = T)
  }
  cat("-- create LWA_glms\n")
  LWdata <-  createLWA_glms(
    data.path1    = data.path,
    #data.pathout1 = data.out,
    data.pathout1 = "data/in/lookup_files",
    reg_list      = c("BS","GOA","AI"),
    spcode_IN     = species_lkup$SPECIES_CODE,
    LWDATA        = "LW_srvy_noObs.Rdata",
    LWDATA_nm     = "LW_srvy",
    flname        = LWname,
    saveit        = T)
  
  cat("-- create LWA_glms\n")
  LWdata <-  createLWA_glms(
    data.path1    = data.path,
    #data.pathout1 = data.out,
    data.pathout1 = "data/in/lookup_files",
    reg_list      = c("BS","GOA","AI"),
    spcode_IN     = species_lkup$SPECIES_CODE,
    LWDATA        = "LWA_srvy_noObs.Rdata",
    LWDATA_nm     = "LWA_srvy",
    flname        = "LWAGlms_srvy_all_regs.Rdata",
    saveit        = T)

# load GLMS
load(file.path("data/in/lookup_files","LWAGlms_srvy_all_regs.Rdata"))
LWA.glm <- LW.glm
# load GLMS
load(file.path("data/in/lookup_files",LWname))

# load GLMS
load(file.path("data/in/lookup_files","LWGlms_srvy_BS.Rdata"))
LW_BS.glm <- LW.glm
# load GLMS
load(file.path("data/in/lookup_files",LWname))


plot(1:1000,exp(predict(LW_BS.glm$WALLEYE_POLLOCK,newdata=data.frame(L=1:1000))),type="l",col="red")
lines(1:1000,exp(predict(LW.glm$WALLEYE_POLLOCK,newdata=data.frame(L=1:1000))))

plot(1:1000,exp(predict(LW_BS.glm$ARROWTOOTH_FLOUNDR,newdata=data.frame(L=1:1000))),type="l",col="red")
lines(1:1000,exp(predict(LW.glm$ARROWTOOTH_FLOUNDR,newdata=data.frame(L=1:1000))))


cat("completed LW lookup table\n")
