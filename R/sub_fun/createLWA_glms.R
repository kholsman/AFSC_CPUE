#________________________________________________
# createLWA_glms
#________________________________________________
# Fits LWA glms for each spp
#
createLWA_glms<-function(
  data.path1    = data.in,
  data.pathout1 = data.out,
  reg_list      = "GOA",  #c("BS","AI","GOA")
  spcode_IN     =  splist_n,
  LWDATA        = "LWA_srvy.Rdata",
  LWDATA_nm     = "LWA_srvy",
  flname        = "LWGlms_srvy.Rdata",
  #flsubname     = "LWGlms_noObs.Rdata",
  saveit=T){
  
  require(dplyr)
  
  load(file.path(data.path1,LWDATA))
  eval(parse(text=paste("tmp <- ",LWDATA_nm,sep="")))
  
  LWdata  <-  data.frame(
    CRUISE     =tmp$CRUISE,
    START_TIME =tmp$START_TIME,
    BT         =tmp$GEAR_TEMPERATURE,
    SST        =tmp$SURFACE_TEMPERATURE,
    REGION     =tmp$REGION,
    STRATA     =tmp$STRATUM,
    SP_CODE    =tmp$SPECIES_CODE,
    SPECIES    =tmp$GOAPOLL_PRED,
    SEX        =tmp$SEX,
    AGE        =tmp$AGE,
    L          =tmp$LENGTH/10,
    W          =tmp$WEIGHT
  )
  rm(tmp)
  eval(parse(text=paste0("rm(",LWDATA_nm,")")))
  
  LWdata          <-  LWdata%>%filter(REGION%in%reg_list,SP_CODE%in%spcode_IN,W>0,L>0) 
  sp.names        <-  unique(LWdata$SPECIES)
  LW.glm          <-  list()
  LW.glm$qrydate  <-  qrydate
  LW.glm$regions  <-  reg_list
  #LW.glm$sp_code  <-  spcode_IN
  np              <-  length(sp.names)

  for(s in 1:np){
    dat2  <- LWdata%>%filter(SPECIES==sp.names[s],is.na(W)==FALSE,is.na(L)==FALSE)
    print(s)
    print(length(dat2[,1]))
    sname  <-""
    sname  <-strsplit(as.character(sp.names[s]),split="")[[1]]
    sname[sname==" "]<-"_"
    sname  <-paste(sname,sep="",collapse="")
    
    if(!is.na(mean(dat2$W))){
      print(paste("running LW glm for ",sname))
      eval(parse(text=paste("LW.glm$",sname,"<-lm(log(W)~log(L),data=dat2)",sep="")))
      
    }else{
      print(paste("skipping",sname))
    }
    
  }

  if(saveit) save("LW.glm",file=file.path(data.pathout1,flname))
  return(list(LWdata=LWdata[LWdata$W>0,])) #,LW.glm.sub=LW.glm.sub))
  
}

