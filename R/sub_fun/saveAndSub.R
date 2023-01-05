#'
#'
#'
#'
#'

saveAndSub    <-function(nm,saveDB=FALSE,subyr=2011,path1){
  
  eval(parse(text=paste("tmp <- ",nm,sep="")))
  if(saveDB==TRUE)  save2DB(nm)
  
  txt  <-  file.path(path1,paste(nm,".Rdata",sep=""))
  eval(parse(text=paste("save(",nm,",file=txt)",sep="")))
  
  
  tt     <- tmp%>%dplyr::filter(CRUISE_TYPE=="Race_Groundfish")
  nmsub  <-  paste(nm,"_noObs",sep="")
  eval(parse(text=paste(nmsub,"<-tt",sep="")))
  
  txt    <-  file.path(path1,paste(nmsub,".Rdata",sep=""))
  eval(parse(text=paste("save(",nmsub,",file=txt)",sep="")))
  
  rm(rr)
  rr1    <-  which(tmp$CRUISE_TYPE=="Race_Groundfish" & tmp$YEAR <= subyr)
  rr2    <-  which(tmp$GOAPOLL_PRED=="ARROWTOOTH FLOUNDR"
                   |tmp$GOAPOLL_PRED=="WALLEYE POLLOCK"
                   |tmp$GOAPOLL_PRED=="PACIFIC COD"
                   |tmp$GOAPOLL_PRED=="PACIFIC HALIBUT")
  rr     <-  intersect(rr1,rr2)
  nmsub  <-  paste(nm,"_noObs",subyr,sep="")
  eval(parse(text=paste(nmsub,"<-tmp[rr,]",sep="")))
  txt    <-  file.path(path1,paste(nmsub,".Rdata",sep=""))
  eval(parse(text=paste("save(",nmsub,",file=txt)",sep="")))
  print(paste(nm,"complete"))
  
}