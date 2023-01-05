# Make species lookup table (key)
#----------------------

# load GLMS
load(file.path("data/in/lookup_files",LWname))

species_lkup     <- SPECIES%>%filter(COMMON_NAME%in%splist)%>%select(COMMON_NAME,SPECIES_CODE,SPECIES_NAME)
species_lkup     <- species_lkup[match(species_lkup$COMMON_NAME,  splist),]
species_lkup$sp  <- names(splist)
species_lkup$num <- 1:length(splist)
nn   <-  match(paste0(species_lkup$SPECIES_NAME," (",species_lkup$COMMON_NAME,")"),NODC$NAME)

species_lkup$ECOPATH_PRED <-   NODC[nn,"ECOPATH_PRED"]
species_lkup$NODC         <-   NODC[nn,"NODC"]
species_lkup$GOAPOLL_PRED <-   NODC[nn, "GOAPOLL_PRED"]

nn <- grep("walleye pollock",species_lkup[,1])
species_lkup$ECOPATH_PRED[nn] <- NODC[grep("pollock",NODC$NAME)[1],"ECOPATH_PRED"]
species_lkup$NODC[nn]         <- NODC[grep("pollock",NODC$NAME)[1],"NODC"]
species_lkup$GOAPOLL_PRED[nn] <- NODC[grep("pollock",NODC$NAME)[1],"GOAPOLL_PRED"]

species_lkup$LW_a  <- NA
species_lkup$LW_b  <- NA
species_lkup$df  <- NA
species_lkup$r2  <- NA
species_lkup$LWqrydate  <- NA
species_lkup$LWreg      <- NA

# for(s in 3:length(LW.glm)){
  
for(s in c("WALLEYE POLLOCK","PACIFIC COD","ARROWTOOTH FLOUNDR","SABLEFISH","PACIFIC HALIBUT")){
  #nm <- names(LW.glm)[s]
  #nm <- s
  nn <- which(species_lkup$GOAPOLL_PRED==s)
 # nn <- match(nm,gsub(" ","_",species_lkup$GOAPOLL_PRED))
  s2 <- gsub(" ","_",species_lkup$GOAPOLL_PRED[nn])

  
  species_lkup$LW_a[nn] <- as.numeric(exp(coef(LW.glm[[s2]])[1]))
  species_lkup$LW_b[nn] <- as.numeric(coef(LW.glm[[s2]])[2])
  species_lkup$df[nn] <-summary(LW.glm[[s2]])[[7]][2]
  species_lkup$r2[nn] <-summary(LW.glm[[s2]])[[8]]
  
  
  species_lkup$LWqrydate[nn]   <- LW.glm[["qrydate"]]
  species_lkup$LWreg[nn]       <- paste0(LW.glm[["regions"]],sep="",collapse=",")
}

species_lkup$LW_a[species_lkup$sp=="halibut"]<-3.139e-03
species_lkup$LW_b[species_lkup$sp=="halibut"]<-3.24
species_lkup$df[species_lkup$sp=="halibut"]   <- NA
species_lkup$r2[species_lkup$sp=="halibut"]   <- NA
species_lkup$LWqrydate[species_lkup$sp=="halibut"]   <- LW.glm[["qrydate"]]
species_lkup$LWreg[species_lkup$sp=="halibut"]       <- paste0(LW.glm[["regions"]],sep="",collapse=",")

# species_lkup$LW_a[species_lkup$sp=="sablefish"]<-4.74e-03
# species_lkup$LW_b[species_lkup$sp=="sablefish"]<-3.19
# species_lkup$LWqrydate[species_lkup$sp=="sablefish"]   <- LW.glm[["qrydate"]]
# species_lkup$LWreg[species_lkup$sp=="sablefish"]       <- paste0(LW.glm[["regions"]],sep="",collapse=",")

save(species_lkup,file=file.path("data/in/lookup_files",
                                 "species_lkup.Rdata"))
write.csv(species_lkup,file=file.path("data/in/lookup_files",
                                      "species_lkup.csv"))
