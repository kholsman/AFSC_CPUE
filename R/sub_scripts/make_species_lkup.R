# Make species lookup table (key)
#----------------------
tmp<-read.csv ("data/not_shared/pollock_num_cpue.csv")
# load GLMS
load(file.path("data/in/lookup_files",LWname))

species_lkup <-data.frame(sp =splist, num=1:length(splist)) 
myfun<-function(x,y = SPECIES){
  # x is a list
  out <- suppressWarnings(cbind(sp = names(x)[1],spnm = x[1],num=1,y%>%filter(grepl(x[1],COMMON_NAME))))
  for(i in 2:length(x)){
   out <-  suppressWarnings(rbind(out,cbind(sp = names(x)[i],spnm = x[i],num=i,y%>%filter(grepl(x[i],COMMON_NAME)))))
  }
  return(out)
}

species_lkup_all            <- myfun(x=splist, y = SPECIES)%>%
  select(sp,spnm,num,COMMON_NAME,SPECIES_CODE,SPECIES_NAME)%>%
  ungroup()
species_lkup_all$NAME       <- paste0(species_lkup_all$SPECIES_NAME," (",species_lkup_all$COMMON_NAME,")")
species_lkup_all$LW_a       <-
species_lkup_all$LW_b       <-
species_lkup_all$df         <-
species_lkup_all$r2         <-
species_lkup_all$LWqrydate  <-
species_lkup_all$LWreg      <- NA

species_lkup     <- species_lkup_all[match(splist,species_lkup_all$COMMON_NAME),]

nn   <-  match(paste0(species_lkup$SPECIES_NAME," (",species_lkup$COMMON_NAME,")"),NODC$NAME)

sub <- data.frame(NODC[nn,]%>%select("NAME","ECOPATH_PRED","NODC","GOAPOLL_PRED"))

species_lkup <- species_lkup%>%left_join(sub)
species_lkup_all <- species_lkup_all%>%left_join(sub)

species_lkup_all <- species_lkup_all%>%left_join(species_lkup%>%
                                                   select("sp","NAME","ECOPATH_PRED","NODC","GOAPOLL_PRED"))



#for(s in c("WALLEYE POLLOCK","PACIFIC COD","ARROWTOOTH FLOUNDR","SABLEFISH","PACIFIC HALIBUT")){
for(s in species_lkup$GOAPOLL_PRED){
  nn <- which(species_lkup$GOAPOLL_PRED==s)
  s2 <- gsub(" ","_",species_lkup$GOAPOLL_PRED[nn])
  species_lkup$LW_a[nn] <- as.numeric(exp(coef(LW.glm[[s2]])[1]))
  species_lkup$LW_b[nn] <- as.numeric(coef(LW.glm[[s2]])[2])
  species_lkup$df[nn]          <- summary(LW.glm[[s2]])[[7]][2]
  species_lkup$r2[nn]          <- summary(LW.glm[[s2]])[[8]]
  species_lkup$LWqrydate[nn]   <- LW.glm[["qrydate"]]
  species_lkup$LWreg[nn]       <- paste0(LW.glm[["regions"]],sep="",collapse=",")
}

species_lkup$LW_a[species_lkup$sp=="halibut"] <- 3.139e-03
species_lkup$LW_b[species_lkup$sp=="halibut"] <- 3.24
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
save(species_lkup_all,file=file.path("data/in/lookup_files",
                                 "species_lkup_all.Rdata"))
write.csv(species_lkup,file=file.path("data/in/lookup_files",
                                      "species_lkup.csv"))
