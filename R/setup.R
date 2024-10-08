#'
#'setup.R
#'
#'Base script for setting up the queries and code
#'Kirstin.holsman@noaa.gov
#'
#'
  # Set switches for this code
  #----------------------------------------
  if(update_qrydate){
   qrydate <-  format(Sys.time(), "%Y_%m_%d")
   save(qrydate, file=file.path("data","in","lookup_files","qrydate.R"))
  }else{
    # qrydate       <- "2024_09_16" #"2024_08_14" #2023_09_26"  #"2023_03_02"
    load(file=file.path("data","in","lookup_files","qrydate.R"))# Set up data folders:
  }
  # Set up SQL stuff:
  #----------------------------------------
  username_path <- "R"   # copy the template from the data/in folder to not_shared
 
  #------------------------------------
  code.path   <- getwd()
  data.path   <- data.in   <- file.path("data","in",qrydate)
  
  data.out    <- file.path("data","out",qrydate)
  if(update_qrydate)
    if(dir.exists(data.out)) 
      file.rename(data.out,paste0(data.out,format(Sys.time(), "%Y_%m_%d_%H%M")))
  outfile.fig <- file.path("figs")
  
  
  if(!dir.exists("data"))dir.create("data")
  if(!dir.exists( data.in))dir.create( data.in)
  if(!dir.exists(file.path( data.in,"lookup_files"))) dir.create(file.path( data.in,"lookup_files"))
  # if(!dir.exists("data/in/newest"))dir.create("data/in/newest")
  if(!dir.exists(data.out))dir.create(data.out)
  
  if(!dir.exists(data.out))dir.create(data.out)
  if(!dir.exists(data.path))dir.create(data.path)
  
  if(!dir.exists(outfile.fig)) dir.create(outfile.fig)
  
 
  #walleye pollock, Pacific cod, arrowtooth flounder, sablefish Pacific ,halibut
  

 
  # if(!dir.exists(file.path(data.out,"propB"))) 
  #   dir.create(file.path(data.out,"propB"))
  

  cat("----------------------\n")
  cat("--     setup.R      --\n")
  cat("----------------------\n")

  cat("update_LWdata   = ",update_LWdata,"\n"  )
  cat("update_lkups    = ",update_lkups,"\n"  )
  cat("update_qrydate  = ",update_qrydate,"\n"  )
  cat("data.path = ",data.path,"\n"  )
  cat("data.out  = ",data.out,"\n"  )
  cat(" \n")
  
  
  # KEY : Specify the size bins & species and order for the query (get_CPUE_srvy.R) 
  #----------------------------------------
  subreg       <-  "GOA"		# can be "BS", "GOA", or "AI" 
  
  splist <- c(
    plk       = "walleye pollock",
    pcod      = "Pacific cod",
    atf       = "arrowtooth flounder",
    sablefish = "sablefish",
    halibut   = "Pacific halibut",
    yfs       = "yellowfin sole",
    nrs       = "northern rock sole"
    )
  # splist_n      <-  c(21740, 21720, 10110, 20510, 10120)
  # sp_bins <- list()
  # for(s in 1:length(splist))
  #   sp_bins[[ names(splist)[s]]] <-c(0,400,800,2000) # juvenile adult designations:
  # 
  sp_bins <- list()
  for(s in 1:length(splist))
    sp_bins[[ names(splist)[s]]] <- seq(0,1500,10)  # 1 cm bins
  
  NEBS_strata <- c(70,71,81)
  srvys <- data.frame(reg=c("sebs","nebs","goa","ai","slope"),num=c(98, 143,47,52,78)  )
  # surveys%>%filter(REGION=="GOA",SURVEY_NAME=="Gulf of Alaska Bottom Trawl Survey")%>%select(SURVEY_DEFINITION_ID)
    # survey = 143 Northern Bering Sea survey
    # survey = 98  Eastern Bering Sea Crab/Groundfish Bottom Trawl Survey
    # survey = 52  Aleutian Islands Bottom Trawl Survey
    # survey = 47  Gulf of Alaska Bottom Trawl Survey
    
    
  # General setup
  #----------------------------------------
 
  thisYr       <- as.numeric(format(Sys.time(), "%Y"))
  yrcutoff     <- thisYr-2
  
  # Define the seasons:
  #----------------------------------------
  summer              <-  6:9
  fall                <-  c(10,11)
  winter              <-  c(12,1,2)
  spring              <-  3:5
  
  # Set file names for the main Rdata files
  #----------------------------------------
  LWname              <-   "LWGlms_srvy_all_regs.Rdata" #"LWGlms_all.Rdata"
  

save.image(file.path(data.out,"setup.Rdata"))

