#'
#'setup.R
#'
#'Base script for setting up the queries and code
#'Kirstin.holsman@noaa.gov
#'
#'

  # Set up SQL stuff:
  #----------------------------------------
  username_path <- "R"   # copy the template from the data/in folder to not_shared
  qrydate       <-  "2023_01_03"
  #walleye pollock, Pacific cod, arrowtooth flounder, sablefish Pacific ,halibut
  
  # Set switches for this code
  #----------------------------------------
  update_LWdata       <-  TRUE  # set to 1 to update LW regressions
  update_lkups        <-  TRUE # if updating the lookuptables
  update_qrydate      <-  TRUE # if updating the lookuptables
  if(update_qrydate)
  qrydate       <-  format(Sys.time(), "%Y_%m_%d")
  
  # Set up data folders:
  #------------------------------------
  code.path   <- getwd()
  data.path   <- file.path("data/in",qrydate)
  data.out    <- file.path("data/out",qrydate)
  outfile.fig <- file.path("figs")
  
  
  if(!dir.exists("data"))dir.create("data")
  if(!dir.exists("data/in"))dir.create("data/in")
  if(!dir.exists("data/in/lookup_files"))dir.create("data/in/lookup_files")
  # if(!dir.exists("data/in/newest"))dir.create("data/in/newest")
  if(!dir.exists("data/out"))dir.create("data/out")
  
  if(!dir.exists(data.out))dir.create(data.out)
  if(!dir.exists(data.path))dir.create(data.path)
  
  if(!dir.exists(outfile.fig)) dir.create(outfile.fig)
  # 
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
  
  srvys <- data.frame(reg=c("ebs","goa","ai"),num=c(98,47,52)  )
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

