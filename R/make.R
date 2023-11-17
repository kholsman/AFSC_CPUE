# make.R

source("R/packages.R")       # loads packages
source("R/setup.R")          # load other switches and controls
cat("----------------------\n")
cat("   --  setup complete \n")
source("R/load_functions.R") # defines the create_plot() function
cat("   --  load functions complete \n")
source("R/load_data.R")      # load other switches and controls
cat("   --  load data complete \n")

# options(clustermq.scheduler = "multicore") # optional parallel computing. Also needs parallelism = "clustermq"
# make(
#   plan, # defined in R/plan.R
#   verbose = 2
# )

# to change the password: password *new password here*