# make.R

source("R/packages.R")       # loads packages
source("R/setup.R")          # load other switches and controls
source("R/load_functions.R") # defines the create_plot() function
source("R/load_data.R")      # load other switches and controls

# options(clustermq.scheduler = "multicore") # optional parallel computing. Also needs parallelism = "clustermq"
# make(
#   plan, # defined in R/plan.R
#   verbose = 2
# )