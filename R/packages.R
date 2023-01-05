# ----------------------------------------
# packages.R
# load or install packages
# kirstin.holsman@noaa.gov
# updated 2021
# ----------------------------------------

lib_list <- c(
  # these for reshaping and manipulating data:
  "ncdf4",
  "RODBC",
  "devtools",
  "svMisc",
  #"magrittr",
  #"httr",
  "reshape",
  "reshape2",
  "dplyr", 
  "purrr",
  "readxl", 
  #"tidyverse",
#  "usethis",
  
  # # these for ggplot mapping:
  #   "raster",
  #   "ggspatial",             # used for N arrow and scale bar 
  #   "sf",                    # used for shapefiles
  #   "rnaturalearth",         # has more shapefiles; used to make the "world" object 
  #   "maps",                  # has some state shapefiles, need to be converted with st_as_sf
  #   "akima",                 # Interpolation of Irregularly and Regularly Spaced Data
  # 
  
  # markdown stuff:
  "knitr",
  "kableExtra",
  
  # These for making plots:
  "RColorBrewer",
  "ggplot2", 
  "mgcv",
  "cowplot",               # 
  "wesanderson",
  #"scales",
  #"ggforce",
  #"grid",
  #"processx",
  "plotly",
  "extrafont"
)

# Install missing libraries:
missing <- setdiff(lib_list, installed.packages()[, 1])
if (length(missing) > 0) install.packages(missing)

# Load libraries:
for(lib in lib_list)
  eval(parse(text=paste("library(",lib,")")))

# ## same for git libraries
# lib_list_git <- c(
#   #  "rnaturalearthhires",
#   "thredds"
# )
# 
# 
# missing <- setdiff(lib_list_git, installed.packages()[, 1])
# if (length(missing) > 0) devtools::install_github("bocinsky/thredds")
