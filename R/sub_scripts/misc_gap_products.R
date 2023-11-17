# get GAP Indices

library(devtools)
devtools::install_github("afsc-gap-products/gapindex")

library(gapindex)

regions <- c("AI" = 52, "GOA" = 47, "EBS" = 98, "BSS" = 78, "NBS" = 143)
data_tables <- c("CPUE", "BIOMASS", "SIZECOMP", "AGECOMP")
support_tables <- c("METADATA_COLUMN", "METADATA_TABLE")

sql_channel <- gapindex::get_connected()
username
password

suppressWarnings(channel <- RODBC::odbcConnect(dsn = paste("AFSC"), 
                                               uid = paste(username), 
                                               pwd = paste(password), believeNRows = FALSE))


tbls <- sqlTables(channel)
tbls_gap <- sqlTables(channel, schema  = "GAP_PRODUCTS")
db_list <- RODBC::sqlQuery(channel = channel, query = "SELECT * name FROM GAP_PRODUCTS")

for (idata in support_tables)
 aa<-RODBC::sqlQuery(channel = channel,
                                query = paste0("SELECT * FROM GAP_PRODUCTS.", 
                                               idata))

AREA <- RODBC::sqlQuery(channel = channel,
                           query = paste0("SELECT * FROM GAP_PRODUCTS.", 
                                          "AREA"))
AGECOMP <- RODBC::sqlQuery(channel = channel,
                           query = paste0("SELECT * FROM GAP_PRODUCTS.", 
                                          "AGECOMP"))
CPUE    <- RODBC::sqlQuery(channel = channel,
                           query = paste0("SELECT * FROM GAP_PRODUCTS.", 
                                          "CPUE"))
BIOMASS    <- RODBC::sqlQuery(channel = channel,
                           query = paste0("SELECT * FROM GAP_PRODUCTS.", 
                                          "BIOMASS"))
SIZECOMP    <- RODBC::sqlQuery(channel = channel,
                              query = paste0("SELECT * FROM GAP_PRODUCTS.", 
                                             "SIZECOMP"))

outdir<-file.path("data/out/2023_09_26/GAP_PRODUCTS")
if(!dir.exists(outdir))
   dir.create(outdir)
save(CPUE,file = file.path(outdir,"CPUE.Rdata"))
save(BIOMASS,file = file.path(outdir,"BIOMASS.Rdata"))
save(SIZECOMP,file = file.path(outdir,"SIZECOMP.Rdata"))
save(AGECOMP,file = file.path(outdir,"AGECOMP.Rdata"))


AGECOMP%>%filter(SPECIES_CODE==species_lkup$SPECIES_CODE[1],SURVEY_DEFINITION_ID %in%c(98,))


gapindex::get_data(year_set = 1982:thisYr,
                                      survey_set = "EBS",
                                      spp_codes = species_lkup$SPECIES_CODE[1],
                                      pull_lengths = TRUE,
                                      haul_type = 3,
                                      abundance_haul = "Y",
                                      sql_channel = channel)



