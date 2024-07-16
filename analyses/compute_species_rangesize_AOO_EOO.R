#' Compute species range sizes (Area of occupancy and extent of occurrence)
#' 
#' @description
#' 
#' The `.rds` file is combined with any existing range size data and exported in the `outputs/` folder as 
#' `species_range_sizes.rds`.

# load libraries
library(GeoRange)
library(ConR)
library(ape)
library(ggplot2)

## Parameters ----

n_cores <- 1

## Import GBIF downloads metadata ----

downloads_info <- read.csv(here::here("data", "gbif", "gbif_requests_keys.csv"))

aoo_list <- list()

for (i in 1:nrow(downloads_info)) {
  
  cat(paste0("Calculate AOOs / EOOs for block - ", i, "\r"))
  
  ## Import species occurrences ----
  
  occ <- readRDS(here::here("data", "gbif", 
                            paste0(downloads_info[i, "download_key"], 
                                   "_clean.rds")))

  # remove species with less than 3 data points
  # low_data <- names(table(occ$gbif_key)[table(occ$gbif_key)<3])
  # occ <- occ[!occ$gbif_key%in%low_data,]
  
  # reorder columns to match input requirements
  
  occ <- occ[,c(3,2,1)]
  
  ## Extract unique GBIF keys ----
  
  species <- unique(occ[ , "gbif_key", drop = TRUE])
  
  for(j in 1:length(species)){
    
    occ_tmp <- occ[occ$gbif_key==species[j],]
    
    # cat(
    #   paste("Calculating EOO and AOO for", species[j]),
    #   "- Number",
    #   j,
    #   "of",
    #   length(species),
    #   "\n"
    # )
    
    # NOTE: MEMORY ISSUE
    # try(
    #   
    #   eoo_tmp <- EOO.computing(XY = occ_tmp, method.range = "convex.hull", export_shp = FALSE, show_progress = FALSE)
    #   
    #   h <- dismo::convHull(occ_tmp[,c(1:2)])
    #   
    # )
    # 
    # eoo <- rbind(eoo,eoo_tmp)
    
    try(
      aoo_tmp <- AOO.computing(occ_tmp,
                               cell_size_AOO = 50, #grid size in km - NOTE: 10 causes memory error
                               #nbe.rep.rast.AOO = 10, #could increase, was 30 before
                               show_progress = FALSE) 
    )
    
    if(i == 1 && j == 1){
      
      aoo <- aoo_tmp
      
    } else {
      
      aoo <- rbind(aoo,aoo_tmp)
      
    }

  }
  
}

## Compute range sizes ----

# range_sizes <- parallel::mclapply(1:length(occs), function(i) {
#   
#   data.frame("gbif_key" = names(occs)[i],
#              "n_cells"  = length(unique(occs[[i]])))
#   
# }, mc.cores = n_cores)
# 
# range_sizes <- do.call(rbind.data.frame, range_sizes)

## Export data ----

saveRDS(aoo, paste(here::here(),"/outputs/","aoo.rds", sep=""))
