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


#TEST ----
occ <- readRDS(here::here("data", "gbif", 
                          paste0(downloads_info[1, "download_key"], 
                                 "_clean.rds")))

library(adehabitatHR)

for (i in 1:nrow(downloads_info)) {
  
  cat(paste0("Calculate MCP for block - ", i, "\r"))
  
  ## Import species occurrences ----
  
  occ <- readRDS(here::here("data", "gbif", 
                            paste0(downloads_info[i, "download_key"], 
                                   "_clean.rds")))
  
  # remove species with less than 3 data points
  
  low_data <- names(table(occ$gbif_key)[table(occ$gbif_key)<5])
  occ <- occ[!occ$gbif_key%in%low_data,]
  
  ## Extract unique GBIF keys ----
  
  species <- unique(occ[ , "gbif_key", drop = TRUE])
  
  for(j in 1:length(species)){
    
    occ_tmp <- occ[occ$gbif_key==species[j],]
    
    occ1 <-occ[occ$gbif_key==unique(occ$gbif_key)[j],]
    # str(occ1)
    
    # adehabitatHR
    
    xysp <- sp::SpatialPoints(occ1[,2:3])
    
    # Set the coordinate reference system (CRS)
    # More information on CRS here: 
    # https://www.nceas.ucsb.edu/~frazier/RSpatialGuides/OverviewCoordinateReferenceSystems.pdf
    # The sample data are UTM points in WGS84 from zone 18N
    proj4string(xysp) <- CRS("+init=epsg:4326")
    
    cp<-adehabitatHR::mcp(xysp, percent = 100, )
    
    # https://cran.r-project.org/web/packages/changeRangeR/vignettes/singleSpeciesMetrics.pdf
    # changeRangeR
    eoo <- changeRangeR::mcp(occ1[,2:3], crs = "+proj=longlat +datum=WGS84 +no_defs")
    
    ## area is measured in metersˆ2
    area <- raster::area(eoo)/1000000
    
    # plot(cp)
    # plot(xysp, add=TRUE)
    
    if(j == 1){
      
      cp_tmp <- data.frame(species[j],cp$area)
      eoo_tmp <- data.frame(species[j], area)
      
    } else {
      
      cp_tmp <- rbind(cp_tmp,data.frame(species[j],cp$area))
      eoo_tmp <- rbind(eoo_tmp,data.frame(species[j],area))
      
    }
    
  }
  
  if(i == 1){
    
    cp_table <- cp_tmp
    eoo_table <- eoo_tmp
    
  } else {
    
    cp_table <- rbind(cp_table,cp_tmp)
    eoo_table <- rbind(eoo_table,eoo_tmp)
    
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
