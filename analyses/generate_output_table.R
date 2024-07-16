#read in range sizes
range_sizes <- readRDS(here::here("outputs", "species_range_sizes.rds"))

#read in gbif ids table
gbif_ids <- readRDS(here::here("outputs", "species_list_w_gbif_id.rds"))

#make output table
out_df <- gbif_ids

# NOTE:CHANGE NAME
out_df$range_size <- rep(NA,length(gbif_ids[,1]))

for(i in 1:length(range_sizes[,1])){
  
  out_df$range_size[grep(range_sizes$gbif_key[i],gbif_ids$gbif_key)] <- range_sizes$n_cells[i]
  
}

hist(out_df$range_size) 

####
# ---- AOO ----
####

aoo <- readRDS(paste(here::here(),"/outputs/","aoo.rds", sep=""))

out_df$aoo <- rep(NA, length(out_df$species_name))

for(i in 1:length(aoo$tax)){
  
  out_df$aoo[grep(aoo$tax[i], gbif_ids$gbif_key)] <- aoo$aoo[i]
  
}

hist(out_df$range_size)

####
# ---- Combined -----
####

#number of species with ranges
out_df_ranges_only <- out_df[!is.na(out_df$range_size) | !is.na(out_df$aoo),]
head(out_df_ranges_only)
str(out_df_ranges_only)

table(is.na(out_df_ranges_only$aoo))
table(is.na(out_df_ranges_only$range_size))

#NOTE: Some AOO appear to have failed when cells-based did not

#compare different range size estimates
plot(out_df$range_size~out_df$aoo)

#read in original species list
species_list <- read.csv(here::here("data", "species_list.csv"),row.names = 1)
str(species_list)
species_list$species_name <- gsub("_"," ",species_list$species_name)

#species with no range size
setdiff(species_list$species_name,out_df_ranges_only$species_name)

# CHECK THIS: 
# remove duplicates
out_df_ranges_only <- unique(out_df_ranges_only)

sort(table(out_df_ranges_only$species_name),decreasing = TRUE)
table(table(out_df_ranges_only$species_name))

#info per family
table(out_df_ranges_only$gbif_family)

#write table with range sizes
saveRDS(out_df_ranges_only, here::here("outputs", "species_with_ranges.rds"))

