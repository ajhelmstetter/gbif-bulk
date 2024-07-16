#read in gbif ids table ----
gbif_ids <- readRDS(here::here("outputs", "species_list_w_gbif_id.rds"))

#some IDs are duplicates, so need to remove those, and then figure out the issues later
unique_ids <- na.omit(unique(gbif_ids$gbif_key))

gbif_download_occ(unique_ids, by = as.integer(1000))
