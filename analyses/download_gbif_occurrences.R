#some IDs are duplicates, so need to remove those, and then figure out the issues later
unique_ids <- unique(gbif_ids$gbif_key)

gbif_download_occ(unique_ids, by = as.integer(1000))
