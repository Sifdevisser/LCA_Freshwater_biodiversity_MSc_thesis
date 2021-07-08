setwd("~/Desktop/R")

library(raster); library(foreach); library(sf); library(dplyr); library(matrixStats); library(ggplot2); library(rgdal); library(data.table)

# Input e.g. 'gfdl'
species_files <- function(mod){
  files <- list.files(path = paste0('proc/',mod,'/modelled_occurrence'), recursive = 'TRUE', pattern = ".rds", full.names = "TRUE")
  return(files)
} 

read <- function(species_files){
  df <- readRDS(species_files) 
  dt <- tibble::rownames_to_column(df, "row_no")
  species_id <- data.frame("species_file" = basename(species_files))
  bind_cols(dt,species_id)
}

compile_table <- function(mod) {
  lapply(species_files(mod), read) %>% 
    bind_rows() 
}

# Input e.g. (gfdl (NO QUOTES!), "hist")
categorize_hist <- function(m, t) {
  columns <- colnames(m)
  a <- grep(t, columns)
  c_hist <- m %>% select('row_no', 'area', 'species_file', all_of(a))
  return(c_hist)
}

categorize_scen <- function(m, t) {
  columns <- colnames(m)
  te <- paste0(format(t, digits = 2),"C")
  a <- grep(te, columns)
  c_scen <- m %>% select('row_no', 'area', 'species_file', all_of(a))
  return(c_scen)
}

count_false <- function(s) {
  count <- rowCounts(s, value = FALSE, na.rm=TRUE, dim. = dim(s))
}

temp <- readRDS('proc/ssp/points_template.rds')
basins <- raster('proc/basins/basins_5min_pcrglobwb_adjusted.tif')
basin_ID <- raster::extract(basins, temp)
basin <- cbind(temp, basin_ID)
# used for function area_sum:
basi <- data.frame("row_no" = as.character(basin$row_no), "basin_ID" = as.character(basin$basin_ID))

rm(temp, basins, basin)

# e.g. s = gfdl_hist or gfdl_1.5
area_sum <- function(s) {
  cf <- data.frame("count_false" = count_false(as.matrix(s)))
  sd <- bind_cols(s, cf) %>%
    filter(count_false == 0) %>% 
    na.omit()
  se <- inner_join(basi,sd,by = 'row_no') %>%
    select('row_no', 'basin_ID', 'species_file', 'area') %>%
    group_by(basin_ID, species_file) %>%
    summarise(total_area = sum(area), .groups = 'keep') 
  return(se)
}

# GFDL
gfdl <- compile_table('gfdl')
gfdl_hist <- categorize_hist(gfdl, "hist") 
gfdl_1.5 <- categorize_scen(gfdl, 1.5)
rm(gfdl)
gfdl_area_hist <- area_sum(gfdl_hist)
gfdl_area_1.5 <- area_sum(gfdl_1.5)
rm(gfdl_hist, gfdl_1.5)

# HADGEM
hadgem <- compile_table('hadgem')
hadgem_hist <- categorize_hist(hadgem, "hist") 
hadgem_1.5 <- categorize_scen(hadgem, 1.5)
hadgem_2.0 <- categorize_scen(hadgem, 2.0)
hadgem_3.2 <- categorize_scen(hadgem, 3.2)
hadgem_4.5 <- categorize_scen(hadgem, 4.5)
rm(hadgem)
hadgem_area_hist <- area_sum(hadgem_hist)
hadgem_area_1.5 <- area_sum(hadgem_1.5)
hadgem_area_2.0 <- area_sum(hadgem_2.0)
hadgem_area_3.2 <- area_sum(hadgem_3.2)
hadgem_area_4.5 <- area_sum(hadgem_4.5)
rm(hadgem_hist, hadgem_1.5, hadgem_2.0, hadgem_3.2, hadgem_4.5)


# IPSL
ipsl <- compile_table('ipsl')
ipsl_hist <- categorize_hist(ipsl, "hist") 
ipsl_1.5 <- categorize_scen(ipsl, 1.5)
ipsl_2.0 <- categorize_scen(ipsl, 2.0)
ipsl_3.2 <- categorize_scen(ipsl, 3.2)
ipsl_4.5 <- categorize_scen(ipsl, 4.5)
rm(ipsl)
ipsl_area_hist <- area_sum(ipsl_hist)
ipsl_area_1.5 <- area_sum(ipsl_1.5)
ipsl_area_2.0 <- area_sum(ipsl_2.0)
ipsl_area_3.2 <- area_sum(ipsl_3.2)
ipsl_area_4.5 <- area_sum(ipsl_4.5)
rm(ipsl_hist, ipsl_1.5, ipsl_2.0, ipsl_3.2, ipsl_4.5)

# MIROC
miroc <- compile_table('miroc')
miroc_hist <- categorize_hist(miroc, "hist") 
miroc_1.5 <- categorize_scen(miroc, 1.5)
miroc_2.0 <- categorize_scen(miroc, 2.0)
miroc_3.2 <- categorize_scen(miroc, 3.2)
miroc_4.5 <- categorize_scen(miroc, 4.5)
rm(miroc)
miroc_area_hist <- area_sum(miroc_hist)
miroc_area_1.5 <- area_sum(miroc_1.5)
miroc_area_2.0 <- area_sum(miroc_2.0)
miroc_area_3.2 <- area_sum(miroc_3.2)
miroc_area_4.5 <- area_sum(miroc_4.5)
rm(miroc_hist, miroc_1.5, miroc_2.0, miroc_3.2, miroc_4.5)

# NORESM
noresm <- compile_table('noresm')
noresm_hist <- categorize_hist(noresm, "hist") 
noresm_1.5 <- categorize_scen(noresm, 1.5)
noresm_2.0 <- categorize_scen(noresm, 2.0)
noresm_3.2 <- categorize_scen(noresm, 3.2)
noresm_area_hist <- area_sum(noresm_hist)
noresm_area_1.5 <- area_sum(noresm_1.5)
noresm_area_2.0 <- area_sum(noresm_2.0)
noresm_area_3.2 <- area_sum(noresm_3.2)
rm(noresm_hist, noresm_1.5, noresm_2.0, noresm_3.2)

### Model average
hist <- Reduce(function(x,y) full_join(x=x, y=y, by = c("basin_ID", "species_file")), list(gfdl_area_hist, hadgem_area_hist, ipsl_area_hist, miroc_area_hist, noresm_area_hist)) %>%
  cbind("area_mean_hist" = rowMeans(.[,3:7], na.rm = TRUE)) %>%
  select('basin_ID', 'species_file', 'area_mean_hist')

t_1.5 <- Reduce(function(x,y) full_join(x=x, y=y, by = c("basin_ID", "species_file")), list(gfdl_area_1.5, hadgem_area_1.5, ipsl_area_1.5, miroc_area_1.5, noresm_area_1.5)) %>%
  cbind("area_mean_1.5" = rowMeans(.[,3:7], na.rm = TRUE)) %>% select('basin_ID', 'species_file', 'area_mean_1.5') %>%
  select('basin_ID', 'species_file', 'area_mean_1.5')

t_2.0 <- Reduce(function(x,y) full_join(x=x, y=y, by = c("basin_ID", "species_file")), list(hadgem_area_2.0, ipsl_area_2.0, miroc_area_2.0, noresm_area_2.0)) %>%
  cbind("area_mean_2.0" = rowMeans(.[,3:6], na.rm = TRUE)) %>% 
  select('basin_ID', 'species_file', 'area_mean_2.0')

t_3.2 <- Reduce(function(x,y) full_join(x=x, y=y, by = c("basin_ID", "species_file")), list(hadgem_area_3.2, ipsl_area_3.2, miroc_area_3.2, noresm_area_3.2)) %>%
  cbind("area_mean_3.2" = rowMeans(.[,3:6], na.rm = TRUE)) %>% 
  select('basin_ID', 'species_file', 'area_mean_3.2')

t_4.5 <- Reduce(function(x,y) full_join(x=x, y=y, by = c("basin_ID", "species_file")), list(hadgem_area_4.5, ipsl_area_4.5, miroc_area_4.5)) %>%
  cbind("area_mean_4.5" = rowMeans(.[,3:5], na.rm = TRUE)) %>% 
  select('basin_ID', 'species_file', 'area_mean_4.5')

F_table <- Reduce(function(x,y) full_join(x=x, y=y, by = c("basin_ID", "species_file")), list(hist, t_1.5, t_2.0, t_3.2, t_4.5))

# Used in Basin_EF script:
saveRDS(F_table, file = "Average_area.rds")

