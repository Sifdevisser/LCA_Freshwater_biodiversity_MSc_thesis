setwd('Desktop/R/')

library(data.table); library(dplyr); library(purrr); library(ggplot2); library(tidyr);
library(sf); library(raster); library(cowplot); library(ggpubr)

temp <- readRDS('proc/ssp/points_template.rds')
basins <- raster('proc/basins/basins_5min_pcrglobwb_adjusted.tif')
basin_ID <- raster::extract(basins, temp)
basin <- cbind(temp, basin_ID)
basin_size <- data.frame("area" = as.numeric(basin$area), "basin_ID" = as.character(basin$basin_ID)) %>%
  group_by(basin_ID) %>%
  summarise(total_area = sum(area), .groups = 'keep') 
basin_polygons <- st_read('proc/basins/basins_5min_pcrglobwb_adjusted.gpkg') %>%
  rename("basin_ID" = "ID") 

average_area <- readRDS('proc/Average_area.rds') 

# EB = (1/n) sum[1 - (Anew/Aoriginal)^z]
# 1) Calculate part EB = [1 - (Anew/Aoriginal)^z], 2) add column with 1 (sums up to number of species later)

z <- 0.21

EB_basin <- average_area %>%
  filter(!is.na(area_mean_hist) & !is.na(basin_ID))  %>%
  transmute(area_mean_hist = area_mean_hist,
            area_mean_1.5 = ifelse(area_mean_1.5 > area_mean_hist, area_mean_hist, area_mean_1.5),
            area_mean_2.0 = ifelse(area_mean_2.0 > area_mean_hist, area_mean_hist, area_mean_2.0),
            area_mean_3.2 = ifelse(area_mean_3.2 > area_mean_hist, area_mean_hist, area_mean_3.2),
            area_mean_4.5 = ifelse(area_mean_4.5 > area_mean_hist, area_mean_hist, area_mean_4.5)) %>%
  mutate(affected = ifelse(area_mean_4.5 < area_mean_hist, 'yes', 'no')) %>%
  mutate(partEB_1.5 = 1-(area_mean_1.5 / area_mean_hist)^z,
         partEB_2.0 = 1-(area_mean_2.0 / area_mean_hist)^z,
         partEB_3.2 = 1-(area_mean_3.2 / area_mean_hist)^z,
         partEB_4.5 = 1-(area_mean_4.5 / area_mean_hist)^z,
         n_1.5 = ifelse(is.na(partEB_1.5) == TRUE, 0, 1),
         n_2.0 = ifelse(is.na(partEB_2.0) == TRUE, 0, 1),
         n_3.2 = ifelse(is.na(partEB_3.2) == TRUE, 0, 1),
         n_4.5 = ifelse(is.na(partEB_4.5) == TRUE, 0, 1))
         
EB_basin2 <- EB_basin %>%
  group_by(basin_ID) %>%
  summarise(sum_partEB_1.5 = sum(partEB_1.5, na.rm = T), 
            sum_partEB_2.0 = sum(partEB_2.0, na.rm = T),
            sum_partEB_3.2 = sum(partEB_3.2, na.rm = T),
            sum_partEB_4.5 = sum(partEB_4.5, na.rm = T),
            sum_n_1.5 = sum(n_1.5),
            sum_n_2.0 = sum(n_2.0),
            sum_n_3.2 = sum(n_3.2),
            sum_n_4.5 = sum(n_4.5),
            .groups = 'keep') 

# Filter out small basins (<500km^2)
EB_basin3 <- inner_join(EB_basin2, basin_size, by = "basin_ID") %>%
  filter(total_area > 500) %>%
  transmute(EB_1.5 = 1/sum_n_1.5 * sum_partEB_1.5,
            EB_2.0 = 1/sum_n_2.0 * sum_partEB_2.0,
            EB_3.2 = 1/sum_n_3.2 * sum_partEB_3.2,
            EB_4.5 = 1/sum_n_4.5 * sum_partEB_4.5)

### Effect factor
# EF_marginal = slope @ current state
# EF_average = average between 1.2 and 4.5

EF_basin <- EB_basin3 %>%
  mutate(Marginal = EB_1.5/1.5,
         Average = (EB_4.5-(EB_1.5/1.5*1.2))/(4.5-1.2)) %>%
  transmute(Marginal = Marginal,
            Average = ifelse(Average>0, Average, 0))
EF_basin$Marginal[is.nan(EF_basin$Marginal)]<- NA

EF_basin$basin_ID <- as.integer(EF_basin$basin_ID)
EF_basin2 <- full_join(EF_basin, basin_polygons, by = 'basin_ID', all.y = TRUE)

# Weighted EF
basin_size$basin_ID <- as.integer(basin_size$basin_ID)
EF_basin3 <- inner_join(EF_basin2, basin_size, by = 'basin_ID')
EF_basin3 <- EF_basin3 %>% mutate(Marginal_W = Marginal * total_area,
                                  Average_W = Average * total_area)
Ma_EF_Wm <- sum(EF_basin3$Marginal_W, na.rm = T)/sum(EF_basin3$total_area)
Av_EF_Wm <- sum(EF_basin3$Average_W, na.rm = T)/sum(EF_basin3$total_area)

#### Results
cat("Reported ranges for the effect factors per river basins are:")
cat(scientific(range(EF_basin$Average, na.rm = T), digits = 3), "PDF/°C for the average approach")
cat(scientific(range(EF_basin$Marginal, na.rm = T), digits = 3), "PDF/°C for the marginal approach")
cat("Mean values are:")
cat(scientific(mean(EF_basin$Average, na.rm = T), digits = 3), "PDF/°C for the average approach")
cat(scientific(mean(EF_basin$Marginal, na.rm = T), digits = 3), "PDF/°C for the marginal approach")
cat("Weighted arithmetic mean values by river basin size are:")
cat(scientific(Av_EF_Wm, digits = 3), "PDF/°C for the average approach")
cat(scientific(Ma_EF_Wm, digits = 3), "PDF/°C for the marginal approach")
cat("Median values are:")
cat(scientific(median(EF_basin$Average, na.rm = T), digits = 3), "PDF/°C for the average approach")
cat(scientific(median(EF_basin$Marginal, na.rm = T), digits = 3), "PDF/°C for the marginal approach")

basins <- nrow(EF_basin2)

NA_A <- sum(is.na(EF_basin2$Average))/basins*100
NA_M <- sum(is.na(EF_basin2$Marginal))/basins*100
Z_A <- length(which(EF_basin2$Average==0))/sum(!is.na(EF_basin2$Average))*100
Z_M <- length(which(EF_basin2$Marginal==0))/sum(!is.na(EF_basin2$Marginal))*100

cat("It was not possible to derive effect factors for each river basin. For the average approach,", round(NA_A,digits=1),
    "% has a missing value, and for the marginal approach this is", round(NA_M,digits=1),"%.")
cat("From the basins which have an effect factor,", round(Z_A, digits = 1), "% and", round(Z_M, digits=1), "% have an 
    effect factor of zero (average and marginal approach respectively).")

### Plots

rm(temp, average_area, basin, basin_polygons, basin_size, basins, EB_basin, EB_basin2, EB_basin3, EF_basin) 

plot_a <- ggplot() + 
  geom_sf(data = EF_basin2, aes(geometry = geom, fill = Average), lwd = 0) +
  scale_fill_viridis_c(option = 'D', na.value = "grey", 
                       begin = 0,
                       end = 0.9,
                       direction = -1,
                       breaks = seq(0, 0.2, by = 0.02),
                       labels = c(0.00, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.80),
                       limits = c(0,0.2)) +
  theme(text = element_text(size = 12),
        panel.background = element_rect(fill = "white"),
        legend.position = "none",
        legend.key.width = unit(2, 'cm'),
        legend.title = element_blank())

plot_m <- ggplot() + 
  geom_sf(data = EF_basin2, aes(geometry = geom, fill = Marginal), lwd = 0) +
  scale_fill_viridis_c(option = 'D', na.value = "grey", 
                       begin = 0,
                       end = 0.9,
                       direction = -1,
                       breaks = seq(0, 0.2, by = 0.02),
                       labels = c(0.00, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.80),
                       limits = c(0,0.2),
                       oob = squish) +
  theme(text = element_text(size = 12),
        panel.background = element_rect(fill = "white"),
        legend.position = "none",
        legend.key.width = unit(2, 'cm'),
        legend.title = element_blank())

Fig4 <- ggarrange(plot_a, plot_m,
                  labels = c("Average", "Marginal"),
                  font.label = list(size = 10),
                  ncol = 1, nrow = 2,
                  common.legend = TRUE, legend="bottom")
  
ggsave(path = "Figures/", filename = "Fig4.png",Fig4, dpi = 800)



