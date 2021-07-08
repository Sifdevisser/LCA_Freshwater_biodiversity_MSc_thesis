setwd("~/Desktop/R")

# Libraries
library(readr); library(dplyr); library(ggplot2); library(tidyr); library(tidyverse); 
library("writexl"); library(raster); library(viridis); library(sf); library(scales)

RT <- read_delim("extinctionrisk/RT_species.csv", ",", escape_double = FALSE, trim_ws = TRUE,
                 col_types = cols_only("binomial" = col_character(),
                                       "1.5C_no dispersal" = col_double(),
                                       "2.0C_no dispersal" = col_double(),
                                       "3.2C_no dispersal" = col_double(),
                                       "4.5C_no dispersal" = col_double()))
Area <- read_delim("extinctionrisk/species_traits.csv", ",", escape_double = FALSE, trim_ws = TRUE,
                   col_types = cols_only("binomial" = col_character(),
                                         "area" = col_double()))
Table <- inner_join(x= RT, y = Area, by = "binomial") %>%
  rename(RT_1.5 = '1.5C_no dispersal', RT_2.0 = '2.0C_no dispersal', RT_3.2 = '3.2C_no dispersal', RT_4.5 = '4.5C_no dispersal', Aoriginal = area) %>%
  mutate(Anew_1.5 = ((100 - RT_1.5)/100) * Aoriginal,
         Anew_2.0 = ((100 - RT_2.0)/100) * Aoriginal,
         Anew_3.2 = ((100 - RT_3.2)/100) * Aoriginal,
         Anew_4.5 = ((100 - RT_4.5)/100) * Aoriginal)

n = 11425
z = 0.21

### Since RT is maximally 100%, the assumption that Anew cannot be larger than Aoriginal is guaranteed. 

### E1 
# E1 = 1- (sumAnew / sumAoriginal)^z
sum_Aoriginal <- sum(Table$Aoriginal)
E1_1.5 <- 1 - (sum(Table$Anew_1.5) / sum_Aoriginal)^z
E1_2.0 <- 1 - (sum(Table$Anew_2.0) / sum_Aoriginal)^z
E1_3.2 <- 1 - (sum(Table$Anew_3.2) / sum_Aoriginal)^z
E1_4.5 <- 1 - (sum(Table$Anew_4.5) / sum_Aoriginal)^z

### E2
# E2 = 1 - {(1/n) [sum (Anew/Aoriginal) ]}^z

E2 <- Table %>% transmute(E2_fraction_1.5 = Anew_1.5 / Aoriginal,
                       E2_fraction_2.0 = Anew_2.0 / Aoriginal,
                       E2_fraction_3.2 = Anew_3.2 / Aoriginal,
                       E2_fraction_4.5 = Anew_4.5 / Aoriginal)

E2_1.5 <- 1 - ((1/n)*sum(E2$E2_fraction_1.5))^z
E2_2.0 <- 1 - ((1/n)*sum(E2$E2_fraction_2.0))^z
E2_3.2 <- 1 - ((1/n)*sum(E2$E2_fraction_3.2))^z
E2_4.5 <- 1 - ((1/n)*sum(E2$E2_fraction_4.5))^z

### EB
# EB = (1/n) sum[1 - (Anew/Aoriginal)^z]

EB <- Table %>% transmute(EB_fraction_1.5 = 1 - (Anew_1.5 / Aoriginal)^z, 
                          EB_fraction_2.0 = 1 - (Anew_2.0 / Aoriginal)^z,
                          EB_fraction_3.2 = 1 - (Anew_3.2 / Aoriginal)^z,
                          EB_fraction_4.5 = 1 - (Anew_4.5 / Aoriginal)^z)

EB_1.5 <- (1/n) * sum(EB$EB_fraction_1.5)
EB_2.0 <- (1/n) * sum(EB$EB_fraction_2.0)
EB_3.2 <- (1/n) * sum(EB$EB_fraction_3.2)
EB_4.5 <- (1/n) * sum(EB$EB_fraction_4.5)


### Combine in data.frame
E1_all <- c(0, E1_1.5, E1_2.0, E1_3.2, E1_4.5)
E2_all <- c(0, E2_1.5, E2_2.0, E2_3.2, E2_4.5)
EB_all <- c(0, EB_1.5, EB_2.0, EB_3.2, EB_4.5)
dT <- c(0, 1.5, 2.0, 3.2, 4.5)
Global_extinction <- data.frame(E1_all, E2_all, EB_all, dT)

#### EFFECT FACTOR ##############################################################

### Marginal effect factor
# EF_marginal = slope @ current state
# current state: temperature increase 1.2 degrees
# Derive slope between O and 1.5 degrees

Ma_EF_E1 <- E1_1.5 / 1.5
Ma_EF_E2 <- E2_1.5 / 1.5
Ma_EF_EB <- EB_1.5 / 1.5

### Average effect factor
# EF_average = (E(4.5) - E(1.2)) / (4.5 - 1.2)

E1_1.2 <- E1_1.5 / 1.5 * 1.2
E2_1.2 <- E2_1.5 / 1.5 * 1.2
EB_1.2 <- EB_1.5 / 1.5 * 1.2

Av_EF_E1 <- (E1_4.5-E1_1.2)/(4.5-1.2)
Av_EF_E2 <- (E2_4.5-E2_1.2)/(4.5-1.2)
Av_EF_EB <- (EB_4.5-EB_1.2)/(4.5-1.2)


###### PAF Effect factors
PAF_wtg <- function(file_wtg) {
  a = area(raster(res=1/12))
  plot(a)
  r <- raster(file_wtg)
  a_m <- mask(x = a,mask = extend(r,a))
  plot(a_m)
  paf = sum((r*a_m)[],na.rm=T)/sum(a_m[],na.rm=T)
  paf
}

PAF_1.5 <- PAF_wtg('PAF/PAF_no_dispersal_1.5.tif')
PAF_2.0 <- PAF_wtg('PAF/PAF_no_dispersal_2.0.tif')
PAF_3.2 <- PAF_wtg('PAF/PAF_no_dispersal_3.2.tif')
PAF_4.5 <- PAF_wtg('PAF/PAF_no_dispersal_4.5.tif')

PAF_global <- data.frame(
  wtg = c(0, 1.5, 2.0, 3.2, 4.5),
  PAF1 = c(0, PAF_1.5, PAF_2.0, PAF_3.2, PAF_4.5),
  PAF0.5 = c(0, PAF_1.5/2, PAF_2.0/2, PAF_3.2/2, PAF_4.5/2)
)

### Marginal effect factor
Ma_EF_PAF1 <- PAF_1.5 / 1.5
Ma_EF_PAF0.5 <- (PAF_1.5/2) / 1.5

### Average effect factor
# EF_average = (E(4.5) - E(1.2)) / (4.5 - 1.2)

PAF1_1.2 <- PAF_1.5 / 1.5 * 1.2
PAF0.5_1.2 <- (PAF_1.5/2) / 1.5 * 1.2

Av_EF_PAF1 <- (PAF_4.5-PAF1_1.2)/(4.5-1.2)
Av_EF_PAF0.5 <- ((PAF_4.5/2)-PAF0.5_1.2)/(4.5-1.2)

EF_g <- data.frame(
  row.names = c('average', 'marginal'),
  E1 = c(Av_EF_E1, Ma_EF_E1),
  E2 = c(Av_EF_E2, Ma_EF_E2),
  EB = c(Av_EF_EB, Ma_EF_EB),
  PAF1 = c(Av_EF_PAF1, Ma_EF_PAF1),
  PAF0.5 = c(Av_EF_PAF0.5, Ma_EF_PAF0.5)
)

### Fate factors
# First conversion: correct for Dutch excel notation (/100)
# Second conversion: change kton to kg (* 1E6)
FF <- read_delim(file = "hanafia/mix_fate.csv", delim = ";") %>%
  dplyr::select(Unit, TF_20, TF_100, TF_inf) %>%
  mutate(TF_20 = TF_20/100*1E6, TF_100 = TF_100/100*1E6, TF_inf = TF_inf/100*1E6) 

### Sets of CF
# Symbols: A = average, M = marginal, number = time horizon of fate factors
CF <- FF %>% mutate(
  A_20_E1 = FF$TF_20 * EF_g[1,1],
  A_20_E2 = FF$TF_20 * EF_g[1,2],
  A_20_EB = FF$TF_20 * EF_g[1,3],
  A_20_PAF1 = FF$TF_20 * EF_g[1,4],
  A_20_PAF0.5 = FF$TF_20 * EF_g[1,5],
  M_20_E1 = FF$TF_20 * EF_g[2,1],
  M_20_E2 = FF$TF_20 * EF_g[2,2],
  M_20_EB = FF$TF_20 * EF_g[2,3],
  M_20_PAF1 = FF$TF_20 * EF_g[2,4],
  M_20_PAF0.5 = FF$TF_20 * EF_g[2,5],
  A_100_E1 = FF$TF_100 * EF_g[1,1],
  A_100_E2 = FF$TF_100 * EF_g[1,2],
  A_100_EB = FF$TF_100 * EF_g[1,3],
  A_100_PAF1 = FF$TF_100 * EF_g[1,4],
  A_100_PAF0.5 = FF$TF_100 * EF_g[1,5],
  M_100_E1 = FF$TF_100 * EF_g[2,1],
  M_100_E2 = FF$TF_100 * EF_g[2,2],
  M_100_EB = FF$TF_100 * EF_g[2,3],
  M_100_PAF1 = FF$TF_100 * EF_g[2,4],
  M_100_PAF0.5 = FF$TF_100 * EF_g[2,5],
  A_inf_E1 = FF$TF_inf * EF_g[1,1],
  A_inf_E2 = FF$TF_inf * EF_g[1,2],
  A_inf_EB = FF$TF_inf * EF_g[1,3],
  A_inf_PAF1 = FF$TF_inf * EF_g[1,4],
  A_inf_PAF0.5 = FF$TF_inf * EF_g[1,5],
  M_inf_E1 = FF$TF_inf * EF_g[2,1],
  M_inf_E2 = FF$TF_inf * EF_g[2,2],
  M_inf_EB = FF$TF_inf * EF_g[2,3],
  M_inf_PAF1 = FF$TF_inf * EF_g[2,4],
  M_inf_PAF0.5 = FF$TF_inf * EF_g[2,5]
  ) %>%
  .[,-c(2:4)]

### Results 

### Reported ranges
# Appendix - all sets of CFs
write_xlsx(CF,"globalCF\\CF.xlsx")

# stats in results section
cat("The global characterization factors range from:")
cat(scientific(range(CF$A_100_EB), digits = 3),"PDF yr kg^-1 for the average approach")
cat(scientific(range(CF$M_100_EB), digits = 3),"PDF yr kg^-1 for the marginal approach")
cat("Mean values:")
cat(scientific(mean(CF$A_100_EB), digits = 3),"PDF yr kg^-1 for the average approach")
cat(scientific(mean(CF$M_100_EB), digits = 3),"PDF yr kg^-1 for the marginal approach")
cat("Median values:")
cat(scientific(median(CF$A_100_EB), digits = 3),"PDF yr kg^-1 for the average approach")
cat(scientific(median(CF$M_100_EB), digits = 3, ),"PDF yr kg^-1 for the marginal approach")


### Fig. 2 Boxplot - Average and marginal CF EB, 100 years
Fig2 <- ggplot(data = CF) +
  geom_boxplot(aes(y = A_100_EB, x=factor(1), fill = 'A')) +
  geom_boxplot(aes(y = M_100_EB, x=factor(2), fill = "B")) +
  scale_fill_viridis_d(begin = 0.4, end = 0.8) +
  labs(y = expression(paste("Characterization factor [PDF yr ", kg^-1, "]")), x = "Type") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  scale_x_discrete(labels = c("Average", "Marginal")) +
  theme_minimal() +
  theme(text = element_text(size = 12),
        panel.grid.major = element_line(color=NA),
        strip.background = element_rect('white'),
        strip.background.x = element_blank(),
        strip.background.y = element_blank(),
        strip.text = element_text(angle = 0, vjust = -1, size = 12),
        legend.position = "none")
Fig2

ggsave(path = "Figures/", filename = "Fig2.png",Fig2, dpi = 600)


### Fig. 3 Barchart - Comparing EF for all metrics
# reorder data.frame for barchart
EF_bar <- data.frame(
  Metrics = rep(c("E1", "E2", "EPAF0.5", "EB", "EPAF"),2),
  EF = c(Av_EF_E1, Av_EF_E2, Av_EF_PAF0.5, Av_EF_EB, Av_EF_PAF1,
         Ma_EF_E1, Ma_EF_E2, Ma_EF_PAF0.5, Ma_EF_EB, Ma_EF_PAF1),
  Type = c((rep("Average",5)),(rep("Marginal", 5)))
)

Fig3 <- ggplot(data=EF_bar, aes(x=reorder(Metrics, EF), y=EF, fill=Type)) +
  geom_bar(stat="identity", position=position_dodge()) +
  labs(y = expression(paste("Effect factor [PDF °",C^-1,"]")), x = "Extinction metric") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  scale_x_discrete(labels = c(expression("E"[1]), expression("E"[2]), expression("E"[PAF0.5]), expression("E"[B]), expression("E"[PAF1]))) +
  scale_fill_viridis_d(begin = 0.4, end = 0.8) +
  theme(legend.position = "bottom") +
  theme_minimal() +
  theme(text = element_text(size = 12),
        panel.grid.major = element_line(color=NA),
        strip.background = element_rect('white'),
        strip.background.x = element_blank(),
        strip.background.y = element_blank(),
        strip.text = element_text(angle = 0, vjust = -1, size = 12))
Fig3

ggsave(path = "Figures/", filename = "Fig3.png",Fig3, dpi = 500)

### Sensitivity analysis - Z component
z_sa = 0.25

E_sa <- Table %>% transmute(E_sa_fraction_1.5 = 1 - (Anew_1.5 / Aoriginal)^z_sa, 
                          E_sa_fraction_2.0 = 1 - (Anew_2.0 / Aoriginal)^z_sa,
                          E_sa_fraction_3.2 = 1 - (Anew_3.2 / Aoriginal)^z_sa,
                          E_sa_fraction_4.5 = 1 - (Anew_4.5 / Aoriginal)^z_sa)

E_sa_1.5 <- (1/n) * sum(E_sa$E_sa_fraction_1.5)
E_sa_2.0 <- (1/n) * sum(E_sa$E_sa_fraction_2.0)
E_sa_3.2 <- (1/n) * sum(E_sa$E_sa_fraction_3.2)
E_sa_4.5 <- (1/n) * sum(E_sa$E_sa_fraction_4.5)

Ma_EF_sa <- E_sa_1.5 / 1.5
Av_EF_sa <- (E_sa_4.5-(E_sa_1.5/1.5*1.2))/(4.5-1.2)

SA <- data.frame(
  row.names = c("Marginal", "Average"),
  z0.21 = c(Ma_EF_EB, Av_EF_EB),
  z0.25 = c(Ma_EF_sa, Av_EF_sa), 
  increase = c(((Ma_EF_sa - Ma_EF_EB)/Ma_EF_EB*100), ((Av_EF_sa - Av_EF_EB)/Av_EF_EB*100)))

cat("The percentual increase of taking a z component of 0.25 instead of 0.21 amounts to", 
    round(SA[1,3],digits=1), "and", round(SA[2,3],digits=1),"% for the marginal and average approach respectively")

### Appendix - regression

Figa <- ggplot(aes(x = dT, y = EB_all), data = Global_extinction) +
  geom_segment(aes(x = 0, xend = 1.5, y = 0, yend = 0.03586883, col = "Marginal"), linetype = "dashed") +
  geom_segment(aes(x = 1.2, xend = 4.5, y = 0.03586883/1.5*1.2, yend = 0.26319496, col = "Average"), linetype = "dashed") +
  geom_point() + 
  geom_point(aes(x = 1.2, y = (0.03586883/1.5)*1.2), shape=23, col = "red", fill = "red") +
  labs(color = 'Type effect factor', x = 'Temperature increase [°C]', y = 'Extinction risk [PDF]') +
  scale_x_continuous(breaks=seq(0,5.5,1)) +
  scale_color_viridis_d(option = "D", begin = 0.4, end = 0.8) +  
  theme_minimal() +
  theme(text = element_text(size = 12),
        panel.grid.major = element_line(color=NA),
        strip.background = element_rect('white'),
        strip.background.x = element_blank(),
        strip.background.y = element_blank(),
        strip.text = element_text(angle = 0, vjust = -1, size = 12),
        legend.position = "bottom",
        legend.title = element_blank())
Figa

ggsave(path = "Figures/", filename = "Figa.png",Figa, dpi = 500)





