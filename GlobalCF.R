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

# Current temperature increase
GMST <- data.frame(
  year = c(2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020),
  GMST = c(0.81, 0.85, 0.89, 0.95, 1.10, 1.23, 1.13, 1.05, 1.19, 1.21)
)
CT <- mean(GMST$GMST) 

### Marginal effect factor
# EF_marginal = slope @ current state
# current state: temperature increase 0.92 degrees (CT)
# Derive slope between O and 1.5 degrees

Ma_EF_E1 <- E1_1.5 / 1.5
Ma_EF_E2 <- E2_1.5 / 1.5
Ma_EF_EB <- EB_1.5 / 1.5

### Average effect factor
# EF_average = (E(4.5) - E(CT)) / (4.5 - CT)

E1_CT <- E1_1.5 / 1.5 * CT
E2_CT <- E2_1.5 / 1.5 * CT
EB_CT <- EB_1.5 / 1.5 * CT

Av_EF_E1 <- (E1_4.5-E1_CT)/(4.5-CT)
Av_EF_E2 <- (E2_4.5-E2_CT)/(4.5-CT)
Av_EF_EB <- (EB_4.5-EB_CT)/(4.5-CT)


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
# EF_average = (E(4.5) - E(CT)) / (4.5 - CT)

PAF1_CT <- PAF_1.5 / 1.5 * CT
PAF0.5_CT <- (PAF_1.5/2) / 1.5 * CT

Av_EF_PAF1 <- (PAF_4.5-PAF1_CT)/(4.5-CT)
Av_EF_PAF0.5 <- ((PAF_4.5/2)-PAF0.5_CT)/(4.5-CT)

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
# Second conversion: change kton to kg (1E6)
FF <- read_delim(file = "hanafia/mix_fate.csv", delim = ";") %>%
  dplyr::select(Unit, TF_20, TF_100, TF_inf) %>%
  mutate(TF_20 = TF_20/100/1E6, TF_100 = TF_100/100/1E6, TF_inf = TF_inf/100/1E6) 

### Sets of CF
# Symbols: A = average, M = marginal, number = time horizon of fate factors

CF_A_20 <- FF %>% transmute(
  Unit = Unit,
  A_20_E1 = FF$TF_20 * EF_g[1,1],
  A_20_E2 = FF$TF_20 * EF_g[1,2],
  A_20_EB = FF$TF_20 * EF_g[1,3],
  A_20_PAF1 = FF$TF_20 * EF_g[1,4],
  A_20_PAF0.5 = FF$TF_20 * EF_g[1,5]
  ) 
CF_A_20[, 2:6] <- signif(CF_A_20[, 2:6], digits = 3)  
write.csv(CF_A_20, 'Final/CF_A_20.csv')

CF_M_20 <- FF %>% transmute(
  Unit = Unit,
  M_20_E1 = FF$TF_20 * EF_g[2,1],
  M_20_E2 = FF$TF_20 * EF_g[2,2],
  M_20_EB = FF$TF_20 * EF_g[2,3],
  M_20_PAF1 = FF$TF_20 * EF_g[2,4],
  M_20_PAF0.5 = FF$TF_20 * EF_g[2,5]
)
CF_M_20[, 2:6] <- signif(CF_M_20[, 2:6], digits = 3)  
write.csv(CF_M_20, 'Final/CF_M_20.csv')

CF_A_100 <- FF %>% transmute(
  Unit = Unit,
  A_100_E1 = FF$TF_100 * EF_g[1,1],
  A_100_E2 = FF$TF_100 * EF_g[1,2],
  A_100_EB = FF$TF_100 * EF_g[1,3],
  A_100_PAF1 = FF$TF_100 * EF_g[1,4],
  A_100_PAF0.5 = FF$TF_100 * EF_g[1,5]
)
CF_A_100[, 2:6] <- signif(CF_A_100[, 2:6], digits = 3)  
write.csv(CF_A_100, 'Final/CF_A_100.csv')

CF_M_100 <- FF %>% transmute(
  Unit = Unit,
  M_100_E1 = FF$TF_100 * EF_g[2,1],
  M_100_E2 = FF$TF_100 * EF_g[2,2],
  M_100_EB = FF$TF_100 * EF_g[2,3],
  M_100_PAF1 = FF$TF_100 * EF_g[2,4],
  M_100_PAF0.5 = FF$TF_100 * EF_g[2,5]
)
CF_M_100[, 2:6] <- signif(CF_M_100[, 2:6], digits = 3)  
write.csv(CF_M_100, 'Final/CF_M_100.csv')

CF_A_inf <- FF %>% transmute(
  Unit = Unit,
  A_inf_E1 = FF$TF_inf * EF_g[1,1],
  A_inf_E2 = FF$TF_inf * EF_g[1,2],
  A_inf_EB = FF$TF_inf * EF_g[1,3],
  A_inf_PAF1 = FF$TF_inf * EF_g[1,4],
  A_inf_PAF0.5 = FF$TF_inf * EF_g[1,5]
)
CF_A_inf[, 2:6] <- signif(CF_A_inf[, 2:6], digits = 3)  
write.csv(CF_A_inf, 'Final/CF_A_inf.csv')

CF_M_inf <- FF %>% transmute(
  Unit = Unit,
  M_inf_E1 = FF$TF_inf * EF_g[2,1],
  M_inf_E2 = FF$TF_inf * EF_g[2,2],
  M_inf_EB = FF$TF_inf * EF_g[2,3],
  M_inf_PAF1 = FF$TF_inf * EF_g[2,4],
  M_inf_PAF0.5 = FF$TF_inf * EF_g[2,5]
) 
CF_M_inf[, 2:6] <- signif(CF_M_inf[, 2:6], digits = 3)  
write.csv(CF_M_inf, 'Final/CF_M_inf.csv')

### Results 

# stats in results section
cat("The global characterization factors range from:\n")
cat(scientific(range(CF_A_100$A_100_EB), digits = 3),"PDF yr kg^-1 for the average approach\n")
cat(scientific(range(CF_M_100$M_100_EB), digits = 3),"PDF yr kg^-1 for the marginal approach\n")
cat("Mean values:\n")
cat(scientific(mean(CF_A_100$A_100_EB), digits = 3),"PDF yr kg^-1 for the average approach\n")
cat(scientific(mean(CF_M_100$M_100_EB), digits = 3),"PDF yr kg^-1 for the marginal approach\n")
cat("Median values:\n")
cat(scientific(median(CF_A_100$A_100_EB), digits = 3),"PDF yr kg^-1 for the average approach\n")
cat(scientific(median(CF_M_100$M_100_EB), digits = 3, ),"PDF yr kg^-1 for the marginal approach\n")


### Fig. 2 Boxplot - Average and marginal CF EB, 100 years
CF <- inner_join(x = CF_M_100, y = CF_A_100, by = 'Unit')
Fig2 <- ggplot(data = CF) +
  geom_boxplot(aes(y = A_100_EB, x=factor(1)), fill = '#2c7c8c') +
  geom_boxplot(aes(y = M_100_EB, x=factor(2)), fill = '#7cd454') +
  labs(y = expression(paste("Characterization factor [PDF·yr·",kg^-1,"]")), x = "Type") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  scale_x_discrete(labels = c("Average", "Marginal")) +
  theme_minimal() +
  theme(text = element_text(size = 10),
        axis.text=element_text(size=10),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.text = element_text(angle = 0, vjust = -1, size = 10))
Fig2

ggsave(path = "Figures/", filename = "Fig2.png",Fig2, width = 86, height = 86, unit = 'mm', dpi = 600)


### Fig. 4 Barchart - Comparing EF for all metrics
# reorder data.frame for barchart
EF_bar <- data.frame(
  Metrics = rep(c("E1", "E2", "EPAF0.5", "EB", "EPAF"),2),
  EF = c(Av_EF_E1, Av_EF_E2, Av_EF_PAF0.5, Av_EF_EB, Av_EF_PAF1,
         Ma_EF_E1, Ma_EF_E2, Ma_EF_PAF0.5, Ma_EF_EB, Ma_EF_PAF1),
  Type = c((rep("Average",5)),(rep("Marginal", 5)))
)

Fig4 <- ggplot(data=EF_bar, aes(x=reorder(Metrics, EF), y=EF, fill=Type)) +
  geom_bar(stat="identity", position=position_dodge()) +
  labs(y = expression(paste("Effect factor [PDF·°",C^-1,"]")), x = "Extinction metric") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  scale_x_discrete(labels = c(expression("E"[1]), expression("E"[2]), expression("E"[PAF0.5]), expression("E"[B]), expression("E"[PAF1]))) +
  scale_fill_viridis_d(begin = 0.4, end = 0.8) +
  theme(legend.position = "bottom") +
  theme_minimal() +
  theme(text = element_text(size = 10),
        legend.text=element_text(size=10),
        axis.text=element_text(size=10),
        legend.position = "bottom",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.text = element_text(angle = 0, vjust = -1, size = 10))
Fig4

ggsave(path = "Figures/", filename = "Fig4.png",Fig4, width = 86, height = 86, unit = 'mm',dpi = 500)

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
Av_EF_sa <- (E_sa_4.5-(E_sa_1.5/1.5*CT))/(4.5-CT)

SA <- data.frame(
  row.names = c("Marginal", "Average"),
  z0.21 = c(Ma_EF_EB, Av_EF_EB),
  z0.25 = c(Ma_EF_sa, Av_EF_sa), 
  increase = c(((Ma_EF_sa - Ma_EF_EB)/Ma_EF_EB*100), ((Av_EF_sa - Av_EF_EB)/Av_EF_EB*100)))

cat("The percentual increase of taking a z component of 0.25 instead of 0.21 amounts to", 
    round(SA[1,3],digits=1), "and", round(SA[2,3],digits=1),"% for the marginal and average approach respectively\n")

### Appendix - regression

Figa <- ggplot(aes(x = dT, y = EB_all), data = Global_extinction) +
  geom_segment(aes(x = 0, xend = 1.5, y = 0, yend = 0.03586883, col = "Marginal"), linetype = "dashed") +
  geom_segment(aes(x = CT, xend = 4.5, y = 0.03586883/1.5*CT, yend = 0.26319496, col = "Average"), linetype = "dashed") +
  geom_point() + 
  geom_point(aes(x = CT, y = (0.03586883/1.5)*CT), shape=8, col = "red", fill = "red") +
  labs(color = 'Type effect factor', x = 'Temperature increase [°C]', y = 'Extinction risk [PDF]') +
  scale_x_continuous(breaks=seq(0,5.5,1)) +
  scale_color_viridis_d(option = "D", begin = 0.4, end = 0.8) +  
  theme_minimal() +
  theme(text = element_text(size = 10),
        legend.text=element_text(size=10),
        strip.text = element_text(angle = 0, vjust = -1, size = 10),
        legend.position = "bottom",
        legend.title = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank())
Figa

ggsave(path = "Figures/", filename = "Figa.png",Figa, dpi = 500, width = 180, height = 150, unit = 'mm')


########PP

Fig22 <- ggplot(data = CF) +
  geom_boxplot(aes(y = A_100_EB, x=factor(1)), fill = '#2c7c8c') +
  geom_boxplot(aes(y = M_100_EB, x=factor(2)), fill = '#7cd454') +
  labs(y = expression(atop("Characterization factor", paste("[PDF·yr·",kg^-1,"]"))),
                       x = "Type") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  scale_x_discrete(labels = c("Average", "Marginal")) +
  theme_minimal() +
  theme(text = element_text(size = 20),
        axis.text=element_text(size=18),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.text = element_text(angle = 0, vjust = -1, size = 18))
Fig22

ggsave(path = "PP/", filename = "Fig22.png",Fig22, height = 12.09, width = 20, unit = 'cm', dpi = 600)


Fig42 <- ggplot(data=EF_bar, aes(x=reorder(Metrics, EF), y=EF, fill=Type)) +
  geom_bar(stat="identity", position=position_dodge()) +
  labs(y = expression(paste("Effect factor [PDF·°",C^-1,"]")), x = "Extinction metric") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  scale_x_discrete(labels = c(expression("E"[1]), expression("E"[2]), expression("E"[PAF0.5]), expression("E"[B]), expression("E"[PAF1]))) +
  scale_fill_viridis_d(begin = 0.4, end = 0.8) +
  theme(legend.position = "bottom") +
  theme_minimal() +
  theme(text = element_text(size = 20),
        axis.text=element_text(size=18),
        legend.text=element_text(size=18),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.text = element_text(angle = 0, vjust = -1, size = 18))
Fig42
ggsave(path = "PP/", filename = "Fig42.png",Fig42, height = 12.09, width = 29.21, unit = 'cm',dpi = 400)

Figa2 <- ggplot(aes(x = dT, y = EB_all), data = Global_extinction) +
  geom_segment(aes(x = 0, xend = 1.5, y = 0, yend = 0.03586883, col = "Marginal"), linetype = "dashed", size = 1.5) +
  geom_segment(aes(x = CT, xend = 4.5, y = 0.03586883/1.5*CT, yend = 0.26319496, col = "Average"), linetype = "dashed", size = 1.5) +
  geom_point(size = 2) + 
  geom_point(aes(x = CT, y = (0.03586883/1.5)*CT), shape=8, size = 3, col = "red", fill = "red") +
  labs(color = 'Type effect factor', x = 'Temperature increase [°C]', y = 'Extinction risk [PDF]') +
  scale_x_continuous(breaks=seq(0,5.5,1)) +
  scale_color_viridis_d(option = "D", begin = 0.4, end = 0.8) +  
  theme_minimal() +
  theme(text = element_text(size = 20),
        legend.text=element_text(size=18),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.text = element_text(angle = 0, vjust = -1, size = 18),
        legend.title = element_blank())
Figa2
ggsave(path = "PP/", filename = "Figa2.png",Figa2, dpi = 500, height = 12.09, width = 29.21, unit = 'cm')



