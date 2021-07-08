setwd("~/Desktop/R")

library(readr); library(dplyr); library(sars); library(ggplot2)
 
Drainage_Basins_Table <- read_delim("datatoFigshare/Drainage_Basins_Table.csv", ";", escape_double = FALSE, trim_ws = TRUE,
                                    col_types = cols_only("1.Basin.Name" = col_character(), "9.Surface.Area" = col_double()))
Occurrence_Table <- read_delim("datatoFigshare/Occurrence_Table.csv", ";", escape_double = FALSE, trim_ws = TRUE, 
                               col_types = cols_only("1.Basin.Name" = col_character(), "3.Native.Exotic.Status" = col_character())) %>%
  subset(., `3.Native.Exotic.Status` != "exotic") %>%
  dplyr::select(.,-2)

SAR_table <- aggregate(Occurrence_Table, by = list(Occurrence_Table$`1.Basin.Name`), FUN = length) %>%
  rename('Basin.Name' = 'Group.1', "Species.Richness" = "1.Basin.Name") %>%
  inner_join(x = ., y = Drainage_Basins_Table, by = c("Basin.Name" = "1.Basin.Name")) %>%
  rename(., "Surface.Area"="9.Surface.Area")

### Regression
sa1 <- lm(log(Species.Richness)~log(Surface.Area), data= SAR_table)

### Results
summary(sa1)

# Fig. 1 SAR plot 
Fig1 <- ggplot(SAR_table, aes(x=Surface.Area, y=Species.Richness)) + 
  geom_point() +
  scale_y_log10(labels = function(x) format(x, scientific = TRUE)) +
  scale_x_log10(labels = function(x) format(x, scientific = TRUE)) +
  labs(x = expression(paste("Surface area [",km^2,"]")), y = 'Species richness [-]') +
  annotation_logticks() +
  geom_smooth(method = 'lm', col = viridis(1)) +
  theme_minimal() +
  theme(text = element_text(size = 12),
        panel.grid.major = element_line(color=NA),
        strip.background = element_rect('white'),
        strip.background.x = element_blank(),
        strip.background.y = element_blank(),
        strip.text = element_text(angle = 0, vjust = -1, size = 12))
Fig1

ggsave(path = "Figures/", filename = "Fig1.png",Fig1, dpi = 600)


