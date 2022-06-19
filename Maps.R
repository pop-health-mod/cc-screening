#Maps
library(ggplot2)
library(rgdal)
library(sf)
library(sp)
library(spData)
library(gridExtra)
library(ggsci)
library(GGally)
library(colorspace)
library(dplyr)
library(ggpattern)
library(purrr)
library(readxl)
library(RColorBrewer)
source("https://raw.githubusercontent.com/imaddowzimet/drawcrosshatch/master/draw_crosshatch.R") 
region_country <- read_excel("surveys.xlsx", sheet = 7)

#load cleaned + modelled data
load("pred_prob_gr_combined74_all_ssa")
##region
scrn_evr_reg <- prob_by_reg_year_ever[prob_by_reg_year_ever$Year == 2020 ,c("Region", "Year", "50%")]
scrn_evr_reg$Region[scrn_evr_reg$Region == "Central/Western"] <- "Western/Central"
names(scrn_evr_reg)[names(scrn_evr_reg) == "Region"] <- "subregion"
##country
scrn_evr_count <- prob_by_count_year_ever[prob_by_count_year_ever$Year == 2020 ,c("Country", "Region", "Year", "50%")]
scrn_evr_count$Region[scrn_evr_count$Region == "Central/Western"] <- "Western/Central"
scrn_evr_count$Country[scrn_evr_count$Country == "Congo (Rep)"] <- "Republic of the Congo"
scrn_evr_count$Country[scrn_evr_count$Country == "Cote d'Ivoire"] <- "Côte d'Ivoire"
scrn_evr_count$Country[scrn_evr_count$Country == "Eswatini"] <- "eSwatini"
names(scrn_evr_count)[names(scrn_evr_count) == "Region"] <- "subregion"

###by HIV status
##country
scrn_evr_hiv_count <- prob_by_hiv_count_year_ever[prob_by_hiv_count_year_ever$Year == 2020 ,c("Country", "Region", "Year", "hivstat", "50%")]
scrn_evr_hiv_count$Region[scrn_evr_hiv_count$Region == "Central/Western"] <- "Western/Central"
scrn_evr_hiv_count$Country[scrn_evr_hiv_count$Country == "Congo (Rep)"] <- "Republic of the Congo"
scrn_evr_hiv_count$Country[scrn_evr_hiv_count$Country == "Cote d'Ivoire"] <- "Côte d'Ivoire"
scrn_evr_hiv_count$Country[scrn_evr_hiv_count$Country == "Eswatini"] <- "eSwatini"
names(scrn_evr_hiv_count)[names(scrn_evr_hiv_count) == "Region"] <- "subregion"

#shape data
africa <- world[world$region_un == "Africa", c("name_long", "continent", "region_un", "subregion", "geom")]
africa$subregion[africa$subregion == "Middle Africa" | africa$subregion == "Western Africa"] <- "Western/Central"
africa$subregion[africa$subregion == "Eastern Africa"] <- "Eastern"
africa$subregion[africa$subregion == "Northern Africa"] <- "Northern"
africa$subregion[africa$subregion == "Southern Africa"] <- "Southern"
names(africa)[names(africa) == "name_long"] <- "Country"
africa$subregion[africa$Country == "Zimbabwe"] <- "Southern"

#combine modelled data with shape data
africa_evr_reg <- left_join(africa, scrn_evr_reg)
africa_evr_count <- left_join(africa, scrn_evr_count)
africa_evr_hiv_count <- left_join(africa, scrn_evr_hiv_count)

#add information about number of surveys
surveys <- read_excel("surveys.xlsx", sheet = 1)
unq_surveys <- unique(surveys[,c("Country", "df")])
num_surveys <- data.frame(table(unq_surveys$Country))
names(num_surveys)[names(num_surveys) == "Var1"] <- "Country"
num_surveys$Country <- as.character(num_surveys$Country)
num_surveys$Country[num_surveys$Country == "Congo (Rep)"] <- "Republic of the Congo"
num_surveys$Country[num_surveys$Country == "Cote d'Ivoire"] <- "Côte d'Ivoire"
num_surveys$Country[num_surveys$Country == "Eswatini"] <- "eSwatini"

africa_evr_reg <- left_join(africa_evr_reg, num_surveys)
africa_evr_count <- left_join(africa_evr_count, num_surveys)
africa_evr_hiv_count <- left_join(africa_evr_hiv_count, num_surveys)

africa_evr_reg$Freq[is.na(africa_evr_reg$Freq)] <- "0"
africa_evr_count$Freq[is.na(africa_evr_count$Freq)] <- "0"
africa_evr_hiv_count$Freq[is.na(africa_evr_hiv_count$Freq)] <- "0"
africa_evr_reg$Freq[africa_evr_reg$Freq >= 2] <- "2+"
africa_evr_count$Freq[africa_evr_count$Freq >= 2] <- "2+"
africa_evr_hiv_count$Freq[africa_evr_hiv_count$Freq >= 2] <- "2+"

#adjust median value to categorical variable
thresholds <- data.frame(Lower = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7),
                         Upper = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 1))

names(africa_evr_count)[names(africa_evr_count) == "50%"] <- "median"
levels = rep(NA, nrow(thresholds))
africa_evr_count$median_cat <- NA
for(i in 1:nrow(thresholds)){
  if(i == 1){
    name <-  paste0("Less than ", thresholds$Upper[i]*100, "%")
    levels[i] <- name
    africa_evr_count$median_cat[africa_evr_count$median > thresholds$Lower[i] & africa_evr_count$median <= thresholds$Upper[i]] <- name
    } else if(i == nrow(thresholds)){
      name <-  paste0("Greater than ", thresholds$Lower[i]*100, "%")
      levels[i] <- name
    africa_evr_count$median_cat[africa_evr_count$median > thresholds$Lower[i] & africa_evr_count$median <= thresholds$Upper[i]] <- name
  } else{
    name <- paste0(thresholds$Lower[i]*100, "-", thresholds$Upper[i]*100, "%")
    levels[i] <- name
    africa_evr_count$median_cat[africa_evr_count$median > thresholds$Lower[i] & africa_evr_count$median <= thresholds$Upper[i]] <- name
  }
}
africa_evr_count$median_cat[africa_evr_count$Freq == "0" ] <- "Not Applicable"
africa_evr_count$median_cat <- factor(africa_evr_count$median_cat, levels = c(levels, "Not Applicable"))

names(africa_evr_hiv_count)[names(africa_evr_hiv_count) == "50%"] <- "median"
africa_evr_hiv_count$median_cat <- NA
for(i in 1:nrow(thresholds)){
  if(i == 1){
    name <-  paste0("Less than ", thresholds$Upper[i]*100, "%")
    levels[i] <- name
    africa_evr_hiv_count$median_cat[africa_evr_hiv_count$median > thresholds$Lower[i] & africa_evr_hiv_count$median <= thresholds$Upper[i]] <- name
  } else if(i == nrow(thresholds)){
    name <-  paste0("Greater than ", thresholds$Lower[i]*100, "%")
    levels[i] <- name
    africa_evr_hiv_count$median_cat[africa_evr_hiv_count$median > thresholds$Lower[i] & africa_evr_hiv_count$median <= thresholds$Upper[i]] <- name
  } else{
    name <- paste0(thresholds$Lower[i]*100, "-", thresholds$Upper[i]*100, "%")
    levels[i] <- name
    africa_evr_hiv_count$median_cat[africa_evr_hiv_count$median > thresholds$Lower[i] & africa_evr_hiv_count$median <= thresholds$Upper[i]] <- name
  }
}
africa_evr_hiv_count$median_cat[africa_evr_hiv_count$Freq == "0" ] <- "Not Applicable"
africa_evr_hiv_count$median_cat <- factor(africa_evr_hiv_count$median_cat, levels = c(levels, "Not Applicable"))

color_range <- colorRampPalette(c("firebrick2", "lightgoldenrod1"))
my_colors <- c(color_range(nrow(thresholds)-1), "seagreen1","grey60")
df_color = data.frame(val = c(1:9), color = my_colors)
ggplot(df_color, aes(x = val, fill = color)) +
  geom_bar() + 
  scale_fill_manual(name = "Ever screened for cervical cancer", 
                    values = df_color$color,
                    labels = c(levels, "Not Applicable"))

#without countries with 0 surveys
ggplot() +
  geom_sf_pattern(data = africa_evr_count, aes(fill = median_cat, pattern = Freq), 
                  pattern_fill = "grey20", pattern_density = 0.05, pattern_spacing = 0.02) +
  scale_pattern_manual(name = "Number of Surveys", values = c("none", "stripe", "none")) +
  scale_fill_manual(name = "Ever screened for cervical cancer", values = my_colors[c(1:3, 6, 9)]) +
  theme_minimal() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank()) +
  guides(pattern = "none", 
         pattern_spacing = "none")

#without countries with 0 surveys and women HIV+
ggplot() +
  geom_sf_pattern(data = africa_evr_hiv_count[africa_evr_hiv_count$hivstat == "HIV+" | is.na(africa_evr_hiv_count$hivstat),], 
                  aes(fill = median_cat, pattern = Freq), 
                  pattern_fill = "grey20", pattern_density = 0.05, pattern_spacing = 0.02) +
  scale_pattern_manual(name = "Number of Surveys", values = c("none", "stripe", "none")) +
  scale_pattern_spacing_manual(values = c(0.01, 0.01, 0.01)) +
  scale_fill_manual(name = "Ever screened for cervical cancer", values = my_colors[c(1:5,7)]) +
  theme_minimal() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank()) +
  guides(pattern = "none",
         pattern_spacing = "none")

#without countries with 0 surveys and women HIV-
ggplot() +
  geom_sf_pattern(data = africa_evr_hiv_count[africa_evr_hiv_count$hivstat == "HIV-" | is.na(africa_evr_hiv_count$hivstat),], 
                  aes(fill = median_cat, pattern = Freq), 
                  pattern_fill = "grey20", pattern_density = 0.05, pattern_spacing = 0.02) +
  scale_pattern_manual(name = "Number of Surveys", values = c("none", "stripe", "none")) +
  scale_pattern_spacing_manual(values = c(0.01, 0.01, 0.01)) +
  scale_fill_manual(name = "Ever screened for cervical cancer", values = my_colors[c(1:5,7)]) +
  theme_minimal() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank()) +
  guides(pattern = "none",
         pattern_spacing = "none")

#with countries with 0 surveys
ggplot() +
  geom_sf_pattern(data = africa_evr_count, aes(fill = median_cat2, pattern = Freq, pattern_spacing = Freq), 
                  pattern_fill = "grey20", pattern_density = 0.1) +
  scale_pattern_manual(name = "Number of Surveys", values = c("stripe", "stripe", "none")) +
  scale_pattern_spacing_manual(values = c(0.01, 0.03, 0.05)) +
  scale_fill_manual(name = "Ever screened for cervical cancer", values = my_colors) +
  theme_minimal() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank()) +
  guides(pattern = "none",
         pattern_spacing = "none")

#map for small african island countries
#load africa shapefile and clean to match
a <- read_sf('./shapefile','afr')
africa_2 <- a[a$ADM0_NAME %in% c("Cape Verde", "Mauritius", "Sao Tome and Principe"),c("ADM0_NAME", "ISO3", "geometry")]
names(africa_2)[names(africa_2) == "ADM0_NAME"] <- "Country"
africa_2 <- left_join(africa_2, region_country)
africa_2$Region[africa_2$Region == "Central" | africa_2$Region == "Western"] <- "Western/Central"

#combine modelled data with shape data
africa_evr_count_2 <- left_join(africa_2, scrn_evr_count)
africa_evr_hiv_count_2 <- left_join(africa_2, scrn_evr_hiv_count)

#add information about number of surveys
africa_evr_count_2 <- left_join(africa_evr_count_2, num_surveys)
africa_evr_hiv_count_2 <- left_join(africa_evr_hiv_count_2, num_surveys)

names(africa_evr_count_2)[names(africa_evr_count_2) == "50%"] <- "median"
levels = rep(NA, nrow(thresholds))
africa_evr_count_2$median_cat <- NA
for(i in 1:nrow(thresholds)){
  if(i == 1){
    name <-  paste0("Less than ", thresholds$Upper[i]*100, "%")
    levels[i] <- name
    africa_evr_count_2$median_cat[africa_evr_count_2$median > thresholds$Lower[i] & africa_evr_count_2$median <= thresholds$Upper[i]] <- name
  } else if(i == nrow(thresholds)){
    name <-  paste0("Greater than ", thresholds$Lower[i]*100, "%")
    levels[i] <- name
    africa_evr_count_2$median_cat[africa_evr_count_2$median > thresholds$Lower[i] & africa_evr_count_2$median <= thresholds$Upper[i]] <- name
  } else{
    name <- paste0(thresholds$Lower[i]*100, "-", thresholds$Upper[i]*100, "%")
    levels[i] <- name
    africa_evr_count_2$median_cat[africa_evr_count_2$median > thresholds$Lower[i] & africa_evr_count_2$median <= thresholds$Upper[i]] <- name
  }
}
africa_evr_count_2$median_cat[africa_evr_count_2$Freq == "0" ] <- "Not Applicable"
africa_evr_count_2$median_cat <- factor(africa_evr_count_2$median_cat, levels = c(levels, "Not Applicable"))

ggplot() +
  geom_sf_pattern(data = africa_evr_count_2[africa_evr_count_2$Country == "Cape Verde",], aes(fill = median_cat), 
                  pattern_fill = "grey20", pattern_density = 0.05, pattern_spacing = 0.02) +
  scale_pattern_manual(name = "Number of Surveys", values = c("none", "stripe", "none")) +
  scale_fill_manual(name = "Ever screened for cervical cancer", values = my_colors[c(1:3, 5, 7)]) +
  theme_minimal() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank()) +
  guides(pattern = "none", 
         pattern_spacing = "none")
