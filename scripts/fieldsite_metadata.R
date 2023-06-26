#* INPUTS
#* csv of all field samples
#* csv summary of sampling locations
#* csv of only rosette field samples (manually curated)
#* 
#* 
#* OUTPUTS
#* map showing field sampling locations
#* table of rosette field sample metadata
#* histograms of rosette elevation, date, and region
#*
#*

library(tidyverse)
library(janitor)
library(here)
library(lubridate)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(flextable)
library(officer)

theme_set(theme_bw())

samples <- read_csv(here("data/field_sample_meta/all_samples_tidy.csv")) # all sample metadata
my_mtns <- read_csv(here("data/field_sample_meta/mtns_sampled.csv")) # mtn locations

rosettes <- read_csv(here("data/field_sample_meta/rosette_field_meta.csv")) # rosette fieldsite metadata

# map data
states <- rnaturalearth::ne_states(c("united states of america", "canada"), returnclass = "sf")



# number of samples, etc ----------------------------------------
distinct(samples, sample_id) %>% nrow() # 762 samples
distinct(samples, date) %>% arrange(date) %>% print(n=84) # 84 days
distinct(samples, mountain) %>% arrange(mountain) %>% print(n=33)  # 33 different mountains

# number of samples containing Cogs
rosettes %>% 
  filter(morphotypes %>% str_detect("cog")) %>% 
  nrow()

# number of mountains where Cogs were detected
rosettes %>% 
  filter(morphotypes %>% str_detect("cog")) %>% 
  mutate(mountain = str_sub(sample_id, 1, 3) %>% 
           str_remove_all("[:digit:]")) %>% 
  distinct(mountain) %>% 
  arrange(mountain)

# n samples containing petals
rosettes %>% 
  filter(morphotypes %>% str_detect("petal")) %>% 
  nrow()

# map of field sampling locations --------------------------------

america_proj <- "+proj=aea +lon_0=-121.3769531 +lat_1=47.3073282 +lat_2=51.999961 +lat_0=49.6536446 +datum=WGS84 +units=m +no_defs"

my_mtns_sf <- my_mtns %>% 
  st_as_sf(coords = c("lon", "lat")) %>% 
  sf::st_set_crs(4326) %>% 
  mutate(habitat =fct_relevel(habitat, c("glacier", "alpine", "alpine, forest", "forest")), 
         mtn_code = str_to_upper(mtn_code)) %>% 
  st_transform(america_proj)



my_map_theme <-   theme(panel.background = element_rect(fill = "gray50", colour = NA), # make the ocean near black, 
                        panel.grid.major = element_line(colour="gray30", size=0.1)) # political boundaries in grey


cities <- tribble(
  ~name, ~lon,~lat,
  "Calgary", -114.0719, 51.0447,
  "Vancouver",-123.1207, 49.2827,
  "Seattle",-122.3321,47.6062
) %>% 
  st_as_sf(coords = c("lon", "lat")) %>% 
  sf::st_set_crs(4326) %>% 
  st_transform(america_proj)


# 6 class accent
my_pal <- c("#beaed4","#fdc086", "red", "#7fc97f")
ggplot() +  
  geom_sf(data = st_transform(states, america_proj),  fill="gray20", color="gray40") + 
  geom_sf(data = my_mtns_sf, aes(color = habitat), size =2) +
  geom_sf(data = cities) +
  coord_sf(xlim = c(-0.5e6, 0.6e6), ylim = c(-0.3e6, 0.25e6), expand = FALSE) + 
  scale_color_manual(values = my_pal) +
  labs(color = "Habitat") +
  my_map_theme +
  ggspatial::annotation_scale(location = 'bl', width_hint = 0.15) +
  geom_sf_text(data = my_mtns_sf, aes(label = mtn_code, color = habitat))
ggsave(here("figs/sample_sites.pdf"))

# replace red dots with half orange-green dots in inkscape






# rosette sample table --------------------
table_s1 <- rosettes %>% 
  arrange(sample_id) %>% 
  select(sample_id, date, elevation, lon, lat, morphotypes, notes) %>% 
  mutate(elevation = round(elevation, 0)) %>% 
  mutate(morphotypes = str_replace(morphotypes, "cogs", "cog") %>% 
           str_replace("raisin", "poppy") %>% 
           str_replace("thickwall", "globe"))
# 
# # export table
# flextable(table_s1) %>% 
#   save_as_docx( path = here("figs/tableS1.docx"))
# 
# table_s1 %>% 
#   write_csv(here("figs/table_s1/table_s1.csv"))

# rosettes histograms --------------------------

# table_s1 %>% 
my_colors <- c("petal" = "#377eb8", "star"="#984ea3", "thickwall"="#ff7f00", "globe"="#ff7f00", "bigthickwall"="#ff7f00","bigThickwall"="#ff7f00", "poppy"="#a65628", "cog"="#3c893a", 
               "sac"= "gray50", "other"="black")

# 
# all_samples <- samples %>% 
#   mutate(sample_id = str_to_upper(sample_id))

# all_samples %>% 
#   left_join(table_s1) %>% 
#   
hist_dat <- table_s1 %>% 
  mutate(yday = yday(date),
         date_dummy = as.Date(yday, "1900-01-01"),
         mtn_code = str_sub(sample_id, 1,3) %>% 
           str_to_lower() %>% 
           str_replace("p[:digit:]{2}", "per"),
         morphotypes = str_remove(morphotypes, ", daughters"),
         morphotypes = str_split(morphotypes, ", "),
         weight = 1/map_vec(morphotypes, length)) %>% 
  left_join(my_mtns %>% select(mtn_code, region)) %>% 
  mutate(region = fct_relevel(region, c("SW BC", "SE BC"))) %>% 
  unnest(cols = c(morphotypes)) %>% 
  mutate(morphotypes = fct_relevel(morphotypes, "cog","petal","star", "globe", "poppy"))
  

p1 <- hist_dat %>% 
  ggplot(aes(x = date_dummy, fill = morphotypes)) +
  scale_fill_manual(values = my_colors) +
  geom_histogram(aes(weight =weight), color = "black", size = 0.3) +
  labs(x = "Date", y = "N samples", tag = "a") +
  theme(legend.position = "none")
p1

p2 <- hist_dat %>% 
  ggplot(aes(x = elevation, fill = morphotypes)) +
  geom_histogram(aes(weight = weight), color = "black", size = 0.3) +
  labs(x = "Elevation (m)", y = "N samples", tag = "b") +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "none")
p2

p3 <- hist_dat %>% 
  ggplot(aes(x = region, fill = morphotypes)) +
  geom_histogram(aes(weight =weight), stat = "count") +
  scale_fill_manual(values = my_colors) +
  labs(x = "Region", y = "N samples", tag = "c")
p3
leg <- cowplot::get_legend(p3)
ggpubr::ggarrange(p1, p2, p3+theme(legend.position = "none"), leg, nrow = 2, ncol = 2)

ggsave(here("figs/s1_date_elev_hists/fig_s1.tiff"), height = 4, width = 6)
