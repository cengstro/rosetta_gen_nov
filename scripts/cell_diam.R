#* Compare cog diameters of the two genotypes/ CBC clades
#* take photos of cells, measure inner (red part) and outer diameters in ImageJ 
#*
#*


library(tidyverse)
library(janitor)
library(here)
library(fs)

raw <- dir_ls(here("data/cell_diam/"), regexp = ".csv") %>% 
  map_df(read_csv, .id = "filename")

theme_set(theme_bw())

# tidy ---------------------------
diams <- raw %>% 
  janitor::clean_names() %>% 
  mutate(morphotype = basename(filename) %>% str_split_i("_", 1),
         sample_id = basename(filename) %>% str_split_i("_", 2),
         sample_id = str_remove(sample_id, ".csv"),
         # I first measured the tip to tip diameter, then the inner red diameter 
         measurement = if_else(x1%%2==0, "inner_diam", "outer"),
         # 1,2 are one cell, 3,4 the next, etc (2 measurements per cell)
         cell_id = if_else(x1%%2!=0, (x1+1)/2, x1/2)) %>%
  select(morphotype, sample_id, cell_id, measurement, length) %>% 
  pivot_wider(names_from = "measurement", values_from = "length") %>% 
  # calculate wall thickness, and wall to inner diameter ratio derived from inner and tip to tip measurements
  mutate(max_outer_wall_thickness = (outer-inner_diam) / 2,
         wall_to_inner_diam_ratio = max_outer_wall_thickness / inner_diam ) %>% 
  # mutate(morphotype = fct_relevel(morphotype, c("star", "cog", "petal"))) %>% 
  # mutate(morphotype =morphotype %>% fct_relevel(c("cog","petal","star","thickwall","raisin"))) %>% 
  mutate(morphotype = case_when(morphotype=="raisin"~ "poppy",
                                morphotype=="thickwall" ~ "globe",
                                .default = morphotype)) %>% 
  mutate(morphotype = fct_reorder(morphotype, max_outer_wall_thickness))

diam_long <- diams %>% 
  pivot_longer(cols = c(inner_diam, max_outer_wall_thickness), names_to = "measurement")

diams %>% 
  distinct(sample_id, morphotype)

diams %>% 
  count(sample_id, morphotype) %>% 
  arrange(morphotype) %>% 
  group_by(morphotype) %>% 
  summarise(n_cells_meas = sum(n), n_field_sample = n() )

# # A tibble: 5 × 3
# morphotype n_cells_meas n_field_sample
# <fct>             <int>          <int>
# 1 star                 72              2
# 2 cog                 140              7
# 3 petal                46              2
# 4 raisin                8              1
# 5 thickwall           143              4

# plots -------------------------------

a <- diams %>% 
  ggplot(aes(x = morphotype, y = inner_diam)) +
  geom_boxplot() +
  labs(x = "Morphotype", y= "Maximum inner red-to-red diameter (µm)", tag="a")

b <- diams %>% 
  ggplot(aes(x = morphotype, y = max_outer_wall_thickness)) +
  geom_boxplot() +
  labs(x = "Morphotype", y= "Maximum projection length (µm)", tag="b")

ggpubr::ggarrange(a,b, nrow = 2)
ggsave(here("figs/0_all_figs/fig_S2.tiff"), width = 5, height = 8, units = "in")


# stats-----------
diams %>% 
  summarise(mean = mean(inner_diam))


# 
# 
# # scratch ------------------------
# 
# diams %>% 
#   ggplot(aes(x = type, y = inner)) +
#   geom_boxplot() 
# 
# # compare inner diameters 
# diams %>% 
#   ggplot(aes(x = sample_id, y = inner, color = type)) +
#   geom_boxplot()
# 
# lm(inner~sample_id, data = diams) %>% summary()
# 
# # outer diameters
# diams %>% 
#   ggplot(aes(x = sample_id, y = outer, color = type)) +
#   geom_boxplot() 
# lm(outer~sample_id, data = diams) %>% summary()
# 
# # according to this criteria, p6618 should be gear1
# 
# # any difference in wall thickness?
# diams %>% 
#   ggplot(aes(x  = sample_id, y = wall_thickness, color = type)) +
#   geom_boxplot()
# 
# 
# # which small cells are "full grown", and which are "in development"?
# # developing cells should have a thinner cell wall relative to their diameter
# diams %>% 
#   ggplot(aes(wall_to_diam_ratio)) +
#   geom_histogram(bins = 45) 
# 
# diams %>% 
#   arrange(wall_to_diam_ratio)
# 
# arrange(diams, wall_thickness)
# 
# 
# 
# # results --------------
# #* 2 out of 3 "gear2" samples have significantly smaller diameters than "gear1"
# #* but P6618 is not significantly different from gear1
# #* 
# #* 
# #* 
# #* types <- tribble(
# # ~sample_id, ~type,
# # "bak1912", "gear1",
# # "gar1912","gear1",
# # "gar1930","gear1",
# # "pan1918","gear1",
# # "p6618", "gear2",
# # "pan1915", "gear2",
# # "tem2003", "gear2"
# # )
# 
