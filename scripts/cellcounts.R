#* INPUTS: cell count data
#*
#* OUTPUTS: median relative abundance of cogs and petals
#* time series plot (fig 4)
#*


library(tidyverse)
library(janitor)
library(here)
library(lubridate)

theme_set(theme_bw())

plot_cc <- read_csv(here("data/cellcts/plot_cell_counts_tidy.csv"))
glacier_counts <- read_csv(here("data/cellcts/glacier_cell_counts.csv")) 

morph18 <- read_csv(here("data/cellcts/raw_cellct_2018.csv")) # conducted cellcts on n=126 out of 308


rosettes <- read_csv(here("data/field_sample_meta/rosette_field_meta.csv")) # rosette fieldsite metadata

# tidy ---------------------------------------------

# cog time series @ Perley Rock plots
plot_sids <- plot_cc %>% 
  group_by(plot, date) %>% 
  summarise(total = sum(n)) %>% 
  ungroup() %>% 
  mutate(sample_id = paste0(plot, str_sub(date, 7,7), str_sub(date, 9,10))) %>% 
  arrange(plot, date)

plot_clean <- left_join(plot_cc, plot_sids) %>%
  mutate(species = case_when(species=="rid"~"cog",
                             species=="rub"~"cog",
                             species=="rue"~"cog",
                             species=="env"~"sac")) %>%
  filter(species %in% c("cog", "sac")) %>% 
  group_by(plot, date, sample_id, species, total) %>% 
  summarise(n=sum(n)) %>% 
  ungroup() %>% 
  mutate(percent = n/total,
         presence = if_else(n>0,1,0), 
         .keep="unused") %>% 
  select(plot, date, sample_id, species, percent, presence) %>% 
  dplyr::rename(morphotype = species)
# print(plot_clean)

# petal time series @ Weart 2022
glacier_counts_2 <- glacier_counts %>% 
  rename(sac = friedegg) %>% 
  pivot_longer(sphere:sac) %>% # keep 0's!!!
  mutate(percent = value/total, .keep="unused") %>% 
  group_by(date, sample_id, name) %>% 
  summarise(percent = mean(percent)) %>% # take the mean of 2 replicate counts
  ungroup() %>% 
  mutate(site = str_sub(sample_id, 1,3), 
         year = year(date))

glacier_counts_2 %>% 
  filter(name=="sphere") %>% 
  arrange(-percent)

morph18_2 <- morph18 %>% 
  replace_na() %>% 
  dplyr::rename(cog = ruby) %>% 
  pivot_longer(cols = c(-sample_id, -total)) %>% 
  filter(value>0 ) %>% 
  mutate(percent = value/total) %>% 
  mutate(sample_id = sample_id %>% str_remove("\\.")) %>% 
  mutate(name = if_else(name=="orb" & str_detect(sample_id, "wed"), "thickwall", name)) 




# n samples with rosettes ----------------

table_s1 <- rosettes %>% 
  arrange(date) %>% 
  select(sample_id, date, elevation, lon, lat, morphotypes, notes) %>% 
  mutate(elevation = round(elevation, 0),
         morphotypes = str_replace(morphotypes, "cogs","cog"))

# number of samples containing rosettes
table_s1 %>% 
  mutate(morphotypes = str_split(morphotypes, ", ")) %>% 
  unnest(cols = c(morphotypes)) %>%
  count(morphotypes)

# sacs are not included in the list here

# rosette relative abundance ------------------------------

# number of glacier samples w rosettes PRESENT
glacier_counts_2 %>% 
  mutate(presence = if_else(percent>0,1,0),.keep="unused") %>% 
  pivot_wider(names_from = name, values_from = presence) %>% 
  group_by(date, site) %>% 
  summarise(n_sac = sum(sac), n_petal = sum(petal), n=n()) %>% 
  ungroup() %>% 
  mutate(pct_petal = n_petal/n)

# which samples have the most rosettes?
glacier_counts_2 %>% 
  filter(name!="sphere") %>% 
  ggplot(aes(x = sample_id, y = percent, color = name)) +
  geom_point() +
  facet_wrap(vars(site, date), scales = "free") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# median relative abundance
cat <- glacier_counts_2 %>% 
  filter(site=="CAT")
cat %>% 
  group_by(name) %>% 
  summarise(median(percent), sd(percent))

cat %>% 
  filter(percent>0) %>% 
  group_by(name) %>% 
  summarise(median(percent), sd(percent), n = n())

wed <- glacier_counts_2 %>% 
  filter(site=="WED", date>"2022-08-18")
wed %>% 
  group_by(name) %>% 
  summarise(median(percent), sd(percent))

# where present
wed %>% 
  filter(name %>% str_detect("petal")) %>% 
  summarise(median(percent), sd(percent))

glacier_counts_2 %>% 
  filter(name=="petal", percent>0) %>% 
  summarise(median(percent), sd(percent), n=n())


# cog relative abundance 
pa_cogs18 <- morph18_2 %>% 
  filter(name=="cog", value>0) 
pa_cogs18 %>% 
  summarise(median(percent), sd(percent), n())




# relative abundance time series ------------
my_colors <- c("Petal" = "#377eb8", "Cog"="#3c893a", 
               "Sac"= "gray50")


pp1 <- glacier_counts_2 %>% 
  mutate(sample_id = str_remove(sample_id, "WED22"),
         name = str_to_title(name),
         percent = percent*100) %>% 
  filter(site=="WED", name!="Sphere") %>% 
  ggplot(aes(x = sample_id, y = percent, fill = name)) +
  geom_bar(stat = "identity", position = "stack")+
  facet_wrap(vars(date), scales = "free_y", ncol =3) +
  scale_fill_manual(values = my_colors) +
  coord_flip()+
  labs(x = "Weart Glacier Sample ID (WED22_)", y = "Relative abundance (%)", tag = "a", fill = "Morphotype")
  
pp1

cog_plots <- plot_clean %>% 
  filter(plot %in% c("P12", "P5", "P6")) %>% 
  mutate(plot = paste0("P", str_pad(parse_number(plot), width = 2, pad = "0"))) #%>% 
  # mutate(plot = paste(plot, "2020"))

pp2 <- cog_plots %>% 
  filter(date < "2020-07-25") %>% 
  mutate(percent = percent*100,
         morphotype = str_to_title(morphotype)) %>% 
  ggplot(aes(x = plot, y = percent, fill = morphotype)) +
  geom_bar(stat = "identity", position = "stack")+
  facet_wrap(vars(date), ncol=5) +
  scale_fill_manual(values = my_colors) +
  coord_flip()+
  labs(x = "Perley Rock Plot ID", y = "Relative abundance (%)", tag = "b", fill = "Morphotype")
pp2
ggpubr::ggarrange(pp1, pp2, nrow = 2)
ggsave(here("figs/5_timeseries_morph/newlayout.pdf"), width = 6, height = 6)






##  the original line plot-----------

glacier_summary <- glacier_counts_2 %>% 
  filter(site=="WED") %>% 
  select(date, name, percent) %>% 
  group_by(date,name) %>% 
  summarise( sd = sd(percent),percent = mean(percent)) %>% 
  ungroup() %>% 
  filter(name %in% c("sac", "petal")) %>% 
  dplyr::rename(morphotype = name) %>% 
  mutate(plot="Weart Gl. 2022")


cog_plots %>% 
  bind_rows(glacier_summary) %>% 
  mutate(percent = percent*100 ,
         sd = sd*100,
         morphotype = str_to_title(morphotype),
         morphotype = fct_relevel(morphotype,c("Sac","Cog","Petal"))) %>% # %>% fct_recode("Red cell in sac"= "Sac")
  ggplot(aes(x = date, y = percent, color = morphotype)) +# shape=morphotype
  geom_line() +
  geom_point(size = 2.5) +
  facet_wrap(vars(plot), nrow =2,scales="free") +
  scale_color_manual(values = my_colors) +
  scale_x_date(date_labels = "%b%d") +
  # scale_shape_manual(values =c(5,19,19))+
  # geom_errorbar(aes(ymin = percent - sd, ymax = percent+sd))+
  labs(color = "Morphotype", y = "% relative abundance", x = "Date")
ggsave(here("figs/5_timeseries_morph/timeseries_cc_v2.pdf"), width = 7, height = 5, units = "in")






# cog chlain assoc in 2018 data ----------
mycog <- morph18_2 %>% 
  select(-total, -value) %>% 
  pivot_wider(names_from = name, values_from = percent, values_fill = 0) %>% 
  filter(cog>0 | balloon>0)

mycog %>% 
  ggplot(aes(x = cog, y = balloon)) +
  geom_point()
