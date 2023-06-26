#* put the data wrangling outside of the function?
#* tidy this script up, remove redundancy 
#*
#*
#*
#*

library(ape)
library(ggtree)
library(tidyverse)
library(here)
library(phangorn)
library(treeio)

# read treefiles stored in the `phylo` directory --------------------
today <- "2023-06-24" 

# use consistent morphotype
update_morphonames <- function(tree){
  tips <- tree$tip.label
  newtips <- tips %>% 
      str_replace("callous", "globe") %>% 
      str_replace("thickWall", "globe") %>% 
      str_replace("thickwall", "globe") %>% 
      str_replace("bigThickWall", "bigGlobe") %>% 
      str_replace("envelope", "sac") %>% 
      str_replace("friedEgg","sac") %>% 
    str_replace("raisin","poppy") %>% 
    str_replace("martian","sac") %>% 
    # str_replace("_OTU_D_Engstrom_etal_2020_rbcL_consens", "otuD_EngstromEtal2020") %>%
    str_remove("_NA")
  tree$tip.label <- newtips
  return(tree)
}

my_read <- function(marker, date, root){
  ape::read.tree(here::here(paste0("data/phylo/",date,"_", marker,"/aligns.fasta.treefile"))) %>%
    root(outgroup = root, resolve.root=T) %>% 
    update_morphonames
}

# ITS2
its2_dat <- my_read("its2", today, "FR865586.1_Chlamydomonas_reinhardtii_g")

# rbcL
rbcl_dat <- my_read("rbcl", today, "AB511845.1_Chlamydomonas_reinhardtii_c")
# rbcL insert
rbcl_insert_dat <- ape::read.tree(here::here(paste0("data/phylo/",today,"_rbcl_insert/aligns.fasta.treefile"))) %>%
  update_morphonames()

# 18S Full
# x18sfull_dat <- my_read("18sfull", today, "FR865586.1_Chlamydomonas_reinhardtii")

# 18S Gappy
x18sgappy_dat <- my_read("18sgappy", today, "FR865586.1_Chlamydomonas_reinhardtii")

# cat
cat_dat <- my_read("cat", today, "Chlamydomonas_reinhardtii")


# # 18s A
x18sa_dat <- my_read("18sa", today, "FR865586.1_Chlamydomonas_reinhardtii_g")
# 
# # 18s B
# x18sb_dat <- my_read("18sb", today, "FR865586.1_Chlamydomonas_reinhardtii")

# ribo_dat <- ape::read.tree(here::here("data/phylo/2023-04-03_ribocat/aligns.fasta.treefile")) %>%
#   root(outgroup = "Chlamydomonas_reinhardtii", resolve.root=T)






# Functions ----------------------------------------
gb_regex <- "[:alpha:]{2}[:digit:]{6}.*"
species_regex <- "Chloromonas|Chlamydomonas|Chlainomonas|Scotiella|Sanguina|Ploetila"
sid_regex <- "[:alpha:]+[:digit:]{2}.*"
boot_threshold <- 50

my_colors <- c("petal" = "#377eb8", "star"="#984ea3", "globe"="#ff7f00", "bigGlobe"="#ff7f00","bigGlobe"="#ff7f00", "poppy"="#a65628", "cog"="#3c893a", 
               "sac"= "gray50", "other"="black")


# adds new columns
my_ggtree <- function(tree, node1 = FALSE, node2 = FALSE){
  # tree <- its2_dat
  if(node1){
    gg <- ggtree(tree) %>% 
      flip(node1, node2)
  }else{
    gg <- ggtree(tree)
  }
  
  dat <- gg$data
  
  suppressWarnings(
    newdat <- dat %>% 
      # identify label type
      mutate(type = case_when(str_detect(label, gb_regex) ~"genbank",
                              str_detect(label, species_regex) ~ "genbank_no_gba",
                              (!str_detect(label, gb_regex) & str_detect(label, sid_regex) ) ~ "mysample",
                              # str_detect(label, "otuD") ~"genbank_no_gba",
                              str_detect(label, "[:digit:]+") ~ "bootstrap",
                              .default = "other"),
             genbank_label = if_else(type %>% str_detect("genbank"), label, NA_character_),
             my_sample_label = if_else(type=="mysample", label, NA_character_),
             bootstrap_num = if_else(type=="bootstrap", as.numeric(label), NA)) %>% 
      # omit bootstrap vals below threshold
      mutate(bootstrap_label_filt = if_else(bootstrap_num > boot_threshold, as.character(bootstrap_num), "") ) %>% 
      # parse data from each type of label
      separate_wider_delim(cols = my_sample_label,delim = "_", 
                           names = c("morphotype","sample_id", "tube_id"),
                           too_many = "drop", too_few = "align_start") %>% 
      mutate(morphotype_color_group = if_else(morphotype %in% names(my_colors),
                                              morphotype, "other")) %>%
      separate_wider_delim(cols = genbank_label, delim = "_", 
                           names = c("gba","genus", "species"),
                           too_many = "drop", too_few = "align_start") %>% 
      mutate(newlabel = case_when(type=="genbank" ~ paste(gba, genus, species, sep=" "),
                                  type == "genbank_no_gba"~ label,
                                  type=="mysample" ~ paste("GBXXXXXX.1", morphotype, sample_id, paste0(tube_id, "*"), sep=" "),
                                  .default = NA))
  )
  # overwrite the old $data
  gg$data <- newdat
  return(gg)
}

find_mrca <- function(tree, regex){
  tips <- tree$tip.label
  node_ids <- str_which(tips, regex) 
  return(mrca.phylo(tree, node_ids))
}

my_cladelab <- function(data, regex, label, offset=0, txtoffset=0, color="black", rotate = FALSE, ff=3, lt = 1){
  if(rotate){
    geom_cladelab(node=find_mrca(data, regex), label=label, offset = offset, barcolor = color, textcolor = color, fontface = ff, linetype = lt,
                  angle=270, offset.text = txtoffset, hjust = 'center')
  }else{
    geom_cladelab(node=find_mrca(data, regex), label=label, offset = offset, barcolor = color, textcolor = color, fontface = ff, linetype = lt)
  }
}




# ITS2 -----------------------------

its2_offset <- 2.7
its2_txtoff <- 0.1

ggtree(its2_dat) +
  geom_text(aes(label=node))

my_ggtree(its2_dat, 90, 153) +
  geom_text(aes(label=bootstrap_label_filt))


its2_treeplot <- my_ggtree(its2_dat,90,153) + # , 88, 151
  geom_text(aes(label=bootstrap_label_filt),  vjust= 0.01, hjust = 1, size=2.5, color = "#666666") +
  geom_tiplab(aes(label = newlabel, color = morphotype_color_group), size=3) +
  xlim(0,8) +
  scale_color_manual(values = my_colors) +
  theme(legend.position = "none") +
  my_cladelab(its2_dat, "petal|sac_WED2206", "R. floraniva sp. nov.", its2_offset, its2_txtoff) +
  # my_cladelab(its2_dat, "cog", "R. castellata sp. nov.", its2_offset-0.4, its2_txtoff) + #, lt=2
  my_cladelab(its2_dat, "star", "R. stellaria sp. nov.", its2_offset, its2_txtoff) +
  my_cladelab(its2_dat, "globe", "R. bilasso sp. nov.", its2_offset, its2_txtoff) +
  # my_cladelab(its2_dat, "poppy", "R. papaverna sp. nov.", its2_offset, its2_txtoff) +
  my_cladelab(its2_dat, "sac|cog|petal|thickwall|raisin", "Rosetta gen. nov.", its2_offset+1.3, its2_txtoff) +
  my_cladelab(its2_dat, "sac|cog|petal|thickwall|raisin|typhlos", "Chloromonadinia", its2_offset+2, its2_txtoff, ff=1) +
  my_cladelab(data=its2_dat, regex="halo|bigSphere|shortSpike|Sanguina", label="Sanguina", offset=its2_offset+3, txtoffset=its2_txtoff, color="black") +
  my_cladelab(its2_dat, "Chlainomonas|football", "Chlainomonas", its2_offset+1, its2_txtoff, "black") +
  # my_cladelab(its2_dat, "chenangoensis|remiasii", "D", its2_offset+0.12, its2_txtoff, rotate=FALSE, ff=1) +
  my_cladelab(its2_dat, "krienitzii|orangeBrevi|squareBrevi|hindakii", "Group B", its2_offset+1, its2_txtoff, rotate=FALSE, ff=1) +
  my_cladelab(its2_dat, "miwae|pichinchae|muramotoi|fukushimae|hohamii|tughilliensis", "Groups A+C", its2_offset+1, its2_txtoff, rotate=FALSE, ff=1) +
  # my_cladelab(its2_dat, "", "C", its2_offset+0.8, its2_txtoff, rotate=FALSE, ff=1) +
  # my_cladelab("Chloromonas|Chlainomonas|cog", "Chloromonadinia", 0, "black") + not working
  labs(title = "ITS2") +
  geom_treescale(x = its2_offset*2)
# scaleClade(p, node=17, scale=.1) 
its2_treeplot


ggsave(paste0(here("figs/8_its2/"), today, "_its2.pdf"), 
       width = 7, height = 10, units = "in")



# rbcL ------------------------------------
rbcl_offset <- 0.3
rbcl_txtoff <- 0.01

collapse(my_ggtree(rbcl_dat),node=72) +
  geom_text(aes(label=node))  
# view(gg_rbcl$data)
# gg_rbcl$data %>% 
#   mutate(newlabel = if_else(newlabel %>% str_detect())


rbcl_treeplot <-
  my_ggtree(rbcl_dat) +
  geom_text(aes(label=bootstrap_label_filt),  vjust= 0.01, hjust = 1, size=2.5, color = "#666666") +
  geom_tiplab(aes(label = newlabel, color = morphotype_color_group), size=3) +
  xlim(0,0.75)+
  scale_color_manual(values = my_colors) +
  theme(legend.position = "none") +
  my_cladelab(rbcl_dat, "cog", "R. castellata sp. nov.", rbcl_offset-0.05, rbcl_txtoff, "black", rotate=FALSE) +
  my_cladelab(rbcl_dat, "star", "R. stellaria sp. nov", rbcl_offset-0.05, rbcl_txtoff, "black", rotate=FALSE) +
  # my_cladelab(rbcl_dat, "thickwall", "R. rotunda sp. nov", rbcl_offset-0.05, rbcl_txtoff, "black", rotate=FALSE) +
  my_cladelab(rbcl_dat, "petal", "R. floraniva sp. nov", rbcl_offset-0.05, rbcl_txtoff, "black", rotate=FALSE) +
  my_cladelab(rbcl_dat, "smallSphere|LaHoya|aurantia", "Sanguina", rbcl_offset+0.1, rbcl_txtoff, "black", rotate=FALSE) +
  my_cladelab(rbcl_dat, "Chlainomonas", "Chlainomonas", rbcl_offset+0.02, rbcl_txtoff, "black", rotate=FALSE) +
  my_cladelab(rbcl_dat, "sac|star|cog|thickwall", "Rosetta gen. nov.", rbcl_offset+0.1, rbcl_txtoff, "black", rotate=FALSE) +
  # my_cladelab(rbcl_dat, "chenangoensis|remiasii", "Group D", rbcl_offset+0.06, rbcl_txtoff, rotate=FALSE, ff=1) +
  my_cladelab(rbcl_dat, "krienitzii|orangeBrevi|hindakii", "Group B", rbcl_offset+0.06, rbcl_txtoff, rotate=FALSE, ff=1) +
  my_cladelab(rbcl_dat, "miwae|pichinchae|fukushimae|hohamii", "Groups A+C", rbcl_offset+0.03, rbcl_txtoff, rotate=FALSE, ff=1) +
  my_cladelab(rbcl_dat, "sac|cog|petal|thickwall|raisin|typhlos", "Chloromonadinia", rbcl_offset+0.5, rbcl_txtoff, ff=1) +
  # my_cladelab(rbcl_dat, "reticulata|augustae|arctica|kasaiae", "Non-snow Chloromonas", rbcl_offset+0.03, rbcl_txtoff, rotate=FALSE, ff=1) +
  # my_cladelab(rbcl_dat, "fukushimae|hohamii", "C", rbcl_offset+0.06, rbcl_txtoff, rotate=FALSE, ff=1) + # muramotoi was in group A in matsuzaki etal 2019
  labs(title = "rbcL w insert")+
  geom_treescale(x = rbcl_offset*2.2)
rbcl_treeplot

ggsave(paste0(here("figs/6_rbcl/"), today, "_rbcl.pdf"), 
       width = 6, height = 9, units = "in")
# ggsave(paste0(here("figs/trees/rbcl/"), today, "_rbcl.pdf"),
#        width = 6, height = 8, units = "in")


# rbcl insert dat---------------------



my_ggtree(rbcl_insert_dat) +
  geom_text(aes(label=bootstrap_label_filt),  vjust= 0.01, hjust = 1, size=2.5) +
  geom_tiplab(aes(label = newlabel, color = morphotype_color_group), size=3) +
  xlim(0,0.75)+
  scale_color_manual(values = my_colors) +
  theme(legend.position = "none") +
  labs(title = "rbcL intron")+
  geom_treescale(x = 0.2)

ggsave(paste0(here("figs/rbcl_introns/"), today, "_rbcl_insert.pdf"), 
       width = 4, height = 3, units = "in")

# x18sgappy_tree --------------

x18s_offset <- 0.06
x18s_txtoff <- 0.01

# my_ggtree(x18sgappy_dat) +
  # geom_text(aes(label=bootstrap_label_filt),  vjust= 0.01, hjust = 1, size=2.5, color = "#666666")
collapse(collapse(collapse(collapse(collapse(ggtree(x18sgappy_dat), node=90), node = 131), node = 138),node=145),node=152) +
  geom_point2(aes(subset=(node==90)), shape=23, size=5, fill='red')+
  geom_text(aes(label=node))

x18sgappy_treeplot <-
  my_ggtree(x18sgappy_dat) +
  geom_text(aes(label=bootstrap_label_filt),  vjust= 0.01, hjust = 1, size=2.5, color = "#666666") +
  geom_tiplab(aes(label = newlabel, color = morphotype_color_group), size=3) +
  xlim(0,0.17)+
  scale_color_manual(values = my_colors) +
  theme(legend.position = "none") +
  my_cladelab(x18sgappy_dat, "globe", "R. rubriterra sp. nov", x18s_offset-0.02, x18s_txtoff) +
  my_cladelab(x18sgappy_dat, "petal|star", "R. floraniva &\nR. stellaria spp. nov.", x18s_offset+0.025, x18s_txtoff, lt=2) +
  my_cladelab(x18sgappy_dat, "cog", "R. castellata sp. nov.", x18s_offset-0.02, x18s_txtoff) +
  my_cladelab(x18sgappy_dat, "poppy", "R. papavera sp. nov.", x18s_offset-0.02, x18s_txtoff) +
  my_cladelab(x18sgappy_dat, "smallSphere|Sanguina|honeycomb|bigSphere", "Sanguina", x18s_offset+0.015, x18s_txtoff) +
  my_cladelab(x18sgappy_dat, "Chlainomonas", "Chlainomonas", x18s_offset+0.005, x18s_txtoff) +
  my_cladelab(x18sgappy_dat, "cog|raisin|star|thickwall", "Rosetta gen. nov.", x18s_offset+0.01, x18s_txtoff) +
  # my_cladelab(x18sgappy_dat, "chenangoensis|remiasii", "D", x18s_offset, x18s_txtoff, rotate=FALSE, ff=1) +
  my_cladelab(x18sgappy_dat, "krienitzii|tatrae|hindakii", "Group B", x18s_offset+0.005, x18s_txtoff, rotate=FALSE, ff=1) +
  my_cladelab(x18sgappy_dat, "miwae|pichinchae|muramotoi", "Group A", x18s_offset+0.005, x18s_txtoff, rotate=FALSE, ff=1) +
  my_cladelab(x18sgappy_dat, "fukushimae|hohamii|tughilliensis", "Group C", x18s_offset+0.007, x18s_txtoff, rotate=FALSE, ff=1) + # muramotoi was in group A in matsuzaki etal 2019
  labs(title = "18S")+
  my_cladelab(x18sgappy_dat, "sac|cog|petal|thickwall|raisin|typhlos", "Chloromonadinia", x18s_offset+0.05, x18s_txtoff, ff=1) +
  geom_treescale(x = x18s_offset*2.3) 
x18sgappy_treeplot



ggsave(paste0(here("figs/7_18s/"), today, "_18s.pdf"), 
       width = 6, height =11, units = "in")
# ggsave(paste0(here("figs/trees/18s/"), today, "_18s_gappy.pdf"), 
#        width = 6, height = 11, units = "in")



# CAT concatenated sequence tree  ----------------------------
cat_offset <- 0.3
cat_txtoff <- 0.001

ggcat <- my_ggtree(cat_dat)
ggcat$data$newlabel <- 
  str_remove( string = ggcat$data$newlabel, 
              pattern = "_NA"
              )
# view(ggcat$data)

cat_treeplot <-
  ggcat +
  geom_text(aes(label=bootstrap_label_filt),  vjust= 0.01, hjust = 1, size=2.5, color = "#666666") +
  geom_tiplab(aes(label = newlabel, color = morphotype_color_group), size=4) +
  xlim(0,1.1)+
  scale_color_manual(values = my_colors) +
  theme(legend.position = "none") +
  my_cladelab(cat_dat, "smallSphere|Sanguina", "Sanguina", cat_offset+0.15, cat_txtoff) +
  my_cladelab(cat_dat, "Chlainomonas", "Chlainomonas", cat_offset+0.07, cat_txtoff) +
  my_cladelab(cat_dat, "cog|raisin|petal|star", "Rosetta gen. nov.", cat_offset+0.2, cat_txtoff) +
  my_cladelab(cat_dat, "cog", "R. castellata sp. nov.", cat_offset+0.05, cat_txtoff) +
  my_cladelab(cat_dat, "petal", "R. floraniva sp. nov.", cat_offset+0.05, cat_txtoff) +
  my_cladelab(cat_dat, "star", "R. stellaria sp. nov.", cat_offset, cat_txtoff) +
  my_cladelab(cat_dat, "globe", "R. bilasso sp. nov.", cat_offset, cat_txtoff) +
  my_cladelab(cat_dat, "pich|mura", "Group A", cat_offset+0.15, cat_txtoff, ff=1) +
  # my_cladelab(cat_dat, "remias|chenangoensis", "Group D", cat_offset+0.18, cat_txtoff, ff=1) +
  my_cladelab(cat_dat, "hoham|tughill", "Group C", cat_offset+0.15, cat_txtoff, ff=1) +
  my_cladelab(cat_dat, "cryophila|krienitzii", "Group B", cat_offset+0.15, cat_txtoff, ff=1) +
  my_cladelab(cat_dat, "sac|cog|petal|thickwall|raisin|reticulata|arctica", "Chloromonadinia", cat_offset+0.06, cat_txtoff, ff=1) +
  labs(title = "its2+18s+rbcL")+
  geom_treescale(x = cat_offset*3)
cat_treeplot


ggsave(paste0(here("figs/9_cat/"), today, "_cat.pdf"), 
       width = 6, height = 8, units = "in")
# ggsave(paste0(here("figs/trees/cat/"), today, "_cat.pdf"), 
#        width = 6, height = 8, units = "in")


# all trees on one panel------------------------

ggpubr::ggarrange(
  its2_treeplot,
  rbcl_treeplot,
  x18sgappy_treeplot,
  cat_treeplot,
  nrow = 2, ncol = 2
)
ggsave(paste0(here("figs/trees/"), today, "_tree_multipanel.pdf"),
       width = 13, height = 24, units = "in")

# 
# # 18s full --------------------
# x18sfull_treeplot <-
#   my_ggtree(x18sfull_dat) +
#   geom_text(aes(label=bootstrap_label_filt),  vjust= 0.01, hjust = 1, size=2.5) +
#   geom_tiplab(aes(label = newlabel, color = morphotype_color_group), size=4) +
#   xlim(0,0.1)+
#   scale_color_manual(values = my_colors) +
#   theme(legend.position = "none") +
#   my_cladelab(x18sfull_dat, "smallSphere|Sanguina", "Sanguina-clade", 0.02, 0.001, "black", ff=1) +
#   my_cladelab(x18sfull_dat, "Chlainomonas", "Chlainomonas-clade", 0.018, 0, "black", rotate = FALSE, ff=1) +
#   my_cladelab(x18sfull_dat, "cog|raisin", "Rosetta", 0.01, 0, "black") +
#   labs(title = "18s full")
# x18sfull_treeplot
# 
# ggsave(paste0(here("figs/trees/18s/"), today, "_18s_full.png"), 
#        width = 6, height =11, units = "in")
# ggsave(paste0(here("figs/trees/18s/"), today, "_18s_full.pdf"), 
#        width = 6, height = 11, units = "in")
# 
# 


# 18S A ----------------------------

x18sa_treeplot <-
  my_ggtree(x18sa_dat) +
  geom_text(aes(label=bootstrap_label_filt),  vjust= 0.01, hjust = 1, size=2.5) +
  geom_tiplab(aes(label = newlabel, color = morphotype_color_group), size=4) +
  xlim(0,0.15)+
  scale_color_manual(values = my_colors) +
  theme(legend.position = "none") +
  my_cladelab(x18sa_dat, "smallSphere|Sanguina", "Sanguina-clade", 0.02, 0.001, "black") +
  my_cladelab(x18sa_dat, "Chlainomonas", "Chlainomonas-clade", 0.02, 0, "black", rotate = FALSE) +
  my_cladelab(x18sa_dat, "cog|raisin", "Roseta", 0.02, 0, "black") +
  labs(title = "18S A")
x18sa_treeplot
# 
# ggsave(paste0(here("figs/trees/18s/"), today, "_18s_a.pdf"),
#        width = 6, height = 11, units = "in")
ggsave(paste0(here("figs/7_18s/"), today, "_18s_a_incl_sac_GAR1930.png"),
       width = 9.5, height =11, units = "in")

# 
# 
# 
# 
# # 18S B ----------------------------
# x18sb_treeplot <-
#   my_ggtree(x18sb_dat) +
#   geom_text(aes(label=bootstrap_label_filt),  vjust= 0.01, hjust = 1, size=2.5) +
#   geom_tiplab(aes(label = newlabel, color = morphotype_color_group), size=4) +
#   xlim(0,0.15)+
#   scale_color_manual(values = my_colors) +
#   theme(legend.position = "none") +
#   my_cladelab(x18sb_dat, "smallSphere|Sanguina", "Sanguina-clade", 0.02, 0.001, "black") +
#   my_cladelab(x18sb_dat, "cog|raisin", "Roseta", 0.02, 0, "black") +
#   labs(title = "18S B")
# x18sb_treeplot
# 
# 
# 
# ggsave(paste0(here("figs/trees/18s/"), today, "_18s_b_tree.pdf"),
#        width = 6, height = 8, units = "in")
# ggsave(paste0(here("figs/trees/18s/"), today, "_18s_b_tree.png"),
#        width = 6, height = 8, units = "in")
# 
# 







# # RIBOSOMAL concatenated sequence tree  ----------------------------
# 
# 
# ribo_ggtree <- make_ggtree(ribo_tree) # expect some error messages
# 
# ribo_ggtree +
#   geom_text(aes(label=conf_label),  vjust= 0.01, hjust = 1, size=2.5) +
#   geom_tiplab(aes(label = tip_label), size=4) +
#   xlim(0,1.5) +
#   ggtitle("its2+18S")
# 
# ggsave(paste0(here("figs/trees/ribo/"), today, "_ribo_tree.pdf"), 
#        width = 6, height = 8, units = "in")
# 
# 
# 







# SCRATCH ------------------------
# add a grouping column to identify sacs from which genus

  
  
  
# function to make a ggtree object
make_ggtree <- function(treedat){
  # convert to tidy object
  gg <- ggtree(treedat)
  # gg$data %>% count(isTip)
  
  # format tip labels
  tips <- gg$data %>% 
    filter(isTip == TRUE) %>% 
    separate(label,into = c("morphospecies", "sample_id", "lysate_tube_id"), sep = "_") %>% 
    mutate(tip_label = paste(morphospecies, sample_id, lysate_tube_id, sep = "_") %>% 
             str_replace("callous", "thickWall") %>% 
             str_replace("envelope", "sac") %>% 
             str_replace("sac","sac"),
           morpho_label = case_when(morphospecies=="petal"~"petal",
                                    (morphospecies=="sac" & str_detect(sample_id, "WED") )~"petal",
                                    morphospecies=="star"~"star",
                                    morphospecies %in% c("thickWall", "bigThickWall")~"thickwall",
                                    morphospecies=="martian"~"thickwall",
                                    morphospecies=="raisin"~"raisin",
                                    morphospecies=="cog"~"cog",
                                    morphospecies=="sac"~"cog",
                                    .default = "other"))
  # nrow(tips)
  
  # format bayesian/bootstrap confidence for each node
  boot_thresh <- 60 
  
  boots <- gg$data %>% 
    filter(isTip == FALSE) %>% 
    mutate(conf_label = as.numeric(label),
           conf_label = if_else(label<boot_thresh | is.na(label) | label=="Root", "", as.character(label))) 
  
  # # for output from online IQTree engine  
  # # split label into boots and bayes vals
  # separate(label, c("conf_boots", "conf_bayes"), sep = "/") %>%
  # # only show bootstrap vals greater than 70
  # mutate(conf_boots = as.numeric(conf_boots) %>% round(0),
  #        conf_bayes = as.numeric(conf_bayes),
  #        conf_label = if_else(conf_boots < boot_thresh | is.na(conf_boots), 
  #                        "",
  #                        conf_boots))#paste(conf_boots, conf_bayes, sep = "/")) )
  #          
  # nrow(boots)
  
  # update the gg$data 
  gg$data <- bind_rows(tips, boots)
  gg
}

# parse_label <- function(tree){
#   boots <- tibble(boots = tree$node.label)
#   tips <- tibble(label = tree$tip.label)
#   
#   suppressWarnings(
#     # omit bootstrap vals below threshold
#     newboots <- boots %>% 
#       mutate(val = as.numeric(boots),
#              display_val = case_when(val>boot_threshold ~ as.character(val),
#                                      val < boot_threshold ~ "",
#                                      .default = boots)) %>% 
#       pull()
#   )
#   suppressWarnings(
#     # tidy the tip labels
#     newlabs <- tips %>% 
#       # identify if node is internal, genbank, or from current study
#       mutate(type = case_when(str_detect(label, gb_regex) ~"genbank",
#                               (!str_detect(label, gb_regex) & str_detect(label, "_") ) ~ "mysample",
#                               str_detect(label, "[:digit:]+") ~ "bootstrap",
#                               .default = "other"),
#              genbank_label = if_else(type=="genbank", label, NA_character_),
#              my_sample_label = if_else(type=="mysample", label, NA_character_),
#              bootstrap_label = if_else(type=="bootstrap", as.numeric(label), NA)) %>% 
#       # omit bootstrap vals below threshold
#       mutate(bootstrap_label = if_else(bootstrap_label > boot_threshold, as.character(bootstrap_label), "")) %>% 
#       # parse data from each type of label
#       separate_wider_delim(cols = my_sample_label,delim = "_", 
#                            names = c("morphotype","sample_id", "tube_id"),
#                            too_many = "drop", too_few = "align_start") %>% 
#       separate_wider_delim(cols = genbank_label, delim = "_", 
#                            names = c("gba","genus", "species"),
#                            too_many = "drop", too_few = "align_start") %>% 
#       # reassemble new labels from selected components
#       mutate(newlabel = case_when(type=="genbank" ~ paste(gba, genus, species, sep="_"),
#                                   type=="mysample" ~ paste(morphotype, sample_id, tube_id, sep="_"),
#                                   type == "bootstrap"~ bootstrap_label,
#                                   .default = label)) %>% 
#       pull(newlabel)
#   )
#   
#   # assign the new components as tree properties
#   tree$node.label <- newboots
#   tree$tip.label <- newlabs
#   return(tree)
# }
# # polish: collapse identical sequences
# tips %>% 
#     count(morphospecies, sample_id, parent) %>%
#     filter(n>1)

# # select nodes to collapse
# tips %>% 
#   count(tip_label) %>% 
#   filter(n>1) %>% 
#   left_join(select(tips, tip_label, parent))
# 
# gg + geom_text(aes(label=node), hjust=-.3, size = 2)
# 
# nodes_to_collapse <-  %>% 
#   pull(parent)

# # remove the label data as well
# tips_to_remove <- tips %>% 
#   filter(parent %in% nodes_to_collapse) %>% 
#   group_by(parent) %>% 
#   filter(node == min(node)) %>% # choose one of the labels to keep (should be identical labels)
#   pull(node)
# 
# # remove the collapsed tip labels
# tips <- tips %>% filter(!node %in% tips_to_remove)
# 

# %>%
  # collapse(node = nodes_to_collapse)

# 
# +
#   guides(size="none") +
#   labs(tag="C") +
#   theme(legend.position = "none") +
#   geom_tippoint(aes(color= insignificant_edge, size=rel_abund), alpha=1, na.rm=TRUE) + 
#   theme(legend.position = c(0.89,0.7),
#         legend.title = element_blank(),
#         legend.key = element_blank()) +
#   scale_color_manual(values= colors) +