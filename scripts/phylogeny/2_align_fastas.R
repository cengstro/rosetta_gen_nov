# align and QC fastas

# this script takes the sorted fasta files from 1_sort_fastas.py
# and aligns them using the DECIPHER aligner. 
# The output alignment files are written as fastas, which are stored singley in a folder within the phylo directory
# --the phylogenetic data using each alignment fasta will be added to that folder

# To do: how does IQTree handle gaps? indels? if a position if missing on one seq (N) can it still use the position for all other seqs?
# problem: need to clean up non-biological gaps before feeding to IQTree

#*
#*
#* TO DO
#* output 18s and cat tree should not delete tube ID/species name (puts NA)
#* use dataframes individually for each marker (apply with a fnctino), then combine at the end
#* format the output to include



library(here)
library(DECIPHER)
library(fs)
library(tidyverse)
library(janitor)
library(geodist)

today <- "2023-06-24" # for input and output

# read in the data --------------------------------

# genbank_seqs <- readDNAStringSet(here("data/seq/genbank/gb.fasta"))

# using the up-to-date fasta for each marker

# helper function to trim the length of each fasta header
shorten_names <- function(dnaStringSet, max_char = 38){
  names(dnaStringSet) <- names(dnaStringSet) %>% 
    str_sub(1,max_char) %>% 
    str_replace_all(" ", "_")
  return(dnaStringSet)
}



its2 <- readDNAStringSet(here(paste0("data/seq/its2_", today, ".fasta"))) %>% shorten_names()
x18s <- readDNAStringSet(here(paste0("data/seq/x18s_", today, ".fasta"))) %>% shorten_names()
rbcl <- readDNAStringSet(here(paste0("data/seq/rbcL_", today, ".fasta"))) %>% shorten_names()
# rbcl_primers <- readDNAStringSet(here("data/seq/rbcL_primers.fasta"))
# x18s_primers <- readDNAStringSet(here("data/seq/18s_primers.fasta"))
# its2_primers <- readDNAStringSet(here("data/seq/its2_primers.fasta"))

engstrom2020 <- readDNAStringSet(here("data/seq/engstrom2020CommonChloroASVs.fasta"))

selectTheseNames <- function(dss, names_to_keep){
  regex <- paste0(names_to_keep, sep="|", collapse='') %>% str_sub(end = -2)
  return( dss[str_detect(names(dss), regex)] )
}

removeTheseNames <- function(dss, names_to_discard){
  regex <- paste0(names_to_discard, sep="|", collapse='') %>% str_sub(end = -2)
  return( dss[!str_detect(names(dss), regex)] )
}

#geographical data
rosettes <- read_csv(here("data/field_sample_meta/rosette_field_meta.csv")) 
ak_coords <- read_csv(here("data/field_sample_meta/genbank_geo_dists.csv"))




# ITS2 ---------------------------------------------------


its2_a <- AlignSeqs(its2)
BrowseSeqs(its2_a, highlight = T)
its2_trim <- subseq(its2_a, 2440, 3420) # manually change these based on the alignment

# second round
its2_realign <- RemoveGaps(its2_trim) %>% 
  AlignSeqs()
BrowseSeqs(its2_realign, highlight = TRUE)

# # try adding Segawa et al data into this
# ss <- readDNAStringSet(here("data/seq/genbank/segawa_etal_2018/98OTUs.fasta"))
# cc <- c(its2_realign, ss)
# cc2 <- RemoveGaps(cc) %>% AlignSeqs()
# its2_realign <- cc2
# Chloro kasaie data is not aligning, discard
# its2_kasaiae <- c(removeTheseNames(its2_trim, "kasaiae"), selectTheseNames(its2, "kasaiae") ) %>% 
#   RemoveGaps() %>% 
#   AlignSeqs()
# BrowseSeqs(its2_kasaiae, highlight = TRUE)

# write out the final aligned fasta
its2_out_dir <- here(paste0("data/phylo/", today, "_its2/"))
dir_create(its2_out_dir) # does NOT clobber the old one by default
subseq(its2_realign, 2, 930) %>% 
  writeXStringSet(paste0(its2_out_dir, "aligns.fasta")) # DOES clobber by default

selectTheseNames(its2_realign, "GAR1930") %>%
  RemoveGaps() %>% 
  AlignSeqs() %>% 
  BrowseSeqs(highlight = T)

selectTheseNames(its2_realign, "GAR1930") %>%
  RemoveGaps() %>% 
  AlignSeqs() %>% 
  writeXStringSet(here("data/phylo/GAR1930alignments.fa"))

# # locate primers
# aa <- its2 %>% 
#   selectTheseNames("DL06") %>% 
#   c(its2_primers %>% selectTheseNames("chlor")) %>% 
#   AlignSeqs() 
# BrowseSeqs(aa)

# chloroITS2 2392 - 2410

## distm -----------------
# gb_regex <- "[:alpha:]{2}[:digit:]{6}.*"
# my_its2 <- removeTheseNames(its2_realign, gb_regex)

# subset rosettes
rosette_regex <- "petal|cog|hickWall|raisin|star|martian|friedEgg"
# rosette_regex <- "petal|friedEgg_WED2206" # all identical
# rosette_regex <- "cog|friedEgg_GAR1930|friedEgg_PAN1915"


# |AB903027|LC371435|LC371405|LC371434|AB903027|AB903026|LC371439"
my_rosettes <- selectTheseNames(its2_realign, rosette_regex)

# arrange alphabetically by sample ID to allow joining
my_rosettes <-
  my_rosettes[order(names(my_rosettes) %>% str_split_i("_", 2)),]

its2_distm <- DistanceMatrix(my_rosettes)
# remove diagonal and upper tri
its2_distm[upper.tri(its2_distm, diag = TRUE)] <- NA

genetic_distances <- its2_distm %>% 
  as_tibble(rownames = "name1") %>% #view()
  pivot_longer(-name1, names_to = "name2", values_drop_na = TRUE, values_to = "gendist") %>% 
  mutate(nn1 = str_split_i(name1, "_", 2),
         nn2 = str_split_i(name2, "_", 2),
         morphotype1 = str_split_i(name1, "_", 1),
         morphotype2 = str_split_i(name2, "_", 1), .keep = "unused") %>% #view() 
  dplyr::rename(name1 = nn1, name2 = nn2) %>% 
  # lump "sac" morphos with their respective types
  mutate(morphotype1 = case_when(
    morphotype1=="friedEgg" & name1 == "GAR1930" ~ "cog",
    morphotype1=="friedEgg" & name1 == "PAN1915" ~ "cog",
    morphotype1=="bigThickWall" ~ "thickWall",
    morphotype1=="martian" ~ "thickWall",
    .default = morphotype1),
    morphotype2 = case_when(
      morphotype2=="friedEgg" & name2 == "GAR1930" ~ "cog",
      morphotype2=="friedEgg" & name2 == "PAN1915" ~ "cog",
      morphotype2=="bigThickWall" ~ "thickWall",
      morphotype2=="martian" ~ "thickWall",
      .default = morphotype2)
  )
genetic_distances %>% arrange(gendist)


my_geodata <- rosettes %>% 
  filter((sample_id %in% genetic_distances$name1 | sample_id %in% genetic_distances$name2)) %>% 
  select(sample_id, lon, lat) %>% 
  # set the order the same, to allow for joining
  arrange(sample_id)

gdists <- my_geodata %>% 
  dplyr::rename(x = lon, y = lat) %>% 
  geodist(measure = "geodesic") # dist in meters
gdists[upper.tri(gdists, diag = TRUE)] <- NA
colnames(gdists) <- my_geodata$sample_id
rownames(gdists) <- my_geodata$sample_id
geo_dists <- gdists %>% 
  as_tibble(rownames = "name1") %>% 
  pivot_longer(-name1, names_to = "name2", values_drop_na = TRUE, values_to = "geodist") 


mydists <- genetic_distances %>% 
  full_join(geo_dists) %>% 
  drop_na() %>% 
  # only compare within the same morphotype
  filter(morphotype1==morphotype2) %>% 
  mutate(morphotype1 = fct_recode(morphotype1, "globe" = "thickWall", "poppy" = "raisin") %>% 
           fct_relevel(c("cog","petal","star","globe","poppy"))) 

my_colors <- c("petal" = "#377eb8", "star"="#984ea3", "globe"="#ff7f00", "poppy"="#a65628", "cog"="#3c893a")


mydists %>% 
  mutate(geodist = geodist/1000) %>% 
  ggplot(aes(x = geodist, y = gendist, color = morphotype1)) +
  geom_point() +
  scale_color_manual(values = my_colors) +
  labs(y = "Genetic distance", x = "Geographic distance (km)", color="Morphotype")
ggsave(here("figs/0_all_figs/fig_S8.tiff"))

mydists %>% 
  filter(morphotype1 == morphotype2) %>% 
  mutate(geodist = geodist/1000) %>% 
  ggplot(aes(x = geodist, y = gendist, color = morphotype1)) +
  geom_point() +
  facet_wrap(vars(morphotype1)) +
  scale_color_manual(values = my_colors) +
  labs(y = "Genetic distance", x = "Geographic distance (km)", color="Morphotype")# +
  # geom_text(aes(label = paste0(name1, name2)))
ggsave(here("figs/v5_all_figs/fig_S8.tiff"))



mydists %>% 
  filter(morphotype1 =="star")


  # geom_smooth(method ="lm")
  # geom_text(aes(label = paste(name1, name2)))
lm(gendist~geodist, data = mydists) %>% summary() 
# 1.2e-7, se = 5.8e-8, t=2, p = 0.06

cor(mydists$geodist, mydists$gendist) # Pearsons r = 0.1



# 18s  --------------------

# initial alignment
x18s_align <- x18s %>% AlignSeqs()
# BrowseSeqs(x18s_align, highlight = T)
x18s_trim <- subseq(x18s_align, 150, 1640) 
x18s_realign <- RemoveGaps(x18s_trim) %>% AlignSeqs()
# BrowseSeqs(x18s_realign, highlight = TRUE)


# just 18s A --------------
x18sA_aa <- selectTheseNames(x18s, "18SA|.1_") %>% #names()
  removeTheseNames("18SB") %>% 
  AlignSeqs()
BrowseSeqs(x18sA_aa)
x18sA_aa2 <- subseq(x18sA_aa, 105, 960) %>% 
  RemoveGaps() %>% 
  AlignSeqs()
BrowseSeqs(x18sA_aa2)

x18sa_out_dir <- here(paste0("data/phylo/", today, "_18sa/"))
dir_create(x18sa_out_dir)

x18sA_aa2 %>% 
  removeTheseNames("AlaskaHarding") %>% #names()
  writeXStringSet(paste0(x18sa_out_dir, "aligns.fasta"))




## concatenate 18S A and B -----------
nn_mat <- str_split_fixed(names(x18s_realign), pattern= "_", n=4)[,1:3] 
newnames <- paste(nn_mat[,1], nn_mat[,2], nn_mat[,3], sep  ="_") %>% 
  str_replace("gear", "cog")

x18s_tbl <- tibble(
  names = newnames,
  seq = as.vector(x18s_realign)
)

# combine like names
nesty <- x18s_tbl %>% 
  group_by(names) %>% 
  mutate(n = n()) %>%
  # filter(n>1) %>% 
  nest(data = c(seq)) %>% 
  ungroup()

nesty2 <- nesty %>%  
  mutate(dss = map(data, ~DNAStringSet(pull(.x))),
         cat = map(dss, ConsensusSequence, threshold = 0.5),
         catv = map(cat, as.vector)) %>% 
  unnest(catv) %>% 
  mutate(catv = catv %>% str_replace_all("\\+", "N"))

final18s <- DNAStringSet(nesty2$catv)
names(final18s) <- nesty2$names #%>% str_remove("_.+[:digit:]\\-[:digit:]+.+")
BrowseSeqs(final18s, highlight = T)
length(final18s)  


# troubleshoot a few outliers
# # Why the "NNNNNNR··NK·K····NNY" -- i thought intermediate bases weren't allowed?
# tt <- filter(nesty2, str_detect(names, "8-17L1"))$dss[[1]] %>% subseq(770,900)
# BrowseSeqs(tt)
# ConsensusSequence(tt, threshold = 0.5) %>% BrowseSeqs()

# xx <- filter(nesty2, str_detect(names, "cog_PAN1918_5-25T3"))$dss[[1]]
# BrowseSeqs(xx)
# ConsensusSequence(xx, threshold = 0.5) %>% BrowseSeqs()
# 
# xx <- filter(nesty2, str_detect(names, "5-13T1"))$dss[[1]]
# BrowseSeqs(xx)
# ConsensusSequence(xx, threshold = 0.5) %>% BrowseSeqs()
# xx <- filter(nesty2, str_detect(names, "thickWall"))$dss
# xx2 <- c(xx[[1]][2], xx[[3]], xx[[3]], xx[[4]], xx[[5]]) # sloppy coding
# 
# BrowseSeqs(xx2, highlight=TRUE)

final18s %>%
  selectTheseNames(names_to_keep = "thickWall|martian|bigThickWall") %>%
  BrowseSeqs(highlight = T)
# 
# xx <- filter(nesty2, str_detect(names, "1-11L1"))$dss[[1]]
# BrowseSeqs(xx, highlight=TRUE)

## write out full, gappy sequence -----------------
x18s_out_dir <- here(paste0("data/phylo/", today, "_18sgappy/"))
dir_create(x18s_out_dir)
writeXStringSet(final18s, paste0(x18s_out_dir, "aligns.fasta"))

# see if Alps are still outliers if I remove the gappy end bit
x18s_out_dir2 <- here(paste0("data/phylo/", today, "_18sgappy_shortalps/"))
dir_create(x18s_out_dir2)
writeXStringSet(subseq(final18s, 1, 880) %>% removeTheseNames("7-5L3"), paste0(x18s_out_dir2, "aligns.fasta"))


# ## split into A, B, complete alignments ---------------------------
# not_in_a <- c("devil_NIL1901", "honeycomb_BDW1906B", "squareBrevi_TEN2207","cog_BLU2203")
# not_in_b <- c("football_PAN1918", "nowall_BAK1918", "spiky_ECH2201", "bigSphere_WED2206")
# not_in_full <- c(not_in_a,not_in_b)
# 
# 
# 
# 
# full <- removeTheseNames(final18s, not_in_full)
# # BrowseSeqs(full)
# length(full)
# aa <- removeTheseNames(final18s, not_in_a) %>% subseq(10,820)
# # BrowseSeqs(aa)
# length(aa)
# 
# bb <- removeTheseNames(final18s, not_in_b) %>% subseq(960, 1480)
# # BrowseSeqs(bb)
# length(bb)
# 
# # write out the final aligned fasta
# 
# 
# x18full_out_dir <- here(paste0("data/phylo/", today, "_18sfull/"))
# dir_create(x18full_out_dir)
# writeXStringSet(full, paste0(x18full_out_dir, "aligns.fasta"))
# 
# x18sa_out_dir <- here(paste0("data/phylo/", today, "_18sa/"))
# dir_create(x18sa_out_dir)
# writeXStringSet(aa, paste0(x18sa_out_dir, "aligns.fasta"))
# 
# x18sb_out_dir <- here(paste0("data/phylo/", today, "_18sb/"))
# dir_create(x18sb_out_dir)
# writeXStringSet(bb, paste0(x18sb_out_dir, "aligns.fasta"))
# 
# 
# x18s_aa <- readDNAStringSet(here("data/phylo/2023-04-23_18sgappy/aligns.fasta"))
# x18s_aa %>% 
#   RemoveGaps() %>% 
#   writeXStringSet(here("data/phylo/2023-04-23_18sgappy/x18s_nogap.fasta"))


## distm -----------------------
# are there BP differences between the P6618 clade and the main Cog clade?
selectTheseNames(final18s, "cog") %>% 
  RemoveGaps() %>% 
  AlignSeqs() %>% 
  BrowseSeqs(highlight=T)
# 
# x18s_mysamples <- selectTheseNames(final18s, rosette_regex)
# 
# x18s_mysamples <-
#   x18s_mysamples[order(names(x18s_mysamples) %>% str_split_i("_", 2)),]
# 
# x18s_distm <- DistanceMatrix(x18s_mysamples)
# # remove diagonal and upper tri
# x18s_distm[upper.tri(x18s_distm, diag = TRUE)] <- NA
# x18s_genetic_distances <- x18s_distm %>% 
#   as_tibble(rownames = "name1") %>% 
#   pivot_longer(-name1, names_to = "name2", values_drop_na = TRUE, values_to = "gendist") %>% 
#   mutate(nn1 = str_split_i(name1, "_", 2),
#          nn2 = str_split_i(name2, "_", 2),
#          morphotype1 = str_split_i(name1, "_", 1),
#          morphotype2 = str_split_i(name1, "_", 1), .keep = "unused") %>% 
#   dplyr::rename(name1 = nn1, name2 = nn2) %>% 
#   # lump "sac" morphos with their respective types
#   mutate(morphotype1 = case_when(
#     morphotype1=="friedEgg" & name1 == "GAR1930" ~ "cog",
#     morphotype1=="friedEgg" & name1 == "PAN1915" ~ "cog",
#     morphotype1=="bigThickWall" ~ "thickWall",
#     morphotype1=="martian" ~ "thickWall",
#     .default = morphotype1),
#     morphotype2 = case_when(
#       morphotype2=="friedEgg" & name2 == "GAR1930" ~ "cog",
#       morphotype2=="friedEgg" & name2 == "PAN1915" ~ "cog",
#       morphotype2=="bigThickWall" ~ "thickWall",
#       morphotype2=="martian" ~ "thickWall",
#       .default = morphotype2)
#   )
# x18s_genetic_distances



# rbcL ---------------------------------------------------
# commented out primer design stuff
# primers <- readDNAStringSet(here("data/seq/genbank/f0r3primers.fasta"))
# rbcl <- c(rbcl,primers)
# otuD <- readDNAStringSet(here("data/seq/genbank/otu_d_all_asvs.fasta"))
# names(otuD) <- paste0("asv", names(otuD))
# rbcl <- c(rbcl,otuD)

rbcl_a <- AlignSeqs(rbcl)
BrowseSeqs(rbcl_a, highlight = T)
# ConsensusSequence(rbcl_a, threshold=0.5) %>% as.character() # 441. For f0r3 seq primer design
# GCAGGTTTTAAAGCTGGTGTAAAAGATTATCGTTTAACATATTACACTCCAGATTACGTTGTAAAAGATACAGACATTCTTGCKGCTTTCCGTATGACTCCWCAAGCAGGTGTTCCAATTGAAGAAGCTGGTGCTGCTGTTGCTGCTGAATCTTCTACAGGTACTTGGACAACTGTATGGACTGATGGTTTAACAAGTCTTGACCGTTACAAAGGTCGTTGTTACGATATCGAACCAGTTGCWGGTGAAGAAAACCAATACATTGCTTACGTTGCTTACCCTATTGACCTTTTTGAAGAAGGTTCTGTAACTAACTTATTTACATCTATTGTAGGTAACGTTTTTGGTTTCAAAGCTCTTCGTGCTCTACGTCTTGAAGATTTACGTATTTCTTGTGCATATGCTAAAACATTCCAAGGACCTCCTCACGGGAT

BrowseSeqs(selectTheseNames(rbcl_a, "thickwall|martian|thickWall|bigThickwall"), highlight = T)


rbcl_out_dir <- here(paste0("data/phylo/", today, "_rbcl/"))
dir_create(rbcl_out_dir)
writeXStringSet(rbcl_a, paste0(rbcl_out_dir, "aligns.fasta"))


## Rbcl introns--------------------
rbcl_insertions <- subseq(rbcl_a, 680, 1040) %>%
  RemoveGaps()
# rbcl_insertions[width(rbcl_insertions)>200] %>% # set to 130 to remove N
#   names()
# rbcl_insertions[width(rbcl_insertions)<200] %>%
#   names()
# all stars, petals have the insertion
# thickwalls, cogs do not have it, except for P6618 and TEM2003

insertion_alignments <- rbcl_insertions[width(rbcl_insertions)>200] %>%
  as.vector() %>% 
  str_replace_all("N", "") %>% 
  DNAStringSet() %>% 
  AlignSeqs()
names(insertion_alignments) <- names(rbcl_insertions[width(rbcl_insertions)>200])
DistanceMatrix(insertion_alignments)
BrowseSeqs(insertion_alignments, highlight = TRUE)
353


BrowseSeqs(insertion_alignments, highlight = T)
rbcl_out_dir2 <- here(paste0("data/phylo/", today, "_rbcl_insert/"))
dir_create(rbcl_out_dir2)
writeXStringSet(insertion_alignments, paste0(rbcl_out_dir2, "aligns.fasta"))

dm <- insertion_alignments %>% DistanceMatrix()


# convert to RNA for alignment in RNAfold web server
insertion_alignments[8] %>% RNAStringSet() %>% RemoveGaps() %>% view()
  

# translate(RemoveGaps(insertion_alignments[1]))

## rbcl with introns removed ---------------------
rbcl_trim1 <- subseq(rbcl_a, 230, 679)
# rbcl_comb <- rbcl_trim1
# 419
# rosettes (but not envelopes) have a insertion in their rbcL?
rbcl_trim2 <- subseq(rbcl_a, 1039, 1300)
# 255
rbcl_comb <- DNAStringSet(paste0(rbcl_trim1, rbcl_trim2))

names(rbcl_comb) <- names(rbcl_trim2)
# #674
# rbcl_realign <- RemoveGaps(rbcl_comb) %>% 
#   AlignSeqs()
# BrowseSeqs(rbcl_realign, highlight = TRUE)



rbcl_out_dir2 <- here(paste0("data/phylo/", today, "_rbcl_noinsert/"))
dir_create(rbcl_out_dir2)
writeXStringSet(rbcl_comb, paste0(rbcl_out_dir2, "aligns.fasta"))


# # are there BP differences between the P6618 clade and the main Cog clade?
# selectTheseNames(rbcl_a, "cog") %>% 
#   RemoveGaps() %>% 
#   AlignSeqs() %>% 
#   BrowseSeqs(highlight=T)


## including Engstrom 2020 Illumina data -----------

# only showing ASVs with >1% relative abundance in at least one sample
engstrom2020 <- engstrom2020 %>% 
  selectTheseNames("ChloromonasD|ChloromonasF")
names(engstrom2020)

engs2020_aligns <- c(rbcl_a, engstrom2020) %>% RemoveGaps() %>% 
  AlignSeqs()


rbcl_out_dir3 <- here(paste0("data/phylo/", today, "_rbcl_plus_2020_asvs/"))
dir_create(rbcl_out_dir3)
writeXStringSet(engs2020_aligns, paste0(rbcl_out_dir3, "aligns.fasta"))

# engs2020_aligns %>% BrowseSeqs(highlight = T)
# engs2020_sub <- subseq(engs2020_aligns, 504, 690)# %>% BrowseSeqs(highlight = T)

# concatenated sequence --------------------------------


# include "big" thickwall in the concatenated sequence
pattern <- "bigThickWall_BLU2203_9-28L1"
replacement <- "thickWall_BLU2203_9-28L1"
ind <- str_which(names(final18s),pattern)
names(final18s)[ind] <- replacement
names(final18s)


# convert to tbl
ll <- list(its2_realign, final18s, rbcl_a)
names(ll) <- c("ITS2", "18S", "rbcL")

dd <- ll %>% 
  map(~enframe(as.character(.x, use.names=TRUE))) %>% 
  enframe() %>%
  mutate(marker = name, .keep="unused") %>% 
  unnest() %>% 
  # identify genbank seqs with regex
  mutate(name = str_extract(name, "[:alpha:]+_[:alnum:]+"))  


dd$name
dd %>% 
  count(marker)

# remove redundant labels (or labels that look identical to start)
dd %>% 
  count(name, marker) %>% 
  arrange(-n)

dd2 <- dd %>% 
  group_by(name, marker) %>% 
  slice_head(n=1) %>% 
  ungroup()


# remove seqs with gaps
all_markers_present <- dd2 %>% 
  count(name) %>% 
  filter(n>=3, !is.na(name)) %>% 
  pull(name)
print(all_markers_present)
length(all_markers_present)

# subset with seq for all 3 markers
finalcat <- dd2 %>% 
  filter(name %in% all_markers_present)

# concatenate
cats <- finalcat %>% 
  select(marker, name, value) %>% 
  pivot_wider(names_from = "marker") %>% 
  unite("seq", c("ITS2", "18S", "rbcL"), sep = "") %>% 
  filter(name!="<NA>") %>% 
  deframe() %>% 
  DNAStringSet()
length(cats)

names(cats) %>% str_count("BLU2203")
names(cats) %>% str_count("raisin")


cat_out_dir <- here(paste0("data/phylo/", today, "_cat/"))
dir_create(cat_out_dir)
writeXStringSet(cats, paste0(cat_out_dir, "aligns.fasta"))



# get tube IDs of all sequences in trees ---------------
gb_regex <- "[:alpha:]{2}[:digit:]{6}.*"


names_df <- tibble(name = names(its2), gene = "ITS2") %>% 
  bind_rows(tibble(name = names(rbcl), gene = "rbcL")) %>% 
  bind_rows(tibble(name = names(x18s), gene = "18S")) %>% 
  filter(!name %>% str_detect(gb_regex)) %>% 
  separate(name, c("morphotype", "sample_id", "tube_id", "marker2"), sep = "_") %>%
  mutate(gene2 = gene) %>% 
  select(-marker2) %>% 
  nest_by(gene2)
names_df$data[[1]] %>% 
  full_join(names_df$data[[2]], by = c("morphotype", "sample_id","tube_id")) %>% 
  full_join(names_df$data[[3]], by =c("morphotype", "sample_id", "tube_id")) %>% 
  dplyr::rename(x18s = gene.x, its2 = gene.y, rbcl = gene) %>% 
  arrange(tube_id) %>% 
  distinct() %>% view()
  # write_csv(here("data/seq/genbank_submissions/tube_ids_in_trees.csv"))

names(rbcl) %>% 
  selectTheseNames("3-12")

# get seq length stats ----------------
library(officer)
library(flextable)

getStats <- function(dss, lab){
  # dss=final18s
  nogaps <- dss %>% 
    DECIPHER::RemoveGaps() %>% 
    as.character() %>% 
    str_remove_all("N") %>% 
    DNAStringSet()
  # Get names of min width seqs
  nn <- tibble(width  = BiocGenerics::width(nogaps),
               names = names(dss)) %>%
    arrange(width)
  
  minname <- nn %>% 
    slice_head(n=1) %>% 
    pull(names)

  statz <- nogaps %>% 
    BiocGenerics::width() %>% 
    enframe() %>% 
    summarise(min_bp = min(value), median_bp = median(value), max_bp = max(value)) %>% 
    tibble::add_column("tree name" = lab, .before = 1) %>% 
    add_column("min seq name" = minname)
  
  return(statz)
}



tt <- getStats(its2_realign, "ITS2") %>% 
  bind_rows(getStats(final18s, "18S")) %>% # gappy
  bind_rows(getStats(rbcl_realign, "rbcL")) %>% 
  bind_rows(getStats(cats, "concatenated")) %>% 
  mutate(across(where(is.numeric), round))
tt


flextable(tt) %>% 
  save_as_docx(path = here("figs/trees/n_bases_per_tree.docx"))

# ribosomal contatendted sequenct -----------------------


# convert to tbl
ll3 <- list(its2_realign, x18s_a_realign)
names(ll3) <- c("ITS2", "18SA")

dd3 <- ll3 %>% 
  map(~enframe(as.character(.x, use.names=TRUE))) %>% 
  enframe() %>%
  mutate(marker = name, .keep="unused") %>% 
  unnest() %>% 
  mutate(gba = name %>% str_detect(gg)) %>% 
  mutate(type = if_else(!gba, 
                        paste(str_split_fixed(name, "_", 4)[,1],
                              str_split_fixed(name, "_", 4)[,2],
                              sep="_"),
                        #  if it's a genbank entry...
                        str_extract(name, spp))) # this is a terrible way to do this, better edit source data manually

dd3$name
dd3 %>% 
  count(marker)

# remove duplicate entries
dd4 <- dd3 %>% 
  group_by(type, marker) %>% 
  slice_head(n=1) %>% 
  ungroup()
nrow(dd4)/2 # 53

all_markers_present4 <- dd4 %>% 
  count(type) %>% 
  filter(n==2, !is.na(type)) %>% 
  pull(type)
length(all_markers_present4) # 42

# subset with seq for all 2 markers
final4 <- dd4 %>% 
  filter(type %in% all_markers_present4)

# concatenate
cats4 <- final4 %>% 
  select(marker, type, value) %>% 
  pivot_wider(names_from = "marker") %>% 
  unite("seq", c("ITS2", "18SA"), sep = "") %>% 
  deframe() %>% 
  DNAStringSet()

# later: add back in annotation of GBAs used to make concat sequence, etc
cat_out_dir4 <- here(paste0("data/phylo/", today, "_ribocat/"))
dir_create(cat_out_dir4)
writeXStringSet(cats4, paste0(cat_out_dir4, "aligns.fasta"))





# 
# 
# 
# ## SCRATCH -------------------------------------------
# 



# just 18SA ------------------------------------

my_18s_a <- x18s[names(x18s) %>% str_detect("A$")]
my_refs <- x18s[names(x18s) %>% str_detect("^[:alpha:]{2}[:digit:]{6}")]
comb_18sa <- c(my_18s_a,my_refs) %>% AlignSeqs()
BrowseSeqs(comb_18sa, highlight = T)
comb_18sa_trim <- subseq(comb_18sa, 102, 968) # manually change these based on the alignment

# second round
x18s_a_realign <- RemoveGaps(comb_18sa_trim) %>% AlignSeqs()
BrowseSeqs(x18s_a_realign, highlight = TRUE)

# write out the final aligned fasta

# 18s ----------------------------------------------------



# THIS  has a BUG -- see below ****************************************************

# First, combine the two paired reads into the final 18S contig

# split by name into separate DNAStringSets
# trim the trailing A or B off each fasta header
names(x18s) <- names(x18s) %>% str_remove("(?<=18S)[A|B]$")

# get a vector of the unique names
unique_names <- names(x18s) %>% unique()

# split the DNAStringSet by name, and align each pair
splitty <- unique_names %>% 
  map(function(n){
    dna <- x18s[names(x18s)==n]
    if(length(dna)>1){ # if there are two sequences, align them
      return(AlignSeqs(dna, verbose=FALSE)) ### BUG: what if the alignments are no good?
    }else{
      return(dna)
    }
  })
i <- 2
BrowseSeqs(splitty[[i]])
# this will try and force an alignment even if there is no alignment BUG*********************^^^

# make a consensus for each pair
consensus_list <- splitty %>% map(function(dss){
  dss %>% 
    ConsensusSequence(threshold=0.4) # lower the threshold to dela with N's
})
# # sanity check:  compare individual vs consensus at position i

BrowseSeqs(splitty[[i]])
BrowseSeqs(consensus_list[[i]])

# combine back into a single DNAStringSet
final_18s <- Reduce(c, consensus_list) # unlist dosen't work here
names(final_18s) <- unique_names # assuming the order has not changed


# Second, now proceed as with the others to align the samples, 

x18s_a <- AlignSeqs(RemoveGaps(final_18s))
length(x18s_a)
BrowseSeqs(x18s_a, highlight = TRUE)

# # since not all entries have both portions of 18S, split into A and B
# x18sA <- subseq(x18s_a[c(8:10, 13:33)], 103, 950) 
# x18sB <- subseq(x18s_a[c(1:7, 11:14, 25:33)], 1186, 1652) 

x18sA <- subseq(x18s_a[c(7, 10:12, 15:39)], 103, 950) # start at 103
x18sB <- subseq(x18s_a[c(1:7, 11:14, 25:97)], 1186, 1652) 

fix_plus_char <- function(dss){
  nm <- names(dss)
  dss2 <- dss %>% 
    str_replace_all("\\+", "G") %>% 
    DNAStringSet()
  names(dss2) <- nm # replace the names
  return(dss2)
}
# c("c++", "c+") %>% str_replace_all("\\+", "-")

x18sA_align <- AlignSeqs(x18sA %>% RemoveGaps()) %>% fix_plus_char()
x18sB_align <- AlignSeqs(x18sB %>% RemoveGaps()) %>% fix_plus_char()
# x18sA_align %>% str_detect("\\+")

BrowseSeqs(x18sA_align, highlight = TRUE)
BrowseSeqs(x18sB_align, highlight = TRUE)

x18sa_out_dir <- here(paste0("data/phylo/", today, "_x18sa/"))
x18sb_out_dir <- here(paste0("data/phylo/", today, "_x18sb/"))
dir_create(x18sa_out_dir)
dir_create(x18sb_out_dir)
writeXStringSet(x18sA_align, paste0(x18sa_out_dir, "x18sA_align.fasta"))
writeXStringSet(x18sB_align, paste0(x18sb_out_dir, "x18sB_align.fasta"))

# the only bases that differ between Jouds gear and cog 18S are downstream of nt=1660
# this is outside the captured region of my snow18S primer. 



# # include an annotated reference: 
# # ITS1 1794..2036
# # 5.8S 2037..2195
# # ITS2: 2196..2401
# 
# sanguina_nivaloides <- genbank_seqs[names(genbank_seqs) %>% str_detect("GU117577.1")]
# x18s_s <- AlignSeqs(c(sanguina_nivaloides, final_18s))
#   # BUG: the longer sequence is shifted off to one side. The shorter seqs are not padded.
#   # x18s_a %>% subseq(1000,1200) %>% RemoveGaps() %>% AlignSeqs() %>% BrowseSeqs()
# # ^still not working
# 
# BrowseSeqs(x18s_s, highlight=T)
# 
# 


# check which samples have data for all markers
ff <- function(str) {
  splitz <- str_split_fixed(names(str), "_", 4)
  # return(splitz[,3])
  return(paste(splitz[,1],splitz[,2], sep="_"))
}
intersect_names <- Reduce(intersect, list(ff(its2_realign),ff(x18s_a_realign),ff(rbcl_realign))) %>% sort()
# inner_join might be better here




# trim the _rbcL, GBA
trim_chr <- function(nn){
  names(nn) %>% 
    str_remove("_[:alnum:]+$") %>% # remove last delim
    str_remove(gba) %>% 
    str_remove("_18S") %>% 
    str_remove("(?<=miwae)_[:graph:]+$")# tidy up miwae trailing delim
}

its2_names_stripped <- trim_chr(its2_realign) %>% print()
x18sa_names_stripped <- trim_chr(x18s_a_realign) %>% print()
rbcl_names_stripped <- trim_chr(rbcl_realign) %>% print()
# Segawa et al lack AK, Antarctica single cell data for rbcL
# lack muramotoi 18s

# get list present in all 3 sets
in_all <- Reduce(intersect, list(its2_names_stripped, x18sa_names_stripped, rbcl_names_stripped)) %>% print()

# for each in list, stitch together 18sa, its2, rbcl
dna_tbl <- tibble(
  name = in_all,
  its2 = its2_realign[its2_names_stripped %in% in_all] %>% as.character(),
  x18sa = x18s_a_realign[x18sa_names_stripped %in% in_all] %>% as.character(),
  rbcl = rbcl_realign[rbcl_names_stripped %in% in_all] %>% as.character()
)
cats <- dna_tbl %>% 
  unite("cat", c("its2", "x18sa", "rbcl"), sep = "") %>% 
  deframe() %>% 
  DNAStringSet()





# # TrimDNA(aligns, "------", "------") %>% DNAStringSet()
# # trim between 1000 and 1400
# # automated way to do this?
# sub <- subseq(aligns, 2359, 3353)
# BrowseSeqs(sub, highlight = TRUE)
# 
# 
# # only seqs with ungapped length>500 bp
# ungapped <- RemoveGaps(sub)
# sub[width(ungapped)>500]  %>% 
#   writeXStringSet(here("data/phylo/its2_aligns_v2_gt500.fasta"))
# 
# # check just callous
# callous <- sub[str_detect(names(sub), "callous")]
# callous %>% 
#   RemoveGaps() %>% 
#   AlignSeqs() %>% 
#   BrowseSeqs(highlight = TRUE)
# 
# # check just cog/gear
# cog <- sub[str_detect(names(sub), "gear|cog")]
# cog %>% 
#   RemoveGaps() %>% 
#   AlignSeqs() %>% 
#   BrowseSeqs(highlight = TRUE)
# 
# # only the segment for which Joud got data
# subsub <- subseq(sub, 105, 305)
# writeXStringSet(subsub, here("data/phylo/its2_aligns_v2_joud.fasta"))

