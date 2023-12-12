## genetic diversity analysis

## library
library(DBI)
library(tidyverse)
library(LEA)
library(tibble)
library(viridis)
library(ape)
library(ggtree)
library(aplot)
library(readxl)
library(cowplot)
library(StAMPP)
library(adegenet)
library(dartR)
library(PopGenReport)
library(hierfstat)
library(ggtree)
library(grid)
library(factoextra)
library(FactoMineR)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)



group9_col <- c(turbo(5), turbo(4, begin = 0.1, end = 0.9)[-3], "gray50")
names(group9_col) <- c("E", "O", "C", "A", "D", "R", "B", "G", "Vietnam")
group9_col <- group9_col[c("D", "C", "G", "A", "B", "O", "R", "E", "Vietnam")]
group6_col <- group9_col[!names(group9_col) %in% c("B", "R")]
names(group6_col) <- c("D", "C", "G", "A", "OB", "ER", "Vietnam")

group9_size <- c(rep(4, 8), 3)
names(group9_size) <- names(group9_col)
group6_size <- group9_size[!names(group9_size) %in% c("B", "R")]
names(group6_size) <- c("D", "C", "G", "A", "OB", "ER", "Vietnam")

group9_alpha <- c(rep(1, 8), 0.8)
names(group9_alpha) <- names(group9_col)
group6_alpha <- group9_alpha[!names(group9_alpha) %in% c("B", "R")]
names(group6_alpha) <- c("D", "C", "G", "A", "OB", "ER", "Vietnam")

group9_shape <- rep(16, 9)
names(group9_shape) <- names(group9_col)
group9_shape["R"] <- 17
group9_shape["B"] <- 18
group6_shape <- group9_shape[!names(group9_shape) %in% c("B", "R")]
names(group6_shape) <- c("D", "C", "G", "A", "OB", "ER", "Vietnam")


#################
## import data ##
#################
acc_info <- read.delim("./raw_genotype/all_accessions_info.tsv")
acc_info <- acc_info %>% arrange(code)
af_info <- acc_info %>% filter(group != "Vietnam")

all_geno <- read.delim("./raw_genotype/all_genotypes.tsv", row.names = 1)
all_geno <- all_geno[sort(rownames(all_geno)),]

plotdir <- "../plot_paper"
# plotdir <- "./plot_paper"
dir.create(plotdir, showWarnings = F)

poppr::missingno(all_genind)


##############
## assign African populations ##
##############
af_geno <- all_geno[!grepl("S-|TR", rownames(all_geno)),]
af_geno[is.na(af_geno)] <- 9
af_geno_file <- file.path(pca_dir, "af.geno")
write.geno(af_geno, af_geno_file)
af_snmf <- snmf(input.file = af_geno_file, K = 1:10, 
                 repetitions = 100, CPU = 20, 
                 entropy = T, seed = 987654321, project = "new")
plot(all_snmf)

af_prop_list <- sapply(5:8, function(k) {
  prop_5 <- Q(af_snmf, K = k, run = which.min(cross.entropy(af_snmf, K = k))) %>% as.data.frame()
  prop_5$code <- af_info$code
  prop_5 <- left_join(prop_5, af_info, by = "code")
  ind_order <- LEA::barchart(af_snmf, K = 8, run = which.min(cross.entropy(af_snmf, K = 8)), plot = F)
  # prop_5 <- prop_5[ind_order$order,]
  # prop <- prop[id_ord$order,]
  prop_5 <- pivot_longer(prop_5, cols = starts_with("V"), 
                         names_to = "ancestry", values_to = "proportion")
  prop_5$k <- k
  prop_5$code <- factor(prop_5$code, levels = af_info$code[ind_order$order])
  return(prop_5)
}, USE.NAMES = T, simplify = F)

af_snmf_prop <- do.call(rbind, af_prop_list)

ggplot(af_snmf_prop) +
  facet_grid(rows = vars(paste("K =", k)), 
             cols = vars(group),
             scales = "free", space = "free") +
  geom_col(aes(x = code, y = proportion, fill = ancestry), width = 1) +
  scale_fill_manual(values = as.character(group9_col)) +
  # theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, h = 1, v = .5),
        panel.background = element_rect(fill = "white"), 
        panel.grid = element_line(color = "gray", size = .1)) +
  scale_y_continuous(expand = c(0,0))

af_snmf_prop %>% 
  filter(k == 7) %>% 
  group_by(code) %>% 
  slice_max(proportion) %>% 
  mutate(assign_group = if_else(proportion > 0.7, ancestry, "hybrid")) 


all_geno[is.na(all_geno)] <- 9
all_geno_file <- file.path(pca_dir, "all.geno")
write.geno(all_geno, all_geno_file)
all_snmf <- snmf(input.file = all_geno_file, K = 1:10, 
                 repetitions = 100, CPU = 56, 
                 entropy = T, seed = 987654321, project = "new")
tiff(filename = file.path(plotdir, "snmf_kasp_cross_entropy.tiff"),
     width = 3, height = 2, units = "in", res = 600)
plot(all_snmf)
dev.off()


## get results for K=5:8
all_prop_list <- sapply(5:8, function(k) {
  prop_5 <- Q(all_snmf, K = k, run = which.min(cross.entropy(all_snmf, K = k))) %>% as.data.frame()
  prop_5$code <- acc_info$code
  prop_5 <- left_join(prop_5, acc_info, by = "code")
  ind_order <- LEA::barchart(all_snmf, K = k, run = which.min(cross.entropy(all_snmf, K = k)), plot = F)
  # prop_5 <- prop_5[ind_order$order,]
  # prop <- prop[id_ord$order,]
  prop_5 <- pivot_longer(prop_5, cols = starts_with("V"), 
                         names_to = "ancestry", values_to = "proportion")
  prop_5$k <- k
  return(prop_5)
}, USE.NAMES = T, simplify = F)

all_snmf_prop <- do.call(rbind, all_prop_list)

ggplot(all_snmf_prop) +
  facet_grid(rows = vars(paste("K =", k)), cols = vars(group),
             scales = "free", space = "free") +
  geom_col(aes(x = code, y = proportion, fill = ancestry), width = 1) +
  scale_fill_manual(values = as.character(group9_col)) +
  theme(axis.text.x = element_text(angle = 90, h = 1, v = .5),
        panel.background = element_rect(fill = "white"), 
        panel.grid = element_line(color = "gray", size = .1))

## -> k = 6 is the most fit

## plot K = 6 

#### reorder references
ref_ord6 <- all_snmf_prop %>% 
  filter(k == 6, group != "Vietnam") %>% 
  select(code, ancestry, proportion, group) %>% 
  pivot_wider(names_from = ancestry, values_from = proportion) %>% 
  rowwise() %>% 
  mutate(snmf_group = which.max(c(V1,V2,V3,V4,V5,V6))) %>% 
  group_by(group) %>% 
  arrange(desc(V3), V6, V1, V2, V4, V5, .by_group = T) %>% pull(code) %>% as.character()

#### reorder vn samples
vn_ord6 <- all_snmf_prop %>% 
  filter(k == 6, group == "Vietnam") %>% 
  select(code, ancestry, proportion, group) %>% 
  group_by(code) %>% 
  arrange(desc(proportion), .by_group = T) %>% 
  slice_head(n = 2) %>% 
  slice_tail(n =1) %>% 
  mutate(ancestry = if_else(proportion < 0.001, "V6", ancestry)) %>% 
  mutate(ancestry = factor(ancestry, levels = c("V1", "V5", "V2", "V3", "V4", "V6"))) %>%
  group_by(ancestry) %>%
  arrange(desc(proportion), .by_group = T) %>% pull(code)



ind_ord6 <- c(ref_ord6, vn_ord6)

#### match snmf ancestral group with prior ancestries
ancestry_ord6 <- all_snmf_prop %>% 
  filter(k == 6) %>% 
  # mutate(group = if_else(source == "sample", "Vietnam", group)) %>% 
  filter(!group %in% c("B", "R")) %>% 
  rowwise() %>% 
  mutate(new_ancestry = which(names(group9_col) %in% group)) %>% 
  group_by(ancestry) %>% 
  arrange(desc(proportion)) %>% 
  slice_head(n = 1) %>% 
  # mutate(new_ancestry = as.numeric(group)) %>% 
  select(ancestry, new_ancestry)


# left_join(all_snmf_prop, ancestry_ord, by = c("ancestry", "k"))
# anc_cols <- c(turbo(8, begin = 0.08), "gray40")
anc_cols <- group9_col
names(anc_cols) <- 1:9

#### barplot all inds
tiff(filename = file.path(plotdir, "snmf_6.tiff"),
     width = 14, height = 7, units = "in", res = 1200)
all_snmf_prop6 <- all_snmf_prop %>% 
  left_join(., ancestry_ord6, by = "ancestry") %>% 
  filter(k == 6) %>%
  mutate(code = factor(code, levels = ind_ord6)) %>% 
  mutate(group = fct_relevel(group, names(group9_col))) %>% 
  mutate(new_ancestry_name = names(group9_col[group])) %>% 
  mutate(new_ancestry_name = case_when(new_ancestry_name %in% c("O", "B") ~ "OB",
                                       new_ancestry_name %in% c("E", "R") ~ "ER",
                                       TRUE ~ new_ancestry_name)) %>% 
  mutate(new_ancestry_name = factor(new_ancestry_name, names(group6_col)))

snmf1 <- ggplot(all_snmf_prop6) +
  facet_grid(cols = vars(new_ancestry_name),
             scales = "free", space = "free") +
  geom_col(aes(x = code, y = proportion, fill = as.factor(new_ancestry)), 
           width = 1, show.legend = F) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_rect(fill = "white"), 
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        panel.border = element_rect(fill = NA, color = "gray", size = 1),
        panel.spacing.y = unit(0.04, "npc"),
        panel.spacing.x = unit(0.01, "npc")) +
  scale_fill_manual(values = anc_cols) +
  xlab("Individual") +
  ylab("Ancestry proportion") +
  scale_y_continuous(expand = c(0, 0))
dev.off()

#### barplot only reference
ref_snmf6 <- ggplot(all_snmf_prop6 %>% filter(group != "Vietnam")) +
  facet_grid(cols = vars(new_ancestry_name),
             scales = "free", space = "free") +
  geom_col(aes(x = code, y = proportion, fill = as.factor(new_ancestry)), 
           width = 1, show.legend = F) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_rect(fill = "white"), 
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 16),
        axis.title = element_text(size = 12),
        panel.border = element_rect(fill = NA, color = "gray", size = 1),
        panel.spacing.y = unit(0.04, "npc"),
        panel.spacing.x = unit(0.01, "npc")) +
  scale_fill_manual(values = anc_cols) +
  xlab("Individual") +
  ylab("Ancestry proportion") +
  scale_y_continuous(expand = c(0, 0))

#### barplot only vn
vn_snmf6 <- ggplot(all_snmf_prop6 %>% filter(group == "Vietnam")) +
  geom_col(aes(x = code, y = proportion, fill = as.factor(new_ancestry)), 
           width = 1, show.legend = F) +
  scale_fill_manual(values = anc_cols) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_rect(fill = "white"), 
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 16),
        axis.title = element_text(size = 12),
        panel.border = element_rect(fill = NA, color = "gray", size = 1),
        panel.spacing.y = unit(0.04, "npc"),
        panel.spacing.x = unit(0.01, "npc"),
        plot.title = element_text(size = 16, hjust = 0.5)) +
  xlab("Individual") +
  ylab("Ancestry proportion") +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("Vietnam")


## plot K = 5:8
all_snmf_prop5to8 <- all_snmf_prop %>% 
  left_join(., ancestry_ord6, by = "ancestry") %>% 
  # filter(k == 6) %>%
  mutate(code = factor(code, levels = ind_ord6)) %>% 
  mutate(group = fct_relevel(group, names(group9_col)))

ggplot(all_snmf_prop5to8) +
  facet_grid(cols = vars(group), rows = vars(k),
             scales = "free", space = "free") +
  geom_col(aes(x = code, y = proportion, fill = as.factor(new_ancestry)), 
           width = 1, show.legend = F) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_rect(fill = "white"), 
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        panel.border = element_rect(fill = NA, color = "gray", size = 1),
        panel.spacing.y = unit(0.04, "npc"),
        panel.spacing.x = unit(0.01, "npc")) +
  # scale_fill_manual(values = turbo(8)) +
  scale_fill_manual(values = as.character(anc_cols)) +
  xlab("Individual") +
  ylab("Ancestry proportion") +
  scale_y_continuous(expand = c(0, 0))


#############
## run PCA ##
#############
all_pca <- PCA(all_geno)

## get PC1 and PC2
all_pca_proj <- all_pca$ind$coord[,1:2]
colnames(all_pca_proj) <- c("PC1", "PC2")
all_pca_proj <- all_pca_proj %>% 
  as.data.frame() %>% 
  rownames_to_column("code") 

all_pca_proj <- left_join(all_pca_proj, acc_info, by = "code")

all_pca_proj <- all_pca_proj %>% 
  mutate(group = ifelse(is.na(group), "Vietnam", group)) %>%
  mutate(label = ifelse(group == "Vietnam", as.character(code), group)) %>% 
  # mutate(group = factor(group, levels = names(group9_col))) %>% 
  mutate(new_ancestry_name = names(group9_col[group])) %>% 
  mutate(new_ancestry_name = case_when(new_ancestry_name %in% c("O", "B") ~ "OB",
                                       new_ancestry_name %in% c("E", "R") ~ "ER",
                                       TRUE ~ new_ancestry_name)) %>% 
  mutate(new_ancestry_name = factor(new_ancestry_name, names(group6_col)))

eigen <- all_pca$eig[,1]


## plot PCA results
# tiff(filename = file.path(plotdir, "pca_all.tiff"),
#      width = 6, height = 5, units = "in", res = 1200)
# pca1 <- all_pca_proj %>% 
#   arrange(desc(group)) %>% 
#   ggplot() +
#   geom_point(aes(x = PC1, y = PC2, color = group, shape = group), size = 5) +
#   # geom_text(aes(x = PC1, y = PC2, label = label, color = group)) +
#   xlab(paste0("PC1 (", format(eigen[1]/sum(eigen)*100, digits = 2), "%)")) +
#   ylab(paste0("PC2 (", format(eigen[2]/sum(eigen)*100, digits = 2), "%)")) +
#   scale_color_manual(values = group9_col) + 
#   scale_shape_manual(values = 0:8) +
#   theme_minimal() + 
#   theme(panel.border = element_rect(color = "gray", fill = NA),
#         axis.title = element_text(size = 12))
# pca1

pca1 <- all_pca_proj %>%
  arrange(desc(group)) %>%
  ggplot() +
  geom_point(aes(x = PC1, y = PC2,
                 fill = new_ancestry_name), shape = 21, color = "black", alpha = 0.7, size = 3) +
  # geom_text(aes(x = PC1, y = PC2, label = label, color = group)) +
  xlab(paste0("PC1 (", format(eigen[1]/sum(eigen)*100, digits = 2), "%)")) +
  ylab(paste0("PC2 (", format(eigen[2]/sum(eigen)*100, digits = 2), "%)")) +
  # scale_color_manual(values = c(as.character(group6_col), "gold3")) +
  # scale_color_manual(values = c(rep("black", 9), "gold3")) +
  scale_fill_manual(values = group6_col) +
  guides(fill = guide_legend(title = "group")) +
  # scale_shape_manual(values = group6_shape) +
  # scale_size_manual(values = c(rep(3, 9), 2)) +
  # scale_alpha_manual(values = c(rep(1, 8), 0.7, 1)) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "gray", fill = NA),
        axis.title = element_text(size = 12))

tiff(filename = file.path(plotdir, "pca_all.tiff"),
     width = 8, height = 7, units = "in", res = 800)
pca1
dev.off()


###################
## NJ clustering ##
###################

## euclidean distance
all_geno[all_geno == 9] <- NA
all_dist <- dist(all_geno, method = "euclidean")

## nj tree
nj_tree <- ape::nj(all_dist)

## plot nj tree
# group6_col <- group9_col
# group6_col["R"] <- group6_col["E"]
# group6_col["B"] <- group6_col["O"]

nj1 <- ggtree(nj_tree, ladderize = F, layout = "ape", color = "gray") %<+% (all_snmf_prop6 %>% distinct(code, .keep_all = T)) +
  # geom_tiplab(aes(subset = grepl("S-", label), color = group), hjust = -0.5,
  #             linesize = .2, size = 5) +
  geom_tippoint(aes(color = new_ancestry_name, size = new_ancestry_name, alpha = new_ancestry_name)) +
  guides(color = guide_legend(title = "group"),
         size = guide_legend(title = "group"),
         alpha = guide_legend(title = "group")) +
  scale_color_manual(values = group6_col) +
  scale_size_manual(values = group6_size) +
  scale_alpha_manual(values = group6_alpha) +
  # coord_flip() +
  scale_x_reverse(expand = c(0.02,0.02)) +
  scale_y_reverse(expand = c(0.02,0.02))
nj1


#################
## african map ##
#################
## african coordinates
africa <- ne_countries(scale = "medium", continent = "africa", returnclass = "sf")

## reference coordiates
rob_countries_names <- c("Ivory Coast", "Cameroon", "Angola", 
                         "Gabon", "Central African Republic",
                         "Uganda", "Democratic Republic of the Congo")
rob_countries <- africa %>% filter(sovereignt %in% rob_countries_names)
rob_countries[8,] <- rob_countries[4,]
rob_countries$sovereignt[8] <- paste("North", rob_countries$sovereignt[4])
rob_groups <- sapply(rob_countries$sovereignt, function(n) {
  g <- case_when(n == "Uganda" ~ "O",
                 n == "Gabon" ~ "A",
                 n == "Ivory Coast" ~"D",
                 n == "Democratic Republic of the Congo" ~ "E",
                 n == "Central African Republic" ~ "B",
                 n == "Cameroon" ~"C",
                 n == "Angola" ~ "G",
                 n == "North Democratic Republic of the Congo" ~ "R")
  return(g)
})
group_coord <- sp::coordinates(as(rob_countries, "Spatial")) %>% as.data.frame()
group_coord[,3] <- rob_groups
colnames(group_coord) <- c("x_mean", "y_mean", "group")
group_coord$x_mean[group_coord$group == "R"] <- group_coord$x_mean[group_coord$group == "E"] + 3
group_coord$y_mean[group_coord$group == "R"] <- group_coord$y_mean[group_coord$group == "E"] - 7

## plot map
map_col <- group9_col[!names(group9_col) %in% c("R", "Vietnam")]
map_col["B"] <- map_col["O"]
names(map_col)[7] <- "ER"

map_col <- group9_col[!names(group9_col) == "Vietnam"]
map_col["B"] <- map_col["O"]
map_col["R"] <- map_col["E"]

afmap <- ggplot(data = africa) +
  geom_sf(fill = "white", color = "grey40", size = 0.3) +
  coord_sf(xlim = c(-15, 35), ylim = c(-18, 13), expand = TRUE) +
  geom_point(data = group_coord,
             aes(x = x_mean, y = y_mean, 
                 fill = group, color = group, size = group),
             alpha = 0.2, shape = 21, show.legend = F) +
  geom_text(data = group_coord,
            aes(x = x_mean, y = y_mean,
                label = group, color = group),
            size = 10, show.legend = F) +
  theme_void() +
  theme(legend.box.margin = margin(l = 20), 
        legend.key.size = unit(14, "pt"),
        legend.key.height = unit(15, "pt")) +
  scale_color_manual(values = map_col) +
  scale_fill_manual(values = map_col) + 
  scale_size_manual(values = c(25, 25, 25, 25, 30, 25, 20, 17)) + 
  scale_x_continuous(expand = c(0,0.02)) + 
  scale_y_continuous(expand = c(0,0))



###################
## combine plots ##
###################
blank <- ggplot() + theme_void()
map_tree <- plot_grid(afmap, nj1, nrow = 1, ncol = 2, 
                      rel_widths = c(0.4, 0.6), 
                      labels = c("A", "D"))
snmf_plots <- plot_grid(ref_snmf6, blank, vn_snmf6, nrow = 3, ncol = 1,
                        rel_heights = c(0.4, 0.05, 0.4),
                        labels = c("B", "", "C"))
# pdf(file.path(plotdir, "map_nj_snmf.pdf"),
#     width = 12, height = 12)
tiff(file.path(plotdir, "map_nj_snmf.tiff"),
    width = 10, height = 10, units = "in", res = 1200)
plot_grid(map_tree, blank, snmf_plots, nrow = 3, ncol = 1, 
          rel_heights = c(0.4, 0.05, 0.7))
dev.off()


## assign groups

all_snmf_prop6 <- all_snmf_prop6 %>% 
  group_by(code) %>% 
  slice_max(proportion) %>% 
  mutate(assign_group = if_else(proportion >= 0.7, ancestry, "hybrid")) %>% 
  mutate(assign_group = case_when(assign_group == "V1" ~ "ER",
                                  assign_group == "V2" ~ "D",
                                  assign_group == "V3" ~ "G",
                                  assign_group == "V4" ~ "OB",
                                  assign_group == "V5" ~ "A",
                                  assign_group == "V6" ~ "C",
                                  TRUE ~ "hybrid")) %>% 
  mutate(assign_group = if_else(group == "Vietnam", "Vietnam", assign_group)) %>% 
  mutate(code = as.character(code)) %>% 
  arrange(code)

all_snmf_prop6 %>% 
  group_by(assign_group) %>%
  summarise(n = n())


################
## statistics ##
################
## redundancy
# vn_all_diss <- poppr::diss.dist(all_genind[all_genind@pop == "Vietnam"])
# hist(vn_all_diss)
# vn_all_diss <- as.matrix(vn_all_diss)
# vn_all_diss[sapply(1:nrow(vn_all_diss), function(r) any(vn_all_diss[r,-r] == 0)),
#             sapply(1:ncol(vn_all_diss), function(c) any(vn_all_diss[-c,c] == 0))]


## 
all_df <- all_geno
all_df[all_df == 0] <- 11
all_df[all_df == 1] <- 12
all_df[all_df == 2] <- 22
all_genind <- df2genind(all_df, ncode = 1, pop = all_snmf_prop6$assign_group)
all_genind@strata <- all_snmf_prop6
# levels(all_genind@pop) <- c("D", "C", "G", "A", "OB", "ER", "Vietnam", "hybrid")


grp_genind <- all_genind[all_genind@pop != "hybrid"]
grp_df <- data.frame(pop = grp_genind@pop, all_df[all_genind@pop != "hybrid",])

all_sum <- summary(grp_genind)

ar <- allelic.richness(grp_df)
# ar_df <- colMeans(ar$Ar)
# ar_df[, c("D", "C", "G", "A", "OB", "ER", "Vietnam")]

all_bs <- basic.stats(grp_df)
colMeans(all_bs$Fis, na.rm = T)

all_geno[all_geno == 9] <- NA
Fis <- fis.dosage(all_geno[all_genind@pop != "hybrid",], pop = as.character(all_genind@pop[all_genind@pop != "hybrid"]))



bs_stats <- as.data.frame(rbind(as.numeric(sprintf("%d", all_sum$n.by.pop)),
                                as.numeric(sprintf("%.2f", colMeans(ar$Ar))),
                                as.numeric(sprintf("%.2f", colMeans(all_bs$Ho))),
                                as.numeric(sprintf("%.2f", colMeans(all_bs$Hs))),
                                as.numeric(sprintf("%.2f", Fis[colnames(ar$Ar)]))))
colnames(bs_stats) <- colnames(all_bs$Ho)
rownames(bs_stats) <- c("N", "AR", "Ho", "He", "Fis")

grp_gl <- gi2gl(grp_genind)
# levels(grp_gl@pop) <- c("D", "C", "G", "A", "OB", "ER", "Vietnam")

grp_fst <- stamppFst(grp_gl, percent = 99, nboots = 1000)

bs_stats <- rbind(bs_stats, t(grp_fst$Fsts))
# bs_stats[c("N", "AR", "Ho", "He", "Fis", "D", "C", "G", "A", "OB", "ER"), c("D", "C", "G", "A", "OB", "ER", "Vietnam")]

write.table(bs_stats, file = file.path(plotdir, "bs_stats.txt"))

### congo-related inds

cg_vn <- all_snmf_prop6 %>% 
  filter(group %in% c("R", "E") | (group == "Vietnam" & new_ancestry == 8 & proportion > 0.9)) 

cg_vn_genind <- all_genind[unique(cg_vn$code)]
cg_vn_dist <- dist(cg_vn_genind, method = "euclidean")
cg_vn_hclust <- hclust(cg_vn_dist, method = "ward.D")
plot(cg_vn_hclust)

cg_vn_tree <- nj(cg_vn_dist)

ggtree(cg_vn_tree, ladderize = T) +
  geom_tiplab(hjust = 1,  align = F, linesize = .2, nudge_x = .2, size = 2)

cg_vn_hclust$labels <- gsub("m$", "", cg_vn_hclust$labels)

tiff(filename = file.path(plotdir, "hclust_ER.tiff"),
     width = 12, height = 6, units = "in", res = 800)
ggtree(cg_vn_hclust, layout = "dendrogram") %<+% (acc_info %>% filter(code %in% cg_vn$code)) +
  geom_tippoint(aes(color = group), shape = 16, size = 2) +
  geom_tiplab(aes(subset = grepl("(S-|TR)", label)), angle = 90, hjust = 1, size = 2, nudge_x = -0.3) +
  guides(color = guide_legend(title = "group")) +
  scale_shape_manual(values = 0:8) +
  scale_color_manual(values = group9_col[unique(cg_vn$group)])
dev.off()


## core set

core6 <- read_xlsx("./core_selection/corehunter_selection_SH_EE_35_acc_cuttings_20211020.xlsx")
core6 <- core6 %>% 
  mutate(core6 = if_else(grepl("TR", code), "core", core6)) 

pca2 <- core6 %>% 
  select(code, core6) %>% 
  right_join(., all_pca_proj, by = "code") %>% 
  mutate(new_ancestry_name = if_else(!is.na(core6), "Vietnam core set", as.character(new_ancestry_name))) %>% 
  mutate(new_ancestry_name = factor(new_ancestry_name, levels = c("Vietnam", "D", "C", "G", "A", "OB", "ER", "Vietnam core set"))) %>% 
  arrange(new_ancestry_name) %>% 
  ggplot() +
  geom_point(aes(x = PC1, y = PC2,
                 fill = new_ancestry_name), shape = 21, color = "black", alpha = 0.7, size = 3) +
  # geom_text(aes(x = PC1, y = PC2, label = label, color = group)) +
  xlab(paste0("PC1 (", format(eigen[1]/sum(eigen)*100, digits = 2), "%)")) +
  ylab(paste0("PC2 (", format(eigen[2]/sum(eigen)*100, digits = 2), "%)")) +
  # scale_color_manual(values = c(as.character(group6_col), "gold3")) +
  # scale_color_manual(values = c(rep("black", 9), "gold3")) +
  scale_fill_manual(values = c(group6_col, "Vietnam core set" = "gold3")) +
  guides(fill = guide_legend(title = "group")) +
  # scale_shape_manual(values = group6_shape) +
  # scale_size_manual(values = c(rep(3, 9), 2)) +
  # scale_alpha_manual(values = c(rep(1, 8), 0.7, 1)) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "gray", fill = NA),
        axis.title = element_text(size = 12))


eval_core_6 <- readRDS("./R_data/eval_core_6.RDS")
tiff(filename = file.path(plotdir, "eval_core_6.tiff"),
     width = 8, height = 8, units = "in", res = 1200)
# pdf(file = file.path(plotdir, "eval_core_6.pdf"),
#     width = 8, height = 8)
eval_plots <- plot_grid(eval_core_6$plot[[5]], 
                        eval_core_6$plot[[7]],
                        eval_core_6$plot[[8]],
                        nrow = 1, ncol = 3, 
                        hjust = 0,
                        labels = c("A", "B", "C"))
plot_grid(eval_plots, pca2, 
          nrow = 2, ncol = 1, rel_heights = c(.7, 1),
          hjust = 0,
          labels = c("", "D"))
dev.off()
