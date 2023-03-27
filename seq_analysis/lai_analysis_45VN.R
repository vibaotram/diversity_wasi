
library(LEA)
library(Rsamtools)
library(tidyverse)
library(viridis)
library(Rcpp)
library(data.table)
library(karyoploteR)
library(gt)
library(gridExtra)
library(ggpubr)
library(aplot)
# library(S4Vectors)


# combine elai sets


all_chr <- GRanges()
all_chr <- do.call(GRangesList, sapply(c(1:11), function(i) {
  chr <- readRDS(paste0("/shared/projects/elai_most/data/snakelai_out_all/elai3runs_45VN_merged/chrom_", sprintf("%02d", i), ".RDS"))
  return(chr)
}, simplify = T)) %>% unlist
seqlevels(all_chr) <- gsub("CC1.8.Chr", "Chr ", seqlevels(all_chr))
all_chr$individual <- gsub(".", "-", all_chr$individual, fixed = T)
all_chr$individual <- gsub("m", "", all_chr$individual, fixed = T)




all_windows_large <- all_chr[width(all_chr) > 1e6,]


group_col <- c("#30123BFF", "#28BBECFF", "#A2FC3CFF", "#FB8022FF", "#7A0403FF")
names(group_col) <- c("ER", "OB", "C", "AG", "D")

# genome
genome_fa <- "/shared/projects/vietcaf/data/reference/CC1.8_v2_pseudomolecule_cat.fa"
genome <- scanFa(genome_fa)
genome <- genome[grepl("Chr", names(genome))] # remove contigs
names(genome) <- gsub(" .+", "", names(genome))
names(genome) <- gsub("CC1.8.Chr", "Chr ", names(genome))
genome_GR <- GRanges(seqnames = names(genome), ranges = IRanges(start = 1, end = width(genome)))



plotdir <- "../plot_paper"

## summary table for overall proportion

global <- do.call(rbind, lapply(unique(all_windows_large$individual), function(i) {
  # cat("> individual:", i, "\n")
  all_windows_ind <- all_windows_large[all_windows_large$individual == i,]
  all_windows_ind <- all_windows_ind[width(all_windows_ind) > 1e6]
  grps <- mcols(all_windows_ind) %>% as.data.frame() %>% dplyr::select(-individual) %>% colnames()
  glob_adm <- do.call(rbind, lapply(grps, function(g) {
    pw <- sum(mcols(all_windows_ind)[,g]*width(all_windows_ind))/sum(width(genome_GR))
    return(pw)
  }))
  glob_adm <- as.data.frame(t(glob_adm))
  colnames(glob_adm) <- grps
  rownames(glob_adm) <- i
  return(glob_adm)
}))

global <- global[stringr::str_sort(rownames(global), numeric = T), stringr::str_sort(colnames(global))]
undetermined <- 1 - rowSums(global)
global$undetermined <- undetermined
colnames(global) <- gsub("group_", "ancestry ", colnames(global))
rownames(global) <- gsub("m", "", rownames(global))
rownames(global) <- gsub("\\.", "-", rownames(global))

t <- global %>% 
  as.data.frame() %>% 
  gt(rownames_to_stub = T) %>% 
  fmt_number(columns = everything(), decimals = 3) %>% 
  data_color(colors = scales::col_numeric(palette = c("white", group_col["AG"]), domain = c(0,1)),
             columns = "AG") %>% 
  data_color(colors = scales::col_numeric(palette = c("white", group_col["ER"]), domain = c(0,1)),
             columns = "ER") %>% 
  data_color(colors = scales::col_numeric(palette = c("white", group_col["OB"]), domain = c(0,1)),
             columns = "OB") %>% 
  data_color(colors = scales::col_numeric(palette = c("white", group_col["D"]), domain = c(0,1)),
             columns = "D") %>% 
  data_color(colors = scales::col_numeric(palette = c("white", group_col["C"]), domain = c(0,1)),
             columns = "C") %>% 
  data_color(colors = scales::col_numeric(palette = c("white", "black"), domain = c(0,1)),
             columns = "undetermined") 
gtsave(t, file.path(plot_dir, "global_ancestry_elai.pdf"))
gt <- as_ggplot(t, vwidth = 2400, vheight = 2000)

## plot all individuals together



cum_min <- function(x) {
  y <- x
  tt <- 0
  for (i in 1:length(x)) {
    if (x[i] == 0 || is.na(x[i])) {
      y[i] <- 0
    } else {
      y[i] <- tt
      tt  <- tt + x[i]
    }
  }
  return(y)
}

cum_min(c(0,0.5,0,0,0.5))

cum_max <- function(x) {
  y <- x
  tt <- 0
  for (i in 1:length(x)) {
    if (x[i] == 0 || is.na(x[i])) {
      y[i] <- 0
    } else {
      tt  <- tt + x[i]
      y[i] <- tt
    }
  }
  return(y)
}

cum_max(c(0,0.5,0.5,0,0))

inds_ER <- rownames(global)[global$AG == 0 & global$D == 0 & global$OB == 0]
inds_AG <- rownames(global)[global$AG > .1 & global$D < 0.1 & global$OB < 0.03]
inds_D <- rownames(global)[global$D > .1 & global$AG == 0]
inds_other <- rownames(global)[!rownames(global) %in% c(inds_ER, inds_D, inds_AG)]
ind_3adm <- c("S-70", "S-134", "S-105", "S-122")
ind_order <- unique(c(inds_ER, inds_other[!inds_other %in% ind_3adm], 
                      inds_AG, inds_D, ind_3adm))

all_blocks <- mcols(all_windows_large) %>% 
  as.data.frame() %>% 
  mutate(chromosome = unlist(mapply(rep, 
                                    all_windows_large@seqnames@values, 
                                    all_windows_large@seqnames@lengths)),
         start = start(all_windows_large),
         end = end(all_windows_large)) %>% 
  relocate(chromosome, start, end, individual)

lai_plot <- all_blocks %>% 
  pivot_longer(cols = OB:ER, names_to = "ancestry", values_to = "proportion") %>% 
  mutate(start = start/1e6,
         end = end/1e6,
         individual = factor(individual, 
                             levels = ind_order)) %>% 
  group_by(chromosome, start, end, individual) %>% 
  mutate(proportion_min = cum_min(proportion),
         proportion_max = cum_max(proportion)) %>% 
  ggplot() + 
  facet_grid(rows = vars(individual), cols = vars(chromosome),
             scale = "free_x", space = "free_x", switch = "x", shrink = F) + 
  geom_rect(aes(xmin = start, xmax = end,
                ymin = proportion_min, ymax = proportion_max,
                fill = ancestry)) +
  scale_fill_manual(values = group_col, na.value = "gray") +
  theme_minimal() + 
  theme(panel.grid = element_blank(), 
        axis.ticks = element_line(),
        axis.text.y = element_blank(),
        strip.background.x = element_rect(fill = "#dadde6", color = "#dadde6"),
        strip.background.y = element_rect(fill = NA, color = NA),
        # strip.background.y = element_rect(fill = "#dadde6", color = "#dadde6"),
        strip.placement = "outside",
        strip.text.y = element_text(angle = 0, h = 0),
        panel.spacing.x = unit(4, "pt"), 
        panel.spacing.y = unit(0.5, "pt"), 
        panel.background = element_rect(fill = "gray60", color = "white")) +
  scale_y_continuous(breaks = 0.5) +
  scale_x_continuous(breaks = seq(0, 80, 10), expand = c(0, 0)) +
  xlab("position (Mb)") + ylab("ancestry assignment")
# ggsave("./plots/all_lai.pdf", width = 20, height = 12)


adm_chr <- all_windows_large[!all_windows_large$individual %in% inds_ER]

avg_local <- do.call(rbind, lapply(seqlevels(adm_chr), function(ch) {
  print(ch)
  chr <- adm_chr[adm_chr@seqnames == ch]
  checkpoints <- sort(unique(c(start(chr), end(chr))))
  chr_score <- do.call(rbind, lapply(checkpoints, function(p) {
    print(paste(ch, ">", p))
    check_gr <- pintersect(chr, 
                           GRanges(ch, 
                                   IRanges(start = p, end = p), 
                                   strand = "+"))
    check_gr <- mcols(check_gr[check_gr$hit])
    check_score <- check_gr %>% 
      as.data.frame() %>% 
      select(-individual, -hit) %>% 
      summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>% 
      mutate(chromosome = ch,
             position = p) %>% 
      relocate(chromosome, position)
    return(check_score)
  }))
  return(chr_score)
}))

avg_plot <- avg_local %>% 
  mutate(position = position/1e6,
         set = "") %>% 
  pivot_longer(cols = OB:ER, names_to = "ancestry", values_to = "proportion") %>% 
  ggplot() + 
  facet_grid(cols = vars(chromosome), rows = vars(set), scale = "free", space = "free", switch = "x") + 
  geom_line(aes(x = position, y = proportion, color = ancestry), show.legend = F) +
  scale_color_manual(values = group_col) +
  theme_minimal() + 
  theme(panel.grid.minor = element_blank(), 
        axis.ticks = element_line(),
        # panel.grid.major.y = element_blank(),
        strip.background.x = element_rect(fill = "#dadde6", color = "#dadde6"),
        strip.background.y = element_rect(fill = NA, color = NA),
        strip.text.y = element_text(angle = 0, h = 0),
        strip.placement = "outside") +
  scale_x_continuous(breaks = seq(0, 80, 10)) +
  xlab("Position (Mb)") + ylab("Average ancestry proportion")
avg_plot

inds_ER <- rownames(global)[global$AG == 0 & global$D == 0 & global$OB == 0]
inds_AG <- rownames(global)[global$AG > .1 & global$D < 0.1 & global$OB < 0.03]
inds_D <- rownames(global)[global$D > .1 & global$AG == 0]
adm1 <- rownames(global)[global$AG > 0 & global$D == 0 & global$OB > 0]
adm2 <- rownames(global)[global$AG > 0 & global$D > 0 & global$OB == 0]
small_adm <- rownames(global)[!rownames(global) %in% c(inds_ER,inds_AG, inds_D, adm1, adm2)]
adm_order <- unique(c(small_adm, inds_AG, inds_D, adm1, adm2))

adm_blocks <- mcols(all_windows_large) %>% 
  as.data.frame() %>% 
  mutate(chromosome = unlist(mapply(rep, 
                                    all_windows_large@seqnames@values, 
                                    all_windows_large@seqnames@lengths)),
         start = start(all_windows_large),
         end = end(all_windows_large)) %>% 
  relocate(chromosome, start, end, individual) %>% 
  filter(individual %in% adm_order)

lai_amd_plot <- adm_blocks %>% 
  pivot_longer(cols = OB:ER, names_to = "ancestry", values_to = "proportion") %>% 
  mutate(start = start/1e6,
         end = end/1e6,
         individual = factor(individual, 
                             levels = adm_order)) %>% 
  group_by(chromosome, start, end, individual) %>% 
  mutate(proportion_min = cum_min(proportion),
         proportion_max = cum_max(proportion)) %>% 
  ggplot() + 
  facet_grid(rows = vars(individual), cols = vars(chromosome),
             scale = "free_x", space = "free_x", switch = "x", shrink = F) + 
  geom_rect(aes(xmin = start, xmax = end,
                ymin = proportion_min, ymax = proportion_max,
                fill = ancestry)) +
  scale_fill_manual(values = group_col, na.value = "gray") +
  theme_minimal() + 
  theme(panel.grid = element_blank(), 
        axis.ticks = element_line(),
        axis.text.y = element_blank(),
        strip.background.x = element_rect(fill = "#dadde6", color = "#dadde6"),
        strip.background.y = element_rect(fill = NA, color = NA),
        # strip.background.y = element_rect(fill = "#dadde6", color = "#dadde6"),
        strip.placement = "outside",
        strip.text.y = element_text(angle = 0, h = 0),
        panel.spacing.x = unit(4, "pt"), 
        panel.spacing.y = unit(0.5, "pt"), 
        panel.background = element_rect(fill = "gray60", color = "white")) +
  scale_y_continuous(breaks = 0.5) +
  scale_x_continuous(breaks = seq(0, 80, 10), expand = c(0, 0)) +
  xlab("Position (Mb)") + ylab("Ancestry assignment")

tiff(filename = file.path(plotdir, "lai_adm.tiff"),
     width = 16, height = 12, units = "in", res = 1200)
# pdf(file = file.path(plotdir, "lai_adm.pdf"),
#     width = 16, height = 12)
plot_list(gglist = list(lai_amd_plot, avg_plot), 
          ncol = 1, heights = c(4,1),
          labels = c("A", "B"))
dev.off()


avg_window <- avg_local %>%
  group_by(chromosome) %>% 
  mutate(window = cut(position, 
                     c(seq(0, max(position), 100000)),
                     dig.lab = 5,
                     include.lowest = T, right = F)) %>% 
  group_by(chromosome, window) %>% 
  summarise(across(OB:ER, mean))


marey_chr_data <- "./Coffea_canephora_Crouzillat2020.txt"
recom <- read.table(marey_chr_data, header = T)
recom$map <- sapply(recom$map, function(i) 
  paste0("Chr ", sprintf("%02d", which(LETTERS == i))), 
  USE.NAMES = F)


recom_window <- recom %>% 
  mutate(rec = gen/(phys/1e6)) %>% 
  rename(position = phys, chromosome = map) %>% 
  group_by(chromosome) %>% 
  mutate(window = cut(position, 
                      c(seq(0, max(position), 100000)),
                      dig.lab = 5,
                      include.lowest = T, right = F)) %>% 
  group_by(chromosome, window) %>% 
  summarise(rec = mean(rec))

merge(avg_window, recom_window) %>% 
  ggplot(aes(x = rec, y = ER)) + 
  geom_point() + 
  geom_smooth(method='lm', se = T)
  # geom_point(aes(x = gen, y = ER))


merge(avg_window, recom_window) %>% 
  ggplot(aes(x = gen, y = AG)) + 
  geom_point() + 
  geom_smooth(method='lm', se = T)

marey_chr_data <- list.files("./snps_data", pattern = "chromosome", full.names = T)
rec_rate <- do.call(rbind, lapply(marey_chr_data, function(f) {
  c <- gsub("(.+chromosome|.txt)", "", f)
  df <- read.table(f, header = T)
  df <- df %>% 
    mutate(chromosome = paste0("Chr ", sprintf("%02d", which(LETTERS == c))))
  return(df)
}
  ))


avg_window <- avg_local %>%
  group_by(chromosome) %>% 
  mutate(window = cut(position, 
                      c(seq(0, max(position), 50000)),
                      dig.lab = 5,
                      include.lowest = T, right = F)) %>% 
  group_by(chromosome, window) %>% 
  summarise(across(OB:ER, mean))


recom_window <- rec_rate %>% 
  mutate(position = phys * 1e6) %>% 
  filter(!is.na(rec.rate)) %>% 
  group_by(chromosome) %>% 
  mutate(window = cut(position, 
                      c(seq(0, max(position), 50000)),
                      dig.lab = 5,
                      include.lowest = T, right = F)) 


merge(avg_window, recom_window) %>% 
  ggplot(aes(x = rec.rate, y = AG)) + 
  geom_point() + 
  geom_smooth(method='lm', se = T)


merge(avg_window, recom_window) %>% 
  ggplot(aes(x = rec.rate, y = D)) + 
  geom_point() + 
  geom_smooth(method='lm', se = T)
format(summary(lm(D ~ rec.rate, merge(avg_window, recom_window)))$r.squared, digits = 3)
