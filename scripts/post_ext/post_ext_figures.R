#!/usr/bin/env Rscript
# post_ext_figures.R -- figures for a `symbiomics post` project, in the style of post_analysis/RNA-seq_re/walkthrough.R
#   Rscript scripts/post_ext/post_ext_figures.R <ext_config.json>
# Reads only tables written by `symbiomics post run` and post_ext_analysis.R (qc, go); no statistics are recomputed
# here except the walkthrough's PCA (FactoMineR) and pvclust, which are descriptive.
#   Fig 1 read composition | Fig 2 sample structure (PCA/fviz, pvclust, pheatmap) | Fig 3 volcanoes
#   Fig 4 GO over-representation (topGO weight01 bubble plot) | Fig 5 candidate genes (box plots)
# Each figure has its own legend (<= 100 words) BELOW the panels. Output: <out_dir>/figures/Fig{1..5}.{pdf,png}
# Environment: symbiomics-post (envs/post.yml).

suppressPackageStartupMessages({
  library(jsonlite); library(ggplot2); library(patchwork); library(ggrepel); library(dplyr); library(tidyr)
  library(FactoMineR); library(factoextra); library(pheatmap); library(pvclust); library(pals); library(stringr)
})
args <- commandArgs(TRUE); CFG_PATH <- normalizePath(args[1]); cfg <- fromJSON(CFG_PATH)
CFG_DIR <- dirname(CFG_PATH); rel <- function(p) normalizePath(file.path(CFG_DIR, p), mustWork = FALSE)
PROJ <- rel(cfg$project_dir); EXT <- rel(cfg$out_dir); FIG <- file.path(EXT, "figures"); dir.create(FIG, showWarnings = FALSE)
RUN <- fromJSON(file.path(PROJ, "output/latest_success.json"))$path
rd <- function(f) read.delim(f, stringsAsFactors = FALSE, check.names = FALSE)
FONT <- "Liberation Sans"; stopifnot(grepl("Liberation", systemfonts::match_fonts(FONT)$path))
set.seed(123)

# ---- style: Kelly palette and red/blue/grey regulation colours as in walkthrough.R ------------------------------
KELLY <- as.vector(kelly(13))
col_grp <- c(control = KELLY[3], high = KELLY[4], treated_low = KELLY[5])
lab_grp <- c(control = "Control", treated_low = "Treated, low matsutake reads", high = "High matsutake reads")
col_org <- c(pine = "#5B8C5A", matsutake = "#C97B2D")
REG_COL <- c("Up-regulated" = "red", "Down-regulated" = "blue", "Unchanged" = "grey50")
theme_pub <- function(bs = 8) theme_minimal(bs, base_family = FONT) +
  theme(plot.title = element_text(face = "bold", size = bs + 1), plot.tag = element_text(face = "bold", size = bs + 4),
        legend.position = "bottom", legend.title = element_text(face = "bold"), panel.grid.minor = element_blank())
theme_set(theme_pub())
short <- function(x) sub("_root$", "", x)
trunc_lab <- function(x, n = 30) ifelse(nchar(x) > n, paste0(substr(x, 1, n - 1), "…"), x)
nwords <- function(x) length(strsplit(trimws(x), "\\s+")[[1]])
save_fig <- function(p, name, w, h, number, title, legend) {            # legend BELOW the panels, <= 100 words
  stopifnot(nwords(legend) + nwords(title) + 3 <= 100)
  cap <- paste(strwrap(sprintf("Figure %d | %s. %s", number, title, legend), width = 165), collapse = "\n")
  p <- p + plot_annotation(tag_levels = "A", caption = cap, theme = theme(plot.caption = element_text(size = 6.2, family = FONT, hjust = 0, lineheight = 1.1, margin = margin(t = 8))))
  ggsave(file.path(FIG, paste0(name, ".pdf")), p, width = w, height = h, units = "in", device = cairo_pdf)
  ggsave(file.path(FIG, paste0(name, ".png")), p, width = w, height = h, units = "in", dpi = 600, device = ragg::agg_png, bg = "white")
  cat("[fig] saved", name, "\n")
}

q <- rd(file.path(EXT, "library_qc.tsv")); q$exposure_group <- factor(q$exposure_group, levels = c("control", "high", "treated_low"))
AN1 <- cfg$main_analysis; AN2 <- cfg$secondary_analysis
analyses <- rd(file.path(PROJ, "analyses.tsv")); fdr_cut <- as.numeric(analyses$fdr_cut[analyses$analysis_id == AN1]); lfc_cut <- as.numeric(analyses$lfc_cut[analyses$analysis_id == AN1])
ex <- cfg$exposure; qc <- cfg$qc; fg <- cfg$figures
lab_tab <- NULL; ann_f <- file.path(PROJ, "annotations.tsv")
if (file.exists(ann_f)) { an <- rd(ann_f); an <- an[an$kind == "gene_labels", ]; if (nrow(an)) { f <- file.path(PROJ, an$path[1]); if (file.exists(f)) lab_tab <- rd(f) } }
gene_label <- function(ids, with_id = FALSE) { out <- ids
  if (!is.null(lab_tab)) { m <- match(ids, lab_tab$gene_id); has <- !is.na(m) & nzchar(lab_tab$label[m]); out[has] <- if (with_id) paste0(lab_tab$label[m][has], " (", ids[has], ")") else lab_tab$label[m][has] }
  out }
deg_of <- function(an) { f <- list.files(file.path(RUN, "analyses", an, "contrasts"), pattern = "^DEG.tsv$", recursive = TRUE, full.names = TRUE)[1]; rd(f) }

# =============================== Figure 1: read composition ==================================================
d1 <- q %>% arrange(factor(age_group, levels = c("5y", "15y", "30y", "1y", "cont")), exposure_group, desc(matsutake_ppm))
d1$age_group <- factor(d1$age_group, levels = c("5y", "15y", "30y", "1y", "cont"))
d1$lib <- factor(short(d1$sample_id), levels = rev(short(d1$sample_id)))
comp <- d1 %>% select(lib, age_group, pine = primary_aligned_pct, matsutake = matsutake_aligned_pct) %>%
  pivot_longer(c(pine, matsutake), names_to = "genome", values_to = "pct") %>% mutate(genome = factor(genome, levels = c("pine", "matsutake")))
d1$txt_col <- ifelse(d1$qc_pass, "black", "#B2182B")
comp <- d1 %>% transmute(lib, age_group, pine = primary_aligned_pct, matsutake = matsutake_aligned_pct, other = pmax(0, 100 - primary_aligned_pct - matsutake_aligned_pct)) %>%
  pivot_longer(c(other, matsutake, pine), names_to = "genome", values_to = "pct") %>% mutate(genome = factor(genome, levels = c("other", "matsutake", "pine")))
pA1 <- ggplot(comp, aes(pct, lib)) + geom_col(aes(fill = genome), width = 0.75) +
  geom_point(data = d1, aes(x = -6, y = lib, colour = exposure_group), shape = 15, size = 2.6, inherit.aes = FALSE) +
  geom_text(data = d1, aes(x = 104, y = lib, label = sprintf("%.1f M", primary_aligned / 1e6), colour = I(txt_col)), size = 1.9, family = FONT, hjust = 0, inherit.aes = FALSE) +
  facet_grid(age_group ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_manual(values = c(other = "grey85", col_org), breaks = c("pine", "matsutake", "other"),
                    labels = c(pine = "P. densiflora genome", matsutake = "T. matsutake genome", other = "Unaligned / other"), name = "Read pairs") +
  scale_colour_manual(values = col_grp, labels = lab_grp, name = "Matsutake-read group") +
  scale_x_continuous(limits = c(-9, 120), breaks = c(0, 25, 50, 75, 100), expand = c(0, 0)) +
  labs(x = "Read pairs (%)", y = NULL, title = "Read composition (HISAT2)") +
  theme(axis.text.y = element_text(size = 5.5), strip.placement = "outside", strip.text.y.left = element_text(angle = 0, face = "bold", size = 7), panel.grid = element_blank())
pB1 <- ggplot(d1, aes(matsutake_ppm + 1, lib)) + geom_vline(xintercept = ex$high_ppm, linetype = "dashed", linewidth = 0.3) +
  geom_segment(aes(x = 1, xend = matsutake_ppm + 1, yend = lib, colour = exposure_group), linewidth = 0.4, show.legend = FALSE) +
  geom_point(aes(colour = exposure_group), size = 1.7, show.legend = FALSE) + facet_grid(age_group ~ ., scales = "free_y", space = "free_y") +
  scale_x_log10(breaks = 10^(0:5), labels = c("1", "10", "100", "1,000", "10,000", "100,000")) + scale_colour_manual(values = col_grp) +
  labs(x = "Matsutake reads per million pairs", y = NULL, title = "Matsutake read abundance") +
  theme(axis.text.y = element_blank(), strip.text = element_blank(), panel.grid = element_blank())
save_fig(pA1 + pB1 + plot_layout(widths = c(1.5, 1), guides = "collect") & theme(legend.position = "bottom", legend.box = "vertical", legend.key.size = unit(0.3, "cm"), legend.text = element_text(size = 6.5), legend.title = element_text(size = 7, face = "bold"), legend.spacing.y = unit(0.05, "cm")),
  "Fig1_read_composition", 7.2, 8.4, 1, "Read composition of the 26 root libraries",
  sprintf("Reads were aligned with HISAT2 to the pine HA genome and to the matsutake genome and counted with HTSeq (pairs). Stacked bars: percentage of read pairs aligned to each genome (unaligned/other = remainder); numbers: pine-aligned pairs (red: below %.0f M, not modelled). Dashed line: %s matsutake reads per million pairs, used only to define the high group; it is not an infection rate. 15-/30-y trees were not modelled.",
          qc$min_aligned_pairs / 1e6, format(ex$high_ppm, big.mark = ",")))

# =============================== Figure 2: sample structure (walkthrough Steps 4, 4.5, 7) ====================
lc <- as.matrix(read.delim(file.path(RUN, "analyses", AN1, "normalized/logCPM.tsv"), row.names = 1, check.names = FALSE))
ny <- as.data.frame(t(lc))                                                     # sample x gene (log2-CPM, TMM, filterByExpr genes)
grp <- factor(q$exposure_group[match(rownames(ny), q$sample_id)], levels = c("control", "high")); np <- nlevels(grp)
res_pca <- PCA(ny, graph = FALSE)
pA2 <- fviz(res_pca, element = "ind", label = "none", habillage = grp, legend.title = "Group", addEllipses = FALSE, geom = "point",
            pointshape = 16, pointsize = 4, invisible = "quali", axes = 1:2) +
  scale_color_manual(values = as.vector(kelly(13)[3:(2 + np)])) + scale_fill_manual(values = as.vector(kelly(13)[3:(2 + np)])) +
  geom_text_repel(data = data.frame(res_pca$ind$coord[, 1:2], id = rownames(ny)), aes(Dim.1, Dim.2, label = short(id)), inherit.aes = FALSE,
                  size = 2.3, family = FONT, seed = 1, box.padding = 0.4, segment.size = 0.2) +
  labs(title = "PCA (FactoMineR, all tested genes)") + theme_pub() + theme(legend.position = "bottom")
n_top <- ceiling(fg$top_variable_frac * ncol(ny)); vg <- apply(t(ny), 1, var); tv <- names(sort(vg, decreasing = TRUE))[seq_len(n_top)]
vdf <- t(ny)[tv, ]; colnames(vdf) <- short(colnames(vdf))
pv <- pvclust(vdf, method.hclust = "average", method.dist = "correlation", nboot = fg$pvclust_nboot, quiet = TRUE)
dendro_gg <- function(pv) {                                                    # pvclust dendrogram with AU (red) p-values, as plot(pv)
  hc <- pv$hclust; n <- length(hc$order); pos <- setNames(match(seq_len(n), hc$order), seq_len(n)); ht <- setNames(rep(0, n), seq_len(n))
  seg <- list(); lab <- list(); x <- c(pos, rep(NA, n - 1)); names(x) <- c(seq_len(n), paste0("m", seq_len(n - 1))); h <- c(rep(0, n), hc$height); names(h) <- names(x)
  for (i in seq_len(n - 1)) { a <- hc$merge[i, ]; key <- function(v) if (v < 0) as.character(-v) else paste0("m", v)
    k1 <- key(a[1]); k2 <- key(a[2]); x[[paste0("m", i)]] <- mean(c(x[[k1]], x[[k2]])); hi <- hc$height[i]
    seg[[length(seg) + 1]] <- data.frame(x = c(x[[k1]], x[[k1]], x[[k2]]), xend = c(x[[k1]], x[[k2]], x[[k2]]), y = c(h[[k1]], hi, h[[k2]]), yend = c(hi, hi, hi))
    lab[[length(lab) + 1]] <- data.frame(x = x[[paste0("m", i)]], y = hi, au = round(100 * pv$edges$au[i])) }
  list(seg = do.call(rbind, seg), lab = do.call(rbind, lab), leaves = data.frame(x = seq_len(n), label = hc$labels[hc$order]))
}
dg <- dendro_gg(pv)
pB2 <- ggplot() + geom_segment(data = dg$seg, aes(x = x, xend = xend, y = y, yend = yend), linewidth = 0.4) +
  geom_text(data = dg$lab, aes(x = x, y = y, label = au), colour = "red", size = 2.2, vjust = -0.4, family = FONT) +
  geom_text(data = dg$leaves, aes(x = x, y = 0, label = label), angle = 90, hjust = 1, vjust = 0.5, size = 2.4, family = FONT, nudge_y = -0.005) +
  scale_y_continuous(expand = expansion(mult = c(0.25, 0.1))) + labs(x = NULL, y = "Height (1 - correlation)", title = sprintf("pvclust (AU %%, %d bootstraps)", fg$pvclust_nboot)) +
  theme(axis.text.x = element_blank(), panel.grid = element_blank())
brk <- c(seq(min(vdf), 0, length.out = 51), seq(0, max(vdf), length.out = 50)[-1]); colr <- colorRampPalette(c("magenta", "black", "yellow"))(length(brk) - 1)
anno <- data.frame(Group = grp, row.names = short(rownames(ny))); anno_col <- list(Group = setNames(as.vector(kelly(13)[3:(2 + np)]), levels(grp)))
ph <- pheatmap(vdf, color = colr, breaks = brk, annotation_col = anno, annotation_colors = anno_col, show_rownames = FALSE, fontsize_col = 8, fontsize = 7,
               main = sprintf("Top %.0f%% most variable genes (n = %d)", 100 * fg$top_variable_frac, n_top), silent = TRUE)
pC2 <- wrap_elements(full = ph$gtable)
save_fig((pA2 | pB2) / pC2 + plot_layout(heights = c(1, 1.25)), "Fig2_sample_structure", 7.2, 8.4, 2, "Sample structure of the 5-y-old root libraries",
  sprintf("log2-CPM from symbiomics post (filterByExpr, TMM). A: PCA of all %s tested genes (FactoMineR, scaled; factoextra fviz). B: pvclust (average linkage, correlation distance, %d bootstraps; red = approximately unbiased p-value) on the top %.0f%% most variable genes (%d). C: heat map of those genes (magenta-black-yellow, log2-CPM). 6 high-read and 3 control libraries.",
          format(ncol(ny), big.mark = ","), fg$pvclust_nboot, 100 * fg$top_variable_frac, n_top))

# =============================== Figure 3: volcano plots (walkthrough Step 6.5) ==============================
n_label <- fg$volcano_labels
vol_plot <- function(df, title, xlim_v, ylim_v) {
  df <- df %>% mutate(FDRp = pmax(FDR, .Machine$double.xmin)) %>% filter(is.finite(logFC), !is.na(FDR))
  df_lab <- df %>% mutate(score = -log10(FDRp) * abs(logFC)) %>% filter(Regulation != "Unchanged") %>% slice_max(score, n = n_label, with_ties = FALSE) %>%
    mutate(label = gene_label(gene_id))
  ggplot(df, aes(logFC, -log10(FDRp))) + geom_point(aes(color = Regulation), size = 1) + ggtitle(title) +
    xlab(expression("log"[2] * "FC")) + ylab(expression("-log"[10] * "FDR")) + scale_color_manual(values = REG_COL) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") + geom_hline(yintercept = -log10(c(0.05, 0.01, 0.001)), linetype = "dashed", color = "black") +
    geom_label_repel(data = df_lab, aes(label = label), size = 2.1, max.overlaps = 30, box.padding = 0.25, point.padding = 0.1, min.segment.length = 0, label.size = 0.15, label.padding = 0.12, seed = 3) +
    coord_cartesian(xlim = xlim_v, ylim = ylim_v)
}
de1 <- deg_of(AN1); de2 <- deg_of(AN2)
xl <- c(-1, 1) * ceiling(max(abs(c(de1$logFC, de2$logFC)), na.rm = TRUE)); yl <- c(0, ceiling(max(-log10(pmax(c(de1$FDR, de2$FDR), .Machine$double.xmin)), na.rm = TRUE)))
cnt <- function(d) sprintf("Up %d, Down %d", sum(d$Regulation == "Up-regulated"), sum(d$Regulation == "Down-regulated"))
p3a <- vol_plot(de1, sprintf("DESeq2 (symbiomics post): %s", cnt(de1)), xl, yl); p3b <- vol_plot(de2, sprintf("edgeR QL (symbiomics post): %s", cnt(de2)), xl, yl)
save_fig(p3a + p3b + plot_layout(guides = "collect") & theme(legend.position = "bottom"), "Fig3_volcano", 7.2, 4.4, 3,
  "Differential expression, high matsutake-read libraries vs control",
  sprintf("Both panels come from symbiomics post on HISAT2 + HTSeq counts (filterByExpr, TMM, same design matrix and contrast), drawn on the same axes. A: effect_test=deseq2 (Wald test, ashr-shrunk log2FC, BH FDR). B: effect_test=qlf (edgeR quasi-likelihood). Colours: FDR < %g and |log2FC| > %g. Dashed lines: |log2FC| = %g; FDR 0.05, 0.01, 0.001. Labels: top %d genes by -log10 FDR x |log2FC|.",
          fdr_cut, lfc_cut, lfc_cut, n_label))

# =============================== Figure 4: GO over-representation (walkthrough Step 9) =======================
go_f <- file.path(EXT, "TopGO_result.tsv")
if (file.exists(go_f) && file.info(go_f)$size > 10) {
  et <- rd(go_f); et <- et %>% filter(Analysis == AN1, Fishers < cfg$go$p_cut, Significant > 1)
  go_thresh <- -log10(c(0.05, 0.01, 0.001)); plots <- list()
  for (reg in c("Up", "Down")) {
    gd <- et %>% filter(Regulation == reg) %>% arrange(Fishers) %>% slice_head(n = cfg$go$top_n)
    if (!nrow(gd)) next
    gd <- gd %>% mutate(Term = str_trunc(Term, 50), lab = paste0(GO.ID, " ", Term))
    plots[[reg]] <- ggplot(gd, aes(x = reorder(lab, -log10(Fishers)), y = -log10(Fishers))) +
      geom_hline(yintercept = go_thresh, linetype = "dashed", color = "grey60") +
      geom_point(aes(size = Significant, fill = Gene_ratio), shape = 21) + geom_text(aes(label = Significant), hjust = 2.5, size = 2.6) +
      scale_size(range = c(3, 10)) + scale_fill_continuous(low = "yellow", high = "red", limits = c(0, 1)) +
      scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) + xlab(NULL) + ylab("Enrichment -log10(P) (weight01)") +
      ggtitle(sprintf("%s-regulated genes (%s)", reg, paste(sort(unique(gd$Ontology)), collapse = "/"))) + coord_flip() + theme_pub(8)
  }
  if (length(plots)) {
    save_fig(wrap_plots(plots, ncol = 1, heights = vapply(plots, function(p) nrow(p$data) + 3, 1)) + plot_layout(guides = "collect") & theme(legend.position = "bottom"),
      "Fig4_GO_enrichment", 7.2, 2.4 + 0.22 * sum(vapply(plots, function(p) nrow(p$data), 1)), 4, "GO over-representation of differentially expressed genes",
      sprintf("DEGs from symbiomics post (DESeq2; FDR < %g, |log2FC| > %g) were tested separately for up- and down-regulation with topGO (weight01, Fisher exact; BP, MF, CC). Background: tested genes with eggNOG-mapper GO terms. Shown: terms with P < %g (top %d). Bubble size and number: DEGs in the term; fill: DEG fraction of the term. Dashed lines: P = 0.05, 0.01, 0.001.",
              fdr_cut, lfc_cut, cfg$go$p_cut, cfg$go$top_n))
  } else cat("[fig] no GO term passes p_cut: Fig4 skipped\n")
} else cat("[fig] TopGO_result.tsv missing: run the `go` stage; Fig4 skipped\n")

# =============================== Figure 5: candidate genes (walkthrough Step 8.11 box plots) =================
cand <- de1 %>% filter(Regulation != "Unchanged") %>% slice_min(FDR, n = 12, with_ties = FALSE) %>% mutate(label = gene_label(gene_id, with_id = TRUE))
dB <- as.data.frame(t(lc[cand$gene_id, , drop = FALSE])) %>% mutate(sample_id = rownames(.)) %>% pivot_longer(-sample_id, names_to = "gene_id", values_to = "logcpm") %>%
  left_join(select(q, sample_id, exposure_group), by = "sample_id") %>% left_join(select(cand, gene_id, label, logFC, FDR), by = "gene_id") %>%
  mutate(panel = factor(sprintf("%s\nlog2FC %.1f, FDR %.2g", trunc_lab(label, 36), logFC, FDR), levels = unique(sprintf("%s\nlog2FC %.1f, FDR %.2g", trunc_lab(cand$label, 36), cand$logFC, cand$FDR))),
         exposure_group = factor(as.character(exposure_group), levels = c("control", "high")))
p5 <- ggplot(dB, aes(exposure_group, logcpm, colour = exposure_group)) + geom_boxplot(outlier.shape = NA, linewidth = 0.3, width = 0.5, fill = NA) +
  geom_point(position = position_jitter(width = 0.12, height = 0, seed = 1), size = 1.6) + facet_wrap(~panel, scales = "free_y", ncol = 4) +
  scale_colour_manual(values = col_grp, guide = "none") + scale_x_discrete(labels = c(control = "Control", high = "High")) +
  labs(x = NULL, y = "log2 CPM (TMM)") + theme(strip.text = element_text(size = 5.8), panel.spacing = unit(2.5, "mm"))
save_fig(p5, "Fig5_candidate_genes", 7.2, 6.2, 5, "Expression of the strongest differentially expressed genes",
  sprintf("The 12 genes with the smallest FDR among the DESeq2 DEGs (FDR < %g, |log2FC| > %g; symbiomics post). Values: TMM-normalised log2 CPM from symbiomics post; points are libraries (6 high-read, 3 control), boxes median and quartiles. Panel titles give the eggNOG-mapper name (gene ID), shrunk log2FC and FDR.", fdr_cut, lfc_cut))
cat("[done]\n")
