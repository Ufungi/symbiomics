#!/usr/bin/env Rscript
# post_ext_analysis.R -- analyses that sit NEXT TO `symbiomics post` (they never replace it).
#
#   Rscript scripts/post_ext/post_ext_analysis.R <ext_config.json> qc      # before `symbiomics post run`
#   Rscript scripts/post_ext/post_ext_analysis.R <ext_config.json> annot   # eggNOG-mapper output -> post annotation tables
#   Rscript scripts/post_ext/post_ext_analysis.R <ext_config.json> go      # after  `symbiomics post run` (topGO, as in walkthrough.R)
#
# qc    multi-organism read partition from HTSeq special counters (primary genome + companion genomes); signal groups
#       from companion-organism reads per million pairs (never from host expression; read abundance only defines the
#       groups, it is not an infection rate); library QC; the `samples.tsv` that `symbiomics post` reads.
# annot eggNOG-mapper `.emapper.annotations` -> `gene_go` and `gene_labels` tables declared in annotations.tsv.
# go    GO over-representation of the DEGs called by `symbiomics post` (Regulation column): topGO weight01 + Fisher,
#       Up and Down separately, BP/MF/CC, background = tested genes with GO annotation (same as walkthrough.R Step 9).
#
# All project-specific values live in <ext_config.json>. Environment: `symbiomics-post` (envs/post.yml).

suppressPackageStartupMessages(library(jsonlite))
args <- commandArgs(TRUE)
stopifnot(length(args) >= 2, args[2] %in% c("qc", "annot", "go"))
CFG_PATH <- normalizePath(args[1]); STAGE <- args[2]
cfg <- fromJSON(CFG_PATH, simplifyVector = TRUE)
CFG_DIR <- dirname(CFG_PATH)
`%||%` <- function(a, b) if (is.null(a)) b else a
rel <- function(p) normalizePath(file.path(CFG_DIR, p), mustWork = FALSE)
PROJ <- rel(cfg$project_dir)
OUT <- rel(cfg$out_dir); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
set.seed(cfg$seed %||% 20261022)
rd <- function(f) read.delim(f, stringsAsFactors = FALSE, check.names = FALSE)
wr <- function(x, f) write.table(x, f, sep = "\t", quote = FALSE, row.names = FALSE)

# ------------------------------------------------------------------------------------------------ qc ---------
if (STAGE == "qc") {
  read_htseq <- function(f) { d <- read.delim(f, header = FALSE, col.names = c("gene", "n"), stringsAsFactors = FALSE)
    list(assigned = sum(d$n[!startsWith(d$gene, "__")]), not_aligned = d$n[d$gene == "__not_aligned"], total = sum(d$n)) }
  meta <- rd(rel(cfg$metadata_table))
  st <- rd(rel(cfg$pine_stats_table))                       # pipeline counters: sample, assigned, ..., not_aligned
  counter_cols <- intersect(c("assigned", "alignment_not_unique", "ambiguous", "no_feature", "not_aligned", "too_low_aQual"), names(st))
  st$fragments <- rowSums(st[, counter_cols])
  q <- merge(meta, st[, c("sample", "fragments", "assigned", "not_aligned")], by.x = "sample_id", by.y = "sample")
  q$primary_aligned <- q$fragments - q$not_aligned
  q$primary_aligned_pct <- 100 * q$primary_aligned / q$fragments
  q$primary_assigned_pct <- 100 * q$assigned / q$fragments
  comp <- rd(rel(cfg$companion_table))
  for (org in unique(comp$organism)) {
    cc <- comp[comp$organism == org, ]
    v <- t(vapply(seq_len(nrow(cc)), function(i) { h <- read_htseq(rel(cc$file[i])); c(h$assigned, h$total - h$not_aligned, h$total) }, numeric(3)))
    d <- data.frame(sample_id = cc$sample_id); d[[paste0(org, "_assigned")]] <- v[, 1]; d[[paste0(org, "_aligned_pct")]] <- 100 * v[, 2] / v[, 3]
    q <- merge(q, d, by = "sample_id"); q[[paste0(org, "_ppm")]] <- 1e6 * q[[paste0(org, "_assigned")]] / q$fragments
  }
  ppm <- q[[paste0(unique(comp$organism)[1], "_ppm")]]
  ex <- cfg$exposure
  q$exposure_group <- ifelse(q[[ex$label_column]] == ex$control_label, "control", ifelse(ppm >= ex$high_ppm, "high", "treated_low"))
  q$qc_pass <- q$primary_aligned >= cfg$qc$min_aligned_pairs
  q$in_model <- q$qc_pass & q$age_group %in% cfg$qc$model_ages
  q$exclusion_reason <- ifelse(q$in_model, "", ifelse(!q$qc_pass, sprintf("primary-genome aligned read pairs < %g", cfg$qc$min_aligned_pairs), "age group not modelled"))
  q <- q[match(meta$sample_id, q$sample_id), ]
  wr(q, file.path(OUT, "library_qc.tsv"))
  ps <- cfg$post_samples
  hdr <- c("sample_id", "input_id", "source_sample_id", "stats_source_path", "biological_unit_id", "library_id", "technical_merge_id", "study_id",
           "sample_sets", "include_analysis", "exclusion_reason", "subject_id", "treatment", "timepoint", "time_h", "batch", "tissue", "stats_sample_id",
           "count_unit_override", "technical_partition_id", "technical_disjoint", "age_group", "biological_unit_status", "metadata_status")
  s <- as.data.frame(matrix("", nrow(q), length(hdr), dimnames = list(NULL, hdr)), stringsAsFactors = FALSE)
  s$sample_id <- q$sample_id; s$input_id <- ps$input_id; s$source_sample_id <- q$sample_id; s$biological_unit_id <- q$sample_id; s$library_id <- q$sample_id
  s$study_id <- ps$study_id; s$sample_sets <- ifelse(q$in_model, ps$sample_set, ""); s$include_analysis <- ifelse(q$in_model, "true", "false")
  s$exclusion_reason <- q$exclusion_reason; s$treatment <- q$exposure_group; s$tissue <- "root"; s$age_group <- q$age_group
  s$biological_unit_status <- "verified"; s$metadata_status <- "supplied"
  wr(s, file.path(PROJ, "samples.tsv"))
  cat(sprintf("[qc] %d libraries; QC-pass %d; in model %d (%s)\n", nrow(q), sum(q$qc_pass), sum(q$in_model),
              paste(names(table(q$exposure_group[q$in_model])), table(q$exposure_group[q$in_model]), collapse = ", ")))
  quit(save = "no")
}

# ---------------------------------------------------------------------------------------------- annot -------
if (STAGE == "annot") {
  a <- cfg$annotation
  em_f <- rel(a$emapper_annotations)
  hd <- strsplit(sub("^#", "", grep("^#query", readLines(em_f), value = TRUE)[1]), "\t")[[1]]
  em <- read.delim(em_f, comment.char = "#", header = FALSE, quote = "", stringsAsFactors = FALSE)
  names(em) <- hd[seq_len(ncol(em))]
  go <- em[em$GOs != "-" & em$GOs != "", c("query", "GOs")]
  gg <- unique(data.frame(gene_id = rep(go$query, lengths(strsplit(go$GOs, ","))), go_id = unlist(strsplit(go$GOs, ",")), stringsAsFactors = FALSE))
  gg <- gg[grepl("^GO:", gg$go_id), ]
  lab <- ifelse(em$Preferred_name != "-" & em$Preferred_name != "", em$Preferred_name, ifelse(em$Description != "-", substr(em$Description, 1, 40), ""))
  gl <- data.frame(gene_id = em$query, label = lab, stringsAsFactors = FALSE); gl <- gl[nzchar(gl$label), ]
  wr(gg, file.path(PROJ, a$gene_go_out)); wr(gl, file.path(PROJ, a$gene_labels_out))
  cat(sprintf("[annot] %d genes with eggNOG hits; %d with GO (%d gene-GO pairs); %d labels\n", nrow(em), length(unique(gg$gene_id)), nrow(gg), nrow(gl)))
  quit(save = "no")
}

# ------------------------------------------------------------------------------------------------- go --------
suppressPackageStartupMessages({ library(topGO); library(dplyr); library(stringr) })
latest <- fromJSON(file.path(PROJ, "output/latest_success.json")); RUN <- latest$path
g2go_tab <- rd(file.path(PROJ, cfg$annotation$gene_go_out)); g2go_full <- lapply(split(g2go_tab$go_id, g2go_tab$gene_id), unique)
run_topgo <- function(interesting, ontology, g2go) {
  data <- new("topGOdata", ontology = ontology, allGenes = interesting, annot = annFUN.gene2GO, gene2GO = g2go)
  res <- runTest(data, algorithm = "weight01", statistic = "fisher")
  tab <- GenTable(data, Fishers = res, useLevels = TRUE, numChar = 1000, topNodes = max(1, min(500, length(score(res)))))
  tab$Ontology <- ontology; tab
}
et <- data.frame()
for (an in cfg$go$analyses) {
  degf <- list.files(file.path(RUN, "analyses", an, "contrasts"), pattern = "^DEG.tsv$", recursive = TRUE, full.names = TRUE)
  for (f in degf) {
    deg <- rd(f); bg <- intersect(names(g2go_full), deg$gene_id); g2go <- g2go_full[bg]
    for (reg in c("Up", "Down")) {
      sel <- intersect(deg$gene_id[deg$Regulation == paste0(reg, "-regulated")], bg)
      cat(sprintf("[go] %s/%s: %d annotated DEGs (%d background genes)\n", an, reg, length(sel), length(bg)))
      if (length(sel) < 2) next
      interesting <- factor(as.integer(bg %in% sel)); names(interesting) <- bg
      tabs <- bind_rows(lapply(c("BP", "MF", "CC"), function(o) tryCatch(run_topgo(interesting, o, g2go), error = function(e) NULL))) %>%
        mutate(Fishers = suppressWarnings(as.numeric(str_replace(Fishers, "^<\\s*", ""))), Gene_ratio = Significant / Annotated,
               Analysis = an, Regulation = reg) %>% arrange(Fishers)
      et <- bind_rows(et, tabs)
    }
  }
}
wr(et, file.path(OUT, "TopGO_result.tsv"))
cat("[go] wrote", file.path(OUT, "TopGO_result.tsv"), "\n")
