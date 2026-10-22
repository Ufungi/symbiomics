#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3L) stop('Usage: quant_qc_dashboard.R sample_metrics.tsv unassigned_long.tsv output_prefix')
metrics_path <- args[[1L]]
unassigned_path <- args[[2L]]
out_prefix <- args[[3L]]
set.seed(20261008L)
metrics <- read.delim(metrics_path, check.names=FALSE, stringsAsFactors=FALSE,
                      quote='', comment.char='', na.strings=c('', 'NA'))
unassigned <- read.delim(unassigned_path, check.names=FALSE, stringsAsFactors=FALSE,
                         quote='', comment.char='', na.strings=c('', 'NA'))
required <- c('sample_id', 'gene_count_sum', 'genes_detected', 'quantified_gene_count',
              'assignment_counter_rate_pct')
missing <- setdiff(required, names(metrics))
if (length(missing)) stop('sample_metrics.tsv is missing: ', paste(missing, collapse=', '))
if (!('sample_id' %in% names(unassigned))) stop('unassigned_long.tsv is missing sample_id')

number <- function(x) suppressWarnings(as.numeric(as.character(x)))
metrics$gene_count_sum <- number(metrics$gene_count_sum)
metrics$genes_detected <- number(metrics$genes_detected)
metrics$quantified_gene_count <- number(metrics$quantified_gene_count)
metrics$assignment_counter_rate_pct <- number(metrics$assignment_counter_rate_pct)
metrics$supplied_gene_detection_pct <- ifelse(
  is.finite(metrics$genes_detected) & is.finite(metrics$quantified_gene_count) & metrics$quantified_gene_count > 0,
  100 * metrics$genes_detected / metrics$quantified_gene_count, NA_real_)
metrics$group[is.na(metrics$group) | !nzchar(metrics$group)] <- 'unknown'
if (!('study_id' %in% names(metrics))) metrics$study_id <- 'all samples'
metrics$study_id[is.na(metrics$study_id) | !nzchar(metrics$study_id)] <- 'unknown study'
metrics$study_id <- factor(metrics$study_id, levels=unique(metrics$study_id))
metrics$group <- factor(metrics$group, levels=unique(metrics$group))
study_levels <- levels(metrics$study_id)
group_levels <- levels(metrics$group)
study_cols <- setNames(grDevices::hcl.colors(max(3L, length(study_levels)), 'Dark 3')[seq_along(study_levels)], study_levels)
group_pch <- setNames(rep(c(16, 17, 15, 18, 3, 4), length.out=length(group_levels)), group_levels)
study_col <- unname(study_cols[as.character(metrics$study_id)])
group_symbol <- unname(group_pch[as.character(metrics$group)])
group_x <- as.integer(metrics$group)

plot_group_points <- function(y, title, ylab, log_y=FALSE) {
  valid <- is.finite(y)
  if (!any(valid)) {
    plot.new(); title(main=title); text(.5, .5, 'No values available')
    return(invisible(NULL))
  }
  if (log_y) {
    positives <- y[valid & y > 0]
    if (length(positives)) {
      lower <- 10^floor(log10(min(positives)))
      upper <- 10^ceiling(log10(max(positives)))
      if (upper <= lower) upper <- lower * 10
      y_plot <- rep(lower/2, length(y))
      y_plot[valid & y > 0] <- y[valid & y > 0]
      plot(group_x[valid] + stats::runif(sum(valid), -.12, .12), y_plot[valid],
           log='y', ylim=c(lower/2, upper*1.7), xaxt='n', xlab='Input group', ylab=ylab,
           pch=group_symbol[valid], col=study_col[valid], main=title, cex=.85)
      axis(1, at=seq_along(group_levels), labels=group_levels)
      tick <- 10^(seq(floor(log10(lower)), ceiling(log10(upper))))
      axis(2, at=c(lower/2, tick), labels=c('0', format(tick, scientific=TRUE, trim=TRUE)), las=1)
      abline(h=tick, col='grey90', lty=1)
      points(group_x[valid & y==0] + stats::runif(sum(valid & y==0), -.12, .12),
             rep(lower/2, sum(valid & y==0)), pch=4, col=study_col[valid & y==0], lwd=1.2)
      mtext('true zero shown at bottom lane', side=3, line=.15, cex=.63, col='grey35')
      return(invisible(NULL))
    }
  }
  plot(group_x[valid] + stats::runif(sum(valid), -.12, .12), y[valid],
       xaxt='n', xlab='Input group', ylab=ylab, pch=group_symbol[valid],
       col=study_col[valid], main=title, cex=.85)
  axis(1, at=seq_along(group_levels), labels=group_levels)
  abline(h=pretty(y[valid]), col='grey90', lty=1)
}

plot_unassigned <- function() {
  if (!nrow(unassigned) || !('category' %in% names(unassigned)) || !('count' %in% names(unassigned))) {
    plot.new(); title(main='HTSeq unassigned counters'); text(.5, .5, 'No category counts available')
    return(invisible(NULL))
  }
  joined <- merge(unassigned, metrics[, c('sample_id', 'group', 'study_id')], by='sample_id', all.x=TRUE, sort=FALSE)
  joined$count <- number(joined$count)
  joined <- joined[is.finite(joined$count) & !is.na(joined$category), , drop=FALSE]
  if (!nrow(joined)) {
    plot.new(); title(main='HTSeq unassigned counters'); text(.5, .5, 'No category counts available')
    return(invisible(NULL))
  }
  cats <- unique(joined$category)
  cat_x <- match(joined$category, cats)
  g_idx <- match(joined$group, group_levels)
  symbols <- group_pch[as.character(joined$group)]
  colors <- study_cols[as.character(joined$study_id)]
  plot(cat_x + stats::runif(length(cat_x), -.16, .16), log10(joined$count + 1),
       xaxt='n', xlab='HTSeq counter category', ylab='log10(count + 1)',
       pch=unname(symbols), col=unname(colors), main='Unassigned categories by sample', cex=.65)
  axis(1, at=seq_along(cats), labels=cats, las=2, cex.axis=.72)
  abline(h=pretty(log10(joined$count + 1)), col='grey90', lty=1)
}

render_dashboard <- function() {
  old <- par(no.readonly=TRUE)
  on.exit(par(old))
  layout(matrix(c(1,2,3,4,5,5), nrow=3, byrow=TRUE), heights=c(1,1,.28))
  par(oma=c(.2,.3,2.0,.3), mar=c(5.2,4.4,3.4,.8))
  plot_group_points(metrics$assignment_counter_rate_pct,
                    'Native counter assignment (%)', 'Assigned / all HTSeq counters (%)', log_y=TRUE)
  plot_group_points(metrics$supplied_gene_detection_pct,
                    'Detected of supplied genes (%)', 'Detected genes / supplied matrix rows (%)')
  plot_group_points(metrics$gene_count_sum,
                    'Gene count sum', 'Sum of unfiltered gene counts (log10 scale)', log_y=TRUE)
  plot_unassigned()
  par(mar=c(0,0,0,0))
  plot.new()
  mtext('Quantification QC - four views', side=3, outer=TRUE, line=.7, cex=1.2, font=2)
  legend('center', legend=c(paste('Study:', study_levels), paste('Group:', group_levels)),
         col=c(unname(study_cols), rep('grey25', length(group_levels))),
         pch=c(rep(16, length(study_levels)), unname(group_pch)),
         ncol=min(5L, length(study_levels) + length(group_levels)),
         bty='n', cex=.78, title='Color = study; symbol = input group')
}

pdf(paste0(out_prefix, '.pdf'), width=13, height=9, onefile=TRUE, useDingbats=FALSE)
render_dashboard()
dev.off()
graphics.off()
svg(paste0(out_prefix, '.svg'), width=13, height=9, pointsize=10, onefile=TRUE)
render_dashboard()
dev.off()
cat('Created ', paste0(out_prefix, '.pdf'), ' and ', paste0(out_prefix, '.svg'), '\n', sep='')