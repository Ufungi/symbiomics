#!/usr/bin/env Rscript
# General bulk RNA-seq engine. The orchestrator supplies resolved, selected inputs.
options(stringsAsFactors = FALSE, warn = 1)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop('Usage: Rscript engine.R job.json')
if (!requireNamespace('jsonlite', quietly = TRUE)) stop('The symbiomics-post environment requires jsonlite')
job <- jsonlite::fromJSON(args[[1]], simplifyVector = FALSE)
scalar <- function(x, default = '') if (is.null(x) || !length(x) || is.na(x[[1]]) || identical(x[[1]], '')) default else as.character(x[[1]])
a <- job$analysis
out <- scalar(job$out_dir)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
status <- list(status = 'running', analysis_id = scalar(a$analysis_id), stages = list(), warnings = character())
write_json <- function(x, p) jsonlite::write_json(x, p, auto_unbox = TRUE, pretty = TRUE, na = 'null', digits = 16)
write_status <- function() write_json(status, file.path(out, 'status.json'))
record <- function(stage, state, reason = '') { status$stages[[stage]] <<- list(status = state, reason = reason); write_status() }
fail <- function(...) stop(paste0(...), call. = FALSE)
assert <- function(ok, ...) if (!isTRUE(ok)) fail(...)
need <- function(pkg) assert(requireNamespace(pkg, quietly = TRUE), 'Required R package is missing: ', pkg)
path <- function(...) { p <- file.path(out, ...); dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE); p }
read_tsv <- function(p) read.delim(p, check.names = FALSE, colClasses = 'character', quote = '"', comment.char = '', na.strings = '', fileEncoding = 'UTF-8')
write_tsv <- function(x, p) { dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE); write.table(x, p, sep = '\t', quote = TRUE, qmethod = 'double', row.names = FALSE, na = 'NA', fileEncoding = 'UTF-8') }
write_matrix <- function(x, p, id = 'gene_id') { z <- data.frame(rownames(x), x, check.names = FALSE); names(z)[1] <- id; write_tsv(z, p) }
parts <- function(x, sep = '|') { x <- scalar(x); if (!nzchar(x)) character() else strsplit(x, sep, fixed = TRUE)[[1]] }
kv <- function(x) { z <- parts(x, ';'); ans <- character(); for (s in z) { at <- regexpr('=', s, fixed = TRUE)[1]; assert(at > 1L, 'Expected key=value: ', s); k <- substr(s, 1, at - 1L); assert(!k %in% names(ans), 'Duplicate key: ', k); value <- substr(s, at + 1L, nchar(s)); assert(!grepl('%(?![0-9A-Fa-f]{2})',value,perl=TRUE),'Invalid percent escape: ',value); value <- utils::URLdecode(value); assert(!grepl('[\t\r\n]',value),'Encoded control characters are not allowed'); ans[k] <- value }; ans }
num <- function(x, default) { z <- suppressWarnings(as.numeric(scalar(x, as.character(default)))); assert(length(z) == 1 && is.finite(z), 'Invalid numeric setting'); z }
safe_id <- function(x) assert(grepl('^[A-Za-z0-9][A-Za-z0-9_.-]*$', x) && !grepl('..', x, fixed = TRUE), 'Unsafe ID: ', x)
encode <- function(x) utils::URLencode(enc2utf8(x), reserved = TRUE)
rows_df <- function(rows, cols) { if (!length(rows)) return(as.data.frame(setNames(rep(list(character()), length(cols)), cols))); do.call(rbind, lapply(rows, function(r) as.data.frame(setNames(lapply(cols, function(k) scalar(r[[k]])), cols)))) }
with_pdf <- function(p, code, width = 9, height = 7) { dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE); if (capabilities('cairo')) grDevices::cairo_pdf(p, width = width, height = height) else grDevices::pdf(p, width = width, height = height); on.exit(grDevices::dev.off()); force(code) }

build_design <- function(samples) {
  factors <- parts(a$factor_cols); numbers <- parts(a$numeric_cols)
  nf <- parts(a$nuisance_factor_cols); nn <- parts(a$nuisance_numeric_cols)
  subject <- scalar(a$subject_col); subj <- if (nzchar(subject)) subject else character()
  typ <- scalar(a$design_type, 'group')
  assert(typ %in% c('group','additive','factorial','paired','trend'), 'Unsupported design_type: ', typ)
  assert(length(c(factors,numbers)) > 0L, 'At least one explanatory variable is required')
  assert(typ != 'paired' || length(subj), 'paired design requires subject_col')
  assert(typ != 'trend' || length(numbers), 'trend design requires numeric_cols')
  cols <- c(factors,numbers,nf,nn,subj)
  assert(!anyDuplicated(cols), 'A variable is declared more than once')
  assert(all(grepl('^[A-Za-z][A-Za-z0-9_]*$', cols)), 'Model column names must be safe identifiers')
  assert(all(cols %in% names(samples)), 'Missing model columns: ', paste(setdiff(cols,names(samples)), collapse=', '))
  refs <- kv(a$reference_levels); units <- kv(a$numeric_units)
  dat <- samples; maps <- list(); component <- character()
  for (v in cols) assert(!any(is.na(dat[[v]]) | !nzchar(dat[[v]])), 'Missing metadata: ', v)
  for (v in c(factors,nf,subj)) {
    assert(v %in% names(refs), 'Missing reference level for ', v)
    levels0 <- c(refs[[v]], sort(setdiff(unique(dat[[v]]), refs[[v]])))
    assert(refs[[v]] %in% dat[[v]], 'Reference level absent: ', v, '=', refs[[v]])
    assert(length(levels0) >= 2L, 'Factor has fewer than two levels: ', v)
    codes <- sprintf('L%04d', seq_along(levels0)); maps[[v]] <- setNames(codes, levels0)
    dat[[v]] <- factor(unname(maps[[v]][dat[[v]]]), levels = codes)
    contrasts(dat[[v]]) <- stats::contr.treatment(codes, base = 1)
    for (i in seq.int(2L,length(codes))) component[paste0(v,codes[[i]])] <- paste0('factor:', v, ':', encode(levels0[[i]]))
  }
  for (v in c(numbers,nn)) {
    dat[[v]] <- suppressWarnings(as.numeric(dat[[v]])); assert(all(is.finite(dat[[v]])), 'Non-numeric/non-finite variable: ', v)
    assert(length(unique(dat[[v]])) > 1L, 'Numeric variable has no variation: ', v)
    assert(v %in% names(units) && nzchar(units[[v]]), 'Missing numeric unit for ', v)
    component[v] <- paste0('numeric:',v)
  }
  ints <- parts(a$interaction_terms)
  assert(typ != 'factorial' || length(ints), 'factorial requires interaction_terms')
  terms <- cols
  for (term in ints) {
    iv <- strsplit(term, ':', fixed = TRUE)[[1]]
    assert(length(iv) >= 2L && !anyDuplicated(iv) && all(iv %in% c(factors,numbers)), 'Interaction must use distinct explanatory variables: ', term)
    for (k in seq_len(length(iv))) terms <- c(terms, apply(combn(iv,k),2,paste,collapse=':'))
  }
  f <- stats::reformulate(unique(terms)); d <- model.matrix(f, dat); rownames(d) <- samples$sample_id
  ids <- vapply(colnames(d), function(v) {
    if (v == '(Intercept)') return('intercept')
    cc <- strsplit(v, ':', fixed = TRUE)[[1]]; mapped <- unname(component[cc]); assert(!anyNA(mapped), 'Unmapped model coefficient: ',v)
    if (length(mapped) == 1L) mapped else paste0('interaction:',paste(mapped,collapse='*'))
  }, '')
  assert(!anyDuplicated(ids), 'Coefficient ID collision')
  coef <- data.frame(coefficient_id = ids, internal_column = colnames(d), formula_term = c('(Intercept)',attr(terms(f),'term.labels'))[attr(d,'assign')+1L], stringsAsFactors=FALSE)
  write_matrix(d,path('design','design_matrix.tsv'),'sample_id'); write_tsv(coef,path('design','coefficients.tsv'))
  mapdf <- do.call(rbind,lapply(names(maps),function(v) data.frame(variable=v,label=names(maps[[v]]),internal_level=unname(maps[[v]]),reference=names(maps[[v]])==refs[[v]])))
  if (!is.null(mapdf)) write_tsv(mapdf,path('design','factor_levels.tsv'))
  rank <- qr(d)$rank; rdf <- nrow(d)-rank
  report <- list(formula = paste(deparse(f),collapse=''), rank=rank, columns=ncol(d), residual_df=rdf, references=as.list(refs), numeric_units=as.list(units), sample_count=nrow(samples), biological_units=length(unique(samples$biological_unit_id)))
  if (rank < ncol(d)) report$aliased_coefficients <- ids[qr(d)$pivot[seq.int(rank+1L,ncol(d))]]
  write_json(report,path('design','design_report.json'))
  assert(rank == ncol(d), 'Rank-deficient design; aliased coefficients: ',paste(report$aliased_coefficients,collapse=', '))
  assert(rdf > 0L, 'No positive residual degrees of freedom')
  assert('biological_unit_id' %in% names(samples) && !anyNA(samples$biological_unit_id), 'biological_unit_id is required for inference')
  repeated <- anyDuplicated(samples$biological_unit_id) > 0L
  assert(!repeated || length(subj), 'Repeated biological units require an explicit subject model')
  if (length(subj)) {
    pairs <- unique(samples[c('biological_unit_id',subject)])
    assert(!anyDuplicated(pairs$biological_unit_id) && !anyDuplicated(pairs[[subject]]), 'subject and biological_unit_id must have a one-to-one mapping')
    assert(length(unique(samples[[subject]])) >= 2L,'At least two independent subjects are required')
    write_tsv(samples[c('sample_id','biological_unit_id',subject,factors,numbers)],path('design','subject_observations.tsv'))
  }
  if (typ == 'group' && !length(subj)) {
    cells <- interaction(samples[factors], drop = TRUE)
    sizes <- tapply(samples$biological_unit_id, cells, function(z) length(unique(z)))
    assert(all(sizes >= 2L), 'Group inference requires at least two biological units per observed cell')
  }
  if (!length(factors)) assert(length(unique(samples$biological_unit_id)) >= 3L, 'Numeric-only models require at least three independent biological units')
  list(d=d,dat=dat,formula=f,ids=ids,factors=factors,numbers=numbers,nf=nf,nn=nn,subject=subj,maps=maps,refs=refs,units=units,report=report)
}

make_contrasts <- function(model,samples) {
  cr <- rows_df(job$contrasts,c('analysis_id','contrast_id','basis','term','weight','test_id'))
  if (nrow(cr)) cr <- cr[cr$analysis_id==scalar(a$analysis_id),,drop=FALSE]
  vectors <- list(); tests <- list(); details <- list(); failures <- list()
  for (id in unique(cr$contrast_id)) {
    safe_id(id); r <- cr[cr$contrast_id==id,,drop=FALSE]
    tryCatch({
      assert(length(unique(r$basis))==1L && length(unique(r$test_id))==1L,'Inconsistent basis/test_id')
      w <- suppressWarnings(as.numeric(r$weight)); assert(all(is.finite(w)), 'Non-finite contrast weight')
      v <- rep(0,ncol(model$d)); cells <- list(); supports <- list()
      if (r$basis[[1]]=='coefficient') {
        assert(all(r$term %in% model$ids),'Unknown coefficient ID: ',paste(setdiff(r$term,model$ids),collapse=', '))
        for(i in seq_len(nrow(r))) v[match(r$term[[i]],model$ids)] <- v[match(r$term[[i]],model$ids)]+w[[i]]
      } else {
        assert(r$basis[[1]]=='cell','Unsupported contrast basis')
        assert(abs(sum(w)) <= 1e-10*max(1,sum(abs(w))),'Cell comparison weights must sum to zero')
        for(i in seq_len(nrow(r))) {
          cell <- kv(r$term[[i]]); expected <- c(model$factors,model$numbers)
          assert(setequal(names(cell),expected),'Cell must specify every explanatory variable exactly: ',paste(expected,collapse=', '))
          row <- model$dat[1,,drop=FALSE]; support <- rep(TRUE,nrow(samples))
          for(vv in c(model$factors,model$nf,model$subject)) {
            value <- if(vv %in% model$factors) cell[[vv]] else model$refs[[vv]]
            assert(value %in% names(model$maps[[vv]]),'Unknown factor level: ',vv,'=',value)
            row[[vv]] <- factor(unname(model$maps[[vv]][value]),levels=levels(model$dat[[vv]]))
            contrasts(row[[vv]]) <- contrasts(model$dat[[vv]])
            if(vv %in% model$factors) support <- support & samples[[vv]]==value
          }
          for(vv in c(model$numbers,model$nn)) {
            value <- if(vv %in% model$numbers) suppressWarnings(as.numeric(cell[[vv]])) else mean(model$dat[[vv]])
            assert(is.finite(value),'Invalid cell numeric value: ',vv)
            row[[vv]] <- value
            if(vv %in% model$numbers) support <- support & abs(model$dat[[vv]]-value)<1e-10
          }
          assert(any(support),'Cell has no observed samples: ',r$term[[i]])
          drow <- model.matrix(model$formula,row); assert(identical(colnames(drow),colnames(model$d)),'Cell model column mismatch')
          v <- v+w[[i]]*drow[1,]; cells[[i]] <- as.list(row[c(model$nf,model$nn,model$subject)])
          supports[[i]] <- support
          if(scalar(a$design_type)=='group' && !length(model$subject) && w[[i]] != 0) assert(length(unique(samples$biological_unit_id[support]))>=2L,'Cell requires at least two biological units: ',r$term[[i]])
        }
      }
      if(length(model$subject)) {
        if(length(supports)) complete <- Reduce(intersect,lapply(supports[w!=0],function(s) unique(samples[[model$subject]][s]))) else {
          support_cols <- c(model$factors, model$numbers)
          informative <- vapply(split(seq_len(nrow(samples)),samples[[model$subject]]), function(ii) nrow(unique(samples[ii,support_cols,drop=FALSE])) >= 2L, TRUE)
          complete <- names(informative)[informative]
        }
        assert(length(complete)>=2L,'Paired contrast requires at least two independently informative subjects')
      }
      assert(any(abs(v)>1e-12),'Contrast vector is zero')
      names(v) <- model$ids; vectors[[id]] <- v; tests[[id]] <- r$test_id[[1]]
      details[[id]] <- list(basis=r$basis[[1]],terms=as.list(r$term),weights=as.list(w),fixed_nuisance=cells,test_id=r$test_id[[1]])
    },error=function(e) { failures[[id]] <<- conditionMessage(e) })
  }
  if(length(vectors)) write_matrix(do.call(cbind,vectors),path('design','contrast_vectors.tsv'),'coefficient_id')
  write_json(list(contrasts=details,failed=failures),path('design','contrast_report.json'))
  list(vectors=vectors,tests=tests,details=details,failures=failures,declared=cr)
}

load_counts <- function(samples) {
  if(length(job$salmon)) {
    need('tximport'); records <- rows_df(job$salmon,c('sample_id','quant_path','gene_map_path'))
    assert(!anyDuplicated(records$sample_id) && setequal(records$sample_id,samples$sample_id),'Salmon mapping must cover selected samples exactly')
    records <- records[match(samples$sample_id,records$sample_id),,drop=FALSE]
    maps <- lapply(unique(records$gene_map_path),read_tsv)
    maps <- lapply(maps,function(x) { assert(ncol(x)>=2L,'tx2gene needs transcript and gene columns'); x <- x[,1:2]; names(x)<-c('tx_id','gene_id'); x })
    txmap <- unique(do.call(rbind,maps)); assert(!anyDuplicated(txmap$tx_id),'Transcript maps to multiple genes')
    txi <- tximport::tximport(setNames(records$quant_path,records$sample_id),type='salmon',tx2gene=txmap,countsFromAbundance='no')
    saveRDS(txi,path('normalized','tximport.rds')); return(list(counts=txi$counts,length=txi$length,source='salmon_tximport_length_offset'))
  }
  z <- read_tsv(scalar(job$counts_path)); assert(ncol(z)>1L && names(z)[1]=='gene_id','Counts TSV must start with gene_id')
  assert(!anyNA(z$gene_id) && !anyDuplicated(z$gene_id),'Missing or duplicate gene_id')
  assert(all(samples$sample_id %in% names(z)),'Selected sample is missing from count matrix')
  m <- as.matrix(z[,samples$sample_id,drop=FALSE]); storage.mode(m)<-'double'; rownames(m)<-z$gene_id
  assert(all(is.finite(m)) && all(m>=0),'Counts must be finite and non-negative')
  list(counts=m,length=NULL,source='raw_gene_counts')
}

normalize_counts <- function(raw,model,exploratory) {
  y <- edgeR::DGEList(raw$counts); lib <- colSums(y$counts)
  assert(all(lib>0),'Zero-count libraries cannot be normalized; leave them in global QC and explicitly exclude them from this analysis')
  if(exploratory) {
    threshold <- 1; required <- min(2L,ncol(y)); keep <- rowSums(edgeR::cpm(y)>=threshold)>=required
    rule <- list(method='unsupervised_cpm',min_cpm=threshold,min_samples=required)
  } else { keep <- edgeR::filterByExpr(y,design=model$d); rule <- list(method='filterByExpr',design='design/design_matrix.tsv',min_count=10,min_total_count=15,large_n=10,min_prop=0.7) }
  write_tsv(data.frame(gene_id=rownames(y),retained=keep),path('normalized','filter.tsv')); assert(any(keep),'No genes remain after filtering')
  y <- y[keep,,keep.lib.sizes=FALSE]
  post_lib <- colSums(y$counts)
  write_tsv(data.frame(sample_id=colnames(y),raw_count_sum=lib,filtered_count_sum=post_lib,normalizable=post_lib>0),path('normalized','sample_filter_status.tsv'))
  assert(all(post_lib>0), 'Samples have zero retained counts after filtering: ',paste(colnames(y)[post_lib==0],collapse=', '),'. Specify an explicit analysis subset; global raw QC remains available.')
  if(!is.null(raw$length)) {
    len <- raw$length[rownames(y),colnames(y),drop=FALSE]
    assert(all(is.finite(len)) && all(len>0),'Invalid effective lengths for retained genes')
    len <- len/exp(rowMeans(log(len)))
    scaled <- y$counts/len; factors <- edgeR::calcNormFactors(scaled,method='TMM'); eff <- colSums(scaled)*factors
    y <- edgeR::scaleOffset(y,log(sweep(len,2,eff,'*')))
    y$samples$norm.factors <- factors
    write_matrix(y$offset,path('normalized','log_offsets.tsv'))
  } else y <- edgeR::calcNormFactors(y,method='TMM')
  saveRDS(y,path('normalized','dge.rds'))
  write_tsv(data.frame(sample_id=colnames(y),y$samples,check.names=FALSE),path('normalized','normalization_factors.tsv'))
  logcpm <- edgeR::cpm(y,log=TRUE,prior.count=2)
  write_matrix(logcpm,path('normalized','logCPM.tsv'))
  write_json(list(source=raw$source,filter=rule,method='TMM',logCPM_prior_count=2,genes_input=nrow(raw$counts),genes_retained=nrow(y),sample_count=ncol(y),length_offset=!is.null(raw$length)),path('normalized','settings.json'))
  list(y=y,logcpm=logcpm)
}

explore <- function(logcpm,samples,folder='exploratory',label='Unadjusted') {
  dir.create(path(folder),recursive=TRUE,showWarnings=FALSE)
  if(ncol(logcpm)<2L) {
    write_matrix(as.matrix(dist(t(logcpm))),path(folder,'sample_distance.tsv'),'sample_id')
    write_json(list(status='skipped',reason='At least two selected samples are required for PCA/MDS/heatmap'),path(folder,'plot_settings.json'))
    return(invisible(NULL))
  }
  vars <- apply(logcpm,1,var); chosen <- head(order(vars,decreasing=TRUE),500L); mat <- logcpm[chosen,,drop=FALSE]
  groupcol <- head(intersect(parts(a$factor_cols),names(samples)),1)
  labels <- if(length(groupcol)) samples[[groupcol]] else rep('samples',nrow(samples))
  labels[is.na(labels)]<-'missing'; levels0<-unique(labels); colors<-setNames(grDevices::hcl.colors(length(levels0),'Dark 3'),levels0)[labels]
  if(ncol(mat)>=2L && any(vars>0)) {
    pc <- prcomp(t(mat),center=TRUE,scale.=FALSE); coords <- pc$x[,seq_len(min(2L,ncol(pc$x))),drop=FALSE]
    if(ncol(coords)==1L) coords<-cbind(coords,PC2=0)
    explained <- pc$sdev^2/sum(pc$sdev^2)
    write_tsv(data.frame(sample_id=rownames(coords),coords,group=labels),path(folder,'pca_scores.tsv'))
    write_tsv(data.frame(component=seq_along(explained),variance_fraction=explained),path(folder,'pca_variance.tsv'))
    write_matrix(pc$rotation,path(folder,'pca_loadings.tsv'))
    with_pdf(path(folder,'PCA.pdf'), { plot(coords[,1],coords[,2],col=colors,pch=19,xlab=sprintf('PC1 (%.1f%%)',100*explained[1]),ylab=sprintf('PC2 (%.1f%%)',100*ifelse(length(explained)>1,explained[2],0)),main=paste(label,'PCA')); text(coords[,1],coords[,2],labels=rownames(coords),pos=3,cex=max(.35,min(.8,20/nrow(coords)))); legend('topright',legend=levels0,col=unique(colors),pch=19,cex=.7) })
    mds <- limma::plotMDS(logcpm,top=min(500L,nrow(logcpm)),dim.plot=c(1L,min(2L,ncol(logcpm)-1L)),plot=FALSE)
    if(ncol(logcpm)==2L) mds$y <- rep(0,2L)
    write_tsv(data.frame(sample_id=samples$sample_id,MDS1=mds$x,MDS2=mds$y,group=labels),path(folder,'mds_scores.tsv'))
    with_pdf(path(folder,'MDS.pdf'),{ plot(mds$x,mds$y,col=colors,pch=19,xlab='MDS1',ylab='MDS2',main=paste(label,'MDS')); text(mds$x,mds$y,labels=samples$sample_id,pos=3,cex=max(.35,min(.8,20/nrow(samples)))) })
  }
  distm <- as.matrix(dist(t(logcpm))); write_matrix(distm,path(folder,'sample_distance.tsv'),'sample_id')
  hm <- logcpm[head(order(vars,decreasing=TRUE),min(50L,nrow(logcpm))),,drop=FALSE]
  hm <- hm[apply(hm,1,sd)>0,,drop=FALSE]; if(nrow(hm)) hm<-t(scale(t(hm)))
  write_matrix(hm,path(folder,'heatmap_matrix.tsv'))
  if(nrow(hm)>=2L && ncol(hm)>=2L) with_pdf(path(folder,'Heatmap_top_variable.pdf'),stats::heatmap(hm,scale='none',col=colorRampPalette(c('#2166ac','white','#b2182b'))(101),ColSideColors=colors,margins=c(8,8),cexRow=.5,cexCol=.6,main=paste(label,'row-scaled logCPM')),width=11,height=10)
  write_json(list(pca_max_variable_genes=500,heatmap_max_variable_genes=50,pca_center=TRUE,pca_scale=FALSE,heatmap_scale='gene_row_zscore',mds_method='limma_plotMDS_leading_logFC',distance='euclidean_logCPM',group_column=groupcol),path(folder,'plot_settings.json'))
}

adjusted_view <- function(raw,model,samples) {
  if(scalar(a$adjusted_view,'none')=='none') return(record('adjusted_view','skipped','not_requested'))
  tryCatch({
    assert(scalar(a$adjusted_view)=='combat_seq','Unsupported adjusted_view')
    batchcol<-scalar(a$adjusted_batch_col); assert(batchcol %in% names(samples),'Missing adjusted_batch_col')
    batch<-factor(samples[[batchcol]]); assert(!anyNA(batch) && nlevels(batch)>=2L && all(table(batch)>=2L),'ComBat requires at least two batches and two samples in each')
    assert(is.null(raw$length) && all(abs(raw$counts-round(raw$counts))<1e-8),'ComBat-Seq requires integer gene counts; imported length-offset data unsupported')
    need('sva')
    protected<-setdiff(c(parts(a$factor_cols),parts(a$numeric_cols),parts(a$nuisance_factor_cols),parts(a$nuisance_numeric_cols),scalar(a$subject_col)),c('',batchcol))
    if(length(parts(a$interaction_terms))) assert(!batchcol %in% unlist(strsplit(parts(a$interaction_terms),':',fixed=TRUE)),'Cannot remove an interacting batch')
    terms0<-c(protected,parts(a$interaction_terms)); pmat<-if(length(terms0)) model.matrix(reformulate(unique(terms0)),model$dat) else matrix(1,nrow(samples),1)
    covar<-pmat[,-1,drop=FALSE]; joint<-cbind(model.matrix(~batch),covar)
    assert(qr(joint)$rank==ncol(joint),'Protected biology and batch are confounded; adjusted view skipped')
    adjusted<-sva::ComBat_seq(raw$counts,batch=batch,covar_mod=if(ncol(covar)) covar else NULL,full_mod=TRUE)
    z<-edgeR::DGEList(adjusted); k<-edgeR::filterByExpr(z,design=if(is.null(model))matrix(1,nrow(samples),1) else model$d); assert(any(k),'No adjusted genes survive filtering'); z<-edgeR::calcNormFactors(z[k,,keep.lib.sizes=FALSE])
    l<-edgeR::cpm(z,log=TRUE,prior.count=2); write_matrix(l,path('exploratory_adjusted','logCPM_visualization_only.tsv')); explore(l,samples,'exploratory_adjusted','Batch-adjusted visualization only'); record('adjusted_view','ok','visualization_only')
  },error=function(e) record('adjusted_view','skipped',conditionMessage(e)))
}

load_annotations <- function() {
  ans<-list(); for(ann in job$annotations) {
    if(scalar(ann$assay_id)!=scalar(a$assay_id)) next
    id<-scalar(ann$annotation_id); safe_id(id)
    z<-read_tsv(scalar(ann$path)); kind<-scalar(ann$kind)
    req<-switch(kind,gene_go=c('gene_id','go_id'),gene_pathway=c('gene_id','pathway_id','pathway_name'),gene_labels=c('gene_id','label'),fail('Unsupported annotation kind: ',kind))
    assert(all(req %in% names(z)),'Annotation columns missing: ',id)
    assert(!anyNA(z[req]),'Annotation values missing: ',id)
    if(nzchar(scalar(ann$id_map_path))) {
      mp<-read_tsv(scalar(ann$id_map_path)); assert(all(c('source_id','gene_id') %in% names(mp)),'ID mapping requires source_id,gene_id')
      assert(!anyDuplicated(mp$source_id),'Ambiguous annotation ID mapping')
      z$gene_id<-mp$gene_id[match(z$gene_id,mp$source_id)]
      write_json(list(rows=nrow(z),unmapped_rows=sum(is.na(z$gene_id))),path('annotations',paste0(id,'_mapping.json')))
      z<-z[!is.na(z$gene_id),,drop=FALSE]
    }
    universe<-NULL
    if(nzchar(scalar(ann$universe_path))) { u<-read_tsv(scalar(ann$universe_path)); assert('gene_id' %in% names(u),'universe_path requires gene_id'); universe<-unique(u$gene_id) }
    ans[[id]]<-list(kind=kind,data=unique(z),universe=universe)
  }; ans
}

run_enrichment <- function(tab,cv,y,design,anns,dest) {
  requested<-parts(a$enrichment); if(!length(requested)||identical(requested,'none')) return(invisible(NULL))
  assert(all(requested %in% c('go','pathway')),'Unsupported enrichment request')
  for(kind in requested) assert(any(vapply(anns,function(z) z$kind==if(kind=='go')'gene_go' else 'gene_pathway',TRUE)),'Requested ',kind,' enrichment has no matching annotation')
  for(id in names(anns)) {
    ann<-anns[[id]]; go<-ann$kind=='gene_go'
    if(!(go && 'go' %in% requested) && !(ann$kind=='gene_pathway' && 'pathway' %in% requested)) next
    z<-ann$data; column<-if(go)'go_id' else 'pathway_id'; tested<-rownames(y)
    background<-intersect(tested,unique(z$gene_id)); if(!is.null(ann$universe)) background<-intersect(background,ann$universe)
    assert(length(background)>0L,'No tested annotated genes overlap: ',id)
    zz<-z[z$gene_id %in% background,,drop=FALSE]; sets<-lapply(split(zz$gene_id,zz[[column]]),unique)
    write_tsv(data.frame(gene_id=background),file.path(dest,paste0(id,'_background.tsv')))
    write_tsv(zz,file.path(dest,paste0(id,'_gene_sets.tsv')))
    write_json(list(tested_genes=length(tested),annotated_tested=length(background),unmapped_tested=sum(!tested %in% z$gene_id)),file.path(dest,paste0(id,'_mapping.json')))
    if(go) {
      assert(scalar(a$go_method,'fisher_ora')=='fisher_ora','Unsupported go_method')
      selected<-tab$gene_id[tab$Regulation!='Unchanged']
      ora<-lapply(c('All','Up-regulated','Down-regulated'),function(direction) {
        chosen<-intersect(if(direction=='All')selected else tab$gene_id[tab$Regulation==direction],background)
        result<-do.call(rbind,lapply(names(sets),function(term) { genes<-sets[[term]]; overlap<-intersect(chosen,genes); n<-length(background); k<-length(genes); s<-length(chosen); x<-length(overlap); p<-if(!s)1 else phyper(x-1,k,n-k,s,lower.tail=FALSE); data.frame(term_id=term,direction=direction,overlap=x,set_size=k,selected=s,background=n,PValue=p,genes=paste(overlap,collapse='|')) }))
        result$FDR<-p.adjust(result$PValue,'BH'); result
      }); write_tsv(do.call(rbind,ora),file.path(dest,paste0(id,'_go_ora.tsv')))
    } else {
      methods<-parts(a$pathway_methods); assert(length(methods)&&all(methods %in% c('camera','fgsea')),'Declare pathway_methods camera and/or fgsea')
      if('camera' %in% methods) {
        yy <- y[background,,keep.lib.sizes=TRUE]
        idx<-lapply(sets,function(s) which(rownames(yy) %in% s))
        result<-edgeR::camera.DGEList(yy,index=idx,design=design,contrast=cv)
        write_tsv(data.frame(pathway_id=rownames(result),result),file.path(dest,paste0(id,'_camera.tsv')))
      }
      if('fgsea' %in% methods) {
        need('fgsea'); ranks<-setNames(sign(tab$logFC)*sqrt(pmax(tab$F,0)),tab$gene_id)
        assert(all(is.finite(ranks)),'Non-finite fgsea ranking'); ranks<-sort(ranks[names(ranks) %in% background],decreasing=TRUE)
        write_tsv(data.frame(gene_id=names(ranks),rank_statistic=ranks),file.path(dest,paste0(id,'_fgsea_ranks.tsv')))
        result<-as.data.frame(fgsea::fgseaMultilevel(pathways=sets,stats=ranks,minSize=1,maxSize=length(ranks),eps=0,nproc=1))
        if('leadingEdge' %in% names(result)) result$leadingEdge<-vapply(result$leadingEdge,paste,collapse='|',FUN.VALUE='')
        write_tsv(result,file.path(dest,paste0(id,'_fgsea.tsv')))
      }
    }
  }
}

# DESeq2 test (effect_test=deseq2): same filtered genes, design matrix and contrast vectors as the edgeR path.
# Size factors are the TMM factors (normalization stays TMM); Wald test; LFC shrunk with ashr; BH over all tested genes
# (independentFiltering=FALSE). Columns are mapped onto the edgeR table layout so volcano/enrichment code is shared:
# logFC = shrunk log2FC, log2FC_mle = unshrunk, F = Wald stat^2 (ranking sign(logFC)*sqrt(F) = |stat|*sign), logCPM = edgeR AveLogCPM.
build_deseq2 <- function(raw,norm,model,y) {
  need('DESeq2'); need('ashr'); genes<-rownames(norm$y); samples<-colnames(norm$y)
  d<-model$d[samples,,drop=FALSE]; colnames(d)<-paste0('b',seq_len(ncol(d))-1L)
  dds<-DESeq2::DESeqDataSetFromMatrix(round(raw$counts[genes,samples,drop=FALSE]),data.frame(row.names=samples,sample_id=samples),design=d)
  eff<-norm$y$samples$lib.size*norm$y$samples$norm.factors; DESeq2::sizeFactors(dds)<-eff/exp(mean(log(eff)))
  DESeq2::DESeq(dds,quiet=TRUE)
}
deseq2_table <- function(dds,cv,y) {
  r<-DESeq2::results(dds,contrast=as.numeric(cv),independentFiltering=FALSE)
  sh<-suppressMessages(DESeq2::lfcShrink(dds,contrast=as.numeric(cv),type='ashr',quiet=TRUE))
  data.frame(logFC=sh$log2FoldChange,log2FC_mle=r$log2FoldChange,logCPM=y$AveLogCPM,F=r$stat^2,PValue=r$pvalue,FDR=r$padj,row.names=rownames(r),check.names=FALSE)
}

main <- function() {
  need('edgeR'); need('limma'); safe_id(scalar(a$analysis_id))
  samples<-read_tsv(scalar(job$samples_path)); assert('sample_id' %in% names(samples) && nrow(samples)>0,'No selected samples'); assert(!anyNA(samples$sample_id)&&!anyDuplicated(samples$sample_id),'Invalid sample IDs')
  write_tsv(samples,path('design','analysis_samples.tsv'))
  set.seed(num(a$seed,123)); mode<-scalar(job$mode,'run'); assert(mode %in% c('validate','run'),'Unsupported engine mode')
  exploratory<-scalar(a$analysis_mode,'inferential')=='exploratory'; assert(scalar(a$analysis_mode,'inferential') %in% c('exploratory','inferential'),'Unsupported analysis_mode')
  fdr<-num(a$fdr_cut,.05); lfc<-num(a$lfc_cut,1); assert(fdr>0&&fdr<=1&&lfc>=0,'Invalid FDR/LFC threshold')
  assert(scalar(a$normalization,'TMM')=='TMM','Only TMM normalization is supported')
  effect<-scalar(a$effect_test,'qlf'); assert(effect %in% c('qlf','treat','deseq2'),'Unsupported effect_test')
  model<-NULL; cs<-NULL
  if(exploratory) {
    write_json(list(analysis_mode='exploratory',inference=FALSE,reason='Experimental independence/design not asserted; unsupervised analysis only'),path('design','design_report.json'))
    record('design','skipped','exploratory_without_inference')
  } else {
    model<-build_design(samples); cs<-make_contrasts(model,samples); record('design','ok')
    for(id in names(cs$failures)) record(paste0('contrast:',id),'failed',cs$failures[[id]])
  }
  if(mode=='validate') {
    status$status<<-if(!exploratory&&length(cs$failures))'partial_failure' else 'ok'; write_status(); return(invisible(NULL))
  }
  raw<-load_counts(samples); raw$counts<-raw$counts[,samples$sample_id,drop=FALSE]
  norm<-normalize_counts(raw,model,exploratory); record('normalization','ok')
  explore(norm$logcpm,samples); record('exploratory','ok')
  if(exploratory) { record('deg','skipped','exploratory_without_inference'); record('enrichment','skipped','exploratory_without_inference'); if(scalar(a$adjusted_view,'none')!='none') record('adjusted_view','skipped','A verified protection model is required'); status$status<<-'ok'; write_status(); return(invisible(NULL)) }
  adjusted_view(raw,model,samples)
  assert(length(cs$vectors)>0L,'No valid contrasts were declared')
  y<-edgeR::estimateDisp(norm$y,model$d,robust=TRUE)
  # edgeR 4.4.0 glmQLFit(legacy=FALSE) SIGBUSes when the preset NB scalar dispersion >= 4 (later edgeR caps preset dispersions; emulate, no-op below 4)
  disp<-NULL
  if(!is.null(y$trended.dispersion)){ if(is.null(y$AveLogCPM)) y$AveLogCPM<-edgeR::aveLogCPM(y); ql<-mean(y$trended.dispersion[order(y$AveLogCPM,decreasing=TRUE)[1:ceiling(0.1*nrow(y))]]); if(is.finite(ql)&&ql>=4){ warning(sprintf('preset NB dispersion %.3f capped to 3.9 for glmQLFit (edgeR 4.4.0 segfault guard)',ql)); disp<-3.9 } }
  fit<-edgeR::glmQLFit(y,model$d,dispersion=disp,robust=TRUE)
  saveRDS(list(y=y,fit=fit,design=model$d),path('normalized','model.rds'))
  dds<-NULL; if(effect=='deseq2') dds<-build_deseq2(raw,norm,model,y)
  anns<-list(); annotation_error<-if(nzchar(scalar(job$annotation_error)))scalar(job$annotation_error) else NULL
  if(!identical(scalar(a$enrichment,'none'),'none') && is.null(annotation_error)) tryCatch({anns<-load_annotations()},error=function(e) annotation_error<<-conditionMessage(e))
  for(id in names(cs$vectors)) {
    tryCatch({
      cv<-cs$vectors[[id]]; dest<-path('contrasts',id); dir.create(dest,recursive=TRUE,showWarnings=FALSE)
      tab<-if(effect=='deseq2') deseq2_table(dds,cv,y) else { result<-if(effect=='treat')edgeR::glmTreat(fit,contrast=cv,lfc=lfc) else edgeR::glmQLFTest(fit,contrast=cv); edgeR::topTags(result,n=Inf,sort.by='none')$table }
      tab<-data.frame(analysis_id=scalar(a$analysis_id),contrast_id=id,gene_id=rownames(tab),tab,check.names=FALSE)
      tab$Regulation<-ifelse(tab$FDR<fdr & tab$logFC>lfc,'Up-regulated',ifelse(tab$FDR<fdr & tab$logFC< -lfc,'Down-regulated','Unchanged'))
      tab$Regulation[is.na(tab$Regulation)]<-'Unchanged'  # genes without a test statistic (e.g. DESeq2 Cook's-distance outliers) are never called
      write_tsv(tab,file.path(dest,'DEG.tsv')); write_json(list(effect_test=effect,deseq2=if(effect=='deseq2')list(size_factors='TMM',test='Wald',lfc_shrinkage='ashr',independentFiltering=FALSE) else NULL,fdr_cut=fdr,lfc_cut=lfc,FDR_family='tested_genes_within_analysis_contrast',tested_genes=nrow(tab),genes_without_FDR=sum(is.na(tab$FDR)),definition=cs$details[[id]],positive_effect='weighted sum of declared model terms on log2 scale'),file.path(dest,'settings.json'))
      if(!'F' %in% names(tab)) { q<-edgeR::glmQLFTest(fit,contrast=cv); tab$F<-q$table$F; write_json(list(ranking='sign(logFC)*sqrt(QL_F)',ranking_test='QLF separate from effect-size test'),file.path(dest,'ranking_definition.json')) }
      write_tsv(data.frame(gene_id=tab$gene_id,logFC=tab$logFC,minus_log10_FDR=-log10(pmax(tab$FDR,.Machine$double.xmin)),Regulation=tab$Regulation),file.path(dest,'volcano_data.tsv'))
      with_pdf(file.path(dest,'Volcano.pdf'), { cols<-c('Unchanged'='grey65','Up-regulated'='#b2182b','Down-regulated'='#2166ac')[tab$Regulation]; plot(tab$logFC,-log10(pmax(tab$FDR,.Machine$double.xmin)),col=cols,pch=16,cex=.5,xlab='Declared contrast effect (log2)',ylab='-log10(BH FDR)',main=id); abline(v=c(-lfc,lfc),h=-log10(fdr),lty=2,col='grey50') })
      record(paste0('contrast:',id),'ok')
      if(!identical(scalar(a$enrichment,'none'),'none')) tryCatch({ assert(is.null(annotation_error),annotation_error); run_enrichment(tab,cv,y,model$d,anns,dest); record(paste0('enrichment:',id),'ok') },error=function(e) record(paste0('enrichment:',id),'failed',conditionMessage(e)))
    },error=function(e) record(paste0('contrast:',id),'failed',conditionMessage(e)))
  }
  alltests<-unique(cs$declared$test_id); alltests<-alltests[nzchar(alltests)]
  for(tid in alltests) tryCatch({
    safe_id(tid); assert(effect=='qlf','Omnibus requires effect_test=qlf'); ids<-unique(cs$declared$contrast_id[cs$declared$test_id==tid])
    assert(all(ids %in% names(cs$vectors)),'Omnibus has an invalid/missing component contrast')
    cv<-do.call(cbind,cs$vectors[ids]); assert(ncol(cv)>=2L && qr(cv)$rank==ncol(cv),'Omnibus requires at least two linearly independent contrast vectors')
    z<-edgeR::glmQLFTest(fit,contrast=cv); tab<-edgeR::topTags(z,n=Inf,sort.by='none')$table
    tab<-tab[,!grepl('^logFC',names(tab)),drop=FALSE]
    write_tsv(data.frame(analysis_id=scalar(a$analysis_id),test_id=tid,gene_id=rownames(tab),tab),path('omnibus',tid,'DEG.tsv'))
    write_json(list(df=qr(cv)$rank,contrasts=ids,FDR_family='tested_genes_within_analysis_omnibus',direction='not_defined'),path('omnibus',tid,'settings.json')); record(paste0('omnibus:',tid),'ok')
  },error=function(e)record(paste0('omnibus:',tid),'failed',conditionMessage(e)))
  status$status<<-if(any(vapply(status$stages,function(z)z$status %in% c('failed','blocked'),TRUE)))'partial_failure' else 'ok'
  write_status()
}

write_status()
tryCatch(withCallingHandlers(main(),warning=function(w) { status$warnings <<- unique(c(status$warnings,conditionMessage(w))); }),error=function(e) { status$status <<- 'failed'; status$error <<- conditionMessage(e); record('engine','failed',conditionMessage(e)); message('Analysis failed: ',conditionMessage(e)) })
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
package_names <- c('jsonlite','edgeR','limma','statmod','tximport','fgsea','sva','DESeq2','ashr')
package_versions <- setNames(lapply(package_names,function(pkg) tryCatch(as.character(utils::packageVersion(pkg)),error=function(e) NA_character_)),package_names)
write_json(list(R=R.version.string,packages=package_versions,RNG_kind=as.list(RNGkind()),seed=scalar(a$seed,'123'),platform=R.version$platform,locale=Sys.getlocale(),BLAS=extSoftVersion()[['BLAS']],thread_environment=as.list(Sys.getenv(c('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'),unset=NA_character_))),file.path(out,'runtime.json'))
write_status()
quit(status=if(status$status=='ok')0L else 1L)
