##################################################################################
############################ FINAL PROJECT #######################################
############################ ELIZABETH ETAFO #####################################
############################### 12/19/25 #########################################
##################################################################################

library(shiny)   #load shiny package

#####################################################################
######################Increasing Shiny Upload Limit##################
######################################################################

#default r shiny upload limit is 5mb. Some sequencing count files could be above this limit, 
#hence it is essential to increase the upload size limit. The below code increases 
#the upload limit to 1gb.

options(shiny.maxRequestSize=1000*1024^2)

##########################################################################
#####################Installing necessary packages#########################
##########################################################################

#Packages listed here so hosting tools (rsconnect / Posit Connect Cloud) can detect
#every dependency. This block never runs.
if(FALSE){
  library(shiny); library(DESeq2); library(AnnotationDbi); library(ComplexHeatmap)
  library(circlize); library(ggplot2); library(ggrepel); library(dplyr); library(forcats)
  library(RColorBrewer); library(DT)
  library(org.Hs.eg.db); library(org.Mm.eg.db); library(org.Rn.eg.db)
  library(org.Dm.eg.db); library(org.Dr.eg.db); library(org.Sc.sgd.db)
  library(fgsea); library(readxl)
}

loadlib <- function(extra_Bioc=NULL, other_extra=NULL){  #incase you need to load extra packages
  if(!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager") 
  }
  #organism annotation databases are large, so they are installed only when
  #first selected in the app (see load_orgdb below)
  Bioc_<- c("DESeq2", "AnnotationDbi", "ComplexHeatmap")
  Bioc_pack<- unique(c(extra_Bioc, Bioc_)) #concatenates extra biocmanager packages to default pacakages
  Bioc_missing<- Bioc_pack[!vapply(Bioc_pack, requireNamespace, 
                                   FUN.VALUE = logical(1), quietly=TRUE)]
  if(length(Bioc_missing)>0){
    BiocManager::install(Bioc_missing, ask = FALSE, update=FALSE) 
  }
  
  other_<- c("dplyr", "ggplot2", "ggrepel", "forcats", "grid", "RColorBrewer", "circlize", "DT")
  other_pack<- unique(c(other_extra, other_))#concatenates extra cran packages to default cran pacakages
  other_missing<- other_pack[!vapply(other_pack, requireNamespace, 
                                     FUN.VALUE = logical(1), quietly=TRUE)]
  if(length(other_missing)>0){
    install.packages(other_missing, ask = FALSE, update=FALSE) 
  }
  
  suppressPackageStartupMessages({
    lapply(c(Bioc_pack, other_pack), library, character.only= TRUE)
  })
  
  message("All required packages have been successfully installed and loaded")
}

########################################################################
############Install (if needed) and load an organism database##########
########################################################################
load_orgdb<- function(pkg_name){
  #pkg_name: e.g. "org.Hs.eg.db". Installed once, the first time it is used.
  if(!requireNamespace(pkg_name, quietly = TRUE)){
    if(!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(pkg_name, ask = FALSE, update = FALSE)
  }
  if(!requireNamespace(pkg_name, quietly = TRUE))
    stop(paste0("Could not install '", pkg_name, "'. Install it with BiocManager::install('", pkg_name, "')"))
  get(pkg_name, envir = asNamespace(pkg_name))
}

########################################################################
#####################Function to read count matrix######################
#######################################################################
#file name without compression and table extensions (used as a sample name)
file_stem<- function(name){
  s<- basename(name)
  s<- sub("\\.(gz|bz2|zip)$", "", s, ignore.case = TRUE)
  s<- sub("\\.(txt|tsv|csv|tab|counts?|out|htseq|featurecounts)$", "", s, ignore.case = TRUE)
  s<- sub("\\.(txt|tsv|csv|tab|counts?|out|htseq|featurecounts)$", "", s, ignore.case = TRUE)
  s
}

#tidies sample names taken from aligner output (e.g. "/data/S1Aligned.sortedByCoord.out.bam" -> "S1")
clean_sample_names<- function(x){
  x<- basename(x)
  x<- sub("\\.(bam|sam|cram)$", "", x, ignore.case = TRUE)
  x<- sub("[._]?Aligned.*$", "", x)
  x<- sub("[._]?(sorted|sort|dedup|markdup)$", "", x, ignore.case = TRUE)
  x
}

#reads one count table: a full matrix, a featureCounts table, or a per-sample HTSeq / STAR file
read_count_table<- function(path, name = basename(path)){
  #Excel files: first sheet, written to a temporary tab-separated file
  if(grepl("\\.xlsx?$", name, ignore.case = TRUE)){
    if(!requireNamespace("readxl", quietly = TRUE)) stop("Reading Excel files needs the 'readxl' package.")
    x<- as.data.frame(readxl::read_excel(path, sheet = 1, col_names = FALSE, col_types = "text", .name_repair = "minimal"))
    x<- x[!(is.na(x[[1]]) & rowSums(!is.na(x)) == 0), , drop = FALSE]
    while(nrow(x) > 0 && grepl("^#", as.character(x[1, 1]))) x<- x[-1, , drop = FALSE]   #comment lines
    x[is.na(x)]<- ""
    tmp<- tempfile(fileext = ".tsv")
    utils::write.table(x, tmp, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    path<- tmp; name<- sub("\\.xlsx?$", ".tsv", name, ignore.case = TRUE)
  }
  is_csv<- grepl("\\.csv(\\.gz)?$", name, ignore.case = TRUE)
  lines<- readLines(path, n = 60, warn = FALSE)            #readLines also reads .gz files
  lines<- lines[nzchar(lines) & !grepl("^#", lines)]
  if(length(lines) == 0) stop(paste0("'", name, "' is empty."))
  sep<- if(is_csv) "," else "\t"
  first<- strsplit(lines[1], sep, fixed = TRUE)[[1]]
  if(!is_csv && length(first) < 2) sep<- ""                  #space-separated file
  if(sep == "") first<- strsplit(trimws(lines[1]), "\\s+")[[1]]
  second<- suppressWarnings(as.numeric(gsub('"', "", first[min(2, length(first))])))
  header<- !(length(first) >= 2 && !is.na(second))          #HTSeq files have no header line
  df<- utils::read.table(path, sep = sep, header = header, comment.char = "#", quote = "\"",
                         check.names = FALSE, stringsAsFactors = FALSE, fill = TRUE)
  if(!identical(rownames(df), as.character(seq_len(nrow(df))))){
    ids<- rownames(df); vals<- df                            #gene IDs were read as row names
  } else {
    ids<- as.character(df[[1]]); vals<- df[, -1, drop = FALSE]
  }
  #featureCounts: drop the annotation columns
  fc<- c("Chr", "Start", "End", "Strand", "Length")
  if(ncol(vals) > 5 && identical(names(vals)[1:5], fc)) vals<- vals[, -(1:5), drop = FALSE]
  #STAR ReadsPerGene.out.tab: keep the unstranded counts column
  if(grepl("ReadsPerGene", name) && ncol(vals) == 3) vals<- vals[, 1, drop = FALSE]
  #drop summary rows (HTSeq "__no_feature", STAR "N_unmapped", ...)
  drop<- grepl("^__", ids) | ids %in% c("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous")
  ids<- ids[!drop]; vals<- vals[!drop, , drop = FALSE]
  #keep numeric columns only (e.g. a gene_name column is dropped)
  num<- vapply(vals, function(v){
    x<- suppressWarnings(as.numeric(v)); all(is.na(x) == is.na(v) | v %in% c("", "NA"))
  }, logical(1))
  vals<- vals[, num, drop = FALSE]
  if(ncol(vals) == 0) stop(paste0("No numeric count columns found in '", name, "'."))
  m<- data.matrix(vals)
  rownames(m)<- ids
  if(!header) colnames(m)<- if(ncol(m) == 1) file_stem(name) else paste0(file_stem(name), "_", seq_len(ncol(m)))
  else colnames(m)<- clean_sample_names(colnames(m))
  if(ncol(m) == 1 && header && grepl("^V[0-9]+$|^count", colnames(m), ignore.case = TRUE)) colnames(m)<- file_stem(name)
  m
}

#unpacks .tar / .tar.gz / .tgz / .zip uploads (e.g. GEO "_RAW.tar") into a list of files
expand_archives<- function(paths, names){
  out_p<- character(0); out_n<- character(0)
  for(i in seq_along(paths)){
    if(grepl("\\.(tar|tar\\.gz|tgz|zip)$", names[i], ignore.case = TRUE)){
      d<- tempfile("archive_"); dir.create(d)
      if(grepl("\\.zip$", names[i], ignore.case = TRUE)) utils::unzip(paths[i], exdir = d)
      else utils::untar(paths[i], exdir = d)
      f<- list.files(d, recursive = TRUE, full.names = TRUE)
      f<- f[!grepl("(^|/)(\\.|__MACOSX)", f)]
      out_p<- c(out_p, f); out_n<- c(out_n, basename(f))
    } else {
      out_p<- c(out_p, paths[i]); out_n<- c(out_n, names[i])
    }
  }
  list(paths = out_p, names = out_n)
}

#checks that values are raw counts (whole numbers >= 0) and converts them to integers
check_count_values<- function(m, round_estimated = FALSE){
  v<- m[!is.na(m)]
  note<- NULL
  if(any(v < 0))
    stop("The count table has negative values, so it looks like log-transformed or normalized data (for example log2 CPM or microarray values). DESeq2 needs raw read counts (whole numbers, 0 or more).")
  nz<- v[v != 0]
  frac<- if(length(nz)) mean(abs(nz - round(nz)) > 1e-8) else 0
  if(frac > 0.01){
    if(!round_estimated)
      stop("The values are not whole numbers, so this looks like normalized data (FPKM, RPKM, TPM or CPM). DESeq2 needs raw read counts. Look for a 'raw counts' file, or use 'Load from GEO' (NCBI raw counts). If the values are estimated counts from RSEM or Salmon, tick 'Round estimated counts' and upload again.")
    note<- "Non-integer estimated counts were rounded to whole numbers."
  }
  m<- round(m)
  storage.mode(m)<- "integer"
  list(counts = m, note = note)
}

#reads one or several count files and returns a single genes x samples integer matrix
read_counts<- function(paths, names = basename(paths), na_to_zero = FALSE, round_estimated = FALSE){
  ex<- expand_archives(paths, names)
  notes<- character(0)
  mats<- list(); failed<- character(0)
  for(i in seq_along(ex$paths)){
    m<- tryCatch(read_count_table(ex$paths[i], ex$names[i]), error = function(e) NULL)
    if(is.null(m)) failed<- c(failed, ex$names[i]) else mats[[length(mats) + 1]]<- m
  }
  if(length(mats) == 0) stop("No count table could be read from the uploaded file(s).")
  if(length(failed) > 0) notes<- c(notes, paste0("Skipped files that are not count tables: ", paste(failed, collapse = ", ")))
  if(length(mats) == 1){
    counts<- mats[[1]]
  } else {
    genes<- unique(unlist(lapply(mats, rownames)))
    counts<- do.call(cbind, lapply(mats, function(m){
      out<- matrix(NA_real_, length(genes), ncol(m), dimnames = list(genes, colnames(m)))
      out[rownames(m), ]<- m
      out
    }))
    if(anyDuplicated(colnames(counts))) colnames(counts)<- make.unique(colnames(counts))
    notes<- c(notes, paste0("Merged ", length(mats), " count files into one table (", ncol(counts), " samples)."))
  }
  
  #handles NAs
  if(anyNA(counts)){
    if(na_to_zero){
      counts[is.na(counts)]<- 0
    }else{
      stop("NA values in counts, set na_to_zero = TRUE")
    }
  }
  chk<- check_count_values(counts, round_estimated)
  counts<- chk$counts
  if(!is.null(chk$note)) notes<- c(notes, chk$note)
  
  #handles duplicate sample names/ gene IDs
  if (anyDuplicated(colnames(counts))){
    stop("Duplicated sample found in count data")
  }
  if (anyDuplicated(rownames(counts))){
    dup<- sum(duplicated(rownames(counts)))
    message("Removing ", dup, " duplicated genes IDs")
    counts<- counts[!duplicated(rownames(counts)), , drop = FALSE]
  }
  list(counts = counts, notes = notes)
}

#reads a sample sheet (metadata) file
read_meta<- function(met){
  #similar to feature counts, metadata may be .csv, .txt
  if(grepl("\\.csv$", met, ignore.case = TRUE)){
    meta<- read.csv(met, check.names = FALSE, as.is=TRUE)
  } else{
    meta<- read.delim(met, check.names = FALSE, as.is=TRUE)
  }
  meta[[1]]<- as.character(meta[[1]])
  if (anyDuplicated(meta[[1]])){
    meta[[1]]<- make.unique(meta[[1]])
  }
  rownames(meta)<- trimws(meta[,1]) #trimws removes whitespace 
  meta<- meta[, -1, drop=FALSE]
  #trim stray spaces in group labels ("Activated " and "Activated" are the same group)
  meta[]<- lapply(meta, function(x) if(is.character(x)) trimws(x) else x)
  meta
}

#reads counts and (optionally) metadata
#example use: read_file("counts.csv", "metadata.csv")
#without a metadata file, a sample sheet with an empty 'Group' column is created to fill in
read_file<- function(cnts, met = NULL, na_to_zero = FALSE, cnts_names = basename(cnts), round_estimated = FALSE){
  rc<- read_counts(cnts, cnts_names, na_to_zero = na_to_zero, round_estimated = round_estimated)
  counts<- rc$counts
  if(is.null(met)){
    meta<- data.frame(Group = rep("", ncol(counts)), row.names = colnames(counts), stringsAsFactors = FALSE)
  } else {
    meta<- read_meta(met)
    #same samples in a different order: put the metadata in the order of the count columns
    if(setequal(rownames(meta), colnames(counts)) && !identical(rownames(meta), colnames(counts)))
      meta<- meta[colnames(counts), , drop = FALSE]
  }
  invisible(list(counts = counts, metadata = meta, notes = rc$notes))
}

########################################################################
##############Functions to load a dataset directly from GEO#############
########################################################################
#downloads a file; FALSE if it failed or NCBI returned a web page instead of data
download_geo_file<- function(url, dest){
  ok<- tryCatch({ utils::download.file(url, dest, mode = "wb", quiet = TRUE); TRUE },
                error = function(e) FALSE, warning = function(w) FALSE)
  if(!ok || !file.exists(dest) || file.size(dest) == 0) return(FALSE)
  con<- file(dest, "rb"); magic<- readBin(con, "raw", 2); close(con)
  if(identical(magic, as.raw(c(0x1f, 0x8b)))) return(TRUE)        #gzip file
  head<- tolower(paste(readLines(dest, n = 30, warn = FALSE), collapse = " "))
  !grepl("<html|<!doctype|recaptcha|checking your browser", head)
}

#parses the sample information from GEO Series Matrix file(s) into a sample sheet
#attributes: "ids" (names each sample may have in count files) and "organism"
parse_series_matrix<- function(paths){
  rows<- list(); ids<- list(); orgs<- character(0)
  for(p in paths){
    ln<- readLines(p, warn = FALSE)
    ln<- ln[grepl("^!Sample_", ln)]
    sp<- lapply(strsplit(ln, "\t", fixed = TRUE), function(x) gsub('^"|"$', "", x))
    key<- vapply(sp, `[`, character(1), 1)
    get1<- function(k){ i<- which(key == k); if(length(i)) sp[[i[1]]][-1] else NULL }
    gsm<- get1("!Sample_geo_accession")
    if(is.null(gsm)) next
    df<- data.frame(row.names = gsm, check.names = FALSE, stringsAsFactors = FALSE)
    ch<- sp[key == "!Sample_characteristics_ch1"]
    for(r in ch){
      vals<- r[-1]
      for(j in seq_along(vals)){
        v<- vals[j]; if(!nzchar(v)) next
        if(grepl(":", v)){ nm<- trimws(sub(":.*$", "", v)); val<- trimws(sub("^[^:]*:", "", v)) }
        else { nm<- "characteristic"; val<- trimws(v) }
        if(!nm %in% names(df)) df[[nm]]<- NA_character_
        df[gsm[j], nm]<- val
      }
    }
    src<- get1("!Sample_source_name_ch1"); if(!is.null(src)) df[["source"]]<- trimws(src)
    tit<- get1("!Sample_title"); if(!is.null(tit)) df[["title"]]<- trimws(tit)
    org<- get1("!Sample_organism_ch1"); if(!is.null(org)) orgs<- c(orgs, org)
    #other names of each sample (title, library names, descriptions) to match count-file columns
    extra<- c(list(tit), lapply(sp[key %in% c("!Sample_description", "!Sample_library_name")], `[`, -1))
    for(j in seq_along(gsm)){
      v<- unlist(lapply(extra, function(x) if(length(x) >= j) x[j] else NULL))
      v<- trimws(sub("^[^:]*:\\s*", "", v[nzchar(v)]))            #"Library name: sample 1" -> "sample 1"
      ids[[gsm[j]]]<- unique(c(gsm[j], v))
    }
    rows[[length(rows) + 1]]<- df
  }
  if(length(rows) == 0) stop("No sample information found in the GEO Series Matrix.")
  cols<- unique(unlist(lapply(rows, names)))
  meta<- do.call(rbind, lapply(rows, function(d){ for(c in setdiff(cols, names(d))) d[[c]]<- NA_character_; d[, cols, drop = FALSE] }))
  #column names that work in a design formula, e.g. "cell line" -> "cell_line"
  nm<- gsub("[^A-Za-z0-9]+", "_", names(meta)); nm<- gsub("^_|_$", "", nm)
  nm<- ifelse(grepl("^[0-9]", nm) | !nzchar(nm), paste0("X", nm), nm)
  names(meta)<- make.unique(nm, sep = "_")
  #columns with one value for every sample cannot define groups, so they are dropped (title is kept)
  keep<- vapply(names(meta), function(c) c == "title" || length(unique(meta[[c]])) > 1, logical(1))
  meta<- meta[, keep, drop = FALSE]
  meta[is.na(meta)]<- ""
  attr(meta, "ids")<- ids
  attr(meta, "organism")<- if(length(orgs)) names(sort(table(orgs), decreasing = TRUE))[1] else NA_character_
  meta
}

#matches count-table columns to GEO samples (GSM IDs) using the sample names GEO lists
#returns a vector of GSM IDs (NA where no unique match)
match_samples_to_geo<- function(cols, ids){
  norm<- function(x) tolower(gsub("[^A-Za-z0-9]", "", x))
  nc<- norm(cols)
  vapply(seq_along(cols), function(i){
    hit<- names(ids)[vapply(ids, function(v){
      v<- norm(v); v<- v[nchar(v) >= 4]
      any(v == nc[i]) || any(vapply(v, function(x) grepl(x, nc[i], fixed = TRUE), logical(1)))
    }, logical(1))]
    if(length(hit) == 1) hit else NA_character_
  }, character(1))
}

#guesses the gene ID type from the IDs (used to pre-set the annotation)
guess_id_type<- function(ids){
  ids<- utils::head(ids[!is.na(ids)], 500)
  if(mean(grepl("^ENS[A-Z]*G[0-9]", ids)) > 0.5) return("ensembl")
  if(mean(grepl("^[0-9]+$", ids)) > 0.5) return("entrezid")
  "symbol"
}

#maps a GEO organism name to the app's organism choices
organism_choice<- function(org){
  if(is.na(org)) return(NA_character_)
  m<- c("Homo sapiens" = "human", "Mus musculus" = "mouse", "Rattus norvegicus" = "rat",
        "Drosophila melanogaster" = "fly", "Danio rerio" = "zfish", "Saccharomyces cerevisiae" = "yeast")
  if(org %in% names(m)) m[[org]] else "other"
}

#lists file names in a GEO folder (Series Matrix or supplementary files)
list_geo_dir<- function(url){
  idx<- tempfile(fileext = ".html")
  ok<- tryCatch({ utils::download.file(url, idx, mode = "wb", quiet = TRUE); TRUE }, error = function(e) FALSE, warning = function(w) FALSE)
  if(!ok || !file.exists(idx)) return(character(0))
  html<- paste(readLines(idx, warn = FALSE), collapse = " ")
  f<- regmatches(html, gregexpr('href="[^"/?][^"]*"', html))[[1]]
  f<- utils::URLdecode(gsub('^href="|"$', "", f))
  unique(basename(f[grepl("^GSE", basename(f))]))
}

#downloads raw counts and sample information for a GEO series (GSE)
#counts: NCBI-generated raw counts if available, otherwise the authors' count file(s)
load_geo<- function(gse,
                    base_url = getOption("rnaseqapp.geo_base", "https://www.ncbi.nlm.nih.gov"),
                    ftp_url  = getOption("rnaseqapp.geo_ftp",  "https://ftp.ncbi.nlm.nih.gov")){
  gse<- toupper(trimws(gse))
  if(!grepl("^GSE[0-9]+$", gse)) stop("Enter a GEO series accession such as GSE106305.")
  notes<- character(0)
  stub<- sub("[0-9]{1,3}$", "nnn", gse)
  sdir<- paste0(ftp_url, "/geo/series/", stub, "/", gse, "/")
  
  #1. sample information from the Series Matrix file(s)
  sm<- tempfile(fileext = ".txt.gz")
  sm_files<- if(download_geo_file(paste0(sdir, "matrix/", gse, "_series_matrix.txt.gz"), sm)) sm else character(0)
  if(length(sm_files) == 0){
    for(x in grep("_series_matrix\\.txt\\.gz$", list_geo_dir(paste0(sdir, "matrix/")), value = TRUE)){
      d<- tempfile(fileext = ".txt.gz"); if(download_geo_file(paste0(sdir, "matrix/", x), d)) sm_files<- c(sm_files, d)
    }
  }
  meta<- if(length(sm_files) > 0) parse_series_matrix(sm_files) else NULL
  
  #2a. NCBI-generated raw counts (human, then mouse)
  builds<- c(human = "GRCh38.p13", mouse = "GRCm39")
  cnt_file<- tempfile(fileext = ".tsv.gz"); organism<- NULL; id_type<- "entrezid"; counts<- NULL; source<- NULL
  for(org in names(builds)){
    url<- paste0(base_url, "/geo/download/?type=rnaseq_counts&acc=", gse, "&format=file&file=",
                 gse, "_raw_counts_", builds[[org]], "_NCBI.tsv.gz")
    if(download_geo_file(url, cnt_file)){
      organism<- org
      counts<- check_count_values(read_count_table(cnt_file, paste0(gse, "_raw_counts.tsv.gz")))$counts
      source<- "NCBI"
      break
    }
  }
  
  #2b. otherwise: the authors' supplementary count file(s)
  tried<- character(0); normalized<- character(0)
  if(is.null(counts)){
    sup<- list_geo_dir(paste0(sdir, "suppl/"))
    table_ext<- "\\.(txt|tsv|csv|tab|xlsx|xls)(\\.gz)?$"
    cand<- sup[grepl(table_ext, sup, ignore.case = TRUE) & !grepl("_RAW\\.tar$", sup)]
    looks_norm<- grepl("fpkm|rpkm|tpm|cpm|normali[sz]ed|vst|rlog|log2|deseq|edger|limma|diff|result|deg", cand, ignore.case = TRUE)
    looks_count<- grepl("count|featurecount|htseq|raw|reads|expected|matrix", cand, ignore.case = TRUE)
    cand<- c(cand[looks_count & !looks_norm], cand[!looks_count & !looks_norm], cand[looks_norm])
    for(x in utils::head(cand, 6)){
      d<- file.path(tempdir(), paste0("geo_", x))
      if(!download_geo_file(paste0(sdir, "suppl/", utils::URLencode(x)), d)) next
      tried<- c(tried, x)
      m<- tryCatch(read_count_table(d, x), error = function(e) NULL)
      if(is.null(m) || ncol(m) < 2) next
      chk<- tryCatch(check_count_values(m), error = function(e) NULL)
      if(is.null(chk)){ normalized<- c(normalized, x); next }
      counts<- chk$counts; source<- x; id_type<- guess_id_type(rownames(counts))
      break
    }
    organism<- if(!is.null(meta)) organism_choice(attr(meta, "organism")) else NA_character_
  }
  
  if(is.null(counts)){
    msg<- paste0("No raw counts found for ", gse, ". NCBI has not made raw counts for this study")
    if(length(normalized) > 0) msg<- paste0(msg, ", and the authors' file(s) contain normalized values, not raw counts (", paste(normalized, collapse = ", "), "). DESeq2 needs raw counts.")
    else if(length(tried) > 0) msg<- paste0(msg, ", and the authors' file(s) could not be read as a count table (", paste(tried, collapse = ", "), "). Download them from GEO and upload them with 'Upload files'.")
    else msg<- paste0(msg, ", and no count table was found among the authors' files. If the study has a '_RAW.tar' of per-sample count files, download it from GEO and upload it with 'Upload files'.")
    stop(msg)
  }
  
  #3. put counts and sample information together
  if(!is.null(meta)){
    ids<- attr(meta, "ids")
    if(!all(colnames(counts) %in% rownames(meta))){
      #the authors' columns use their own sample names: match them to GEO samples
      gsm<- match_samples_to_geo(colnames(counts), ids)
      if(!anyNA(gsm) && !anyDuplicated(gsm)){
        colnames(counts)<- gsm
        notes<- c(notes, "Matched the count-file columns to the GEO samples by sample name.")
      }
    }
    both<- rownames(meta)[rownames(meta) %in% colnames(counts)]   #GEO sample order
    if(length(both) >= 2){
      if(length(both) < ncol(counts)) notes<- c(notes, paste0(ncol(counts) - length(both), " samples without sample information were left out."))
      if(length(both) < nrow(meta)) notes<- c(notes, paste0(nrow(meta) - length(both), " samples in the series have no counts in this file (for example non-RNA-seq samples) and were left out."))
      counts<- counts[, both, drop = FALSE]; meta<- meta[both, , drop = FALSE]
    } else {
      #could not match: keep the count columns and let the user fill in the groups
      meta<- data.frame(Group = rep("", ncol(counts)), row.names = colnames(counts), stringsAsFactors = FALSE)
      notes<- c(notes, "The sample names in the count file could not be matched to the GEO samples. Fill in the Group column in the Sample Sheet tab (the GEO page lists which sample is which).")
    }
  } else {
    meta<- data.frame(Group = rep("", ncol(counts)), row.names = colnames(counts), stringsAsFactors = FALSE)
    notes<- c(notes, "The sample information could not be downloaded. Fill in the Group column in the Sample Sheet tab.")
  }
  attr(meta, "ids")<- NULL; attr(meta, "organism")<- NULL
  
  if(identical(source, "NCBI")){
    notes<- c(notes, paste0("Loaded NCBI raw counts for ", gse, " (", organism, ", ", ncol(counts), " samples, gene IDs are Entrez IDs). ",
                            "NCBI counts are made with one standard pipeline, so they may differ slightly from the original study."))
  } else {
    notes<- c(notes, paste0("NCBI raw counts are not available for ", gse, ", so the authors' count file was used: ", source,
                            " (", ncol(counts), " samples). Check the organism and ID type settings."))
  }
  list(counts = counts, metadata = meta, organism = organism, id_type = id_type, notes = notes, gse = gse)
}


########################################################################
################Function to quality check data##########################
#######################################################################
qc_check<- function(counts, meta){
  if(!all(colnames(counts)%in% rownames(meta)))
    stop("Some samples in counts data missing from metadata")
  if(!all(colnames(counts) == rownames(meta)))
    stop("Mismatch in sample order between counts and metadata")
  if(ncol(counts) != nrow(meta))
    stop("Difference in number of samples in counts and metadata")
  message("Quality checks passed")
}

#######################################################################
#######################Function to prepare metadata####################
#######################################################################
prepare_meta<- function(meta, Exp_groups=NULL, ref_levels = list()){
  #meta: metadata dataframe
  #Exp_groups: character vector of columns to use as experimental factors, e.g. c("Activation", "Infection")
  #ref_levels: list of reference levels, e.g.list(Infection = "Uninfected", Activation = "Resting")
  
  #Some quality check of metadata 
  if(is.null(Exp_groups)){ #warning message if Exp_group is not provided
    stop("Provide experimental grouping column names (Exp_groups).")
  }
  missing_cols<- setdiff(Exp_groups, names(meta))
  if (length(missing_cols) > 0){  #warning if experiment groups are not found in metadata
    stop(sprintf("These Exp_groups are missing from metadata: %s", paste(missing_cols, collapse = ",")))
  }
  
  #convert all experiment groups columns to factors
  for(col in Exp_groups){
    meta[[col]]<- as.factor(meta[[col]])
  }
  
  #set reference levels
  
  for (nm in names(ref_levels)){
    if (nm %in% names(meta)){
      meta[[nm]]<- relevel(meta[[nm]], ref=ref_levels[[nm]])
    }else{
      warning(sprintf("Column '%s' not found in metadata, skipping reference level.", nm))
    }
  }
  invisible(meta)
}

#######################################################################
#########################Function to run DESeq2#######################
#######################################################################
run_deseq<- function(counts, meta, Exp_design=NULL, min_count=10){
  #counts and meta defined from previous function
  #Exp_design: to be used later in function (e.g. Exp_design= ~type)
  
  #Check if DESeq2 is installed and load the package
  if(!requireNamespace("DESeq2", quietly = TRUE)){
    if(!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager") 
    BiocManager::install("DESeq2", ask = FALSE, update=FALSE)
  }
  suppressPackageStartupMessages(library(DESeq2))
  
  #validate function
  if(is.null(Exp_design))
    stop("Provide a valid DESeq2 design formula (Exp_design)")
  if(!inherits(Exp_design, "formula"))
    stop("Exp_design must be a formula object, e.g. ~Condition or ~Activation + Infection")
  
  #Create DESeq2 dataset
  dds<- DESeqDataSetFromMatrix(countData = counts, 
                               colData = meta,
                               design = Exp_design )
  #Prefilter low count genes
  keep<- rowSums(counts(dds)) >= min_count
  dds <- dds[keep,]
  
  #Run DESeq normalization and model fitting
  deg <- DESeq(dds)
  
  message("DESeq2 analysis completed")
  deg
}

#########################################################################
##############Function to extract results and normalization##############
#########################################################################
extract_results<- function(deg, contrast=NULL, alpha=NULL, normalized= TRUE){
  #Example use: contrast= c("type", "normal", "cancer")----in metadata we can 2 controls in the "type" column, normal and cancer
  #some basic checks
  if(!inherits(deg, "DESeqDataSet"))
    stop("'deg' must be a DESeqDataSet (output of DESeq())")
  if (is.null(alpha)) alpha <- 0.05
  if (!is.numeric(alpha) || length(alpha) !=1 || alpha <=0 || alpha >= 1)
    stop("'alpha' must be a single number in (0, 1)")
  
  #result with and without contrast
  if(is.null(contrast)){
    res<- results(deg, alpha = alpha)
  } else if (is.character(contrast) && length(contrast) ==1){
    res <-results(deg, name = contrast, alpha = alpha)
  } else{
    res<- results(deg, contrast = contrast, alpha= alpha)
  }
  
  #Normalized counts matrix
  norm<- counts(deg, normalized= normalized)
  
  list(count_results= res, normalized_counts= norm)
}

########################################################################
#################Function to annotate gene names#######################
########################################################################
annotate_results<- function(res, orgdb, id_type=NULL, 
                            to_symbol="SYMBOL"){
  
  #basic checks
  if(missing(orgdb)) stop("Please pass an OrgDb object, e.g. org.Hs.eg.db or org.Mm.eg.db")
  if(!inherits(orgdb, "OrgDb")) stop("'orgdb' must be an OrgDb (e.g. org.Hs.eg.db) ")
  
  if(!requireNamespace("AnnotationDbi", quietly = TRUE)) 
    stop("Install package 'AnnotationDbi'. Install via BiocManager::install('AnnotationDbi').")
  
  resdf<- as.data.frame(res)
  genekeys<- rownames(resdf)
  if(is.null(genekeys)) stop("Row names of 'res' must contain gene IDs") 
  
  
  #guess ID type if not provided
  guess_idtype <- function(ids){
    if (all(grepl("^ENSG\\d+$", sub("\\.\\d+$", "", ids)))) return ("ENSEMBL")
    if (all(grepl("^ENSMUSG\\d+$", sub("\\.\\d+$", "", ids)))) return("ENSEMBL")
    if (all(grepl("^\\d+$", ids))) return("ENTREZID")
    if (all(grepl("^[A-Za-z0-9._-]+$", ids))) return("SYMBOL")
    return(NA_character_)
  }
  if(is.null(id_type)){
    id_type<- guess_idtype(genekeys)
  }
  
  #validate ID types
  ktypes <- AnnotationDbi::keytypes(orgdb)
  if (is.na(id_type) || !(id_type%in% ktypes)) {
    stop(sprintf("Could not determine valid id_type. Please specify one explicitly.
                  Available types: %s", paste(ktypes, collapse= ",")))
  }
  
  #perform gene mapping  
  resdf$GeneName<- AnnotationDbi::mapIds(orgdb,
                                         keys= genekeys,
                                         keytype=id_type,
                                         column= to_symbol,
                                         multiVals = "first")
  
  attr(resdf, "id_type_used")<- id_type
  resdf
}

#######################################################################
######################Function to get top genes#######################
#######################################################################
top_genes<- function(resdf, padj="padj", n=15){ #can change n to number of best genes you want to extract
  if(!requireNamespace("dplyr", quietly= TRUE)) #warning message if you dont have package
    stop("Install 'dplyr' package first")
  
  if(!padj %in% names(resdf))  #warning message is padj is not in your result dataframe
    stop("Column 'padj' not found in results dataframe")
  
  best_genes<- dplyr::arrange(resdf, !!as.name(padj))
  head(best_genes, n)
}

#######################################################################
##################Function to generate volcano plot####################
########################################################################
plot_volcano<- function(resdf, best_genes= NULL, out_png= NULL, 
                        log2FoldChange="log2FoldChange",
                        padj= "padj", 
                        lfc_thr= 1, padj_thr= 0.05,
                        label_col= "GeneName", 
                        main= "Differential Expression Volcano Plot") {
  
  #Some checks
  if(!all(c(log2FoldChange, padj) %in% names(resdf))) #checks if columns required for plotting are in data frame
    stop("Missing some columns for volcano plot")
  if(!requireNamespace("ggplot2", quietly = TRUE)) #checks if ggplot2 is installed
    stop("'ggplot2' package required for plotting volcano plot")
  
  #removes NAs
  keep<- !is.na(resdf[[padj]]) & !is.na(resdf[[log2FoldChange]]) 
  vol<- resdf[keep, , drop=FALSE]
  
  # .sig marks significat genes (pvalue <= threshold and log2foldchange >= threshold)
  vol$neglog10<- -log10(vol[[padj]])
  
  vol$status<- "NS"
  vol$status [vol[[padj]]<= padj_thr&vol[[log2FoldChange]]>= lfc_thr]<- "Up"
  vol$status [vol[[padj]]<= padj_thr&vol[[log2FoldChange]]<= -lfc_thr]<- "Down"
  vol$status<- factor(vol$status, levels=c("Up","Down", "NS"))
  
  #building the ggplot using the ggplot2 package
  vol_plot<- ggplot2::ggplot(
    vol,
    ggplot2::aes(x=.data[[log2FoldChange]], y=.data[["neglog10"]])) + 
    ggplot2::geom_point(ggplot2::aes(color=.data[["status"]]), size= 1.3, alpha=0.7) +
    ggplot2::geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype= "dashed") +
    ggplot2::geom_hline(yintercept = -log10(padj_thr), linetype ="dashed") +
    ggplot2::scale_color_manual(
      values=c(Up="#F72585", Down="#3B4CC0", NS="grey50"),
      name="Differential Expression"
    )+
    ggplot2::labs(x="log2 fold change", y="-log10(adjusted p-value)", title= main)+
    ggplot2::theme_minimal(base_size = 14)+
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face="bold"),
                   panel.grid = ggplot2::element_blank(),
                   axis.line = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   legend.position = "right",
                   legend.title = ggplot2::element_text(face="bold", size = 12),
                   legend.text = ggplot2::element_text(size=10)
    )
  
  if(!is.null(best_genes) && nrow(best_genes)>0){
    lab_rows<- intersect(row.names(best_genes), row.names(vol))
    if (length(lab_rows)) {
      label_df <- vol[lab_rows, , drop = FALSE]
      # choose label text: use label_col if present, else fallback to rownames
      label_df$..label.. <- if (label_col %in% names(label_df)) label_df[[label_col]] else rownames(label_df)
      
      #adding labels in lot for top genes
      vol_plot<- vol_plot +
        ggrepel::geom_text_repel(data= label_df, 
                                 ggplot2::aes(x=.data[[log2FoldChange]],
                                              y=.data[["neglog10"]],
                                              label= .data[[label_col]]), 
                                 size=3,
                                 max.overlaps = 20,
                                 box.padding = 0.4,
                                 point.padding = 0.3,
                                 segment.color="grey30",
                                 force = 2,
                                 min.segment.length = 0)  
      
    }
  }
  if (!is.null(out_png)){
    ggplot2::ggsave(filename = out_png, 
                    plot=vol_plot, width = 7,
                    height = 6, dpi = 150)
  } else {
    print(vol_plot)
  }
  invisible(vol_plot)
}
########################################################################
#######################Function to plot heat map########################
########################################################################
plot_Heatmap<- function(deg, resdf, out_png= NULL, n=30, main="Top Differentially Expressed Genes", mat=NULL){
  #check for required package
  if(!requireNamespace("ComplexHeatmap", quietly = TRUE))
    stop("'ComplexHeatmap' package required")
  if(!requireNamespace("circlize", quietly = TRUE))
    stop("'circlize' package required")
  if(!requireNamespace("grid", quietly = TRUE))
    stop("'grid' package required")
  
  use_genes<- top_genes(resdf, n=n) #extracts topgenes from the top genes function
  
  if(is.null(mat)){
    mat<- counts(deg, normalized=TRUE) #extract normalized counts from DESeq object 
  }
  
  genes<- intersect(rownames(use_genes), row.names(mat)) #only keep genes presne ton both top_genes and expression matrix
  mat<- mat[genes, , drop=FALSE]
  mat.z<- t(scale(t(mat))) #zscore normalize each gene across samples
  mat.z[is.na(mat.z)] <- 0 #repace NA with zero
  labs<- use_genes[genes, "GeneName"]
  labs <- ifelse(is.na(labs) | labs == "", genes, labs) 
  rownames(mat.z)<- make.unique(labs)
  
  design_var<- attr(terms(design(deg)), "term.labels")[1] #extracts first design variable from DESeq2 design
  
  group_vec<- NULL
  top_annot<- NULL
  my_cols <- c(    "#4D90BD","#E06C9F", "#8BC77B","#CFA6E8","#F2C14E")
  if(!is.null(design_var)){
    cd<- as.data.frame(SummarizedExperiment::colData(deg))
    cd<- cd[colnames(mat.z), , drop=FALSE]  #match sample metadata to mat.z column
    if(design_var%in% names(cd)){
      group_vec<- factor(cd[[design_var]])
      levels(group_vec) <- tools::toTitleCase(levels(group_vec)) #get groouping factor for samples
      cols<- setNames(my_cols[seq_len(nlevels(group_vec))],
                      levels(group_vec))
      
      #building the heatmap annotaion bar showing sample groups
      top_annot<- ComplexHeatmap::HeatmapAnnotation(
        Group=group_vec,
        col = list(Group=cols),
        annotation_legend_param = list(
          title="Groups",
          title_gp=grid::gpar(fontsize=12, fontface="bold"),
          labels_gp=grid::gpar(fontsize=10)
        )
      )
    }
  }
  #generating colo scale for zscores
  heat_colors<- circlize::colorRamp2(
    c(-3,0,3), c("#3B4CC0", "white","#F72585")
  )
  #generating heatmap using the complexheatmap package
  Heat_plot<- ComplexHeatmap::Heatmap(mat.z, 
                                      name = "Z-score",
                                      heatmap_legend_param = list(
                                        title_gp=grid::gpar(fontsize=12, fontface="bold"),
                                        labels_gp=grid::gpar(fontsize=10)
                                      ),
                                      col=heat_colors,
                                      top_annotation = top_annot,
                                      cluster_rows=TRUE, 
                                      cluster_columns=TRUE,
                                      column_split = group_vec,
                                      show_row_names = TRUE,
                                      show_column_names = FALSE,
                                      column_title = main,
                                      column_title_gp = grid::gpar(fontsize=16, fontface="bold"),
                                      rect_gp = grid::gpar(col="grey70", lwd=0.3))
  if (!is.null(out_png)) {
    png(out_png, width = 1200, height = 1200, res = 150)
    on.exit(dev.off(), add = TRUE)
  }
  suppressMessages(ComplexHeatmap::draw(Heat_plot))
}
######################################################################
####################Function for Barplot##############################
#####################################################################
plot_bar<- function(resdf, n=30,log2FoldChange="log2FoldChange",
                    padj= "padj",label_col="GeneName", x_limit=NULL,
                    direction= c("up", "down", "all"),
                    bar_color="#F72585", color_by_direction=TRUE,
                    up_color="#F72585",
                    down_color="#3B4CC0",
                    neutral_color="grey70",
                    title="Differential Expression Bar Plot",
                    out_png=NULL){
  
  #Some checks
  if(!all(c(log2FoldChange, padj) %in% names(resdf))) #checks if columns required for plotting are in data frame
    stop("Missing some columns for bar plot")
  if(!requireNamespace("ggplot2", quietly = TRUE)) #checks if ggplot2 is installed
    stop("'ggplot2' package required")
  direction<- match.arg(direction)
  #removes NAs
  keep<- !is.na(resdf[[padj]]) & !is.na(resdf[[log2FoldChange]]) 
  bardf<- resdf[keep, , drop=FALSE]
  
  use_genes<- top_genes(bardf, n=n)
  
  #specifies direction of the bar plot
  if(direction == "up") use_genes<- dplyr::filter(use_genes, log2FoldChange>0)
  if(direction == "down") use_genes<- dplyr::filter(use_genes, log2FoldChange<0)
  use_genes<- dplyr::filter(use_genes, !is.na(.data[[label_col]]))
  
  #reorder factor levels by Log2FC so bars are correctly sorted
  use_genes[[label_col]] <- forcats::fct_reorder(
    use_genes[[label_col]], use_genes[[log2FoldChange]]
  )
  
  if(color_by_direction&&direction=="all"){
    use_genes<- dplyr::mutate(
      use_genes,
      Direction=dplyr::case_when(
        .data[[log2FoldChange]]>0 ~"Up",
        .data[[log2FoldChange]]<0 ~"Down",
        TRUE  ~"No change"
      )
    )
    #building bar plot with color matching the direction
    barpl<- ggplot2::ggplot(
      use_genes,
      ggplot2::aes(x=.data[[label_col]],
                   y=.data[[log2FoldChange]],
                   fill=.data$Direction)
    )+
      ggplot2::geom_col(width=0.7)+
      ggplot2::scale_fill_manual(
        values = c(Up=up_color, Down=down_color, `No change`=neutral_color),
        name="Expression Change"
      )
  }else{
    barpl<- ggplot2::ggplot(
      use_genes,
      ggplot2::aes(x=.data[[label_col]],y=.data[[log2FoldChange]])
    )+
      ggplot2::geom_col(width=0.7, fill=bar_color)
  }
  #generates barplot
  barpl<- barpl+ 
    ggplot2::coord_flip()+
    ggplot2::geom_hline(yintercept = 0, color="grey30")+
    ggplot2::labs(x=NULL, y="log2(fold change)", title = title)+
    ggplot2::theme_classic(base_size = 12)+
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust=0.5, face="bold"),
      axis.ticks.y = ggplot2::element_line(),
      axis.line = ggplot2::element_line(color="black"),
      axis.text.y = ggplot2::element_text(size=8, face = "bold"),
      axis.text.x = ggplot2::element_text(size=8, face = "bold"),
      axis.title.x = ggplot2::element_text(size = 10, face = "bold"),
      axis.title.y = ggplot2::element_text(size = 10, face = "bold"),
      legend.text=ggplot2::element_text(size = 8),
      legend.title=ggplot2::element_text(size=9,face = "bold")
    )
  if (!is.null(x_limit)) {
    barpl <- barpl + ggplot2::scale_y_continuous(limits = x_limit)
  }
  if (!is.null(out_png)) {
    ggplot2::ggsave(filename = out_png, plot = barpl, width = 6, height = 4, dpi = 300)
  }
  print(barpl)
  invisible(barpl)
}

########################################################################
##################Quality-control and extra plot functions##############
########################################################################
#All functions return a ggplot object (or draw a base/ComplexHeatmap plot)
#and, when out_png is given, save the plot as a PNG file instead.

save_gg<- function(plot, out_png, width=9, height=7){
  ggplot2::ggsave(out_png, plot, width=width, height=height, dpi=150, bg="white")
  invisible(out_png)
}

##################Variance-stabilised data with a fallback##############
safe_vst<- function(dseq){
  #vst() needs at least ~1000 genes with reasonable counts; for very small
  #datasets fall back to varianceStabilizingTransformation()
  tryCatch(DESeq2::vst(dseq, blind=TRUE),
           error=function(e) tryCatch(DESeq2::varianceStabilizingTransformation(dseq, blind=TRUE),
                                      error=function(e2) NULL))
}

#######################PCA plot of samples#############################
plot_pca<- function(vsd_mat, meta, color_by, ntop=500, label_samples=TRUE,
                    title="PCA of samples (VST, top variable genes)", out_png=NULL){
  #uses the ntop most variable genes, like DESeq2::plotPCA
  rv<- apply(vsd_mat, 1, stats::var)
  sel<- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  pca<- stats::prcomp(t(vsd_mat[sel, , drop=FALSE]))
  pct<- round(100*pca$sdev^2/sum(pca$sdev^2), 1)
  df<- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], sample=colnames(vsd_mat))
  df$group<- if(!is.null(color_by) && color_by %in% names(meta)) as.factor(meta[colnames(vsd_mat), color_by]) else "all"
  pl<- ggplot2::ggplot(df, ggplot2::aes(x=.data$PC1, y=.data$PC2, color=.data$group))+
    ggplot2::geom_point(size=4, alpha=0.85)+
    ggplot2::labs(x=paste0("PC1 (", pct[1], "% variance)"),
                  y=paste0("PC2 (", pct[2], "% variance)"),
                  color=color_by, title=title)+
    ggplot2::theme_minimal(base_size=14)+
    ggplot2::theme(plot.title=ggplot2::element_text(face="bold", hjust=0.5))
  if(label_samples && requireNamespace("ggrepel", quietly=TRUE))
    pl<- pl+ggrepel::geom_text_repel(ggplot2::aes(label=.data$sample), size=3.2, show.legend=FALSE, max.overlaps=30)
  if(!is.null(out_png)) return(save_gg(pl, out_png))
  pl
}

#################Sample-to-sample distance heatmap######################
plot_sample_distances<- function(vsd_mat, meta, annotate_by=NULL, out_png=NULL,
                                 title="Sample-to-sample distances (VST)"){
  d<- as.matrix(stats::dist(t(vsd_mat)))
  top_ann<- NULL
  if(!is.null(annotate_by) && annotate_by %in% names(meta)){
    grp<- as.character(meta[colnames(d), annotate_by])
    lv<- unique(grp)
    cols<- stats::setNames(grDevices::hcl.colors(max(3, length(lv)), "Dark 3")[seq_along(lv)], lv)
    top_ann<- ComplexHeatmap::HeatmapAnnotation(df=stats::setNames(data.frame(grp), annotate_by),
                                                col=stats::setNames(list(cols), annotate_by))
  }
  ht<- ComplexHeatmap::Heatmap(d, name="Distance", column_title=title,
                               col=circlize::colorRamp2(c(0, max(d)), c("#08306B", "white")),
                               top_annotation=top_ann,
                               clustering_distance_rows=stats::as.dist(d),
                               clustering_distance_columns=stats::as.dist(d),
                               row_names_gp=grid::gpar(fontsize=9), column_names_gp=grid::gpar(fontsize=9))
  if(!is.null(out_png)){
    grDevices::png(out_png, width=1500, height=1300, res=150)
    on.exit(grDevices::dev.off())
  }
  ComplexHeatmap::draw(ht)
  invisible(NULL)
}

####################Library size per sample#############################
plot_library_size<- function(dseq, meta, color_by=NULL, out_png=NULL,
                             title="Library size (total counts) per sample"){
  df<- data.frame(sample=colnames(dseq),
                  millions=colSums(DESeq2::counts(dseq, normalized=FALSE))/1e6,
                  size_factor=DESeq2::sizeFactors(dseq))
  df$group<- if(!is.null(color_by) && color_by %in% names(meta)) as.factor(meta[df$sample, color_by]) else "all"
  df$sample<- factor(df$sample, levels=df$sample)
  pl<- ggplot2::ggplot(df, ggplot2::aes(x=.data$sample, y=.data$millions, fill=.data$group))+
    ggplot2::geom_col(width=0.7)+
    ggplot2::geom_text(ggplot2::aes(label=sprintf("SF %.2f", .data$size_factor)), vjust=-0.4, size=3)+
    ggplot2::labs(x=NULL, y="Total counts (millions)", fill=color_by, title=title,
                  subtitle="SF = DESeq2 size factor used for normalization")+
    ggplot2::theme_minimal(base_size=13)+
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45, hjust=1),
                   plot.title=ggplot2::element_text(face="bold", hjust=0.5),
                   plot.subtitle=ggplot2::element_text(hjust=0.5))
  if(!is.null(out_png)) return(save_gg(pl, out_png, width=10, height=6))
  pl
}

#########################Dispersion plot###############################
plot_dispersion<- function(dseq, out_png=NULL){
  #DESeq2's dispersion estimates: gene-wise (black), fitted trend (red), final (blue)
  if(!is.null(out_png)){
    grDevices::png(out_png, width=1400, height=1100, res=150)
    on.exit(grDevices::dev.off())
  }
  DESeq2::plotDispEsts(dseq, main="Dispersion estimates")
  invisible(NULL)
}

########################P-value histogram##############################
plot_pvalue_hist<- function(resdf, out_png=NULL){
  df<- data.frame(pvalue=resdf$pvalue[!is.na(resdf$pvalue)])
  pl<- ggplot2::ggplot(df, ggplot2::aes(x=.data$pvalue))+
    ggplot2::geom_histogram(breaks=seq(0, 1, by=0.025), fill="#2A9D8F", color="white")+
    ggplot2::labs(x="Raw p-value", y="Number of genes", title="P-value distribution",
                  subtitle="Expected shape: flat, with a peak near 0 if there are real differences")+
    ggplot2::theme_minimal(base_size=13)+
    ggplot2::theme(plot.title=ggplot2::element_text(face="bold", hjust=0.5),
                   plot.subtitle=ggplot2::element_text(hjust=0.5))
  if(!is.null(out_png)) return(save_gg(pl, out_png, width=9, height=6))
  pl
}

############################MA plot####################################
plot_ma<- function(resdf, best_genes=NULL, lfc_thr=1, padj_thr=0.05, label_col="GeneName",
                   title="MA plot", out_png=NULL){
  keep<- !is.na(resdf$padj) & !is.na(resdf$log2FoldChange) & resdf$baseMean > 0
  df<- resdf[keep, , drop=FALSE]
  df$status<- "NS"
  df$status[df$padj<= padj_thr & df$log2FoldChange>= lfc_thr]<- "Up"
  df$status[df$padj<= padj_thr & df$log2FoldChange<= -lfc_thr]<- "Down"
  df$status<- factor(df$status, levels=c("Up", "Down", "NS"))
  pl<- ggplot2::ggplot(df, ggplot2::aes(x=log10(.data$baseMean), y=.data$log2FoldChange))+
    ggplot2::geom_point(ggplot2::aes(color=.data$status), size=1.1, alpha=0.6)+
    ggplot2::geom_hline(yintercept=0, color="grey30")+
    ggplot2::geom_hline(yintercept=c(-lfc_thr, lfc_thr), linetype="dashed")+
    ggplot2::scale_color_manual(values=c(Up="#F72585", Down="#3B4CC0", NS="grey60"),
                                name="Differential Expression")+
    ggplot2::labs(x="log10 mean normalized count", y="log2 fold change", title=title)+
    ggplot2::theme_minimal(base_size=14)+
    ggplot2::theme(plot.title=ggplot2::element_text(face="bold", hjust=0.5))
  if(!is.null(best_genes) && requireNamespace("ggrepel", quietly=TRUE)){
    lab<- df[intersect(rownames(best_genes), rownames(df)), , drop=FALSE]
    lab<- lab[!is.na(lab[[label_col]]) & lab[[label_col]] != "", , drop=FALSE]
    lab<- utils::head(lab, 15)
    pl<- pl+ggrepel::geom_text_repel(data=lab, ggplot2::aes(label=.data[[label_col]]),
                                     size=3, max.overlaps=30)
  }
  if(!is.null(out_png)) return(save_gg(pl, out_png))
  pl
}

#####################Counts for a single gene##########################
plot_gene_counts<- function(dseq, gene_id, meta, group_by, gene_label=NULL, out_png=NULL){
  nc<- DESeq2::counts(dseq, normalized=TRUE)
  if(!gene_id %in% rownames(nc)) stop("Gene not found in the analysed data")
  df<- data.frame(sample=colnames(nc), count=as.numeric(nc[gene_id, ]))
  df$group<- if(!is.null(group_by) && group_by %in% names(meta)) as.factor(meta[df$sample, group_by]) else "all"
  lbl<- if(!is.null(gene_label) && nzchar(gene_label)) paste0(gene_label, " (", gene_id, ")") else gene_id
  pl<- ggplot2::ggplot(df, ggplot2::aes(x=.data$group, y=.data$count + 0.5, color=.data$group))+
    ggplot2::geom_boxplot(outlier.shape=NA, width=0.5, color="grey40")+
    ggplot2::geom_jitter(width=0.12, height=0, size=3.5, alpha=0.85)+
    ggplot2::scale_y_log10()+
    ggplot2::labs(x=group_by, y="Normalized count (+0.5, log scale)", title=lbl)+
    ggplot2::theme_minimal(base_size=14)+
    ggplot2::theme(legend.position="none", plot.title=ggplot2::element_text(face="bold", hjust=0.5))
  if(!is.null(out_png)) return(save_gg(pl, out_png, width=7, height=6))
  pl
}



#########################################################################
###########Functions for sample subsets, replicates and shrinkage########
#########################################################################
#keeps only samples whose metadata column 'col' has one of the chosen 'values'
subset_samples<- function(counts, meta, col=NULL, values=NULL){
  if(is.null(col) || !nzchar(col) || length(values)==0) return(list(counts=counts, meta=meta))
  if(!col %in% names(meta)) stop(paste0("Subset column '", col, "' is not in the metadata."))
  keep<- as.character(meta[[col]]) %in% values
  if(sum(keep) < 2) stop("Fewer than 2 samples are left after subsetting. Choose more groups.")
  list(counts= counts[, keep, drop=FALSE], meta= meta[keep, , drop=FALSE])
}

#removes design terms that have a single value left (e.g. CellLine after keeping only PC3)
drop_single_level_terms<- function(design, meta){
  vars<- all.vars(design)
  single<- vars[vapply(vars, function(v) v %in% names(meta) && length(unique(meta[[v]])) < 2, logical(1))]
  if(length(single)==0) return(list(design=design, dropped=character(0)))
  tl<- attr(stats::terms(design), "term.labels")
  keep<- tl[!vapply(tl, function(t) any(all.vars(stats::as.formula(paste("~", t))) %in% single), logical(1))]
  if(length(keep)==0) stop("No design factor has more than one group in the selected samples.")
  list(design= stats::as.formula(paste("~", paste(keep, collapse=" + "))), dropped= single)
}

#returns a warning message when any group of 'col' has fewer than 'min_n' samples
replicate_warning<- function(meta, col, min_n=3){
  if(is.null(col) || !col %in% names(meta)) return(NULL)
  tab<- table(as.character(meta[[col]]))
  low<- tab[tab < min_n]
  if(length(low)==0) return(NULL)
  paste0("Low replication: ", paste0(names(low), " (n=", low, ")", collapse=", "),
         ". With fewer than ", min_n, " samples per group, results are less reliable and small changes can be missed.")
}

#shrinks log2 fold changes with DESeq2's built-in 'normal' method (p-values are unchanged)
shrink_results<- function(deg, res, contrast=NULL, name=NULL){
  sh<- tryCatch(suppressMessages(
    if(!is.null(name)) DESeq2::lfcShrink(deg, coef=name, res=res, type="normal")
    else DESeq2::lfcShrink(deg, contrast=contrast, res=res, type="normal")),
    error=function(e) e)
  if(inherits(sh, "error")) return(list(res=res, shrunk=FALSE, msg=conditionMessage(sh)))
  list(res=sh, shrunk=TRUE, msg=NULL)
}

#########################################################################
###################Functions for pathway analysis (GSEA)#################
#########################################################################
#reads the bundled MSigDB Hallmark gene sets (human gene symbols)
load_gene_sets<- function(path=file.path("gene_sets", "hallmark_human.csv"), uppercase=FALSE){
  if(!file.exists(path)) stop("Gene set file not found: ", path)
  gs<- utils::read.csv(path, stringsAsFactors=FALSE)
  gs<- unique(gs[, c("gs_name", "gene_symbol")])
  sets<- split(gs$gene_symbol, gs$gs_name)
  names(sets)<- sub("^HALLMARK_", "", names(sets))
  sets
}

#ranks genes by the DESeq2 Wald statistic and runs fgsea
run_gsea<- function(resdf, sets, symbol_col="GeneName", stat_col="stat", uppercase=FALSE,
                    min_size=15, max_size=500, seed=1){
  if(!requireNamespace("fgsea", quietly=TRUE)) stop("Package 'fgsea' is needed for pathway analysis.")
  df<- as.data.frame(resdf)
  df<- df[!is.na(df[[stat_col]]) & !is.na(df[[symbol_col]]) & nzchar(df[[symbol_col]]), ]
  sym<- as.character(df[[symbol_col]])
  if(uppercase) sym<- toupper(sym)
  st<- df[[stat_col]]
  ord<- order(abs(st), decreasing=TRUE)    #one value per gene symbol: keep the strongest
  sym<- sym[ord]; st<- st[ord]
  keep<- !duplicated(sym)
  ranks<- stats::setNames(st[keep], sym[keep])
  ranks<- sort(ranks, decreasing=TRUE)
  set.seed(seed)
  res<- suppressWarnings(fgsea::fgsea(pathways=sets, stats=ranks, minSize=min_size, maxSize=max_size))
  res<- res[order(res$pval), ]
  list(table=res, ranks=ranks, sets=sets)
}

#bar plot of normalized enrichment scores
plot_gsea_bar<- function(gtab, padj_thr=0.05, n=20, out_png=NULL, title="Hallmark pathways (GSEA)"){
  df<- as.data.frame(gtab)
  df<- df[!is.na(df$NES), ]
  df<- utils::head(df[order(df$padj), ], n)
  df$Direction<- ifelse(df$NES > 0, "Up in case", "Down in case")
  df$Significant<- ifelse(df$padj < padj_thr, paste0("padj < ", padj_thr), "not significant")
  df$pathway<- factor(df$pathway, levels= df$pathway[order(df$NES)])
  p<- ggplot2::ggplot(df, ggplot2::aes(x=NES, y=pathway, fill=Direction, alpha=Significant))+
    ggplot2::geom_col()+
    ggplot2::scale_fill_manual(values=c("Up in case"="#d6336c", "Down in case"="#3b6fd4"))+
    ggplot2::scale_alpha_manual(values=stats::setNames(c(1, 0.35), c(paste0("padj < ", padj_thr), "not significant")))+
    ggplot2::geom_vline(xintercept=0, colour="grey40")+
    ggplot2::labs(title=title, subtitle=paste0("Top ", nrow(df), " pathways by adjusted p-value"),
                  x="Normalized enrichment score (NES)", y=NULL, fill=NULL, alpha=NULL)+
    ggplot2::theme_bw(base_size=13)
  if(!is.null(out_png)) save_gg(p, out_png, width=10, height=8)
  p
}

#running-score (enrichment) plot for one pathway
plot_gsea_enrichment<- function(gsea, pathway, out_png=NULL){
  p<- fgsea::plotEnrichment(gsea$sets[[pathway]], gsea$ranks)+
    ggplot2::labs(title=pathway, x="Gene rank (most up in case to most down)", y="Enrichment score")+
    ggplot2::theme_bw(base_size=13)
  if(!is.null(out_png)) save_gg(p, out_png, width=9, height=5)
  p
}

#####################################################################
###########################Defining UI##############################
######################################################################
ui<- fluidPage(
  #App title:
  titlePanel("Bulk RNAseq Differential Gene Expression Pipeline"),
  
  tags$head(tags$style(HTML("
                             .mySidebar{
                             font-size: 14px; #font size of the side bar text
                            }
                            "))
  ),
#######################################################################  
#############Sidebar Layouts with Inputs and Output Definitions:#########
##########################################################################
  sidebarLayout(
    sidebarPanel(
      class= "mySidebar",
      
      #data notice shown at the top of the sidebar for visitors of the online version
      div(style="background:#fff8e1; border-left:4px solid #f0ad4e; padding:10px 12px; margin-bottom:15px; font-size:13px;",
          tags$b("Data notice"), tags$br(),
          "The online version of this app runs on Posit Connect Cloud, so files you upload are sent to Posit's servers and processed there during your session. The app does not save your files.",
          tags$br(), tags$br(),
          tags$b("Do not upload confidential, unpublished or patient-identifiable data."),
          " For sensitive data, download the app from ",
          tags$a(href="https://github.com/lizetafo/RNAseq-web", target="_blank", "GitHub"),
          " and run it on your own computer."),
      
      #choose where the data come from
      radioButtons("input_mode", "Data source", inline = TRUE,
                   choices = c("Upload files" = "upload", "Load from GEO" = "geo")),
      conditionalPanel(condition = "input.input_mode == 'upload'",
        fileInput(inputId = "Count", "Upload Counts File(s)", multiple = TRUE,
                  accept = c(".csv", ".txt", ".tsv", ".tab", ".gz", ".tar", ".tgz", ".zip", ".xlsx", ".xls")),
        helpText("One count table (genes x samples), a featureCounts output, several per-sample count files (HTSeq or STAR), or a GEO '_RAW.tar' bundle. CSV, TSV, TXT or Excel, also .gz. Must be raw counts, not FPKM or TPM."),
        fileInput(inputId = "Metadata", "Upload Metadata File (optional)",
                  accept = c(".csv", ".txt", ".tsv")),
        helpText("One row per sample. No sample sheet? Leave this empty and fill one in on the Sample Sheet tab."),
        checkboxInput(inputId = "NA_to_zero", "Convert NA in counts to zero", value = TRUE),
        checkboxInput(inputId = "round_estimated", "Round estimated counts (RSEM / Salmon)", value = FALSE),
        actionButton("upload_data", "Upload", class="btn-info"), #using action button to upload  data
        helpText("Upload data and automatically detect experimental design")
      ),
      conditionalPanel(condition = "input.input_mode == 'geo'",
        textInput("geo_acc", "GEO series accession", placeholder = "e.g. GSE106305"),
        actionButton("load_geo", "Load from GEO", class="btn-info"),
        helpText("Uses NCBI's raw counts for the series when available (many human and mouse studies), otherwise the authors' own count file, and builds the sample sheet from the GEO sample information. Check the Sample Sheet tab and the organism / ID type afterwards.")
      ),
      br(),
      #lets visitors try the app without their own data
      actionButton("load_example", "Or load example data (prostate hypoxia)", class="btn-default btn-sm"),
      
      tags$hr(),
      
      conditionalPanel( #this will show additonal sidebar panel inputs only after user clicks on upload 
        condition = "input.upload_data>0 || input.load_example>0 || input.load_geo>0",
      
      #organism database selection
      selectInput(inputId = "Orgdb", 
                  "Organism Database", 
                  choices = c("Human (org.Hs.eg.db)"= "human", 
                              "Mouse (org.Mm.eg.db)"="mouse",
                              "Rat (org.Rn.eg.db)" ="rat",
                              "Fly (org.Dm.eg.db)"= "fly",
                              "Zebrafish (org.Dr.eg.db)"="zfish",
                              "Yeast (org.Sc.sgd.db)"= "yeast",
                              "Other (not listed)"= "other")
      ),
      #Custom Organism database, incase user has other's that are not listed above
      conditionalPanel(
        condition = "input.Orgdb=='other'",
        textInput("OtherOrg", "Enter organism name or annotation package:")
      ),
      helpText("Choose appropriate organism reference database for gene annotation."),
      
      #gene ID type selction
      selectInput(inputId = "idtype", 
                  "ID Type for Gene Annotation", 
                  choices = c("ENSEMBL" ="ensembl",
                              "ENTREZID"= "entrezid",
                              "SYMBOL"= "symbol",
                              "Other (not listed)"= "other")
      ),
      #similarly, custom ID type can be entered if not provided in the select options
      conditionalPanel(
        condition = "input.idtype=='other'",
        textInput("OtherIDtype", "Enter custom ID type (e.g., REFSEQ, UNIPROT):")
      ),
      helpText("Choose appropriate ID Type for gene names annotation (example: ENSEMBL)."),
      
      #optional: analyze only some samples (e.g. one cell line)
      selectInput("subset_col", "Analyze a subset of samples (optional)", choices = c("All samples" = "")),
      conditionalPanel(
        condition = "input.subset_col != ''",
        selectizeInput("subset_vals", "Keep samples where this column is:", choices = NULL, multiple = TRUE)
      ),
      helpText("Example: choose CellLine, then PC3, to analyze one cell line on its own."),
      
      #text input to enter experiment design types and conditions
      textInput(inputId = "Exp_groups", "Experimental Grouping Factor(s)"),
      helpText("Enter experimental grouping factor(s) (comma-seperated) (example: Activation, Infection)"),
      
      textInput(inputId = "ref_col", "Reference Column"),
      helpText("Enter reference column name from metadata (example: Activation)"),
      selectInput(inputId = "ref_levels", "Reference level (control group)", choices = NULL),
      helpText("Choose the control group that fold changes are compared against (example: Resting)"),
      
      textInput(inputId = "Exp_design", "Experimental Design Formula"),
      helpText("Enter experiment design formula to be used for DESeq2 (example ~Activation)"),
      
      #text input to select for contrast
      textInput(inputId = "contrast", "Contrast (Factor, Case, Control)"),
      helpText("Enter contrast (from metadata)( example: Activation, Activated, Resting)"),
      
      #user can select thresholds
      numericInput("alpha", "Adjusted p-value Threshold (alpha)", value = 0.05, 
                   min = 0.0001, max= 0.5, step = 0.01),
      numericInput("lfc_thr", "log2FC Threshold", value=1, step = 0.1),
      numericInput("padj_thr_plot", "Adjusted p-value Threshold", value=0.05, min=1e-6, max=0.5),
      checkboxInput("shrink", "Shrink log2 fold changes (DESeq2 'normal' method)", value = FALSE),
      helpText("Shrinkage makes fold changes of low-count genes more reliable. P-values do not change."),
      
      #Options for heatmap transformation types
      selectInput("heatmap_transform", "Transformation for Heatmap",
                  choices = c("Normalized counts"="norm",
                              "VST"="vst",
                              "rlog"="rlog"),
                  selected = "norm"),
      #User can determine the number of tope genes they want displayed in their volcano plot and heatmap
      numericInput("n", "Number of Top genes",
                   value = 30, min=10, max = 200, step=5),
      helpText("Enter number of top genes to use for plots"),
      
      #action button to run DESeq2 pipeline
      actionButton("run", "Run", class= "btn-primary"),
      helpText("Run RNA DESeq Pipeline"),
      
      #advanced: show any DESeq2 coefficient (e.g. an interaction term) without re-running
      conditionalPanel(condition = "input.run>0",
        tags$hr(),
        selectInput("result_choice", "Result to show (advanced)", choices = NULL),
        helpText("Default: the contrast above. You can also pick any DESeq2 coefficient, for example an interaction term such as CellLinePC3.ConditionHypoxia (design ~CellLine + Condition + CellLine:Condition).")
      )
      
    )
    ),
########################################################################################
###############################Main Panel##########################################
##################################################################################
    mainPanel(
      uiOutput("run_notes"),
      tabsetPanel(id = "main_tabs",
        #Sample sheet tab: view and edit the metadata
        tabPanel("Sample Sheet",
                 br(),
                 h4("Sample sheet (metadata)"),
                 helpText("One row per sample. Double-click a cell to edit it. Add a column for each grouping (for example Condition), fill in each sample's group, then click 'Use this sample sheet'. Group names are case-sensitive."),
                 conditionalPanel(condition = "input.upload_data>0 || input.load_example>0 || input.load_geo>0",
                   fluidRow(
                     column(4, textInput("new_col", NULL, placeholder = "New column name, e.g. Condition")),
                     column(8, actionButton("add_col", "Add column"),
                            actionButton("apply_meta", "Use this sample sheet", class = "btn-primary"),
                            downloadButton("download_meta", "Download sample sheet")))),
                 DT::DTOutput("meta_table")
        ),
        #Quality-control tab
        tabPanel("Quality Control",
                 br(),
                 h4("Sample and model quality checks"),
                 conditionalPanel(condition="input.run>0",
                   selectInput("qc_color", "Color / group samples by", choices=NULL)),
                 tabsetPanel(
                   tabPanel("PCA",
                            br(),
                            numericInput("pca_ntop", "Number of most variable genes", value=500, min=50, max=5000, step=50),
                            plotOutput("pca_plot", height="650px"),
                            downloadButton("download_pca", "Download PCA Plot"),
                            helpText("Replicates of the same group should cluster together. Outlier samples stand apart.")),
                   tabPanel("Sample Distances",
                            br(),
                            plotOutput("dist_plot", height="700px"),
                            downloadButton("download_dist", "Download Sample Distance Heatmap"),
                            helpText("Darker = more similar samples. Samples should group by condition.")),
                   tabPanel("Library Size",
                            br(),
                            plotOutput("libsize_plot", height="600px"),
                            downloadButton("download_libsize", "Download Library Size Plot"),
                            helpText("Total reads counted per sample, with the DESeq2 size factor used for normalization.")),
                   tabPanel("Dispersion",
                            br(),
                            plotOutput("disp_plot", height="650px"),
                            downloadButton("download_disp", "Download Dispersion Plot"),
                            helpText("Gene-wise (black), fitted (red) and final (blue) dispersion estimates. Final values should follow the red trend.")),
                   tabPanel("P-value Histogram",
                            br(),
                            plotOutput("pval_plot", height="550px"),
                            downloadButton("download_pval", "Download P-value Histogram"),
                            helpText("A flat histogram with a peak near 0 is expected when there are real differences."))
                 )
        ),
        #Graphical plots tab
        tabPanel("Graphical Plots",
                 br(),
                 h4("Differential Gene Expression Plots"),
                 tabsetPanel(
                   #tab for volcano plot
                   tabPanel(
                     "Volcano Plot",
                     br(),
                     fluidRow(
                       column(
                         width = 8,
                         plotOutput("volcano_plot", height = "1000px"),
                         downloadButton("download_volcano", "Download Volcano Plot")
                       ),
                       column(
                         width = 4,
                         conditionalPanel(
                           condition="input.run>0",
                           h4("Top 15 genes (by padj)"),
                           DT::DTOutput("volcano_top15")
                         )
                       )
                     )
                   ),
                   #tab for heatmap
                   tabPanel("Heatmap Plot",
                            br(),
                            plotOutput("heatmap", height = "1000px"),
                            downloadButton("download_heatmap", "Download Heatmap")),
                   #tab for bar plot
                   tabPanel("Bar Plot",
                            br(),
                            plotOutput("bar_plot", height = "1000px"),
                            downloadButton("download_bar", "Download Bar Plot")
                   ),
                   #tab for MA plot
                   tabPanel("MA Plot",
                            br(),
                            plotOutput("ma_plot", height = "750px"),
                            downloadButton("download_ma", "Download MA Plot"),
                            helpText("Fold change against average expression for every gene."))
                 )
        ),

        #Gene explorer tab
        tabPanel("Gene Explorer",
                 br(),
                 h4("Expression of a single gene"),
                 selectizeInput("gene_pick", "Choose a gene (type a gene name or ID)", choices=NULL,
                                options=list(placeholder="Run the analysis first")),
                 plotOutput("gene_plot", height="550px"),
                 downloadButton("download_gene", "Download Gene Plot"),
                 br(), br(),
                 tableOutput("gene_stats")
        ),
        
        #Pathway analysis tab
        tabPanel("Pathways (GSEA)",
                 br(),
                 h4("Gene set enrichment analysis: MSigDB Hallmark pathways"),
                 helpText("All genes are ranked by the DESeq2 test statistic (most up in the case group first) and tested against the 50 Hallmark gene sets with fgsea. Positive NES = pathway higher in the case group."),
                 uiOutput("gsea_note"),
                 plotOutput("gsea_bar", height="650px"),
                 downloadButton("download_gsea_bar", "Download Pathway Plot"),
                 br(), br(),
                 selectInput("gsea_pick", "Enrichment plot for pathway", choices=NULL, width="400px"),
                 plotOutput("gsea_enrich", height="400px"),
                 downloadButton("download_gsea_enrich", "Download Enrichment Plot"),
                 br(), br(),
                 DT::DTOutput("gsea_table"),
                 downloadButton("download_gsea_table", "Download Pathway Results"),
                 helpText("Gene sets: MSigDB Hallmark collection (Liberzon et al., Cell Systems 2015), CC BY 4.0.")
        ),
        
        #Resutls table tab
        tabPanel("Results",
                 br(),
                 h4("Differential Gene Expression Tabular Results"),
                 tabsetPanel(
                   tabPanel("Annotated",
                            br(),
                            DT::DTOutput("results_ann")),
                   tabPanel("Raw DESeq2 Results",
                            br(),
                            DT::DTOutput("results_raw")),
                   tabPanel("Normalized Counts",
                            br(),
                            DT::DTOutput("results_norm"))
                 ),
                 br(),
                 downloadButton("download_ann", "Download Annotated Results"),
                 downloadButton("download_counts", "Download DESeq2 Results"),
                 downloadButton("download_norm", "Download Normalized Counts")
        )
      )
      
    )
    
    
  )
  
)



########################################################################
###########################Defining Server##############################
########################################################################

server<- function(input, output, session){
  loadlib()  #load required libraries
  
  #makes a reactive to sotre uploaded files
  data_store<- reactiveValues(
    dat=NULL
  )
  
  #Ensure that files are uploaded before user can click upload button
  #stores the data and pre-fills the design settings (used by Upload and by the example button)
  #sample sheet being edited in the Sample Sheet tab (meta_ver changes force a redraw)
  meta_draft<- reactiveVal(NULL)
  meta_ver<- reactiveVal(0)
  set_draft<- function(m){ meta_draft(m); meta_ver(isolate(meta_ver()) + 1) }
  
  #turns an analysis error into a readable message in the app
  #(the hosting server hides ordinary error messages as "An error has occurred")
  show_errors<- function(expr){
    tryCatch(expr, error = function(e){
      if(inherits(e, "shiny.silent.error")) stop(e)   #validation messages pass through unchanged
      validate(need(FALSE, paste("Problem:", conditionMessage(e))))
    })
  }
  
  #value of a reactive, or NULL if it failed (so observers never stop the app on an analysis error)
  safe_value<- function(expr) tryCatch(expr, error = function(e) NULL)
  
  #shows an error as a notification instead of closing the app
  report_error<- function(e) showNotification(conditionMessage(e), type = "error", duration = NULL)
  
  #stores the data and pre-fills the design settings (used by Upload, GEO, the example button and the sample sheet)
  setup_design<- function(dat, update_draft = TRUE){
    #uploaded data are stored in the reactive container initially created to be unsed in pipeline
    data_store$dat<- dat
    if(update_draft) set_draft(dat$metadata)
    if(length(dat$notes) > 0) showNotification(paste(dat$notes, collapse = " "), type = "message", duration = 15)
    
    #Autodetect experimental desing from metadata (User has option to change these)
    #a grouping column has 2-10 groups, no empty values, and replicates (not one value per sample)
     meta<- dat$metadata
      factor_cols<- names(meta)
      #columns with typical names for the comparison are tried first (e.g. treatment before cell_line)
      pref<- grepl("condition|treatment|group|genotype|disease|status|protocol|agent|time|stimul|infect", factor_cols, ignore.case = TRUE)
      factor_cols<- c(factor_cols[pref], factor_cols[!pref])
      
      candidate<-NULL
      for (col in factor_cols){
        v<- as.character(meta[[col]])
        u<- unique(v)
        if (length(u)>= 2 && length(u)<=10 && length(u) < length(v) && all(nzchar(v) & !is.na(v))){
          candidate <- col
          break
        }
      }
      
      updateSelectInput(session, "subset_col", choices = c("All samples" = "", names(meta)), selected = "")
      updateSelectizeInput(session, "subset_vals", choices = character(0), selected = character(0))
      
      if (is.null(candidate)){
        #no usable grouping yet: send the user to the sample sheet
        for(id in c("Exp_groups", "ref_col", "Exp_design", "contrast")) updateTextInput(session, id, value = "")
        updateSelectInput(session, "ref_levels", choices = c("-- choose control group --" = ""), selected = "")
        updateTabsetPanel(session, "main_tabs", selected = "Sample Sheet")
        showNotification("Data loaded. Now fill in the sample sheet: add a column (for example Condition), type each sample's group, then click 'Use this sample sheet'.",
                         type = "warning", duration = 20)
        return(invisible(NULL))
      }
      f<- factor(meta[[candidate]])
      lv<- levels(f)
      
      #Prefils the UI fiels with what was guessed, may not be want the user wants so user can change to desired settings
      updateTextInput(session, "Exp_groups", value= candidate)
      updateTextInput(session, "ref_col", value= candidate)
      #the control group is NOT guessed: the user must choose it
      updateSelectInput(session, "ref_levels",
                        choices = c("-- choose control group --" = "", lv), selected = "")
      updateTextInput(session, "Exp_design", value= paste0("~", candidate))
      updateTextInput(session, "contrast", value= "")
      
      showNotification("Data loaded. Check the sample sheet, choose the control (reference) group and the settings, then click Run.",
                       type="message", duration = 10)
  }

  observeEvent(input$upload_data,{
    if(is.null(input$Count)){
      showNotification("Please choose a counts file first.", type = "warning")
      return(NULL)
    }
    tryCatch({
      #reads counts (one or several files) and the optional metadata
      dat<- read_file(
        cnts= input$Count$datapath,
        met= if(is.null(input$Metadata)) NULL else input$Metadata$datapath,
        na_to_zero = input$NA_to_zero,
        cnts_names = input$Count$name,
        round_estimated = isTRUE(input$round_estimated))
      
      #basic quality checks
      qc_check(dat$counts, dat$metadata)
      setup_design(dat)
    }, error = report_error)
  })

  #loads the bundled example dataset (prostate cancer cells, hypoxia vs normoxia; GEO GSE106305)
  observeEvent(input$load_example,{
    ex_counts<- file.path("example_data", "prostate_hypoxia", "counts.csv")
    ex_meta<- file.path("example_data", "prostate_hypoxia", "metadata.csv")
    if(!file.exists(ex_counts) || !file.exists(ex_meta)){
      showNotification("Example data not found next to app.R.", type="error")
      return(NULL)
    }
    tryCatch({
      dat<- read_file(cnts= ex_counts, met= ex_meta, na_to_zero= TRUE)
      qc_check(dat$counts, dat$metadata)
      updateSelectInput(session, "Orgdb", selected = "human")
      updateSelectInput(session, "idtype", selected = "ensembl")
      setup_design(dat)
    }, error = report_error)
  })
  
  #downloads a dataset from GEO (NCBI raw counts + sample information)
  observeEvent(input$load_geo,{
    tryCatch({
      g<- withProgress(message = "Downloading from GEO (can take a minute)...", value = 0.3, {
        load_geo(input$geo_acc)
      })
      if(!is.na(g$organism)) updateSelectInput(session, "Orgdb", selected = g$organism)
      updateSelectInput(session, "idtype", selected = g$id_type)
      setup_design(list(counts = g$counts, metadata = g$metadata, notes = g$notes))
      updateTabsetPanel(session, "main_tabs", selected = "Sample Sheet")
    }, error = report_error)
  })
  
  #### SAMPLE SHEET EDITOR
  output$meta_table<- DT::renderDT({
    meta_ver()
    m<- isolate(meta_draft())
    validate(need(!is.null(m), "Upload data, load a GEO series, or load the example data to see the sample sheet."))
    DT::datatable(cbind(Sample = rownames(m), m), rownames = FALSE, selection = "none",
                  editable = list(target = "cell", disable = list(columns = 0)),
                  options = list(pageLength = 50, scrollX = TRUE, dom = "tip"))
  })
  
  #saves edited cells
  observeEvent(input$meta_table_cell_edit, {
    info<- input$meta_table_cell_edit
    m<- meta_draft(); req(m)
    for(k in seq_len(nrow(info))){
      i<- info$row[k]; j<- info$col[k]   #column 0 is the sample name
      if(j >= 1 && j <= ncol(m)) m[i, j]<- trimws(as.character(info$value[k]))
    }
    meta_draft(m)
  })
  
  #adds an empty column to fill in
  observeEvent(input$add_col, {
    m<- meta_draft(); req(m)
    nm<- gsub("^_|_$", "", gsub("[^A-Za-z0-9]+", "_", trimws(input$new_col)))
    if(!nzchar(nm)){ showNotification("Type a column name first, for example Condition.", type = "warning"); return(NULL) }
    if(grepl("^[0-9]", nm)) nm<- paste0("X", nm)
    if(nm %in% names(m)){ showNotification(paste0("Column '", nm, "' already exists."), type = "warning"); return(NULL) }
    m[[nm]]<- ""
    set_draft(m)
    updateTextInput(session, "new_col", value = "")
  })
  
  #uses the edited sample sheet for the analysis
  observeEvent(input$apply_meta, {
    m<- meta_draft(); req(m, data_store$dat)
    blank<- names(m)[vapply(m, function(v) any(!nzchar(as.character(v)) | is.na(v)) && any(nzchar(as.character(v))), logical(1))]
    if(length(blank) > 0)
      showNotification(paste0("Some samples have no value in: ", paste(blank, collapse = ", "), ". Fill them in, or that column cannot be used as a group."),
                       type = "warning", duration = 15)
    dat<- data_store$dat
    dat$metadata<- m
    dat$notes<- NULL
    tryCatch({
      qc_check(dat$counts, dat$metadata)
      setup_design(dat, update_draft = FALSE)
    }, error = report_error)
  })
  
  output$download_meta<- downloadHandler(
    filename = function() "sample_sheet.csv",
    content = function(file){
      m<- meta_draft()
      utils::write.csv(cbind(ID = rownames(m), m), file, row.names = FALSE)
    }
  )

  #when a subset column is chosen, offer its values
  observeEvent(input$subset_col, {
    req(data_store$dat)
    meta<- data_store$dat$metadata
    if(nzchar(input$subset_col) && input$subset_col %in% names(meta)){
      v<- sort(unique(as.character(meta[[input$subset_col]])))
      updateSelectizeInput(session, "subset_vals", choices = v, selected = v)
    } else {
      updateSelectizeInput(session, "subset_vals", choices = character(0), selected = character(0))
    }
  }, ignoreInit = TRUE)

  #when the reference column changes, offer its groups as control-group choices
  #metadata of the samples that will be analyzed (after the optional subset)
  meta_in_use<- function(){
    meta<- data_store$dat$metadata
    sc<- isolate(input$subset_col); sv<- isolate(input$subset_vals)
    if(!is.null(sc) && nzchar(sc) && sc %in% names(meta) && length(sv) > 0)
      meta<- meta[as.character(meta[[sc]]) %in% sv, , drop = FALSE]
    meta
  }
  
  #fills the contrast from the control group: case = another group present in the analyzed samples
  fill_contrast<- function(){
    col<- trimws(isolate(input$ref_col)); ref<- isolate(input$ref_levels)
    meta<- meta_in_use()
    if(!nzchar(ref) || !col %in% names(meta)) return(invisible(NULL))
    others<- setdiff(sort(unique(as.character(meta[[col]]))), ref)
    updateTextInput(session, "contrast",
                    value = if(length(others) >= 1) paste(col, others[1], ref, sep=",") else "")
  }
  
  #offers the groups of the reference column (only those in the analyzed samples) as control-group choices
  refresh_ref_choices<- function(){
    req(data_store$dat)
    col<- trimws(isolate(input$ref_col))
    meta<- meta_in_use()
    if(col %in% names(meta)){
      lv<- sort(unique(as.character(meta[[col]])))
      keep<- if(isTRUE(isolate(input$ref_levels) %in% lv)) isolate(input$ref_levels) else ""
      updateSelectInput(session, "ref_levels",
                        choices = c("-- choose control group --" = "", lv), selected = keep)
      if(nzchar(keep)) fill_contrast()
    }
  }
  
  #when the reference column or the sample subset changes, update the control-group choices and the contrast
  observeEvent(input$ref_col, refresh_ref_choices(), ignoreInit = TRUE)
  observeEvent(input$subset_vals, refresh_ref_choices(), ignoreInit = TRUE, ignoreNULL = FALSE)

  #when a control group is chosen, fill in the contrast (case vs control)
  observeEvent(input$ref_levels, {
    req(data_store$dat, nzchar(input$ref_levels))
    fill_contrast()
  }, ignoreInit = TRUE)
   
  
  ###building ref levels from UI
  ref_levels<- reactive({
    ref_list<- list()
    if(nzchar(input$ref_col)&&nzchar(input$ref_levels)){
      ref_list[[input$ref_col]] <- input$ref_levels
    }
    ref_list
  })
  
  ###chooisng organism data base
  orgdb <- reactive({
    if(input$Orgdb != "other"){
      pkg<- switch(input$Orgdb,
                   human = "org.Hs.eg.db",
                   mouse = "org.Mm.eg.db",
                   rat = "org.Rn.eg.db",
                   fly = "org.Dm.eg.db",
                   zfish = "org.Dr.eg.db",
                   yeast = "org.Sc.sgd.db"
      )
      #installs the database the first time it is used (one-time download)
      if(!requireNamespace(pkg, quietly = TRUE))
        showNotification(paste0("Installing ", pkg, " (one-time download, may take a few minutes)..."),
                         type = "warning", duration = NULL, id = "orgdb_install")
      db<- load_orgdb(pkg)
      removeNotification("orgdb_install")
      return(db)
    }
    
    #if the user selected "other", the user is required to provide the organism annotation package name
    validate(need(nzchar(input$OtherOrg),
                  "Please enter the organism annotation package name.")) 
    pkg_name<- input$OtherOrg
    
    #informs user if the package is not installed
    if(!requireNamespace(pkg_name, quietly=TRUE)){
      validate(
        need(FALSE,
             paste0("Package '", pkg_name,
                    "'is not installed. Install it with BiocManager::install('", pkg_name, "')"))
      )
    }
    
    #loads package and return the annotation object
    ns<- loadNamespace(pkg_name)
    get(pkg_name, envir = ns)
  })
  
  #choose gene ID type
  id_type <- reactive({
    if(input$idtype != "other"){
      id<- switch(input$idtype,
                  ensembl = "ENSEMBL",
                  entrezid= "ENTREZID",
                  symbol= "SYMBOL"
      )
      return(id)
    }
    
    #similarly, user need to enter the custom ID type if the select other
    validate(need(nzchar(input$OtherIDtype),
                  "Please enter custom ID type.")) 
    toupper(trimws(input$OtherIDtype))
  })
  
  
  ##reactive to run RNA pipeline (DESeq2 model; settings are read when Run is clicked)
  core<- eventReactive(input$run, { show_errors({
    validate(   #ensures count and metadata have been uploaded
      need(!is.null(data_store$dat),
           "Please upload data and detect design first (click 'Upload')"),
      need(nzchar(input$ref_levels),
           "Please choose the reference (control) group, then click Run.")
    )
    
    #Retrieve stored data
    dat<- data_store$dat
    notes<- as.character(dat$notes)
    
    #optional: keep only some samples (e.g. one cell line)
    sub<- subset_samples(dat$counts, dat$metadata, col= input$subset_col, values= input$subset_vals)
    if(ncol(sub$counts) < ncol(dat$counts))
      notes<- c(notes, paste0("Analyzing ", ncol(sub$counts), " of ", ncol(dat$counts), " samples (",
                              input$subset_col, ": ", paste(input$subset_vals, collapse=", "), ")."))
    dat_use<- list(counts= sub$counts, metadata= sub$meta)
    
    #Prepare experiemtnal groping vector 
    exp_groups<- strsplit(input$Exp_groups, ",")[[1]]
    exp_groups<- trimws(exp_groups)
    exp_groups<- exp_groups[nzchar(exp_groups)]
    
    #the control group must still be present after subsetting
    validate(need(input$ref_levels %in% as.character(dat_use$metadata[[trimws(input$ref_col)]]),
                  "The chosen control group has no samples in the selected subset."))
    
    #prepare metadata
    meta<-prepare_meta(dat_use$metadata, Exp_groups = exp_groups, ref_levels = ref_levels())
    
    #give error if desing formula is invalid
    Design_formula<- tryCatch(as.formula(input$Exp_design),
                              error=function(e){
                                stop("Design formula is not valid. Example: ~Activation or ~Activation + Infection")
                              })
    #factors with a single value left after subsetting are removed from the design
    dd<- drop_single_level_terms(Design_formula, meta)
    if(length(dd$dropped) > 0)
      notes<- c(notes, paste0("Removed ", paste(dd$dropped, collapse=", "),
                              " from the design because only one value is left in the selected samples. Design used: ",
                              paste(deparse(dd$design), collapse=""), "."))
    
    #warn when groups have few replicates
    rw<- replicate_warning(meta, trimws(input$ref_col))
    if(!is.null(rw)) notes<- c(notes, rw)
    
    #Building Contrast Vector, gives error if contrast length is equal to 3
    contrast_vec<- strsplit(input$contrast, ",")[[1]]
    contrast_vec<- trimws(contrast_vec)
    #empty contrast: compare the other group with the control group (only when there are exactly 2 groups)
    rc<- trimws(input$ref_col)
    if(!nzchar(trimws(input$contrast)) && rc %in% names(meta)){
      others<- setdiff(levels(factor(meta[[rc]])), input$ref_levels)
      if(length(others) == 1){
        contrast_vec<- c(rc, others, input$ref_levels)
        notes<- c(notes, paste0("Contrast filled in automatically: ", paste(contrast_vec, collapse = ", "), "."))
      }
    }
    validate(need(length(contrast_vec) == 3,
                  "Invalid contrast. Enter three items: factor, case, control (for example Condition, Hypoxia, Normoxia)."))
    validate(need(contrast_vec[1] %in% names(meta),
                  paste0("The contrast column '", contrast_vec[1], "' is not in the sample sheet.")))
    missing_lv<- setdiff(contrast_vec[2:3], unique(as.character(meta[[contrast_vec[1]]])))
    validate(need(length(missing_lv) == 0,
                  paste0("The contrast group(s) ", paste0("'", missing_lv, "'", collapse = ", "),
                         " have no samples in the analysis. Check the contrast and the sample subset.")))
    
    #Run DESeq2 pipeline
    dseq<- run_deseq(dat_use$counts, meta, Exp_design = dd$design)
    
    #variance-stabilised data for PCA, sample distances and heatmaps
    #(rlog is slow, so it is only computed when chosen for the heatmap; see rld_obj)
    vsd<- safe_vst(dseq)
    
    
    list(dat=dat_use, meta=meta, dseq=dseq, vsd=vsd, notes=notes,
         contrast=contrast_vec, alpha=input$alpha, shrink=isTRUE(input$shrink), n=input$n,
         db=orgdb(), id=id_type(), org=input$Orgdb)
  })})
  
  #after a run: list every DESeq2 coefficient so advanced users can pick one (e.g. an interaction)
  observeEvent(safe_value(core()), {
    rn<- setdiff(DESeq2::resultsNames(core()$dseq), "Intercept")
    lab<- paste(core()$contrast, collapse=", ")
    updateSelectInput(session, "result_choice",
                      choices= c(stats::setNames("contrast", paste0("Contrast: ", lab)), stats::setNames(rn, rn)),
                      selected= "contrast")
  })
  
  ##results for the chosen comparison (changing 'Result to show' does not re-run DESeq2)
  pipeline<- reactive({ show_errors({
    cr<- core()
    choice<- input$result_choice
    use_name<- !is.null(choice) && nzchar(choice) && choice != "contrast" && choice %in% DESeq2::resultsNames(cr$dseq)
    notes<- cr$notes
    
    #extract DESeq2 results
    res_list<- extract_results(cr$dseq, contrast = if(use_name) choice else cr$contrast,
                               alpha = cr$alpha, normalized = TRUE)
    if(use_name) notes<- c(notes, paste0("Showing DESeq2 result '", choice, "' instead of the sidebar contrast."))
    
    #optional log2 fold-change shrinkage (p-values stay the same)
    if(cr$shrink){
      sh<- shrink_results(cr$dseq, res_list$count_results,
                          contrast = if(use_name) NULL else cr$contrast, name = if(use_name) choice else NULL)
      res_list$count_results<- sh$res
      notes<- c(notes, if(sh$shrunk) "Log2 fold changes are shrunk (DESeq2 'normal' method); p-values are unchanged."
                       else paste0("Fold-change shrinkage was not applied: ", sh$msg))
    }
    
    ##annoate gene names
    ann<- annotate_results(res_list$count_results, orgdb = cr$db, id_type = cr$id, to_symbol = "SYMBOL")
    
    ##top genese
    top__genes<- top_genes(ann, padj= "padj", n= cr$n)
    
    list(dat=cr$dat,
         meta= cr$meta,
         dseq=cr$dseq,
         vsd= cr$vsd,
         res_list= res_list,
         ann= ann,
         top__genes= top__genes,
         notes= notes,
         org= cr$org)
  })})
  
  #notes shown above the tabs (subset used, low replication, shrinkage, ...)
  output$run_notes<- renderUI({
    req(input$run > 0)
    nt<- pipeline()$notes
    if(length(nt)==0) return(NULL)
    div(style="background:#eef4fb; border-left:4px solid #3b6fd4; padding:8px 12px; margin:10px 0; font-size:13px;",
        tags$b("Notes on this analysis"),
        tags$ul(style="margin:4px 0 0 0; padding-left:18px;", lapply(nt, tags$li)))
  })

  #rlog-transformed data, computed only when the heatmap uses it
  rld_obj<- reactive({
    p<- core()
    withProgress(message="Computing rlog transformation (can take a while)...", {
      tryCatch(DESeq2::rlog(p$dseq, blind = TRUE), error=function(e) NULL)
    })
  })

  #matrix used by the heatmap
  heatmap_matrix<- function(p){
    if(input$heatmap_transform== "vst" && !is.null(p$vsd)) return(SummarizedExperiment::assay(p$vsd))
    if(input$heatmap_transform== "rlog"){
      rld<- rld_obj()
      if(!is.null(rld)) return(SummarizedExperiment::assay(rld))
    }
    p$res_list$normalized_counts
  }

  #after a run: offer metadata columns for coloring and genes for the explorer
  observeEvent(safe_value(pipeline()), {
    p<- pipeline()
    cols<- names(p$meta)
    first<- trimws(strsplit(input$Exp_groups, ",")[[1]])[1]
    updateSelectInput(session, "qc_color", choices = cols,
                      selected = if(isTRUE(first %in% cols)) first else cols[1])
    ids<- rownames(p$ann)
    labs<- ifelse(is.na(p$ann$GeneName) | p$ann$GeneName == "", ids, paste0(p$ann$GeneName, " (", ids, ")"))
    best<- rownames(p$top__genes)[1]
    updateSelectizeInput(session, "gene_pick", choices = stats::setNames(ids, labs),
                         selected = best, server = TRUE)
  })

  qc_vsd_matrix<- function(p){
    validate(need(!is.null(p$vsd), "Variance-stabilised data could not be computed for this dataset."))
    SummarizedExperiment::assay(p$vsd)
  }
  
  #### PLOTS OUTPUTS
  ###Rendering plots and including download options
  output$volcano_plot<- renderPlot({
    p<- pipeline()
    plot_volcano(
      p$ann, p$top__genes, out_png = NULL,
      lfc_thr = input$lfc_thr,
      padj_thr = input$padj_thr_plot
    )
  })  
  
  output$volcano_top15<- DT::renderDT({
    req(input$run>0)
    p<- pipeline()
    df<- p$ann
    keep<- !is.na(df$padj)& !is.na(df$log2FoldChange) &
      df$padj<= input$padj_thr_plot &
      abs(df$log2FoldChange)>= input$lfc_thr
    df<- df[keep, , drop=FALSE]
    top15<- top_genes(df, padj="padj", n=15 )
    top15[, c("GeneName", "log2FoldChange", "padj"), drop=FALSE]
    }, options= list(pageLength=15, searching=FALSE, ordering=TRUE))
  

  
  output$download_volcano<- downloadHandler(
    filename = function() "Volcano_plot.png",
    content = function(file){
      p<- pipeline()
      plot_volcano(p$ann, p$top__genes, out_png=file, 
                   lfc_thr =input$lfc_thr, padj_thr = input$padj_thr_plot )
    }
  )
  output$heatmap <- renderPlot({
    p<- pipeline()
    plot_Heatmap(p$dseq,p$ann,out_png = NULL, n=input$n, main = paste("Top DE Genes (", input$heatmap_transform,")", sep=""),
                 mat=heatmap_matrix(p))
  })
  
  
  output$download_heatmap<- downloadHandler(
    filename = function() "Heatmap.png",
    content = function(file){
      p<- pipeline()
      plot_Heatmap(p$dseq, p$ann, out_png=file, n=input$n, 
                   main = paste("Top DE Genes (", input$heatmap_transform,")", sep=""),
                   mat=heatmap_matrix(p))
    }
  )
  output$bar_plot <- renderPlot({
    p<- pipeline()
    plot_bar(p$ann, n=input$n, direction="all",out_png = NULL)
  })
  output$download_bar<- downloadHandler(
    filename = function() "Barplot.png",
    content = function(file){
      p<- pipeline()
      plot_bar(p$ann, direction="all", out_png=file, n=input$n)
    }
  )
  
  
  #####Quality-control plots#########
  output$pca_plot<- renderPlot({
    p<- pipeline(); req(input$qc_color)
    plot_pca(qc_vsd_matrix(p), p$meta, color_by=input$qc_color, ntop=input$pca_ntop)
  })
  output$download_pca<- downloadHandler(
    filename = function() "PCA_plot.png",
    content = function(file){
      p<- pipeline()
      plot_pca(qc_vsd_matrix(p), p$meta, color_by=input$qc_color, ntop=input$pca_ntop, out_png=file)
    }
  )
  output$dist_plot<- renderPlot({
    p<- pipeline(); req(input$qc_color)
    plot_sample_distances(qc_vsd_matrix(p), p$meta, annotate_by=input$qc_color)
  })
  output$download_dist<- downloadHandler(
    filename = function() "Sample_distance_heatmap.png",
    content = function(file){
      p<- pipeline()
      plot_sample_distances(qc_vsd_matrix(p), p$meta, annotate_by=input$qc_color, out_png=file)
    }
  )
  output$libsize_plot<- renderPlot({
    p<- pipeline(); req(input$qc_color)
    plot_library_size(p$dseq, p$meta, color_by=input$qc_color)
  })
  output$download_libsize<- downloadHandler(
    filename = function() "Library_size.png",
    content = function(file){
      p<- pipeline()
      plot_library_size(p$dseq, p$meta, color_by=input$qc_color, out_png=file)
    }
  )
  output$disp_plot<- renderPlot({
    p<- pipeline()
    plot_dispersion(p$dseq)
  })
  output$download_disp<- downloadHandler(
    filename = function() "Dispersion_plot.png",
    content = function(file){
      p<- pipeline()
      plot_dispersion(p$dseq, out_png=file)
    }
  )
  output$pval_plot<- renderPlot({
    p<- pipeline()
    plot_pvalue_hist(as.data.frame(p$res_list$count_results))
  })
  output$download_pval<- downloadHandler(
    filename = function() "Pvalue_histogram.png",
    content = function(file){
      p<- pipeline()
      plot_pvalue_hist(as.data.frame(p$res_list$count_results), out_png=file)
    }
  )

  #####MA plot#########
  output$ma_plot<- renderPlot({
    p<- pipeline()
    plot_ma(p$ann, p$top__genes, lfc_thr=input$lfc_thr, padj_thr=input$padj_thr_plot)
  })
  output$download_ma<- downloadHandler(
    filename = function() "MA_plot.png",
    content = function(file){
      p<- pipeline()
      plot_ma(p$ann, p$top__genes, lfc_thr=input$lfc_thr, padj_thr=input$padj_thr_plot, out_png=file)
    }
  )

  #####Gene explorer#########
  gene_label<- function(p, id){
    g<- p$ann[id, "GeneName"]
    if(length(g) == 0 || is.na(g)) "" else g
  }
  output$gene_plot<- renderPlot({
    p<- pipeline()
    req(input$gene_pick, input$gene_pick %in% rownames(p$ann))
    grp<- trimws(strsplit(input$Exp_groups, ",")[[1]])[1]
    plot_gene_counts(p$dseq, input$gene_pick, p$meta, group_by=grp, gene_label=gene_label(p, input$gene_pick))
  })
  output$download_gene<- downloadHandler(
    filename = function() paste0("Gene_", gsub("[^A-Za-z0-9._-]", "_", input$gene_pick), ".png"),
    content = function(file){
      p<- pipeline()
      grp<- trimws(strsplit(input$Exp_groups, ",")[[1]])[1]
      plot_gene_counts(p$dseq, input$gene_pick, p$meta, group_by=grp,
                       gene_label=gene_label(p, input$gene_pick), out_png=file)
    }
  )
  output$gene_stats<- renderTable({
    p<- pipeline()
    req(input$gene_pick, input$gene_pick %in% rownames(p$ann))
    r<- p$ann[input$gene_pick, c("GeneName", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")]
    sci<- function(x) ifelse(is.na(x), "NA", formatC(x, format="e", digits=2))
    data.frame(Gene=ifelse(is.na(r$GeneName), input$gene_pick, r$GeneName),
               baseMean=format(round(r$baseMean, 1), nsmall=1), log2FoldChange=format(round(r$log2FoldChange, 3), nsmall=3),
               lfcSE=format(round(r$lfcSE, 3), nsmall=3), pvalue=sci(r$pvalue), padj=sci(r$padj))
  })

  #####Tables and downloads#########
  
  #annotated results preview
  
  #### PATHWAY ANALYSIS (GSEA)
  gsea_obj<- reactive({
    p<- pipeline()
    validate(need(requireNamespace("fgsea", quietly=TRUE),
                  "Pathway analysis needs the 'fgsea' package (run setup/install_packages.R)."))
    validate(need(p$org %in% c("human", "mouse", "rat"),
                  "Pathway analysis is available for human, mouse and rat data."))
    validate(need("stat" %in% names(p$ann), "No test statistic available for ranking genes."))
    withProgress(message="Running pathway analysis (GSEA)...", {
      run_gsea(p$ann, load_gene_sets(), uppercase = p$org != "human")
    })
  })
  
  observeEvent(safe_value(gsea_obj()), {
    tab<- gsea_obj()$table
    updateSelectInput(session, "gsea_pick", choices= tab$pathway, selected= tab$pathway[1])
  })
  
  output$gsea_note<- renderUI({
    req(input$run > 0)
    org<- pipeline()$org
    if(org %in% c("mouse", "rat"))
      helpText("Mouse and rat genes are matched to the human Hallmark sets by gene symbol (approximate).")
  })
  
  output$gsea_bar<- renderPlot({
    g<- gsea_obj()
    plot_gsea_bar(g$table, padj_thr = input$padj_thr_plot)
  })
  output$download_gsea_bar<- downloadHandler(
    filename = function() "gsea_hallmark_pathways.png",
    content = function(file){ g<- gsea_obj(); plot_gsea_bar(g$table, padj_thr = input$padj_thr_plot, out_png = file) }
  )
  
  output$gsea_enrich<- renderPlot({
    req(input$gsea_pick)
    g<- gsea_obj()
    validate(need(input$gsea_pick %in% names(g$sets), "Choose a pathway."))
    plot_gsea_enrichment(g, input$gsea_pick)
  })
  output$download_gsea_enrich<- downloadHandler(
    filename = function() paste0("gsea_", input$gsea_pick, ".png"),
    content = function(file){ plot_gsea_enrichment(gsea_obj(), input$gsea_pick, out_png = file) }
  )
  
  gsea_table_df<- function(){
    g<- gsea_obj()
    tab<- as.data.frame(g$table)
    data.frame(Pathway= tab$pathway, NES= round(tab$NES, 3), pvalue= signif(tab$pval, 3),
               padj= signif(tab$padj, 3), Size= tab$size,
               LeadingEdge= vapply(tab$leadingEdge, function(g) paste(utils::head(g, 15), collapse=", "), character(1)))
  }
  output$gsea_table<- DT::renderDT({
    DT::datatable(gsea_table_df(), rownames=FALSE, options=list(pageLength=10, scrollX=TRUE))
  })
  output$download_gsea_table<- downloadHandler(
    filename = function() "gsea_hallmark_results.csv",
    content = function(file){
      g<- gsea_obj()
      tab<- as.data.frame(g$table)
      tab$leadingEdge<- vapply(tab$leadingEdge, paste, character(1), collapse=";")
      utils::write.csv(tab, file, row.names=FALSE)
    }
  )

  output$results_ann<- DT::renderDT({
    p<- pipeline()
    head(p$ann, 100)
  }, options =list(pageLength=25, scrollX=TRUE))
  #raw DESeq2 results preview
  output$results_raw<- DT::renderDT({
    p<- pipeline()
    head(as.data.frame(p$res_list$count_results), 100)
  }, options =list(pageLength=25, scrollX=TRUE))
  
  #normalized cunt preview
  output$results_norm<- DT::renderDT({
    p<- pipeline()
    head(as.data.frame(p$res_list$normalized_counts), 100)
  }, options =list(pageLength=25, scrollX=TRUE))
  
  #download annotated results
  output$download_ann<- downloadHandler(
    filename = function() "Annotated_results.csv",
    content = function(file){
      p<- pipeline()
      write.csv(p$ann, file)
    }
  )
  
  #download DESeq2 results
  output$download_counts<- downloadHandler(
    filename = function() "DESeq2_results.csv",
    content = function(file){
      p<- pipeline()
      write.csv(as.data.frame(p$res_list$count_results), file)
    }
  )
  
  #download norm counts
  output$download_norm<- downloadHandler(
    filename = function() "Normalized_counts.csv",
    content = function(file){
      p<- pipeline()
      write.csv(as.data.frame(p$res_list$normalized_counts), file)
    }
  )
   
}

##########################Run SHiny APP##############

shinyApp(ui=ui, server= server)

