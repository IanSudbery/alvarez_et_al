library(dplyr)

filter_and_process_genesets <- function(dbs_df, geneset_categories=NULL, id_type="ensembl_gene_id",  min_size=0, max_size=Inf) {
  if (!(id_type %in% c("ensembl_gene_id", "gene_symbol"))) {
    stop(paste0("ID type ", id_type, "is not 'ensembl_gene_id' or 'gene_symbol'"))
  }

  if (is.null(geneset_categories)) {
    geneset_categories <- unique(dbs_df$gs_cat)
  }
  obj <- dbs_df %>%
    filter(gs_cat %in% geneset_categories) %>%
    group_by(gene_set) %>%
    filter(n() > min_size, n() <=max_size) %>%
    ungroup() %>%
    na.omit() 
    
  obj <- split(obj$geneset_id, obj[,id_type], drop=TRUE)
  return(obj)
    
}

do_enrichment_on_geneset <- function(pwf, genesets) {
  
  go_seq_result <- suppressMessages(goseq::goseq(pwf, gene2cat = genesets))
  
  all_genes <- intersect(rownames(pwf), names(genesets))
  numDEgenes <- sum(pwf[all_genes,]$DEgenes)
  fDEgenes <- numDEgenes/length(all_genes)
  fInCat <- go_seq_result$numDEInCat/go_seq_result$numInCat
  
  n11 = go_seq_result$numDEInCat
  n01 = go_seq_result$numInCat - go_seq_result$numDEInCat
  n10 = numDEgenes - go_seq_result$numDEInCat
  n00 = length(all_genes) - n11 - n01 - n10
  
  go_seq_result$OR = (n00*n11)/(n01*n10)
  go_seq_result$enrichment <- fInCat/fDEgenes
  
  return(go_seq_result)
  
}

run_goseq_on_cat <- function(category, pwf, database, gene_id) {
  geneset_df <- dbGetQuery(database,
                           paste0("SELECT gene_set, ", 
                                  gene_id, 
                                  " as gene_id FROM genecats WHERE gs_cat='",
                                  category, "'"))
  geneset_list <- split(geneset_df$gene_set, geneset_df$gene_id)
  return(do_enrichment_on_geneset(pwf, geneset_list))
}


run_goseq_on_enrichr_cats <- function(pwf, database, gene_id="ensembl_gene_id", cats=NA,cores=1) {
  
  db <- RSQLite::dbConnect(RSQLite::SQLite(), database)
  geneset_info <- DBI::dbGetQuery(db, "SELECT * FROM geneset_info")
  
 
  if (is.na(cats)) {
    cats <- unique(geneset_info$gs_cat)
  }
  
  if (gene_id == "ensembl") {
    gene_id <- "ensembl_gene_id"
  } else if (gene_id=="symbol") {
    gene_id <- "gene_symbol"
  }
  
  if (!(gene_id %in% c("ensembl_gene_id", "gene_symbol"))) {
    stop("Geneset not recognised")
  }
  if (cores==1) {
    
    enrichement_results <- sapply(cats,
                                  run_goseq_on_cat,
                                  pwf=pwf,
                                  database=db,
                                  gene_id=gene_id,
                                  simplify=FALSE)
  } else {
    cl = makeCluster(cores)
    enrichment_results <- parSapply(cl,
                                    run_goseq_on_cat,
                                    pwf=pwf,
                                    database=db,
                                    gene_id=gene_id,
                                    simplify=FALSE)
  }
 
  enrichment_results <- dplyr::bind_rows(enrichement_results, .id="gs_cat")
  enrichment_results <- dplyr::inner_join(enrichment_results, geneset_info, by=c("gs_cat", "category"="gene_set"))
  enrichment_results$FDR <- p.adjust(enrichment_results$over_represented_pvalue, method="BH")
  return(enrichment_results)
}

get_default_genesets <- function(database) {
  db <- RSQLite::dbConnect(RSQLite::SQLite(), database)
  geneset_info <- DBI::dbGetQuery(db, "SELECT * FROM geneset_info")
  data.frame(gs_cat=unique(geneset_info$gs_cat)) %>%
    tidyr::extract(gs_cat, "(.+)(20[0-9]{2}.*)", into=c("cat", "year"), remove=FALSE ) %>%
    arrange(cat, year) %>%
    group_by(cat) %>%
    filter(is.na(year) | year==max(year)) %>%
    ungroup() %>%
    filter(!grepl("[mM]ouse", libraryName),
           !grepl("NIH", libraryName),
           !grepl("Enrichr", libraryName)) -> filtered_libraries
  
  return(filtered_libraries)
}

