library(Gviz)
library(GenomicRanges)
library(rtracklayer)
library(EnsDb.Hsapiens.v86)
make_granges_from_hic <- function(contig, 
                            start, 
                            end, 
                            anchor, 
                            resolution = 5,
                            norm_method = "VC",
                            base_path = ".",
                            path_to_liftOver=NULL){
  if(!is.null(path_to_liftOver)) {
     ch = import.chain(path_to_liftOver[1])
     window_lift = GRanges(seqnames = contig, ranges=IRanges(start, end))
     hg19_window = unlist(liftOver(window_lift, ch))
     start = start(hg19_window)
     end = end(hg19_window)
  
     anchor = GRanges(seqnames = contig, ranges=IRanges(anchor, anchor+1))
     hg19_anchor <- unlist(liftOver(anchor, ch))
     anchor <- start(hg19_anchor)
  }
  
  anchor_bin <- (anchor %/% (resolution*1000)) * (resolution*1000)
  
  path_to_files = file.path(base_path, paste0(resolution, "kb_resolution_intrachromosomal"),
               contig,
               "MAPQGE30")
  RAW_file <- file.path(path_to_files,
                        paste0(contig, "_", resolution, "kb.RAWobserved"))
  raw_data <- read.delim(RAW_file, header=FALSE, col.names=c("from", "to", "count"))
  to_subset <- raw_data[raw_data$to == anchor_bin & 
                          raw_data$from >= start & 
                          raw_data$from < end, 
                        c("from", "count")]
  from_subset <- raw_data[raw_data$from == anchor_bin & 
                            raw_data$to >= start & 
                            raw_data$to < end, 
                          c("to", "count")]
  rm(raw_data)
  colnames(to_subset) <- c("start", "count")
  colnames(from_subset) <- c("start", "count")
  full_data = rbind(to_subset, from_subset)
  full_data = full_data[order(full_data$start),]
  full_data$end = full_data$start+(resolution*1000)
  full_data$count[full_data$start == anchor_bin] <- 0
  
  norm_file <- file.path(path_to_files,
                         paste0(contig, "_", resolution, "kb.", norm_method, "norm"))
  norm_data = read.delim(norm_file, header=FALSE, col.names=c("norm_factor"))
  norm_data$start <- (1:nrow(norm_data) -1 ) * resolution * 1000
  full_data <- merge(full_data, norm_data)
  rm(norm_data)
  full_data$normalized_count <- full_data$count/full_data$norm_factor
  full_data$chromosome = contig
  final_results <- makeGRangesFromDataFrame(full_data, keep.extra.columns = TRUE)
  if (!is.null(path_to_liftOver)) {
    final_results <- unlist(liftOver(final_results, import.chain( path_to_liftOver[2])))
    
  }
  return(final_results)
  
}

gene_plot <- function(chr, start, end, 
                      atac_mm, atac_pc, atac_ylim=4,
                      h3k27ac=NA,h3k27ac_ylim=15,
                      maf=NA, maf_ylim=60,
                      hic=FALSE, hic_ylim=400,
                      footprints=NA, fp_ylim=75,  
                      promoter=NA) {
  displayPars(atac_pc) <- list(name = "PC ATAC", ylim=c(0,atac_ylim))
  if (length(atac_mm) > 0){
    for (i in 1:length(atac_mm)) {
       displayPars(atac_mm[[i]]) <- list(alpha = 1, alpha.title=1,  ylim=c(0,atac_ylim))
      print(getPar(atac_mm[[i]], "alpha"))
    }
    #atac_mm <- OverlayTrack(atac_mm)
  }
  
  track_list = c(atac_mm,
                 atac_pc)
  track_heights = c(rep(1, times=length(atac_mm)),1)
  
  if (!is.na(h3k27ac)) {
    if (length(h3k27ac) > 0) {
      for( i in 1:length(h3k27ac)) {
        displayPars(h3k27ac[[i]]) <- list(name="H3K27ac", ylim=c(0,h3k27ac_ylim))
      }
      print(getPar(h3k27ac[[2]], "alpha"))
      h3k27ac <- OverlayTrack((h3k27ac))
    }
    track_list = c(h3k27ac, track_list)
    track_heights = c(1, track_heights)
  }
  
  if (!is.na(footprints) & length(footprints) > 0) {
    for( i in 1:length(footprints)) {
      displayPars(footprints[[i]]) <- list(alpha = max(1/i, 0.35), name="Footprints", ylim=c(0,fp_ylim))
    }
    footprints <- OverlayTrack(footprints)
    track_list = c(track_list, footprints)
    track_heights = c(track_heights, 1)
  }
  
  if (!is.na(maf)) {
    displayPars(maf) <- list(ylim=c(0, maf_ylim))
    track_list= c(maf, track_list)
    track_heights = c(0.7, track_heights)
  }

  if (hic) {
    hic = make_granges_from_hic(chr, start, end, promoter, resolution=5, base_path="../data/HiC/",
                                path_to_liftOver = c("/shared/sudlab1/General/mirror/ucsc/hg38/hg38ToHg19.over.chain",
                                                     "/shared/sudlab1/General/mirror/ucsc/hg19/hg19ToHg38.over.chain"))
    hic <- DataTrack(hic, data="count", name="HiC", type="polygon", col.mountain="grey75", fill.mountain=c("grey75", "grey75"),
                     ylim=c(0,hic_ylim))
    track_list=c(hic, track_list)
    track_heights=c(0.7,track_heights)
  }

   seqlevelsStyle(EnsDb.Hsapiens.v86) <- "UCSC"
   genes_track <- GeneRegionTrack(getGeneRegionTrackForGviz(EnsDb.Hsapiens.v86, chr=chr, start=start, end=end),
                                          transcriptAnnotation="symbol", collapseTranscripts="meta", name=NULL, fontsize.group=8,
                                          cex.group=1,
                                          fill="grey20",
                                          col="grey20")
   track_list <- c(track_list, genes_track, GenomeAxisTrack())
   
   track_heights = c(track_heights, 2,2)
   plotTracks(track_list,
              chromosome=chr, from=start, to=end, 
              background.title="white", col.title="black", fontcolor.title="black",
              cex.title=0.8,
              showAxis=FALSE, 
              fontface.title=1,
              fontsize=9,
              rotation.title=0,
              sizes=track_heights
   ) 
  
}
