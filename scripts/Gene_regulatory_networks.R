library(ChIPpeakAnno)

#granges function
bed_to_granges <- function(file){
   df <- read.table(file,
                    header=F,
                    stringsAsFactors=F)
 
   if(length(df) > 6){
      df <- df[,-c(7:length(df))]
   }
 
   if(length(df)<3){
      stop("File has less than 3 columns")
   }
 
   header <- c('chr','start','end','id','score','strand')
   names(df) <- header[1:length(names(df))]
 
   if('strand' %in% colnames(df)){
      df$strand <- gsub(pattern="[^+-]+", replacement = '*', x = df$strand)
   }
 
   library("GenomicRanges")
 
   if(length(df)==3){
      gr <- with(df, GRanges(chr, IRanges(start, end)))
   } else if (length(df)==4){
      gr <- with(df, GRanges(chr, IRanges(start, end), id=id))
   } else if (length(df)==5){
      gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score))
   } else if (length(df)==6){
      gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score, strand=strand))
   }
   return(gr) }
   
IDs_to_names<-read.table("IDs_to_TF_names.txt")
TFs_ranges<-bed_to_granges("TF_genes_loci.bed")
Res1_TFs<-read.table("R1_Originals_MAX_express_TFs.txt", header=TRUE)

annotatedNDs_20kb<-annotatedNDs[(annotatedNDs$distancetoFeature <20000) & (annotatedNDs$distancetoFeature >-20000)]
write.table(annotatedNDs_100kb,file="AssignmentNDs.txt", sep="\t")


annotatedHD_20kb_ValidTFs <- annotatedHD_20kb[annotatedHD_20kb$id %in% Res1_TFs$HD,]
write.table(annotatedHD_20kb_ValidTFs,file="AssignmentHD_20kb_ValidTFs.txt", sep="\t")

res1_converted <- apply(Res1_TFs, 2, function(x) sub("^(.+)\\.txt", "\\1", x))

### Example for MMSET subgroup
MMSET_ranges<-bed_to_granges("ALL_FTs_MMSET_sorted.bed")
annotatedMMSET<-annotatePeakInBatch(MMSET_ranges, AnnotationData=TFs_ranges)
annotatedMMSET_20kb<-annotatedMMSET[(annotatedMMSET$distancetoFeature <20000) & (annotatedMMSET$distancetoFeature >-20000)]
write.table(annotatedMMSET_20kb,file="AssignmentMMSET_20kb.txt", sep="\t")
annotatedMMSET_20kb_ValidTFs <- annotatedMMSET_20kb[annotatedMMSET_20kb$id %in% as.character(Res1_TFs$MMSET),]
write.table(annotatedMMSET_20kb_ValidTFs,file="AssignmentMMSET_20kb_ValidTFs.txt", sep="\t")
