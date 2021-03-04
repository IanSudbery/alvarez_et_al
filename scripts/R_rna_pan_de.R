library(DESeq2)
library(data.table) # To convert row names into first column
library(limma)
library(tools)
library(gplots)

# Inputs a library of raw counts with format specified and a column data matrix specifying
# condition and batch.
# 1) Performs a Wald Test for condition accounting for batch. All results are outputted with outlier detection (genes
#	where an outlier is found they are outputted with their p-value/padj).
# 2) Outputs the result from the Wald test to a table (de_res_batch_condition_output_file).
# 3) Gets genes which are significant for MM vs. PC condition 
#   [padj<pdaj_threshold and padj not NA and log2FoldChange >= abs(very_significant_threshold)]
# 4) Outputs 
#   A table with the Wald results:
#     -basemean
#     -pvalue
#     -padj
#     -log2FoldChange_MM_vs_ND
#     -lfcSE_MM_vs_ND: Standard error for the log2FoldChange
#     -DE_subgroup_vs_ND: Whether MM vs. PC is Differentially expressed with respect to ND in terms of the log2FoldChange:
#       -1 if log2FoldChange <= -1
#       +1 if log2FoldChange >= +1
#       0 if abs(log2FoldChange) < 1
#
# 6) Starting from the raw counts, performs Rlog with blind = TRUE. On the Rlog results, 
#    removes Batch Effects while taking into account condition 
#    (see /home/mbp15ja/dev/AuxiliaryPrograms/Atacseq_Aux_Mains/atac_seq_all_quality_MM_ND_samples_by_16_04_2018_sign_dif_atac_seq_rlog_batch_correction.R)
# 7) Add the cytogenetics condition groups to the sample names on the table (6) 
#    and plots heatmaps (/home/mbp15ja/dev/AuxiliaryPrograms/Atacseq_Aux_Mains/atac_seq_all_quality_MM_ND_samples_by_16_04_2018_sign_dif_atac_seq_rlog_batch_correction_heatmap_biclustering_all_metrics_linkages.R) 
#    where the samples and rows are structured on the correlation metric on the inputted data and the euclidean metric on the standardized data.
#    Note, that the heatmaps have as rows all the regions which are significant in the Log Ratio Test testing for condition and accounting for batch + significant abs log2FoldChanges. 
#    This means that the regions can have overexpression of one or more subgroups and NDs (we not only have underexpressed ND).
#
# 
# 
# Inputs:
#   -very_significant_threshold: Only genes where the log2FoldChange of MM vs ND >= than
#       abs(this) will be outputted.
#   -pdaj_threshold: Only genes where the padj of the Wald test < pdaj_threshold
#       will be outputted.
#   -counts_matrix_file: A file with format:
#         id      RS_1.10 RS_1.11 RS_1.12
#         ENSG00000000003 36.0    745.0   361.0
#         ENSG00000000005 0.0     0.0     1.0
# 
#   -coldata_file: A matrix specifying each column of the counts_matrix_file
#         and other attributes.
#         NOTE: The samples are row names (no header).
#         NOTE: The order of the rows is in this matrix is the same as the order of the sample columns
#               in the counts_matrix_file
#         NOTE: "condition" and "batch" columns are compulsory in whatever order
#         NOTE: The base level for "condition" is "ND"
# 
#         A file with format:
#         condition       cd19    batch   donor_id
#       RS_1.10 MM      other   pool1   s1
#       RS_1.11 MM      other   pool1   s2
#
#
#   -de_res_batch_condition_output_file: A matrix with the results from Wald test.
#
#   -output_heatmaps: A png file where to use as basename to output all the heatmaps.
# 
# Outputs:
#   -de_res_batch_condition_output_file: A matrix with the results from Wald test (all).
# 
# Exception:
#   -If the samples in the coldata_file don't coincide with the samples in counts_matrix_file
#     whether in order, number of samples or name.


very_significant_threshold = 1.5
pdaj_threshold = 0.1
counts_matrix_file = "salmon.dir/genes.donors_collapsed.tsv.gz"
coldata_file = "salmon.dir/genes.donors_collapsed.col_data.tsv"
de_res_batch_condition_output_file = "DE.dir/all_pan_rna.tsv.gz"
sign_batch_MM_vs_PC_file = "DE.dir/sign_pan_rna.tsv.gz"
output_heatmaps = "DE.dir/sign_pan_rna.png"



# Check all parameters are not empty
if(very_significant_threshold == '' || is.null(very_significant_threshold)){
  stop("ERROR: very_significant_threshold is empty")
}
# Check all parameters are not empty
if(counts_matrix_file == '' || is.null(counts_matrix_file)){
  stop("ERROR: counts_matrix_file is empty")
}
# Check all parameters are not empty
if(coldata_file == '' || is.null(coldata_file)){
  stop("ERROR: coldata_file is empty")
}
# Check all parameters are not empty
if(de_res_batch_condition_output_file == '' || is.null(de_res_batch_condition_output_file)){
  stop("ERROR: de_res_batch_condition_output_file is empty")
}
# Check all parameters are not empty
if(sign_batch_MM_vs_PC_file == '' || is.null(sign_batch_MM_vs_PC_file)){
  stop("ERROR: sign_batch_MM_vs_PC_file is empty")
}





# Store the column name of the first column to output it later
input_matrix <- read.table(counts_matrix_file,
                           sep="\t",
                           check.names=FALSE, # To remove the trailing X from number           columns
                           header=FALSE)

# Prevent conversion of strings to factors
column <- sapply(input_matrix, is.factor)
input_matrix[column] <- lapply(input_matrix[column], as.character)

first_col_name = input_matrix[1,1]


# Read in the dataframe with the counts
input_matrix <- as.matrix(read.table(counts_matrix_file,
                                     sep="\t",
                                     check.names=FALSE, # To remove the trailing X from number columns
                                     header=TRUE,
                                     row.names=1)) # Get the gene name as index

# Read in the dataframe with the column data
coldata_df <- read.table(coldata_file,
                             sep="\t",
                             check.names=FALSE, # To remove the trailing X from number columns
                             header=TRUE,
                             row.names=1)

# Remove cell lines
coldata_df <- coldata_df[coldata_df$MM.ND %in% c("MM", "ND"),]
input_matrix <- input_matrix[, colnames(input_matrix) %in% rownames(coldata_df)]

# batch has to be a factor
coldata_df$batch <- as.factor(coldata_df$rna_batch)

# condition has to be a factor
coldata_df$condition <- as.factor(coldata_df$MM.ND)



# Put the healthy (ND) condition as first level so the comparison results (only "log2 fold change") is all levels vs ND.
coldata_df$condition = relevel(coldata_df$condition, "ND")

# Put the pool1 as the reference level for batch
coldata_df$batch = relevel(coldata_df$batch, "batch1")

print(coldata_df)
# Get the data into DESeqDataSet format
dds <- DESeqDataSetFromMatrix(countData = input_matrix,
                              colData = coldata_df,
                              design = ~ batch + condition)


# Perform a Wald test since we only have 2 conditions (MM and PC) accounting for batch.
# Instead of a Likelihood Ratio Test (LRT) 
dds = DESeq(dds)

# dds.LRT <- DESeq(dds,
#                  betaPrior=FALSE,
#                  test="LRT",
#                  full=~ batch + condition,
#                  reduced=~ batch)


# Get the result with outlier detection
res_batch_condition=as.data.frame(results(dds))

# setDT modifies the input, make a copy of the dataframe
res_batch_condition_ori <- data.frame(res_batch_condition)


# Output the different differential results
# Convert row names into first column
output_df_res_batch_condition = setDT(res_batch_condition_ori, keep.rownames = TRUE)[]

# Remove the variable
remove(res_batch_condition_ori)


# Rename the id column
colnames(output_df_res_batch_condition)[1] <- first_col_name


# Output the file
write.table(output_df_res_batch_condition, gzfile(de_res_batch_condition_output_file),
            append = FALSE,
            quote = FALSE,
            sep = "\t",
            eol = "\n",
            na = "NA",
            dec = ".",
            row.names = FALSE, # The row names are already included
            col.names = TRUE)


# We get anything which is significant for MM vs. PC
sign_batch_condition = res_batch_condition[
  (!(is.na(res_batch_condition$padj))&(res_batch_condition$padj<pdaj_threshold)&(abs(res_batch_condition$log2FoldChange)>=very_significant_threshold)),]


# Now add a column for MM vs. PC which specifies the log2FoldChange
# 0 -> abs(log2FoldChange) < very_significant_threshold
# 1 -> log2FoldChange >= very_significant_threshold
# -1 -> log2FoldChange <= very_significant_threshold

# Start by assigning 0 to everything, then assign the corresponding columns the other values
sign_batch_condition[,"DE_MM_vs_ND"] = 0
sign_batch_condition[sign_batch_condition$log2FoldChange<=-very_significant_threshold,"DE_MM_vs_ND"] = -1
sign_batch_condition[sign_batch_condition$log2FoldChange>=very_significant_threshold,"DE_MM_vs_ND"] = 1


# Output table
# setDT modifies the input, make a copy of the dataframe
sign_batch_condition_ori <- data.frame(sign_batch_condition)


# Output the different differential results
# Convert row names into first column
output_sign_batch_condition = setDT(sign_batch_condition_ori, keep.rownames = TRUE)[]

# Remove the variable
remove(sign_batch_condition_ori)


# Rename the id column
colnames(output_sign_batch_condition)[1] <- first_col_name


# Output the file
write.table(output_sign_batch_condition, 
            gzfile(sign_batch_MM_vs_PC_file),
            append = FALSE,
            quote = FALSE,
            sep = "\t",
            eol = "\n",
            na = "NA",
            dec = ".",
            row.names = FALSE, # The row names are already included
            col.names = TRUE)


# Get the unique very significant regions
unique_very_sign_regions = rownames(sign_batch_condition)


############################### Heatmap of the sample size normalized input matrix ###############################
############################### only for (significant (LRT condition + batch) and 
############################### abs(log2FoldChange) >= very_sign_threshold) features ###############################

# From /home/mbp15ja/dev/AuxiliaryPrograms/Atacseq_Aux_Mains/atac_seq_all_quality_MM_ND_samples_by_16_04_2018_sign_dif_atac_seq_rlog_batch_correction.R

# Perform the RLog transformation.
# NOTE: the use of blind=TRUE
# We normalize the raw sample reads in peaks
rld <- rlog(dds, blind=TRUE)


# Perform removeBatchEffect


# batch has to be a factor (already done before)
# condition has to be a factor (already done before)

# Put the healthy (ND) condition as first level so the comparison results (only "log2 fold change") is all levels vs ND.
# Already done

# Create the model matrix for condition
model_matrix = model.matrix(~condition, coldata_df)


# Separate the coldata_df into a matrix with batch (code independent of column order)
batch_coldata_df = as.data.frame(coldata_df$batch)
rownames(batch_coldata_df) = rownames(coldata_df)
colnames(batch_coldata_df) = c("batch")



# According to this https://support.bioconductor.org/p/60879/
# We specify the original ~condition + batch
# use design=model.matrix(~condition), and batch=batch
matrix_without_batch <- removeBatchEffect(x=as.matrix(assay(rld)),
                                          batch=as.matrix(batch_coldata_df),
                                          design=model_matrix)

# Once the data matrix with all the features is normalized and batch effects are removed
# subselect on it only the differential ATAC-seq regions (LRT) and significant log2FoldChanges.
matrix_without_batch_significant_LRT = as.data.frame(matrix_without_batch[unique_very_sign_regions,])



# Add the cytogenetics groups to the sample names

new_matrix_without_batch_significant_LRT_samples = list()

for(i in 1:ncol(matrix_without_batch_significant_LRT)){

  # Get the column name
  sample_name = colnames(matrix_without_batch_significant_LRT)[i]

  # Get the corresponding cytogenetic group
  condition_sample = coldata_df[sample_name, "condition"]

  new_sample_name = paste(condition_sample, sample_name, sep="_")

  new_matrix_without_batch_significant_LRT_samples <- c(new_matrix_without_batch_significant_LRT_samples, new_sample_name)

}

# Assign the new sample names
colnames(matrix_without_batch_significant_LRT) = as.list(new_matrix_without_batch_significant_LRT_samples)


# Produce the heatmaps based on /home/mbp15ja/dev/AuxiliaryPrograms/Atacseq_Aux_Mains/atac_seq_all_quality_MM_ND_samples_by_16_04_2018_sign_dif_atac_seq_rlog_batch_correction_heatmap_biclustering_all_metrics_linkages.R

metrics=c("euclidean", "cor")
linkage_methods=c("ward.D2", "complete", "average")


# Prevent conversion of strings to factors
column <- sapply(matrix_without_batch_significant_LRT, is.factor)
matrix_without_batch_significant_LRT[column] <- lapply(matrix_without_batch_significant_LRT[column], as.character)



# Initialize the variable for default names to FALSE, if we don't find any metrics specified
# we will use default names
default_names = FALSE


# If multiple use default names based on the output provided which change from run to run
if (length(metrics) > 1 || length(linkage_methods) > 1){
  default_names = TRUE
}



output_file_extension = file_ext(output_heatmaps)

output_filename = file_path_sans_ext(output_heatmaps, compression = FALSE)

# Create minimum log file to store the progression
log_file = paste(output_filename, ".log", sep="")

# Delete the log file if it exists
unlink(log_file)



# This is done independently of the metric and linkage method so do it only once
# Get the matrix to cluster (coordinates are the rownames)
# The matrix contains features on rows and samples on columns
x = as.matrix(matrix_without_batch_significant_LRT)

#Go through all the possibilities
for(metric in metrics){
  
  for(linkage_method in linkage_methods){
    
    # Create the variables for the names empty
    new_filename = ""
    
    
    # If multiple metrics and linkages are to be used
    # the names to use must vary for every round
    if (default_names == TRUE){
      
      new_filename = paste(output_filename, linkage_method, metric, "heatmap_dendro", sep="_")
      
      new_filename = paste(new_filename, output_file_extension, sep=".")
      
    } else{     # If only one type of output produced use specified names
      
      new_filename = output_heatmaps
      
    }
    
    # Delete the files
    unlink(new_filename)
    
    
    # Output name to log file
    cat(new_filename, file=log_file, sep="\n", append=TRUE)
    
    
    
    # Compute distance matrix for samples and features
    # For correlation and all linkage methods, we will produce euclidean metric
    if(metric == "cor"){
      
      # First we are going to calculate the distance matrix between samples (cell lines)
      # (using the metric specified).
      
      # Compute correlation matrix
      # The function cor() computes pairwise correlation coefficients
      # between the columns of the data. The original matrix has samples
      # on columns.
      samples_distance_matrix.cor <- cor(x, method = "pearson")
      
      # Compute distance matrix for samples
      # Completely disimilar samples will have -1
      # Completely similar samples will have 1
      # This means that the distance between two samples after applying this
      # Goes from 0 (if they have 1 correlation) to 2 if they have correlation -1
      samples_distance_matrix <- as.dist(1 - samples_distance_matrix.cor)
      
      
      # Cleanup variables
      rm(samples_distance_matrix.cor)
      
      
      # Now compute distance matrix between features (using the metric specified).
      
      
      # Compute correlation matrix
      # The function cor() computes pairwise correlation coefficients
      # between the columns of the data. The original matrix has features
      # on rows, so we need to transpose it
      features_distance_matrix.cor <- cor(t(x), method = "pearson")
      
      
      # Compute distance matrix for features
      features_distance_matrix <- as.dist(1 - features_distance_matrix.cor)
      
      
      # Cleanup variables
      rm(features_distance_matrix.cor)
      
      
    } else{ # For non cor metrics
      
      # If we are using euclidean metric, we put the features of each dendrogram
      # in the same footing by scaling them
      if(metric == "euclidean"){
        
        # Scale the rows
        # Scale function scales columns
        row_scale_data = t(scale(t(x), center = TRUE, scale = TRUE))
        
        # Now that the rows are scaled, all the rows are on the same footing
        # And we can hierarchical cluster columns (samples) since all dimensions
        # count the same
        samples_distance_matrix = dist(t(row_scale_data), method=metric)
        
        
        # Scale columns
        col_scale_data = scale(x, center = TRUE, scale = TRUE)
        
        # Now that the samples are scaled, all samples are on equal footing
        # And we can hierarchical cluster rows (features) since all samples
        # count the same
        features_distance_matrix = dist(col_scale_data, method=metric)
        
        
        # Cleanup variables
        rm(row_scale_data)
        
        # Cleanup variables
        rm(col_scale_data)
        
      } else{
        
        # First calculate the samples distance
        
        # dist() calculates distances between rows, samples are on columns, therefore
        # transpose the matrix
        # Calculate the distance matrix between samples
        samples_distance_matrix = dist(t(x), method=metric)
        
        
        # Then calculate the distance between features, features are on rows
        features_distance_matrix = dist(x, method=metric)
        
      }
      
      
    }
    
    # Now produce the heatmap and dendrograms
    # The rows and right will contain features
    # Apply hclust with the linkage method specified
    samples_dend = as.dendrogram(hclust(samples_distance_matrix, method=linkage_method))
    features_dend = as.dendrogram(hclust(features_distance_matrix, method=linkage_method))
    
    
    # Create the file to store the plot
    png(file=new_filename, width=750, height=700)
    
    # Plot it
    # Reduce the margins
    # Get rid of the row dendrogram (regions)
    # Scale the rows so that all the features are in the same range
    # Put a monochromatic pallette to see the relationships as 1 colour scale from white to blue
    # breaks=seq(from=-3.0,to=3.5, length.out = 20)
    heatmap.2(x,
              Rowv = features_dend,
              Colv = samples_dend,
              trace="none", margins=c(8,5),
              labRow = NA,
              dendrogram="column",
              scale="row",
              col=colorRampPalette(c("white","darkblue")))
    
    dev.off()
    
    # Output execution
    cat("Plotting executed", file=log_file, sep="\n", append=TRUE)
    
  }
}




