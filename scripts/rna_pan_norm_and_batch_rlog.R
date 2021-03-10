library(DESeq2)
library(data.table) # To convert row names into first column
library(limma)
library(tools)


# Inputs a library of raw counts with format specified and a column data matrix specifying
# condition and batch.
# 1) Starting from the raw counts specified (all genes), performs Rlog with blind = TRUE. On the Rlog results, 
#    removes Batch Effects while taking into account condition 
#    (see /home/mbp15ja/dev/AuxiliaryPrograms/Atacseq_Aux_Mains/atac_seq_all_quality_MM_ND_samples_by_16_04_2018_sign_dif_atac_seq_rlog_batch_correction.R)
# 2) Add the cytogenetics groups to the sample names on the table (1)
# 
# 
# Inputs:
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
#         NOTE: The base level for "condition" is "PC"
# 
#         A file with format:
#         condition       cd19    batch   donor_id
#       RS_1.10 MM      other   pool1   s1
#       RS_1.11 MM      other   pool1   s2
#
#
# Outputs:
#   -The rlog, batch corrected counts
# 
# Exception:
#   -If the samples in the coldata_file don't coincide with the samples in counts_matrix_file
#     whether in order, number of samples or name.


counts_matrix_file = "salmon.dir/genes.tsv.gz"
coldata_file = "samples.tsv"
output_file = "rlog.dir/rna_pan_rlogs.tsv.gz"

# Check all parameters are not empty
if(counts_matrix_file == '' || is.null(counts_matrix_file)){
  stop("ERROR: counts_matrix_file is empty")
}
# Check all parameters are not empty
if(coldata_file == '' || is.null(coldata_file)){
  stop("ERROR: coldata_file is empty")
}
# Check all parameters are not empty
if(output_file == '' || is.null(output_file)){
  stop("ERROR: output_file is empty")
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
                             row.names=1,
                             stringsAsFactors=FALSE)

# remove cell lines
coldata_df <- subset(coldata_df, Subgroup!="Cell_line")
input_matrix <- input_matrix[,colnames(input_matrix) %in% rownames(coldata_df)]
coldata_df <- coldata_df[rownames(coldata_df) %in% colnames(input_matrix), ]
input_matrix <- input_matrix[, rownames(coldata_df)]

# rename ND to PC for backwards compatitbility
coldata_df$MM.ND[coldata_df$MM.ND=="ND"] <- "PC"

# batch has to be a factor
coldata_df$batch <- as.factor(coldata_df$rna_batch)

# condition has to be a factor
coldata_df$condition <- as.factor(coldata_df$MM.ND)



# Put the healthy (PC) condition as first level so the comparison results (only "log2 fold change") is all levels vs PC.
coldata_df$condition = relevel(coldata_df$condition, "PC")



# Get the data into DESeqDataSet format
dds <- DESeqDataSetFromMatrix(countData = input_matrix,
                              colData = coldata_df,
                              design = ~ batch + condition)


# From /home/mbp15ja/dev/AuxiliaryPrograms/Atacseq_Aux_Mains/atac_seq_all_quality_MM_ND_samples_by_16_04_2018_sign_dif_atac_seq_rlog_batch_correction.R

# Perform the RLog transformation.
# NOTE: the use of blind=TRUE
# We normalize the raw sample reads in peaks
rld <- rlog(dds, blind=TRUE)


# Perform removeBatchEffect


# batch has to be a factor (already done before)
# condition has to be a factor (already done before)

# Put the healthy (PC) condition as first level so the comparison results (only "log2 fold change") is all levels vs PC.
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


# Convert to data.frame
matrix_without_batch = as.data.frame(matrix_without_batch)


# Add the cytogenetics groups to the sample names

new_matrix_without_batch = list()

for(i in 1:ncol(matrix_without_batch)){

  # Get the column name
  sample_name = colnames(matrix_without_batch)[i]

  # Get the corresponding cytogenetic group
  condition_sample = coldata_df[sample_name, "condition"]

  new_sample_name = paste(condition_sample, sample_name, sep="_")

  new_matrix_without_batch <- c(new_matrix_without_batch, new_sample_name)

}

# Assign the new sample names
colnames(matrix_without_batch) = as.list(new_matrix_without_batch)


# Add the starting table feature name
# setDT modifies the input, make a copy of the dataframe
matrix_without_batch_df_ori <- data.frame(matrix_without_batch)


# Output the different differential results
# Convert row names into first column
output_df_matrix_without_batch = setDT(matrix_without_batch_df_ori, keep.rownames = TRUE)[]

# Remove the variable
remove(matrix_without_batch_df_ori)


# Rename the id column
colnames(output_df_matrix_without_batch)[1] <- first_col_name


# Output the file
write.table(output_df_matrix_without_batch, gzfile(output_file),
            append = FALSE,
            quote = FALSE,
            sep = "\t",
            eol = "\n",
            na = "NA",
            dec = ".",
            row.names = FALSE, # The row names are already included
            col.names = TRUE)



