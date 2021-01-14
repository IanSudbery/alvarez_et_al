library(DESeq2)
library(data.table) # To convert row names into first column
library(limma)
library(tools)
library(gplots)

# Inputs a library of raw counts with format specified and a column data matrix specifying
# condition and batch.
# 1) Performs a Log Ratio Test for condition accounting for batch. All results are outputted without outlier detection (regions
#	where an outlier is found are still outputted with p-value/padj).
# 2) Outputs the details for the LRT comparison and all subgroups vs ND.
#   A table with the LRT results:
#     -basemean
#     -pvalue
#     -padj
#
#   And for each subgroup vs ND the individual results:
#     -log2FoldChange_subgroup_vs_ND
#     -lfcSE_subgroup_vs_ND: Standard error for the log2FoldChange
#
# 
# 
# Inputs:
#   -very_significant_threshold: Threshold for each subgroup abs(log2FoldChange) vs ND
#       abs(this) to also have padj < pdaj_threshold and be marked as DE
#   -pdaj_threshold: Only regions where the padj of the LRT < pdaj_threshold
#       will be outputted to the significant subgroup regions (For example, sign_LRT_batch_condition_HD_vs_ND_unique_file)
#
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
#   -output_df_all_details_file: A file with the details for the LRT comparison and
#       all subgroups vs ND. Format:
#                                baseMean pvalue_LRT  padj_LRT log2FoldChange_HD_vs_ND  lfcSE_HD_vs_ND log2FoldChange_CCND1_vs_ND lfcSE_CCND1_vs_ND log2FoldChange_MAF_vs_ND lfcSE_MAF_vs_ND log2FoldChange_MMSET_vs_ND lfcSE_MMSET_vs_ND DE_HD_vs_ND DE_CCND1_vs_ND DE_MAF_vs_ND  DE_MMSET_vs_ND
#   chr1:100028407-100029205 19.39137  0.2165052 0.5307479              -0.409905 
#
#
# 
# Outputs:
#   -output_df_all_details_file
# 
# Exception:
#   -If the samples in the coldata_file don't coincide with the samples in counts_matrix_file
#     whether in order, number of samples or name.
#
very_significant_threshold = 1
pdaj_threshold = 0.1
counts_matrix_file = "tag_counts.dir/subtype_raw_tag_counts.donors_collapsed.tsv.gz"
coldata_file = "tag_counts.dir/subtype_raw_tag_counts.donors_collapsed.col_data.tsv"
output_df_all_details_file = "DE.dir/all_subtype_atac_regions.tsv.gz"


# Check all parameters are not empty
if(counts_matrix_file == '' || is.null(counts_matrix_file)){
  stop("ERROR: counts_matrix_file is empty")
}
# Check all parameters are not empty
if(coldata_file == '' || is.null(coldata_file)){
  stop("ERROR: coldata_file is empty")
}
# Check all parameters are not empty
if(output_df_all_details_file == '' || is.null(output_df_all_details_file)){
  stop("ERROR: output_df_all_details_file is empty")
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


# batch has to be a factor
coldata_df$batch <- as.factor(coldata_df$atac_batch)

# condition has to be a factor
coldata_df$condition <- as.factor(coldata_df$Subgroup)

coldata_df <- subset(coldata_df, !(Subgroup %in% c("Cell_line", "UNKNOWN")))

# Put the healthy (ND) condition as first level so the comparison results (only "log2 fold change") is all levels vs ND.
coldata_df$condition = relevel(coldata_df$condition, "ND")

# Remove from the count matrix the samples that have been removed from the colData
input_matrix <- input_matrix[, rownames(coldata_df)]

# Get the data into DESeqDataSet format
dds <- DESeqDataSetFromMatrix(countData = input_matrix,
                              colData = coldata_df,
                              design = ~ batch + condition)


# Perform a Likelihood Ratio Test (LRT) to compare the good of fit of two statistical models with the null model ("reduced") only taking into account the batch being a special case of the alternative model ("full") taking into account batch and condition.
dds.LRT <- DESeq(dds,
                 betaPrior=FALSE,
                 test="LRT",
                 full=~ batch + condition,
                 reduced=~ batch)



# Get the result without outlier detection
res_LRT_batch_condition=as.data.frame(results(dds.LRT, cooksCutoff=FALSE))





res_HD_vs_ND=as.data.frame(results(dds.LRT, name="condition_HD_vs_ND", cooksCutoff=FALSE))

# The p-values are the same for each comparison (because the p-values mean the
# full=~ batch + condition
# vs
# reduced=~ batch
# is significant.


res_CCND1_vs_ND=as.data.frame(results(dds.LRT, name="condition_CCND1_vs_ND", cooksCutoff=FALSE))

res_MAF_vs_ND=as.data.frame(results(dds.LRT, name="condition_MAF_vs_ND", cooksCutoff=FALSE))

res_MMSET_vs_ND=as.data.frame(results(dds.LRT, name="condition_MMSET_vs_ND", cooksCutoff=FALSE))


# Get a table with the general details (LRT)
output_df_general_details = res_LRT_batch_condition[,c("baseMean","pvalue", "padj")]

# Remove the results object (no longer used)
remove(res_LRT_batch_condition)

# Rename some columns
colnames(output_df_general_details)=c("baseMean", "pvalue_LRT", "padj_LRT")


# Now we have to add the individual details for each subgroup comparison in the table
# When merging by rownames, a new column is created, assign this as row names and delete it
output_df_general_details = transform(merge(output_df_general_details,
                                            res_HD_vs_ND[,c("log2FoldChange","lfcSE"), drop=FALSE], # So it keeps it as a dataframe
                                            by=0,
                                            all=TRUE),
                                      row.names=Row.names,
                                      Row.names=NULL)

# Rename the newly added columns
col_id_log2FoldChange = grep("^log2FoldChange$", colnames(output_df_general_details))
colnames(output_df_general_details)[col_id_log2FoldChange] <- "log2FoldChange_HD_vs_ND"

col_id_lfcSE = grep("^lfcSE$", colnames(output_df_general_details))
colnames(output_df_general_details)[col_id_lfcSE] <- "lfcSE_HD_vs_ND"



# Now we have to add the individual details for each subgroup comparison in the table
# When merging by rownames, a new column is created, assign this as row names and delete it
output_df_general_details = transform(merge(output_df_general_details,
                                            res_CCND1_vs_ND[,c("log2FoldChange","lfcSE"), drop=FALSE], # So it keeps it as a dataframe
                                            by=0,
                                            all=TRUE),
                                      row.names=Row.names,
                                      Row.names=NULL)

# Rename the newly added columns
col_id_log2FoldChange = grep("^log2FoldChange$", colnames(output_df_general_details))
colnames(output_df_general_details)[col_id_log2FoldChange] <- "log2FoldChange_CCND1_vs_ND"

col_id_lfcSE = grep("^lfcSE$", colnames(output_df_general_details))
colnames(output_df_general_details)[col_id_lfcSE] <- "lfcSE_CCND1_vs_ND"





# Now we have to add the individual details for each subgroup comparison in the table
# When merging by rownames, a new column is created, assign this as row names and delete it
output_df_general_details = transform(merge(output_df_general_details,
                                            res_MAF_vs_ND[,c("log2FoldChange","lfcSE"), drop=FALSE], # So it keeps it as a dataframe
                                            by=0,
                                            all=TRUE),
                                      row.names=Row.names,
                                      Row.names=NULL)

# Rename the newly added columns
col_id_log2FoldChange = grep("^log2FoldChange$", colnames(output_df_general_details))
colnames(output_df_general_details)[col_id_log2FoldChange] <- "log2FoldChange_MAF_vs_ND"

col_id_lfcSE = grep("^lfcSE$", colnames(output_df_general_details))
colnames(output_df_general_details)[col_id_lfcSE] <- "lfcSE_MAF_vs_ND"



# Now we have to add the individual details for each subgroup comparison in the table
# When merging by rownames, a new column is created, assign this as row names and delete it
output_df_general_details = transform(merge(output_df_general_details,
                                            res_MMSET_vs_ND[,c("log2FoldChange","lfcSE"), drop=FALSE], # So it keeps it as a dataframe
                                            by=0,
                                            all=TRUE),
                                      row.names=Row.names,
                                      Row.names=NULL)

# Rename the newly added columns
col_id_log2FoldChange = grep("^log2FoldChange$", colnames(output_df_general_details))
colnames(output_df_general_details)[col_id_log2FoldChange] <- "log2FoldChange_MMSET_vs_ND"

col_id_lfcSE = grep("^lfcSE$", colnames(output_df_general_details))
colnames(output_df_general_details)[col_id_lfcSE] <- "lfcSE_MMSET_vs_ND"




# Now add a column for each subgroup which specifies the log2FoldChange of that subgroup
# while taking into account padj significant
# 0 -> abs(log2FoldChange) < very_significant_threshold or padj >=0.1
# 1 -> log2FoldChange >= very_significant_threshold
# -1 -> log2FoldChange <= very_significant_threshold

# Start by assigning 0 to everything, then assign the corresponding columns the other values
output_df_general_details[,"DE_HD_vs_ND"] = 0
output_df_general_details[
  ((output_df_general_details$log2FoldChange_HD_vs_ND<=-very_significant_threshold)
   &(!(is.na(output_df_general_details$padj_LRT))
       &(output_df_general_details$padj_LRT<pdaj_threshold))),"DE_HD_vs_ND"] = -1

output_df_general_details[
  ((output_df_general_details$log2FoldChange_HD_vs_ND>=very_significant_threshold)
   &(!(is.na(output_df_general_details$padj_LRT))
     &(output_df_general_details$padj_LRT<pdaj_threshold))),"DE_HD_vs_ND"] = 1



output_df_general_details[,"DE_CCND1_vs_ND"] = 0
output_df_general_details[
  ((output_df_general_details$log2FoldChange_CCND1_vs_ND<=-very_significant_threshold)
   &(!(is.na(output_df_general_details$padj_LRT))
     &(output_df_general_details$padj_LRT<pdaj_threshold))),"DE_CCND1_vs_ND"] = -1

output_df_general_details[
  ((output_df_general_details$log2FoldChange_CCND1_vs_ND>=very_significant_threshold)
   &(!(is.na(output_df_general_details$padj_LRT))
     &(output_df_general_details$padj_LRT<pdaj_threshold))),"DE_CCND1_vs_ND"] = 1


output_df_general_details[,"DE_MAF_vs_ND"] = 0
output_df_general_details[
  ((output_df_general_details$log2FoldChange_MAF_vs_ND<=-very_significant_threshold)
   &(!(is.na(output_df_general_details$padj_LRT))
     &(output_df_general_details$padj_LRT<pdaj_threshold))),"DE_MAF_vs_ND"] = -1

output_df_general_details[
  ((output_df_general_details$log2FoldChange_MAF_vs_ND>=very_significant_threshold)
   &(!(is.na(output_df_general_details$padj_LRT))
     &(output_df_general_details$padj_LRT<pdaj_threshold))),"DE_MAF_vs_ND"] = 1


output_df_general_details[,"DE_MMSET_vs_ND"] = 0
output_df_general_details[
  ((output_df_general_details$log2FoldChange_MMSET_vs_ND<=-very_significant_threshold)
   &(!(is.na(output_df_general_details$padj_LRT))
     &(output_df_general_details$padj_LRT<pdaj_threshold))),"DE_MMSET_vs_ND"] = -1

output_df_general_details[
  ((output_df_general_details$log2FoldChange_MMSET_vs_ND>=very_significant_threshold)
   &(!(is.na(output_df_general_details$padj_LRT))
     &(output_df_general_details$padj_LRT<pdaj_threshold))),"DE_MMSET_vs_ND"] = 1







# Output each of the tables
# setDT modifies the input, make a copy of the dataframe
output_df_all_details_ori <- data.frame(output_df_general_details)

# Remove the variable
remove(output_df_general_details)

# Output the different differential results
# Convert row names into first column
output_df_all_details = setDT(output_df_all_details_ori, keep.rownames = TRUE)[]

# Remove the variable
remove(output_df_all_details_ori)


# Rename the id column
colnames(output_df_all_details)[1] <- first_col_name


# Output the file with all details
write.table(output_df_all_details, 
            gzfile(output_df_all_details_file),
            append = FALSE,
            quote = FALSE,
            sep = "\t",
            eol = "\n",
            na = "NA",
            dec = ".",
            row.names = FALSE, # The row names are already included
            col.names = TRUE)

# Output file with just the significat values
sign_regions <- subset(output_df_all_details,
                       padj_LRT < pdaj_threshold &
                       (DE_HD_vs_ND != 0 
                       | DE_CCND1_vs_ND != 0
                       | DE_MAF_vs_ND != 0
                       | DE_MMSET_vs_ND != 0))

write.table(sign_regions,
            gzfile(output_df_sign_regions_file),
            append = FALSE,
            quote = FALSE,
            sep = "\t",
            eol = "\t",
            na = "NA",
            dec = ".",
            row.names = FALSE,
            col.names = TRUE)


