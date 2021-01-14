library(DESeq2)
library(data.table) # To convert row names into first column

# Requirements:
# -coldata_file can't have a column name called "final_sample_name"
# 
# Inputs:
#   -counts_matrix_file: A file with format:
#         id      RS_1.10 RS_1.11 RS_1.12
#         ENSG00000000003 36.0    745.0   361.0
#         ENSG00000000005 0.0     0.0     1.0
#
#   -coldata_file: A matrix specifying each column of the counts_matrix_file
#         and a patient and CD19 for each sample.
#         NOTE: The samples are row names (no header).
#         NOTE: The order of the rows is in this matrix is the same as the order of the sample columns
#               in the counts_matrix_file
# 
#         A file with format:
#                   patient	CD19
#         Sample1     ND1	plus
#         Sample2     ND2   plus
#         Sample3     ND1   minus
#         
#
#   -group_by_col_id: The column name in coldata_file specifying the replicate ids (to group replicates)
#
#   -output_file: A matrix with results after collapsing replicates from the column group_by_col_id
#                 Note that the replicates will be summed.
#                          17.4 24.11 24.7 26.18,26.20
#     chr1:778512-778762   85    10   40   20
#     chr1:817236-817440   20     7    7   10
#     chr1:869773-870279  141    12   36   19

args <- commandArgs(trailingOnly=TRUE); # Read Arguments from command line
counts_matrix_file = args[1]
coldata_file = args[2]
group_by_col_id = args[3]
output_file = args[4]


####################################### PARAMETERS CHECKING ###################################
# Check all parameters are not empty
if(counts_matrix_file == '' || is.null(counts_matrix_file)){
  stop("ERROR: counts_matrix_file is empty")
}

if(coldata_file == '' || is.null(coldata_file)){
  stop("ERROR: coldata_file is empty")
}

if(group_by_col_id == '' || is.null(group_by_col_id)){
  stop("ERROR: group_by_col_id is empty")
}

if(output_file == '' || is.null(output_file)){
  stop("ERROR: output_file is empty")
}
####################################### PARAMETERS CHECKING ###################################


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


if(!(group_by_col_id %in% colnames(coldata_df)))
{
  stop(paste("ERROR: the column", group_by_col_id, "doesn't exist in coldata"))
}

# Check the column content is the same
if(!setequal(rownames(coldata_df), colnames(input_matrix))) {
     stop("ERROR: Row names of coldata and counts data are not the same\ncounts_matrix colnames")
   }

coldata_df <- coldata_df[colnames(input_matrix),]

# Get the data into DESeqDataSet format, Using design ~ 1 and overriding it later
dds <- DESeqDataSetFromMatrix(countData = input_matrix,
                              colData = coldata_df,
                              design = ~ 1)

# Create a column "final_sample_name" with the sample name
coldata_df$final_sample_name = rownames(coldata_df)

# Collapse any replicates while storing the combined sample name
# The sample names will be stored in a new default column called "runsCollapsed"
# In rows where multiple samples have been collapsed, comma separated sample ids will be stored.
ddsColl <- collapseReplicates(dds, 
                              groupby = dds[[group_by_col_id]],
                              run = coldata_df$final_sample_name)


# Get the dataframe of the new column descriptions
col_desc_df = as.data.frame(colData(ddsColl))


# Create an empty character array of length the new number of samples after collapsing reps
new_colnames = character(length(colnames(assay(ddsColl))))

# Counter to assign each new sample name
counter_pos = 1

# Substitute the replicate names by the new sample names ("runsCollapsed")
for (row in colnames(assay(ddsColl))) {
  new_colnames[counter_pos] = col_desc_df[(col_desc_df[[group_by_col_id]] == row), "runsCollapsed"]
  
  counter_pos = counter_pos + 1
  
}

# Get the dataframe to output
output_df = assay(ddsColl)

# Assign the new colnames
colnames(output_df) <- new_colnames


# Note that if the output_df is to be outputted now, the first row of the column of the features
# will be empty and when loading the table the header with the sample will be shifted one to the left
# We have to correct this

# Convert row names into first column
output_df = setDT(as.data.frame(output_df), keep.rownames = TRUE)[]


# Rename the id column
colnames(output_df)[1] <- first_col_name


# Output the count file
write.table(output_df, gzfile(output_file),
            append = FALSE,
            quote = FALSE,
            sep = "\t",
            eol = "\n",
            na = "NA",
            dec = ".",
            row.names = FALSE, # The row names are already included
            col.names = TRUE)

# Output the col_data file
col_data_file = sub(".tsv.gz", ".col_data.tsv", output_file)
colnames(col_desc_df)[ncol(col_desc_df)] <- "sample_name"
col_desc_df <- col_desc_df[,c(ncol(col_desc_df), 1:(ncol(col_desc_df)-1))]

write.table(col_desc_df, col_data_file,
            append = FALSE,
            quote = FALSE,
            sep = "\t",
            eol = "\n",
            na = "NA",
            dec = ".",
            row.names = FALSE, # The row names are already included
            col.names = TRUE)