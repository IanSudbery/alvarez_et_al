library("org.Hs.eg.db")


# NOTE: If the file has more rows than the chunk size, allow_multiple_gene_ids will only check duplicates
#     inside chunks but not across chunks!!

# Changelog (rerun anything before this date):
# -03/10/2018: Added support for reading large files by chunk.
# -20/07/2018: Corrected including duplicate results (Different Gene symbols for one gene). Until then it is possible
#              that duplicate rows are outputted (everything the same but different gene symbols).
#
# -19/04/2018: Implemented maintaining row order, also implemented option to filter out
#              duplicate genes in the output (depending on if the input genuinely has duplicate genes)
#
#
# -07/02/2018: Eliminate duplicate same rows with same ENSEMBL id but which contain different SYMBOL.
#
#
# Requirements:
#   -ensembl_id_table_file can't contain a column called "unique_row_identifier"
#
# Gets a tsv table with a specified number indicating the column number of Ensembl ID
# (columns specified begin at 0), looks for each Ensembl ID of the row
# to find the corresponding SYMBOL and GENENAME. Adds the columns to the file and outputs it.
# The output table maintains the row order of the input table.
# The output table only contains input ENSEMBL_ID:output ENSEMBL_ID 1:1 mappings
# If the SYMBOL can't be found on a row, it stores the ENSEMBL Id instead
# If various SYMBOLs are found for a row, it stores only the first one.
# Use: EnsemblTable2Symbols ensembl_id_table_file ensembl_col_id header output_table_file
# Inputs:
#    -ensembl_id_table_file: tsv table file with a column of Ensembl ID
#       Note: The "#" character is considered a normal character, not a
#       comment.
#    -ensembl_col_id: number of column containing the ENSEMBL ID (columns specified begin at 0)
#       to be maintained in the output table
#    -header: TRUE if it contains a header, FALSE if it doesn't
#    -allow_multiple_gene_ids: NOTE: If the file has more rows than the chunk size, allow_multiple_gene_ids will only check duplicates
#     inside chunks but not across chunks!! TRUE if the input table has multiple gene ids in the specified column
#       to start with (and they should be allowed in the final table). For example,
#       an input table which contains different bed regions with multiple genes associated.
#       FALSE if the input table doesn't have gene duplicates, and they shouldn't be allowed
#       in the output file.
#
#    -output_table_file: input_table + SYMBOL and GENENAME
#
#    -outputs:
#       writes the output_table_file
#
#


# Get the arguments
args <- commandArgs(trailingOnly=TRUE); # Read Arguments from command line
ensembl_id_table_file = args[1]
ensembl_col_id = args[2]
header = args[3]
allow_multiple_gene_ids = args[4]
output_table_file = gzfile(args[5])

# The number of lines to read at once
chunkSize <- 100000


header_boolean = FALSE

if (toupper(header) == "TRUE"){
  header_boolean = TRUE
}

allow_multiple_gene_ids_boolean = FALSE

if (toupper(allow_multiple_gene_ids) == "TRUE"){
  allow_multiple_gene_ids_boolean = TRUE
}





# Start a loop for each chunk of lines
index <- 0
con <- file(description=ensembl_id_table_file,open="r")
input_table <- read.table(con,
                                       nrows=chunkSize,
                                       sep="\t",
                                       header=header_boolean,
                                       stringsAsFactors=F,
                                       comment.char = "",
                                       check.names=F, # So that "#" in header are not converted to "."
                                       na.strings=c("", "NA")) # Change blanks to NA


repeat {
  index <- index + 1
  print(paste('Processing rows:', index * chunkSize))



  # Process

  # Get the column + 1 (in R columns start from 1)
  ensembl_details_table = AnnotationDbi::select( org.Hs.eg.db, keys=as.character((input_table[[(strtoi(ensembl_col_id)+1)]])), keytype='ENSEMBL', columns=c('ENSEMBL', 'SYMBOL', 'GENENAME'))


  # Get the column name of the Ensembl ID for the input table
  ensembl_id = colnames(input_table)[strtoi(ensembl_col_id)+1]

  # Build a unique row identifier to maintain row order
  input_table$unique_row_identifier <- 1:nrow(input_table)


  # Merge the tables by Ensembl ID
  output_table = merge(x = input_table,
                       y = ensembl_details_table,
                       by.y = "ENSEMBL",
                       by.x = ensembl_id,
                       all.x = TRUE)


  # Get unique rows, we have to do it in terms of the input_table columns and not the output table
  output_table = output_table[!duplicated(output_table[,c('unique_row_identifier')]),]

  # # For rows where the SYMBOL is not found, fill them with the ensembl id
  output_table$SYMBOL[is.na(output_table$SYMBOL)] <- as.character(output_table[,ensembl_id][is.na(output_table$SYMBOL)])

  if(!allow_multiple_gene_ids_boolean){
    # Remove duplicate rows (only for the ensembl_id column)
    output_table = output_table[!duplicated(output_table[ensembl_id]),]
  }

  # Keep the initial ordering
  output_table = output_table[order(output_table$unique_row_identifier), ]

  # Remove the unique_row_identifier column
  drops <- c("unique_row_identifier")
  output_table = output_table[ , !(names(output_table) %in% drops)]


  # If it's the first chunk, output header (if the file contains) overwritting
  if(index == 1){
    # Write the conversion table
    write.table(output_table,
                append = FALSE,
                file=output_table_file,
                sep="\t",
                col.names = header_boolean,
                row.names = F,
                quote = FALSE)

    }  else{ # If it's not the first chunk, don't output header and append
    # Write the conversion table
    write.table(output_table,
                append = TRUE,
                file=output_table_file,
                sep="\t",
                col.names = FALSE,
                row.names = F,
                quote = FALSE)

  }

  # If the chunk is not complete, it means the last chunk is read
  if (nrow(input_table) != chunkSize){
    print('Processed all files!')
    break}


  # If is not the last chunk read the next
  input_table <- read.table(con,
                          nrows=chunkSize,
                          skip=0,
                          sep="\t",
                          header=FALSE,
                          stringsAsFactors=F,
                          comment.char = "",
                          check.names=F, # So that "#" in header are not converted to "."
                          na.strings=c("", "NA")) # Change blanks to NA



}
