library(SummarizedExperiment)
library(MultiAssayExperiment)
library(MOFA)

# For creating the MultiAssay object using:
# https://bioconductor.org/packages/devel/bioc/vignettes/MultiAssayExperiment/inst/doc/MultiAssayExperiment.html#creating-a-multiassayexperiment-object-a-rich-example
# (4 Creating a MultiAssayExperiment object: a rich example)
# For running MOFA using:
# http://htmlpreview.github.io/?https://github.com/bioFAM/MOFA/blob/master/MOFAtools/vignettes/MOFA_example_scMT.html
#
# Starting from an RNA-seq and ATAC-seq raw counts, RLog normalized and batch corrected accounting for subgroup tables
# with corresponding sample names, performs MOFA analysis with:
# Removes RNA-seq counts with 0 variance
# -Not scaling of variables
# -Allow samples that are not profiled in all omics. The model can cope with missing assays.
# -Allow sparse (incomplete data).
# -Gaussian modeled ATAC and RNA data.
# -Default number of factors being half the number of samples taking into account the minimum variance explained.
# -Maximum number of iterations, increased to 10000
# -Tolerance convergence threshold based on change in the evidence lower bound. For an exploratory run you can use a value between 1.0 and 0.1, but for a "final" model we recommend a value of 0.01.
# -Automatically learn the number of factors based on a minimum variance explained criteria. Factors explaining less than 'DropFactorThreshold' fraction of variation in all views will be removed. For example, a value of 0.01 means that factors that explain less than 1% of variance in all views will be discarded. By default this it zero, meaning that all factors are kept unless they explain no variance at all. Here we use a threshold of 1%.
# -Establish a seed for reproducibility.
#
# It outputs a schematic diagram with the input data, number of dimensions and variables.
# Performs MOFA analysis
#
# Inputs:
#   -ATAC_file: ATAC-seq raw counts, RLog normalized and batch corrected accounting for subgroup tables.
#     Samples coinciding with RNA-seq
#     Format:
#                        MM_OTHER_17.5  HD_24.4  HD_24.7
#     chr1:778514-778760      4.652302 3.617844 4.122264
#     chr1:817230-817447      3.456102 3.991655 2.695754
#     chr1:869781-870211      4.567128 3.941549 4.073439
#
#   -RNA_file: RNA-seq raw counts, RLog normalized and batch corrected accounting for subgroup tables.
#     Samples coinciding with ATAC-seq
#     Format:
#                        MM_OTHER_17.5  HD_24.4  HD_24.7
#     chr1:778514-778760      4.652302 3.617844 4.122264
#     chr1:817230-817447      3.456102 3.991655 2.695754
#     chr1:869781-870211      4.567128 3.941549 4.073439
#
#   -data_scheme_file: File with a schematic of the data
#   -seed_param: A seed to train the model

#ATAC_file="/fastdata/mbp15ja/top_5k_var_peaks_genes_atac_rna_MOFA_peaks_pan_MM_vs_ND_rlog_batch_correct_acc_subgroups_no_gender_no_TSS/rlog_batch_correct_read_counts_in_cons_peaks_no_TSS.tsv.gz"
#RNA_file="/fastdata/mbp15ja/top_5k_var_peaks_genes_atac_rna_MOFA_peaks_pan_MM_vs_ND_rlog_batch_correct_acc_subgroups_no_gender/rna_seq_rlog_batch_correct_accounting_for_subgroup_atac_ids_gene_names_only_no_gender.tsv.gz"
#data_scheme_file="/fastdata/mbp15ja/top_5k_var_peaks_genes_atac_rna_MOFA_peaks_pan_MM_vs_ND_rlog_batch_correct_acc_subgroups_no_gender_no_TSS/data_scheme.png"
#output_MOFA_trained_model_file="/fastdata/mbp15ja/top_5k_var_peaks_genes_atac_rna_MOFA_peaks_pan_MM_vs_ND_rlog_batch_correct_acc_subgroups_no_gender_no_TSS/MOFA.hdf5"
#output_MOFA_object_file="/fastdata/mbp15ja/top_5k_var_peaks_genes_atac_rna_MOFA_peaks_pan_MM_vs_ND_rlog_batch_correct_acc_subgroups_no_gender_no_TSS/MOFA.RData"

args <- commandArgs(trailingOnly=TRUE); # Read Arguments from command line
ATAC_file = args[1]
RNA_file = args[2]
data_scheme_file = args[3]
output_MOFA_trained_model_file = args[4]
output_MOFA_object_file = args[5]
seed_param = args[6]


ATAC_table = read.table(ATAC_file,
                           header = TRUE,
                           sep = "\t",
                           quote = "\"'",
                           row.names = 1,
                           dec = ".")


RNA_table = read.table(RNA_file,
                        header = TRUE,
                        sep = "\t",
                        quote = "\"'",
                        row.names = 1,
                        dec = ".")



# Remove the gender chromosomes
ATAC_table = ATAC_table[- grep("chrX", row.names(ATAC_table)),]
ATAC_table = ATAC_table[- grep("chrY", row.names(ATAC_table)),]


# Get the top 5K features by variance in RNA-seq and ATAC-seq

# Get the variance in RNA-seq
RNA_table["var_RNA"] = apply(RNA_table, 1, var)

# Retain only the top 5k rows by variance
RNA_table=RNA_table[order(RNA_table$var_RNA, decreasing = TRUE),]

# Retain only rows with variance different to 0
RNA_table=RNA_table[1:5000,]

# Remove the var_RNA column
RNA_table$var_RNA <- NULL



# Get the variance in ATAC-seq
ATAC_table["var_ATAC"] = apply(ATAC_table, 1, var)

# Retain only the top 5k rows by variance
ATAC_table=ATAC_table[order(ATAC_table$var_ATAC, decreasing = TRUE),]

# Retain only rows with variance different to 0
ATAC_table=ATAC_table[1:5000,]

# Remove the var_ATAC column
ATAC_table$var_ATAC <- NULL






# Prevent conversion of strings to factors
column <- sapply(ATAC_table, is.factor)
ATAC_table[column] <- lapply(ATAC_table[column], as.character)


column <- sapply(RNA_table, is.factor)
RNA_table[column] <- lapply(RNA_table[column], as.character)


# Convert to matrix
ATAC_table_matrix = as.matrix(ATAC_table)

RNA_table_matrix = as.matrix(RNA_table)


# Create a coldata dataframe with rownames the samples
# and a column called condition (PC and MM) and another called subgroup
sample_names = colnames(ATAC_table)

conditions_list = c()
subgroups_list = c()

for (sample_name in sample_names){
  
  # Split by "_"
  subgroup = strsplit(sample_name, "_")
  
  # Get all but the last element and concatenate them with "_"
  subgroup = paste0(head(subgroup[[1]], -1),sep = "_", collapse ="")
  
  # If the last character is "_" remove it
  if (substring(subgroup, nchar(subgroup)) == "_"){
    subgroup =  substring(subgroup, 1, nchar(subgroup)-1) 
  }
  
  # Append to list
  subgroups_list=c(subgroups_list, subgroup)
  
  # Decide on MM or PC depending on the value of subgroup
  if( subgroup == "PC"){
    conditions_list = c(conditions_list, "PC")
  }else{
    conditions_list = c(conditions_list, "MM")
  }
}

# Create the coldata dataframe
patient_data <- data.frame(condition=conditions_list,
                           subgroup=subgroups_list,
                           row.names=sample_names)





# Create coldata files with the sample names without subgroups
sample_names_ATAC = colnames(ATAC_table)

sample_names_processed_ATAC = c()

for (sample_name_ATAC in sample_names_ATAC){
  
  # Split by "_"
  sample_name_processed = strsplit(sample_name_ATAC, "_")
  
  # Get the last element
  sample_name_processed = tail(sample_name_processed[[1]], n=1)
  
  # Append to list
  sample_names_processed_ATAC=c(sample_names_processed_ATAC, sample_name_processed)
}

# Create the coldata dataframe
coldata_ATAC <- data.frame(sample=sample_names_processed_ATAC,
                     row.names=sample_names_ATAC)


summarized_exp_ATAC = SummarizedExperiment(ATAC_table_matrix, colData = coldata_ATAC)






# Create coldata files with the sample names without subgroups
sample_names_RNA = colnames(RNA_table)

sample_names_processed_RNA = c()

for (sample_name_RNA in sample_names_RNA){
  
  # Split by "_"
  sample_name_processed = strsplit(sample_name_RNA, "_")
  
  # Get the last element
  sample_name_processed = tail(sample_name_processed[[1]], n=1)
  
  # Append to list
  sample_names_processed_RNA=c(sample_names_processed_RNA, sample_name_processed)
}

# Create the coldata dataframe
coldata_RNA <- data.frame(sample=sample_names_processed_RNA,
                           row.names=sample_names_RNA)


summarized_exp_RNA = SummarizedExperiment(RNA_table_matrix, colData = coldata_RNA)




# Creation of the MultiAssayExperiment class object, containing various assays
objlist <- list("ATAC" = summarized_exp_ATAC, "RNA" = summarized_exp_RNA)

myMultiAssay <- MultiAssayExperiment(objlist,
                                     patient_data)


# For debugging purposes
#
# View the different assays
# experiments(myMultiAssay)

# View the coldata (additional sample info not included in the assays)
# colData(myMultiAssay)

# View all the samples and the assays they correspond to (from objlist above)
# sampleMap(myMultiAssay)
# The samples will be concatenated in the order they appear in objlist)




# Create the MOFA object
MOFAobject <- createMOFAobject(myMultiAssay)

# Create the file to store the plot
png(file=data_scheme_file, width=750, height=700)

# Obtain an overview of the data
plotTilesData(MOFAobject, colors=c("#31A354","#377EB8"))

dev.off()



# Specify options for the input data:
# -scaleViews: logical indicating whether to scale views to have unit variance. As long as the scale of the different data sets is not too high, this is not required. Default is FALSE.
#
# -removeIncompleteSamples: logical indicating whether to remove samples that are not profiled in all omics. The model can cope with missing assays, so this option is not required. Default is FALSE
DataOptions <- getDefaultDataOptions()
#DataOptions$scaleViews=TRUE # We want all data in the same footing


# Define model options
# Next, we define model options. The most important are:
#
# -numFactors: number of factors (default is 0.5 times the number of samples). By default, the model will only remove a factor if it explains exactly zero variance in the data. You can increase this threshold on minimum variance explained by setting  TrainOptions$dropFactorThreshold to a value higher than zero.
# -likelihoods: likelihood for each view. Usually we recommend gaussian for continuous data, bernoulli for binary data and poisson for count data. By default, the model tries to guess it from the data.
# -sparsity: do you want to use sparsity? This makes the interpretation easier so it is recommended (Default is TRUE).
# Make sure the data type selected for each assay is Gaussian
ModelOptions <- getDefaultModelOptions(MOFAobject)

data_type_assay=c("gaussian", "gaussian")
names(data_type_assay)=c("ATAC", "RNA")

ModelOptions$likelihood=data_type_assay



#Training options, most important are:
# -maxiter: maximum number of iterations. Ideally set it large enough and use the convergence criteria TrainOptions$tolerance.
# -tolerance: convergence threshold based on change in the evidence lower bound. For an exploratory run you can use a value between 1.0 and 0.1, but for a "final" model we recommend a value of 0.01.
# -DropFactorThreshold: hyperparameter to automatically learn the number of factors based on a minimum variance explained criteria. Factors explaining less than 'DropFactorThreshold' fraction of variation in all views will be removed. For example, a value of 0.01 means that factors that explain less than 1% of variance in all views will be discarded. By default this it zero, meaning that all factors are kept unless they explain no variance at all. Here we use a threshold of 1%.
TrainOptions <- getDefaultTrainOptions()

# Establish a seed from the parameters
TrainOptions$seed <- seed_param

# Automatically drop factors that explain less than 1% of variance in all omics
TrainOptions$DropFactorThreshold <- 0.01

TrainOptions$maxiter=10000

TrainOptions$tolerance=0.01


TrainOptions$verbose=TRUE




# Prepare MOFA
MOFAobject <- prepareMOFA(
  MOFAobject, 
  DataOptions = DataOptions,
  ModelOptions = ModelOptions,
  TrainOptions = TrainOptions
)

MOFAobject <- runMOFA(MOFAobject, outfile=output_MOFA_trained_model_file)

# Save the MOFA object
save(MOFAobject,
     file = output_MOFA_object_file,
     compress = TRUE)
