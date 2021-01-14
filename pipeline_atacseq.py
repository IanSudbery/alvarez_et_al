from ruffus import follows, transform, add_inputs, mkdir, regex, formatter, merge, subdivide, files, collate, suffix
from ruffus.combinatorics import product
#from ruffus.combinatorics import *

import sys
import os
import sqlite3
import cgatcore.experiment as E
import cgatcore.pipeline as P
import pipelineAtacseq
import tempfile
import re
import cgatcore.iotools as IOTools
import cgatpipelines.tasks.rnaseq as rnaseq
import pandas
#sys.path.insert(0, "/home/mbp15ja/dev/AuxiliaryPrograms/Segmentations/")
#import compareSegmentations as compseg
#sys.path.insert(0, "/home/mbp15ja/dev/AuxiliaryPrograms/logParsers/")
#import logParser
#sys.path.insert(0, "/home/mbp15ja/dev/AuxiliaryPrograms/StringOperations/")
#import StringOperations
#sys.path.insert(0, '/home/mbp15ja/dev/AuxiliaryPrograms/Clusterings/')
#import clusterings
#sys.path.insert(0, '/home/mbp15ja/dev/AuxiliaryPrograms/File_operations/')
#import file_operations


import matplotlib as mpl
mpl.use('Agg') # So that matplotlib.pyplot doesn't try to find the X-server and give an error
import matplotlib.pyplot as plt

import math

PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.yml"])


#------------------------------------------------------------------------------
@follows(mkdir("stats.dir"))
@transform("*.bam",
           regex("(.+).bam"),
           r"stats.dir/\1.after_mapping.tsv")
def getInitialMappingStats(infile, outfile):
    ''' Gets the initial mapping rate in terms of total reads '''
    
    # Get the temporal dir specified
    tmp_dir = PARAMS["general_temporal_dir"]
    
    
    # The mapping file is sorted by coordinate, first sort by readname
    temp_file = P.snip(outfile, ".tsv") + ".bam"
    
    log_file = P.snip(outfile, ".tsv") + ".log"
    
    # Samtools creates temporary files with a certain prefix, create a temporal directory name
    samtools_temp_dir = tempfile.mkdtemp(dir=tmp_dir)
    
    samtools_temp_file = os.path.join(samtools_temp_dir, "temp")
    
    disk_space_file = P.snip(outfile, ".tsv") + ".txt"
        
    statement = ''' mkdir -p %(samtools_temp_dir)s &&
                    samtools sort -n -o %(temp_file)s 
                                    -T %(samtools_temp_file)s %(infile)s 
                                    2> %(log_file)s 
    '''
    
    
    # Execute the statement
    P.run(statement)
    
    # Get the stats
    pipelineAtacseq.getMappedUnmappedReads(temp_file, 
                       outfile, 
                       submit=True)
    
    # Remove the temporal file
    statement = '''rm %(temp_file)s; 
    '''
    
    # Execute the statement
    P.run(statement)


#-----------------------------------------------------------------------------
@follows(mkdir("first_filtering.dir"))
@transform("*.bam",
           regex("(.+).bam"),
           r"first_filtering.dir/\1.bam")
def filterOutIncorrectPairsAndExcessiveMultimappers(infile, outfile):
    
    '''Assuming a starting compressed coordinate sorted bam file. 
    Remove  unmapped, mate unmapped and reads failing platform. Keep only
    properly mapped pairs. Sort by name and filter out proper mapped pairs with
    more than the defined number of proper pair alignments'''
     
    allowed_multimappers = PARAMS["filtering_allowed_multimapper_proper_pairs"]
    
    
    log_file = P.snip(outfile, ".bam") + ".log"
    
    first_filtering_bam_output = P.snip(outfile, ".bam") + "_proper_pairs.bam"
       
    # Get the temporal dir specified
    tmp_dir = PARAMS["general_temporal_dir"]
    
    # Temp file: We create a temp file to make sure the whole process goes well
    # before the actual outfile is created
    temp_file = P.snip(outfile, ".bam") + "_temp.bam"
      
    # Samtools creates temporary files with a certain prefix
    samtools_temp_file = (tempfile.NamedTemporaryFile(dir=tmp_dir, delete=False)).name
    
    scripts_dir = os.path.join(os.path.dirname(__file__), "scripts")


    statement = '''samtools view -F 524 -f 2 -u %(infile)s | 
                   samtools sort -n - -o %(first_filtering_bam_output)s 
                                  -T %(samtools_temp_file)s 2> %(log_file)s &&
                   samtools view -h %(first_filtering_bam_output)s | 
                %(scripts_dir)s/assign_multimappers.py -k %(allowed_multimappers)s --paired-end | 
    samtools view -bS - -o %(temp_file)s 2>> %(log_file)s &&    
    mv %(temp_file)s %(outfile)s &&    
    rm %(first_filtering_bam_output)s
    '''

    P.run(statement)


#-------------------------------------------------------------------------
@follows(mkdir("stats.dir"))
@transform(filterOutIncorrectPairsAndExcessiveMultimappers,
           regex(".+/(.+).bam"),
           r"stats.dir/\1.after_first_filter.tsv")
def getFirstFilteringStats(infile, outfile):
    ''' Gets the mapping rate in terms of total reads after the first filtering '''
    
    # Get the temporal dir specified
    tmp_dir = PARAMS["general_temporal_dir"]
    
    
    # The mapping file is sorted by coordinate, first sort by readname
    temp_file = P.snip(outfile, ".tsv") + ".bam"
    
    log_file = P.snip(outfile, ".tsv") + ".log"
    
    # Samtools creates temporary files with a certain prefix, create a temporal directory name
    samtools_temp_dir = tempfile.mkdtemp(dir=tmp_dir)
    
    samtools_temp_file = os.path.join(samtools_temp_dir, "temp")
    
    disk_space_file = P.snip(outfile, ".tsv") + ".txt"
        
    statement = ''' mkdir -p %(samtools_temp_dir)s &&
      samtools sort -n -o %(temp_file)s -T %(samtools_temp_file)s %(infile)s 2> %(log_file)s;
    '''
    # Execute the statement
    P.run(statement)

    
    # Get the stats
    pipelineAtacseq.getMappedUnmappedReads(temp_file, 
                       outfile, 
                       submit=True,
                       job_memory="6G")
    
    # Remove the temporal file
    statement = '''rm %(temp_file)s; 
    '''
    
    # Execute the statement
    P.run(statement)


#----------------------------------------------------------------------------------------
@follows(mkdir("second_filtering.dir"))
@transform(filterOutIncorrectPairsAndExcessiveMultimappers,
           regex(".+/(.+).bam"),
           r"second_filtering.dir/\1.bam")
def filterOutOrphanReadsAndDifferentChrPairs(infile, outfile):
    
    ''' Remove orphan reads (pair was removed) and read pairs mapping to different 
    chromosomes and read pairs which are "facing against one another" with no overlap. 
    Obtain position sorted BAM. Assumes a starting read name sorted BAM file.
    '''
    
    # Get the sample name
    sample_name , _ = os.path.splitext(os.path.basename(outfile))
    
    # Get samples details table
    sample_details = PARAMS["samples_details_table"]
    
    # Get trimmings in the 5' ends done previously (for example in qc).
    five_prime_trim = 0 # pipelineAtacseq.getSampleQCShift(sample_name, sample_details)
    
    integer_five_prime_correction = 0
    
    # To avoid putting "--"
    # Correction is going to be -correction on the start of the + strand
    # Correction is going to be +correction on the end of the - strand
    try:
        integer_five_prime_correction = int(five_prime_trim)  
    except ValueError:
        raise Exception("Five prime trimming argument needs to be an integer.") 
    
    # String with the correction to apply (Eg. "- 2", "+ 5")
    positive_strand_correction = ""
    negative_strand_correction = ""
    
    if integer_five_prime_correction < 0:
        positive_strand_correction = "+ "+str(abs(integer_five_prime_correction))
        negative_strand_correction = "- "+str(abs(integer_five_prime_correction))
    elif integer_five_prime_correction > 0:
        positive_strand_correction = "- "+str(abs(integer_five_prime_correction))
        negative_strand_correction = "+ "+str(abs(integer_five_prime_correction))
    
    # 0 Case: no correction, empty string
    log_file = P.snip(outfile, ".bam") + ".log"  
    
    # Get the temporal dir specified
    tmp_dir = PARAMS["general_temporal_dir"]
    
    # Temp file: We create a temp file to make sure the whole process goes well
    # before the actual outfile is created
    temp_file = P.snip(outfile, ".bam") + "_temp.bam"
      
    # Samtools creates temporary files with a certain prefix
    samtools_temp_file = (tempfile.NamedTemporaryFile(dir=tmp_dir, delete=False)).name 
    first_filtering_bam_output = P.snip(outfile, ".bam") + "_proper_pairs.bam"
    
    # Fixmate seems to fix read pairs which are in different chromosomes but not
    # read pairs "facing away" from each other 
    # Intermediate file outputted to process these cases
    read_name_processing_input_bam = (
        tempfile.NamedTemporaryFile(dir=tmp_dir, delete=False)).name + ".bam"

    read_name_processing_problematic_reads = (
        tempfile.NamedTemporaryFile(dir=tmp_dir, delete=False)).name

    read_name_processing_output_bam = (
        tempfile.NamedTemporaryFile(dir=tmp_dir, delete=False)).name + ".bam"
    
    # Reads are already read name sorted, so fixmate can be run
    # Next filter the "wrong pairs" and get the read names of read pairs with
    # each pair on a different chror "facing away" from each other with no
    # overlap:
    # <----------                    <--------
    #           --------> Allowed             --------> Not allowed
    # Because the enzyme introduces 4 bp in the start of + strand and 5bp in the
    # start of the - strand(end coordinate in bedpe file) a minumum overlap of 
    # 10
    # The 5' ends of the reads is previously extended by any cuts performed in qc
    # which are indicated.

    # Get a bed with both ends of a read pair in each line
    # and extract the problematic read names from there.
    # Then filter them out with picard FilterSamReads
    statement = '''samtools fixmate -r %(infile)s %(first_filtering_bam_output)s 
                                    2> %(log_file)s &&
                   
                   samtools view -F 1804 -f 2 -u %(first_filtering_bam_output)s 
                                 -o %(read_name_processing_input_bam)s 
                                2>> %(log_file)s &&    
                   
                   bedtools bamtobed -bedpe -i %(read_name_processing_input_bam)s 
                   | awk '($1!=$4 || 
                           ($10=="+" && $9=="-" && 
                            ($3-1%(negative_strand_correction)s-5)<($5%(positive_strand_correction)s+4))) 
                          {printf ("%%s\\n", $7)}'
                           > %(read_name_processing_problematic_reads)s 
                           2>> %(log_file)s &&
    
    if [ -s %(read_name_processing_problematic_reads)s ]; 
        then 
        FilterSamReads I=%(read_name_processing_input_bam)s 
                       O=%(read_name_processing_output_bam)s 
                       READ_LIST_FILE=%(read_name_processing_problematic_reads)s 
                            FILTER=excludeReadList 2>> %(log_file)s;
    else 
        ln -s %(read_name_processing_input_bam)s %(read_name_processing_output_bam)s;
    fi &&
     
    samtools sort %(read_name_processing_output_bam)s 
                  -o %(temp_file)s -T %(samtools_temp_file)s 2>> %(log_file)s &&

    mv %(temp_file)s %(outfile)s &&
    
    rm %(first_filtering_bam_output)s 
       %(read_name_processing_input_bam)s 
       %(read_name_processing_problematic_reads)s 
       %(read_name_processing_output_bam)s;
    '''
    
    job_memory="4G"

    P.run(statement)


#------------------------------------------------------------------------------
# Assumes the files are coordinate sorted
@follows(mkdir("dedupped.dir"))
@transform(filterOutOrphanReadsAndDifferentChrPairs,
           regex(".+/(.+).bam"),
           r"dedupped.dir/\1.bam")
def markDuplicates(infile, outfile):
    
    ''' Use picard to mark duplicates in BAM files (not deleted).
    The files are assumed to be coordinate sorted'''
    
    # Used to be 5G
    job_memory = "8G"
    
    # Get the temporal dir specified
    tmp_dir = PARAMS["general_temporal_dir"]
    
    # Temporal dir, to prevent the "No space left on device" error
    temp_dir = tempfile.mkdtemp(dir=tmp_dir)
    
    # Temp file: We create a temp file to make sure the whole process goes well
    # before the actual outfile is created
    temp_file = P.snip(outfile, ".bam") + "_temp.bam"

    statement = ''' MarkDuplicates
                     ASSUME_SORTED=True 
                     INPUT=%(infile)s
                     OUTPUT=%(temp_file)s
                     VALIDATION_STRINGENCY=LENIENT
                     METRICS_FILE=%(outfile)s.metrics
                     REMOVE_DUPLICATES=false
                     TMP_DIR=%(temp_dir)s
                   > %(outfile)s.log &&
                   
                   mv %(temp_file)s %(outfile)s '''

    P.run(statement)


#------------------------------------------------------------------------------
@follows(mkdir("stats.dir"))
@transform(markDuplicates,
           regex(".+/(.+).bam"),
           r"stats.dir/\1.after_marking_dups.tsv")
def getPostDuplicationStats(infile, outfile):
     
    ''' Assuming multimapping is allowed (multiple best alignments can occur)
    Sort the reads by readname, make filterings and get the number
    unique pair mappings:
    1) Correctly mapped pairs and primary alignments only.
    2) Correctly mapped pairs and primary or secondary alignments.
    3) Correctly mapped pairs and secondary alignments only.
    get initial statistics on the reads '''
    
    
    # Get the temporal dir specified
    tmp_dir = PARAMS["general_temporal_dir"] 
    sorted_bam = P.snip(outfile, ".tsv") + "_sorted.bam"
    log_file = P.snip(outfile, ".tsv") + ".log"
    bam_outfile_sec = P.snip(outfile, ".tsv") + "_sec.bam"
    bam_outfile_primary = P.snip(outfile, ".tsv") + "_prim.bam"
    
    
    # Samtools creates temporary files with a certain prefix
    samtools_temp_file = (
        tempfile.NamedTemporaryFile(dir=tmp_dir, delete=False)).name
     
    # First sort the bamfile
    # Then get only the primary alignments 
    statement = '''samtools sort -n 
                            -o %(sorted_bam)s 
                            -T %(samtools_temp_file)s 
                            %(infile)s 
                            2> %(log_file)s &&
    
                   samtools view -h -F 256 
                            %(sorted_bam)s 
                            -o %(bam_outfile_primary)s 
                            -U %(bam_outfile_sec)s 
                            2>> %(log_file)s;
    '''
 
    P.run(statement)
     
    # Now see the mapped pairs and PCR duplicates in each of the 3 files
    primary_stats_file = P.snip(outfile, ".tsv") + "_primary.tsv"
    secondary_stats_file = P.snip(outfile, ".tsv") + "_secondary.tsv"
      
    # Where only primary alignments exist (1 read = 1 alignment) 
    pipelineAtacseq.getUniquelyMappedPairsNoMultimapping(bam_outfile_primary, 
                                                  primary_stats_file, 
                                                  submit=True, 
                                                  job_memory="4G")
    
    # Where multiple alignments can exist (primary + secondary)
    pipelineAtacseq.getCorrectReadPairs(sorted_bam,
                                        outfile, 
                                        submit=True,
                                        job_memory="4G")
      
    # Where multiple alignments can exist (secondary)
    pipelineAtacseq.getCorrectReadPairs(bam_outfile_sec,
                                        secondary_stats_file, 
                                        submit=True,
                                        job_memory="4G")


#------------------------------------------------------------------------------
@subdivide(markDuplicates,
           regex("(.+)/(.+).bam"),
           [(r"\1/\2_pos_sorted.bam"), 
            (r"\1/\2_read_name_sorted.bam")],
            r"\1/\2")
def deduplicate(infile, outfiles, sample):
    '''Remove duplicates, create final name sorted BAM. Assumes a starting position sorted BAM'''
    
    # Get both outfiles
    position_sorted_bam = outfiles[0]
    read_name_sorted_bam = outfiles[1]
    
    log_file = sample + ".log"
    
    # Get the temporal dir specified
    tmp_dir = PARAMS["general_temporal_dir"]
    
    # Create end temp file and intermediate temp file for position sorted and name sorted
    
    # Temp file: We create a temp file to make sure the whole process goes well
    # before the actual outfile is created
    temp_file_pos_sorted_bam = P.snip(position_sorted_bam, ".bam") + "_temp.bam"
      
    # Samtools creates temporary files with a certain prefix
    samtools_pos_sorted_temp_file = (tempfile.NamedTemporaryFile(dir=tmp_dir, delete=False)).name
    
    temp_file_name_sorted_bam = P.snip(read_name_sorted_bam, ".bam") + "_temp.bam"
      
    # Samtools creates temporary files with a certain prefix
    samtools_name_sorted_temp_file = (tempfile.NamedTemporaryFile(dir=tmp_dir, delete=False)).name
    
    statement = '''samtools view -F 1804 
                                 -f 2 
                                 -b %(infile)s 
                                 -o %(temp_file_pos_sorted_bam)s 
                                 2> %(log_file)s &&
        
                    samtools sort -n %(temp_file_pos_sorted_bam)s 
                                 -o %(temp_file_name_sorted_bam)s 
                                 -T %(samtools_name_sorted_temp_file)s 
                                 2>> %(log_file)s &&
    
                   mv %(temp_file_pos_sorted_bam)s %(position_sorted_bam)s &&
    
                   mv %(temp_file_name_sorted_bam)s %(read_name_sorted_bam)s; 
    
    '''
    
    P.run(statement)


#------------------------------------------------------------------------------
@follows(mkdir("library_complexity.dir"))    
@transform(markDuplicates,
           regex(".+/(.+).bam"),
           r"library_complexity.dir/\1.pbc.qc")
def calculateLibrarycomplexity(infile, outfile):
    '''Calculates library complexity'''
      
    # outfile temp file to ensure complete execution before writing outfile
    temp_outfile = P.snip(outfile, ".pbc.qc") + "_temp.pbc.qc"
    
    # outfile temp file to ensure complete execution before writing outfile
    temp_header_outfile = P.snip(outfile, ".pbc.qc") + ".header"
      
    # Get the temporal dir specified
    tmp_dir = PARAMS["general_temporal_dir"]
      
    # Samtools creates temporary files with a certain prefix
    samtools_name_sorted_temp_file = (
        tempfile.NamedTemporaryFile(dir=tmp_dir, delete=False)).name
      
    temp_file_name_sorted_bam = P.snip(infile, ".bam") + "_temp.bam"
      
    log_file = P.snip(outfile, ".pbc.qc") + ".log"
      
      
    # 1) Turns the read name sorted file to bed with both mapped segments in the pair
    # 2) Gets the fields:
    #    -beginning of the most upstream segment
    #    -end of most downstream segment
    #    -mapping strand of each segment.
    # 3) Removes the mitochondrial chromosome regions
    # 4) Sees any repeated regions from 2), counting the times each region appears.
    # 5) Performs calculations for distinct reads, total reads and ratios.
    # 6) Creates a header file and appends the figures calculated
    header_file = "TotalReadPairs\\tDistinctReadPairs\\tOneReadPair\\tTwoReadPairs\\tNRF=Distinct/Total\\tPBC1=OnePair/Distinct\\tPBC2=OnePair/TwoPair" 
    statement = '''samtools sort -n %(infile)s 
                                 -o %(temp_file_name_sorted_bam)s 
                                 -T %(samtools_name_sorted_temp_file)s 
                                 2>> %(log_file)s &&
          
                    bedtools bamtobed -bedpe -i %(temp_file_name_sorted_bam)s
                    | awk 'BEGIN{OFS="\\t"} 
                           (($1==$4) && ($2==$5) && ($3==$6))
                           {$9="+";$10="-"} 
                           {print $0}' 
                    | awk 'BEGIN{OFS="\\t"}{print $1,$2,$4,$6,$9,$10}'
                    | grep -v 'chrM' 
                    | sort 
                    | uniq -c 
                    | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} 
                           ($1==1)
                           {m1=m1+1} 
                           ($1==2){m2=m2+1} 
                           {m0=m0+1} 
                           {mt=mt+$1} 
                           END{printf "%%d\\t%%d\\t%%d\\t%%d\\t%%f\\t%%f\\t%%f\\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' 
                           > %(temp_outfile)s &&
       
                    rm %(temp_file_name_sorted_bam)s %(samtools_name_sorted_temp_file)s* &&
                    echo -e '%(header_file)s' > %(temp_header_outfile)s &&      
                    cat %(temp_outfile)s >> %(temp_header_outfile)s &&
                    rm %(temp_outfile)s &&      
                    mv %(temp_header_outfile)s %(outfile)s;
      
    '''
  
    P.run(statement)


#------------------------------------------------------------------------------
@follows(mkdir("flagstats.dir"), deduplicate)
@transform(deduplicate,
           formatter(".+/(?P<SAMPLE>.+)_pos_sorted\.bam"),
           "flagstats.dir/{SAMPLE[0]}.flagstats")
def index(infile, outfile):    
    '''Index final position sorted BAM, get flag stats.'''
    
    log_file = P.snip(outfile, ".flagstats") + ".log"
    
    # Temp file: We create a temp file to make sure the whole process goes well
    # before the actual outfile is created
    temp_file = P.snip(outfile, ".flagstats") + "_temp.flagstats"
      
    
    # Index Final BAM file
    statement = '''samtools index %(infile)s 2> %(log_file)s &&
                   samtools flagstat %(infile)s > %(temp_file)s &&    
                   mv %(temp_file)s %(outfile)s;
    
    '''
    
    P.run(statement)    
    
    
#------------------------------------------------------------------------------   
@follows(mkdir("tag_align.dir"), index, deduplicate)
@transform(deduplicate,
           formatter(".+/(?P<SAMPLE>.+)_pos_sorted\.bam"),
           "tag_align.dir/{SAMPLE[0]}.PE2SE.tagAlign.gz")
def createTagAlign(infile, outfile):
    '''creates tagAlign file (virtual single end) with (BED 3+3 format)'''
    
    log_file = P.snip(outfile, ".PE2SE.tagAlign.gz") + ".log"
    
    # Temp file: We create a temp file to make sure the whole process goes well
    # before the actual outfile is created
    temp_file = P.snip(outfile, ".PE2SE.tagAlign.gz") + "_temp.PE2SE.tagAlign.gz"
    
    # Create virtual SE file containing both read pairs
    statement = '''bedtools bamtobed -i %(infile)s | 
    awk 'BEGIN{OFS="\\t"}{$4="N";$5="1000";print $0}' | 
    gzip -c > %(temp_file)s 2> %(log_file)s &&
        
    mv %(temp_file)s %(outfile)s;
    '''
    
    P.run(statement)


#------------------------------------------------------------------
@follows(mkdir("final_tag_align.dir"), index)
@transform(createTagAlign,
           regex(".+/(.+?).PE2SE.tagAlign.gz"),
           r"final_tag_align.dir/\1.PE2SE.tagAlign.gz")    
def excludeUnwantedContigsPE2SE(infile, outfile):
    '''Exclude the contigs indicated, performs partial matching for each'''
    
    excluded_chrs = PARAMS["filtering_contigs_to_remove"]
    
    # Temp file: We create a temp file to make sure the whole process goes well
    # before the actual outfile is created
    temp_file = P.snip(outfile, ".PE2SE.tagAlign.gz") + "_temp.PE2SE.tagAlign.gz"
    
    statement = pipelineAtacseq.createExcludingChrFromBedStatement(infile, 
                                                                   excluded_chrs,
                                                                   temp_file)
    
    statement += '''&&
                    mv %(temp_file)s %(outfile)s'''
    
    P.run(statement)


#----------------------------------------------------------------------
@follows(mkdir("final_tag_align.dir"))
@transform(excludeUnwantedContigsPE2SE,
           regex(".+/(.+?).PE2SE.tagAlign.gz"),
           r"final_tag_align.dir/\1.PE2SE.tn5_shifted.tagAlign.gz")
def shiftTagAlign(infile, outfile):
    '''Shifts tag aligns by the TN5 sites and any 5' trimming from qc'''
    
    # Eliminate .PE2SE.tn5_shifted.tagAlign.gz from the sample name
    sample_name = re.sub('\.PE2SE\.tagAlign\.gz$', '', os.path.basename(infile))
    
    # Get samples details table
    sample_details = PARAMS["samples_details_table"]
    
    # Get trimmings in the 5' ends done previously (for example in qc).
    five_prime_trim = 0 # pipelineAtacseq.getSampleQCShift(sample_name, sample_details)
    
    integer_five_prime_correction = 0
    
    # To avoid putting "--"
    # Correction is going to be -correction on the start of the + strand
    # Correction is going to be +correction on the end of the - strand
    try:
        integer_five_prime_correction = int(five_prime_trim)
    except ValueError:   
        raise Exception("Five prime trimming argument needs to be an integer.") 
    
    # String with the correction to apply (Eg. "- 2", "+ 5")
    positive_strand_correction = ""
    negative_strand_correction = ""
    
    if integer_five_prime_correction < 0:
        positive_strand_correction = "+ "+str(abs(integer_five_prime_correction))
        negative_strand_correction = "- "+str(abs(integer_five_prime_correction))
    elif integer_five_prime_correction > 0:
        positive_strand_correction = "- "+str(abs(integer_five_prime_correction))
        negative_strand_correction = "+ "+str(abs(integer_five_prime_correction))
    
    # 0 Case: no correction, empty string
    
    # Get the contigs
    contigs = PARAMS["contigs"]
    log_file = P.snip(outfile, ".PE2SE.tn5_shifted.tagAlign.gz") + ".tn5_shifted.log"
    
    # Temp file: We create a temp file to make sure the whole process goes well
    # before the actual outfile is created
    temp_file = P.snip(outfile, ".PE2SE.tn5_shifted.tagAlign.gz") + "_temp.tn5_shifted.tagAlign.gz"
    
    # Shift the beginning of the elements in the + strand (where the enzume cuts) by +4
    # Shift the end of the elements in the - strand (where the enzume cuts) by -5
    # Apply qc corrections too.
    statement = '''zcat %(infile)s 
                   | awk -F $'\\t' 
                     'BEGIN {OFS = FS}
                     { if ($6 == "+") {
                         $2 = $2 + 4 %(positive_strand_correction)s
                       } else if ($6 == "-") {
                         $3 = $3 - 5 %(negative_strand_correction)s
                       } 
                       print $0}'
                   | gzip -c 
                   > %(temp_file)s 
                   2> %(log_file)s; 
    '''
        
    P.run(statement)
    
    log_file_correction = P.snip(outfile, ".PE2SE.tn5_shifted.tagAlign.gz") + \
                           ".tn5_shifted_slop_correction.log"
    
    # Temp file: We create a temp file to make sure the whole process goes well
    # before the actual outfile is created
    temp_file2 = P.snip(outfile, ".PE2SE.tn5_shifted.tagAlign.gz") + "_temp_correction.tn5_shifted.tagAlign.gz"
    
    # Check that the slop does not surpass the chromosome edges and that the starts are < than ends
    pipelineAtacseq.correctSlopChromosomeEdges(temp_file, 
                                                contigs,
                                                temp_file2,
                                                log_file_correction,
                                                submit=True,
                                                job_memory="4G")
    
    # Remove the temp file
    statement = ''' rm %(temp_file)s &&
    mv %(temp_file2)s %(outfile)s '''
    
    P.run(statement)


#------------------------------------------------------------------------------
@follows(mkdir("filtered_tag_align.dir"))
@transform(shiftTagAlign,
           formatter(".+\.dir/(?P<SAMPLE>(?!pooled_[control|treatment]).+)\.PE2SE\.tn5_shifted\.tagAlign\.gz"),
           "filtered_tag_align.dir/{SAMPLE[0]}.single.end.shifted.filtered.tagAlign.gz")
def filterShiftTagAlign(infile, outfile):
    ''' Filters out regions of low mappability and excessive mappability in the shifted single ends'''
    
    temp_file = P.snip(outfile, ".gz") + "_temp.gz"
    
    excluded_beds = PARAMS["filtering_bed_exclusions"]
    
    statement = pipelineAtacseq.createExcludingBedsFromBedStatement(infile, 
                                                                    excluded_beds, 
                                                                    temp_file)
    
    statement += ''' && mv %(temp_file)s %(outfile)s'''

    P.run(statement)


#--------------------------------------------------------------------------------------
@follows(mkdir("filtered_tag_align_count_balanced.dir"))
@transform(filterShiftTagAlign,
           regex("filtered_tag_align.dir/(.+).single.end.shifted.filtered.tagAlign.gz"),
           r"filtered_tag_align_count_balanced.dir/\1.tsv",
           r"\1")
def calculateNumberOfSingleEnds(infile, outfile, sample):
    '''Get the number of single end reads entering the peak calling'''
        
    statement = '''echo %(sample)s,`zcat %(infile)s | wc -l` > %(outfile)s
                    '''
    
    P.run(statement)


#----------------------------------------------------------------------------------------
@merge(calculateNumberOfSingleEnds, "filtered_tag_align_count_balanced.dir/reads_per_sample.csv")
def mergeSingleEndsCount(infiles, outfile):
    ''' Merge together all the SE counts into a single table '''

    infiles = " ".join(infiles)
    statement = ''' echo sample,n_se > %(outfile)s &&
                    cat %(infiles)s >> %(outfile)s '''
    P.run(statement)


#----------------------------------------------------------------------------------------
@follows(mergeSingleEndsCount)
def get_single_ends():
    pass

##########################################################################################
# Generate Peak sets

@follows(mkdir("pan_balanced_samples.dir"),
         mkdir("subtype_balanced_samples.dir"))
@subdivide(filterShiftTagAlign,
           regex("filtered_tag_align.dir/(.+).bowtie2.single.end.shifted.filtered.tagAlign.gz"),
           add_inputs(mergeSingleEndsCount, "samples.tsv"),
           [r"pan_balanced_samples.dir/\1.single.end.shifted.filtered.tagAlign.gz",
            r"subtype_balanced_samples.dir/\1.single.end.shifted.filtered.tagAlign.gz"],
            r"\1")
def get_balanced_sample_filtered_shifted_SE(infiles, outfiles, sample_name):
    ''' Samples a random sample of reads from each sample. Each sample contains the same
    number of reads: the minimum sample single ends. Returns a file with the same sorting
    as the input'''
    
    # Get the temporal dir specified
    tmp_dir = PARAMS["general_temporal_dir"]
    
    # Separate the infiles
    bed_tags, read_count_file, sample_info = infiles
    
    sample_info = pandas.read_csv(sample_info, sep="\t")
    read_counts = pandas.read_csv(read_count_file, sep=",")

    read_counts["sample"] = read_counts["sample"].str.replace(".bowtie2","")
    sample_info = sample_info.merge(read_counts, left_on="ATAC.sample.code", right_on="sample")
    sample_info = sample_info.set_index("sample")

    #min_per_status = sample_info.groupby("MM.ND").n_se.min()
    #min_per_subtype = sample_info.groupby("Subgroup").n_se.min()

    #get_reads_status = min_per_status[sample_info["MM.ND"][sample_name]]
    #get_reads_subtype = min_per_subtype[sample_info["Subgroup"][sample_name]]

    output_reads = sample_info[sample_info["MM.ND"].isin(["MM","ND"])].n_se.min()
    if sample_name in ["A26.20", "A26.18"]:
        output_reads=output_reads/2
           
    # Extract the file and generate the sample
    statement_template = '''zcat %(bed_tags)s > %(temp_extract_bed)s &&
                    
                    sample -o --preserve-order -d 51 -k %(output_reads)s %(temp_extract_bed)s | gzip > %(outfile)s &&
                    
                    rm %(temp_extract_bed)s 
                    
                    '''
    
    statements = [] 

    for outfile in outfiles:
        temp_extract_bed = (tempfile.NamedTemporaryFile(dir=tmp_dir, delete=False)).name

        statements.append(statement_template % locals())

    # Samples with number of single ends > 100M
    job_memory = "8G"
    
    P.run(statements)

#----------------------------------------------------------------------------------------
def get_files_for_pooling():
    '''This generator makes the tuples for @files that will be used to pool the balanced samples
    into pan-MM/ND and subtype read pools for peak calling'''

    sample_info = pandas.read_csv("samples.tsv", sep="\t")

    for pan in ["MM", "ND"]:
        samples = sample_info[sample_info["MM.ND"] == pan]["ATAC.sample.code"]
        dirname = "pan_balanced_samples.dir"
        outfile = os.path.join(dirname, "%s_pooled_pan.single.end.shifted.filtered.tagAlign.gz" 
                                        % pan)
        infiles = [os.path.join(dirname, "%s.single.end.shifted.filtered.tagAlign.gz") % s 
                   for s in samples]
        yield (infiles, outfile)

    for subtype in set(sample_info["Subgroup"]):
        if subtype == "UNKNOWN":
            continue 

        if subtype =="Cell_line":
            continue

        samples = sample_info[sample_info["Subgroup"] == subtype]["ATAC.sample.code"]
        dirname = "subtype_balanced_samples.dir"
        outfile = os.path.join(dirname, "%s_pooled_subtype.single.end.shifted.filtered.tagAlign.gz" 
                                        % subtype)
        infiles = [os.path.join(dirname, "%s.single.end.shifted.filtered.tagAlign.gz") % s 
                   for s in samples]
        yield (infiles, outfile)

@follows(get_balanced_sample_filtered_shifted_SE)
@files(get_files_for_pooling)
def pool_balanced_single_ends(infiles, outfile):
    '''Pool the downsampled balanced single-end aligned tags for subgroup
    and also for MM and ND'''

    infiles = " ".join(infiles)
    statement = '''zcat %(infiles)s | gzip > %(outfile)s'''
    P.run(statement)

#----------------------------------------------------------------------------------------
@follows(mkdir("filtered_peaks_broad.dir"),
         mkdir("pan_peaks_broad.dir"),
         mkdir("subtype_peaks_broad.dir"))
@subdivide((filterShiftTagAlign,
            pool_balanced_single_ends),
           regex("([^_]+)_.+.dir/(.+).single.end.shifted.filtered.tagAlign.gz"),
           [r"\1_peaks_broad.dir/\2_peaks.broadPeak.gz",
            r"\1_peaks_broad.dir/\2_peaks.gappedPeak.gz"],
           r"\2")
def call_peaks_broad(infile, outfiles, sample):
    ''' Use MACS2 to calculate broad peaks and gapped peaks.
    Sorts the output by -log10pvalue,
    formats the name of the broad and gapped peaks '''
    
    # Get the thresholding values for MACS2
    threshold_method = PARAMS["macs2_threshold_method"].lower()
    
    # If nothing specified default to p (p-value)
    if threshold_method == "":
        threshold_method = "p"
    elif threshold_method not in ["p", "q"]:
        raise Exception("threshold method specified not valid")
    
    
    threshold_quantity = PARAMS["macs2_threshold_quantity"]
    
    # If nothing specified default to 0.01
    if threshold_quantity == "":
        threshold_quantity = "0.01"
        
    
    # Get the read extending and shift values
    shift_parameter = PARAMS["end_extending_shift"]
    
    # If nothing specified default to -100
    if shift_parameter == "":
        shift_parameter = "-100"
      
    extsize_parameter = PARAMS["end_extending_extsize"]
    
    # If nothing specified default to 200
    if extsize_parameter == "":
        extsize_parameter = "200"
    
    # Get the temporal dir specified
    tmp_dir = PARAMS["general_temporal_dir"]
    
    # Get the directory for the outfile
    outdir = os.path.dirname(outfiles[0])
    
    # Create a temporal directory name for the run
    peaks_temp_dir = tempfile.mkdtemp(dir=tmp_dir)
    
    # Get the name of the peak file (_peaks.broadPeak.gz)
    outfile_basename = os.path.basename(outfiles[0])
    outfile_basename_prefix = P.snip(outfile_basename, "_peaks.broadPeak.gz")
    
    # Create the prefix to create the files in the temporal directory
    # and with the name of the run
    outfile_prefix = os.path.join(peaks_temp_dir, outfile_basename_prefix)
    
    # Final outfile names for the post-processed files
    broad_peaks_final_outfile = outfiles[0]
    gappedPeak_peaks_final_outfile = outfiles[1]
    
    
    # 1) Calculate peaks, will create .xls, .broadPeak and .gappedPeak in temp dir
    # 2) Process the .broadPeak and .gappedPeak and output in outdir
    # 3) Delete the .broadPeak and .gappedPeak temp files
    # 4) Copy the remaining created files from temp dir to outdir
    statement = '''macs2 callpeak 
                      -t %(infile)s 
                      -f BED 
                      -n %(outfile_prefix)s 
                      -g hs 
                      -%(threshold_method)s %(threshold_quantity)s 
                      --nomodel 
                      --shift %(shift_parameter)s 
                      --extsize %(extsize_parameter)s 
                      --broad 
                      --broad-cutoff %(threshold_quantity)s
                      --keep-dup all
                      --tempdir %(tmp_dir)s &&
        
    sort -k 8gr,8gr %(outfile_prefix)s_peaks.broadPeak 
    | awk 'BEGIN{OFS="\\t"}{$4="Peak_"NR ; print $0}'
    | gzip -c > %(broad_peaks_final_outfile)s &&
      
    rm %(outfile_prefix)s_peaks.broadPeak &&
      
    sort -k 14gr,14gr %(outfile_prefix)s_peaks.gappedPeak 
    | awk 'BEGIN{OFS="\\t"}{$4="Peak_"NR ; print $0}' | 
    gzip -c > %(gappedPeak_peaks_final_outfile)s &&
        
    rm %(outfile_prefix)s_peaks.gappedPeak &&
    
    mv %(outfile_prefix)s* %(outdir)s
    '''
    
    # The pooled datasets contain a lot of reads, it is advisable to up
    # the memory in this case
    job_memory = "4G"
    
    # Pooled samples for this pipeline
    if ("pooled" in sample):
        
        job_memory = "10G"
    
    P.run(statement, job_condaenv="macs2")


#----------------------------------------------------------------------------------------
@follows(mkdir("filtered_peaks_narrow.dir"),
         mkdir("pan_peaks_narrow.dir"),
         mkdir("subtype_peaks_narrow.dir"))
@subdivide((filterShiftTagAlign,
            pool_balanced_single_ends),
           regex("([^_]+)_.+.dir/(.+).single.end.shifted.filtered.tagAlign.gz"),
           [r"\1_peaks_narrow.dir/\2_peaks.narrowPeak.gz",
             r"\1_peaks_narrow.dir/\2_summits.bed"],
           r"\2")
def call_peaks_narrow(infile, outfiles, sample):
    ''' Use MACS2 to calculate peaks '''
    
    # Get the thresholding values for MACS2
    threshold_method = PARAMS["macs2_threshold_method"].lower()
    
    # If nothing specified default to p (p-value)
    if threshold_method == "":
        threshold_method = "p"
    elif threshold_method not in ["p", "q"]:
        raise Exception("threshold method specified not valid")
    
    
    threshold_quantity = PARAMS["macs2_threshold_quantity"]
    
    # If nothing specified default to p (p-value)
    if threshold_quantity == "":    
        threshold_quantity = "0.01"
    
    # Get the read extending and shift values
    shift_parameter = PARAMS["end_extending_shift"]
    
    # If nothing specified default to -100
    if shift_parameter == "": 
        shift_parameter = "-100"
    
    extsize_parameter = PARAMS["end_extending_extsize"]
    
    # If nothing specified default to 200
    if extsize_parameter == "":
        extsize_parameter = "200"
    
    # Get the temporal dir specified
    tmp_dir = PARAMS["general_temporal_dir"]
    
    # Get the directory for the outfile
    outdir = os.path.dirname(outfiles[0])
    
    # Create a temporal directory name for the run
    peaks_temp_dir = tempfile.mkdtemp(dir=tmp_dir)
    
    # Get the name of the peak file (_peaks.narrowPeak.gz)
    outfile_basename = os.path.basename(outfiles[0])
    
    outfile_basename_prefix = P.snip(outfile_basename, "_peaks.narrowPeak.gz")
    
    # Create the prefix to create the files in the temporal directory
    # and with the name of the run
    outfile_prefix = os.path.join(peaks_temp_dir, outfile_basename_prefix)
    
    # Final outfile names for the post-processed files
    narrow_peaks_final_outfile = outfiles[0]
    
    # 1) Calculate peaks, will create .xls, .narrowPeak in temp dir
    # 2) Process the .narrowPeak and output in outdir: Note that MACS2 can create multiple narrow peaks (Doesn't happen with broad peaks)
    #    with the same coordinate and different significance. According to http://seqanswers.com/forums/showthread.php?t=50394
    #    these are summits within the same peak. To sort this out, whenever there are multiple rows with the same coordinates we
    #    get the row with the highest -log10pvalue (column 8), using sort as in https://unix.stackexchange.com/questions/230040/keeping-first-instance-of-duplicates
    # 3) Delete the .narrowPeak temp files
    # 4) Copy the remaining created files from temp dir to outdir
    statement = '''macs2 callpeak 
                     -t %(infile)s 
                     -f BED 
                      -n %(outfile_prefix)s 
                      -g hs 
                      -%(threshold_method)s %(threshold_quantity)s 
                      --nomodel 
                      --shift %(shift_parameter)s 
                      --extsize %(extsize_parameter)s 
                      -B 
                      --SPMR 
                      --keep-dup all
                      --call-summits
                      --tempdir %(tmp_dir)s &&

    cat %(outfile_prefix)s_peaks.narrowPeak 
    | sort -k1,3 -k8gr 
    | sort -k1,3 -u 
    | sort -k8gr  
    | awk 'BEGIN{OFS="\\t"}{$4="Peak_"NR ; print $0}' 
    | gzip -c > %(narrow_peaks_final_outfile)s &&
    
    rm %(outfile_prefix)s_peaks.narrowPeak &&
    
    mv %(outfile_prefix)s* %(outdir)s
    '''
    
    # The pooled datasets contain a lot of reads, it is advisable to up
    # the memory in this case
    job_memory = "4G"
    
    if ("pooled" in sample):
        job_memory = "10G"        
    
    P.run(statement, job_condaenv="macs2")


#------------------------------------------------------------------------------
@follows(mkdir("filtered_peaks.dir"))
@transform([call_peaks_broad, call_peaks_narrow],
           regex(".+/(.+).gz"),
           r"filtered_peaks.dir/\1.gz")
def filter_peaks(infile, outfile):
    ''' Filters out regions of low mappability and excessive mappability '''
    
    # Temp file: We create a temp file to make sure the whole process goes well
    # before the actual outfile is created
    temp_file = P.snip(outfile, ".gz") + "_temp.gz"
    
    excluded_beds = PARAMS["filtering_bed_exclusions"]
    
    statement = pipelineAtacseq.createExcludingBedsFromBedStatement(infile, 
                                                                    excluded_beds, 
                                                                    temp_file)
    
    statement += ''' &&
                    mv %(temp_file)s %(outfile)s'''

    P.run(statement)
#----------------------------------------------------------------------------------------

@collate(filter_peaks,
        regex(".+/(.+).(broadPeak|narrowPeak).gz"),
        r"filtered_peaks.dir/\1.mergedpeaks.gz")
def merge_broad_narrow_peaks(infiles, outfile):
    '''Merge the broad and narrow peak files for each sample, merging any peaks less than
    100 nt away from another peak'''

    infiles = " ".join(infiles)
    statement = '''zcat %(infiles)s
                   | sort -k1,1 -k2,2n 
                   | awk 'FS="\\t" {printf ("%%s\\t%%s\\t%%s\\t%%s\\n", $1, $2, $3, $6)'} 
                   | bedtools merge -i stdin -d 200
                   | sort -k1,1 -k2,2n 
                   | gzip -c > %(outfile)s '''

    P.run(statement)


@collate(merge_broad_narrow_peaks,
        regex(".+/(.+)_pooled_(.+)_peaks.mergedpeaks.gz"),
        r"filtered_peaks.dir/\2_merged_peaks.bed.gz")
def merge_pooled_peaks(infiles, outfile):
    '''Merge MM and ND peaks to get pan merged peaks, and the peaks for each subtype to get
    balanced subtype peaks'''

    infiles = " ".join(infiles)
    statement = '''zcat %(infiles)s
                   | sort -k1,1 -k2,2n 
                   | awk 'FS="\\t" {printf ("%%s\\t%%s\\t%%s\\t%%s\\n", $1, $2, $3, $6)'} 
                   | bedtools merge -i stdin -d 200
                   | sort -k1,1 -k2,2n 
                   | gzip -c > %(outfile)s '''

    P.run(statement)


@follows(merge_pooled_peaks,
         merge_broad_narrow_peaks)
def call_peaks():
    ''' Dummy task to sync the call of peaks '''
    pass


##########################################################################################
# Peak Quantitation
##########################################################################################
@transform(shiftTagAlign,
           regex(".+\.dir/(.+).PE2SE.tn5_shifted.tagAlign.gz"),
           r"final_tag_align.dir/\1.five.prime.only.single.end.processed.tagAlign.gz")
def get_five_prime_only_single_ends(infile, outfile):
    ''' It creates a 1bp region for each shifted SE with the 5' end only (strand-aware).
    Sorts the output '''
    
    # Get the temporal dir specified
    tmp_dir = PARAMS["general_temporal_dir"]
    
    # To make sure the run is complete before copying the result to the outfile
    outfile_temp = (tempfile.NamedTemporaryFile(dir=tmp_dir, delete=False)).name
    
    
    
    # Gets a 1bp region with the 5' end (strand-aware) and sorts the tagAlignFile
    statement = '''zcat %(infile)s 
                   | awk 'BEGIN {
                            FS ="\\t"; OFS = FS
                          } 

                          { if ($6 == "+") {
                                $3 = $2+1
                          } else if ($6 == "-") {
                                $2 = $3-1
                          } 
                          print $0}'
                   | sort -k 1,1 -k2,2n 
                   | gzip -c > %(outfile_temp)s  &&
    
                   mv %(outfile_temp)s %(outfile)s
    '''
    
    P.run(statement)


@transform(get_five_prime_only_single_ends,
           regex(".+\.dir/(.+).five.prime.only.single.end.processed.tagAlign.gz"),
           r"filtered_tag_align.dir/\1.five.prime.only.single.end.processed.filtered.tagAlign.gz")
def filter_five_prime_only_single_ends(infile, outfile):
    ''' Filters out regions of low mappability and excessive mappability in the 1bp five prime single ends'''
    
    # Reuse the filter_peaks function
    filter_peaks(infile, outfile)


#------------------------------------------------------------------------------
@follows(mkdir("tag_counts.fir"))
@product(filter_five_prime_only_single_ends,
         formatter(".+\.dir/(?P<SAMPLE>.+).five.prime.only.single.end.processed.filtered.tagAlign.gz"),
         merge_pooled_peaks,
         formatter(".+/(?P<PEAKS>.+)_merged_peaks.bed.gz"),
         "tag_counts.dir/{SAMPLE[0][0]}_vs_{PEAKS[1][0]}.tag_counts.gz")
def overlap_sample_with_common_peaks(infiles, outfile):
    ''' For each common peak to all the samples, it calculates the number of filtered
    extended 5' 1bp single ends of reads from the sample which fall in that region.
    An overlap is considered valid if all the SE (1bp) is overlapped by a feature
    (common peak). Produces one line for each common peak, regardless of whether
     there is or not intersection'''
    extended_tag_align = infiles[0]
    merged_common_peaks = infiles[1]
     
    # Get the temporal dir specified
    tmp_dir = PARAMS["general_temporal_dir"]
    
    # To make sure the run is complete before copying the result to the outfile
    outfile_temp = (tempfile.NamedTemporaryFile(dir=tmp_dir, delete=False)).name
     
    # The common_peaks file is already sorted
    # Intersect and count the number of reads in each common_peak
    statement = '''bedtools intersect 
                              -sorted -c 
                              -a %(merged_common_peaks)s 
                              -b %(extended_tag_align)s 
                   | gzip -c > %(outfile_temp)s &&
                   
                   mv %(outfile_temp)s %(outfile)s
    '''
     
    P.run(statement)


#------------------------------------------------------------------------------
@collate(overlap_sample_with_common_peaks,
         regex(".+/(.+)_vs_(.+).tag_counts.gz"),
         r"tag_counts.dir/\2_raw_tag_counts.bed_header.gz")
def group_tag_counts_per_peakset(infiles, outfile):
    '''Take the counts of each sample against each peakset and generate a matrix
    with peak positions and the counts of every sample across that peakset. 
    Sample titles will be shortened. Note that they will have no "A" on the
    start'''

    col_names = " ".join([re.match("tag_counts.dir/(.+).bowtie2_vs_", infile).groups()[0]
                         for infile in infiles] )
    infiles = " ".join(infiles)

    statement = '''bedtools unionbedg 
                             -header
                             -names %(col_names)s
                             -i %(infiles)s
                    | gzip > %(outfile)s'''

    P.run(statement)


#------------------------------------------------------------------------------
@transform(group_tag_counts_per_peakset,
           suffix("bed_header.gz"),
           "tsv.gz")
def combine_bed_cols_to_id(infile, outfile):
    '''Currently the peaks are formatted as bed files - first three columns are contig,
    start and end. These need to be merged into a single column - id - for DE analysis'''

    statement = '''zcat %(infile)s
                 | awk 'BEGIN{FS="\\t";OFS=FS}
                        {$1=$1":"$2"-"$3; 
                         $2=""; $3=""}
                        {printf("%%s\\n",$0)}'
                 | cut --complement -f 2,3
                 | awk 'BEGIN{FS="\\t";OFS=FS} NR==1 {$1="id"} {print $0}'
                 | gzip > %(outfile)s'''

    P.run(statement)


#------------------------------------------------------------------------------
@transform(combine_bed_cols_to_id,
           suffix(".tsv.gz"),
           add_inputs("samples.tsv"),
           ".tech_reps_collapsed.tsv.gz")
def collapse_tech_reps(infiles, outfile):
    '''Technical replicates are seperate libraries created from the sample underlying
    sample. Here we sum counts from technical replicates. Also outputs col_data'''

    script = os.path.join(os.path.dirname(__file__), "scripts/R_collapseReps.R")
    counts, coldata = infiles
    statement = '''Rscript %(script)s %(counts)s %(coldata)s biosample_id %(outfile)s'''
    P.run(statement)


#----------------------------------------------------------------------------------------
@transform(collapse_tech_reps,
           suffix(".tech_reps_collapsed.tsv.gz"),
           add_inputs(r"\1.tech_reps_collapsed.col_data.tsv"),
           ".donors_collapsed.tsv.gz")
def collapse_cd19(infiles, outfile):
    '''Some of the ND libraries were selected for CD19+ and some for CD19-. MOFA analysis
    reveals little to no difference in these, and each were sequenced inthe same batch. 
    These will be merged as doing so will remove a variable from the model. Also outputs
    colldata '''

    script = os.path.join(os.path.dirname(__file__), "scripts/R_collapseReps.R")
    counts, coldata = infiles
    statement = '''Rscript %(script)s %(counts)s %(coldata)s donor_id %(outfile)s'''
    P.run(statement)

#-----------------------------------------------------------------------------------------
@follows(mkdir("DE.dir"))
@subdivide(collapse_cd19,
           regex("(.+)/(.+)_raw_tag_counts.donors_collapsed.tsv.gz"),
           add_inputs(r"\1/\2_raw_tag_counts.donors_collapsed.col_data.tsv"),
           [r"DE.dir/all_\2_atac_regions.tsv.gz",
            r"DE.dir/sign_\2_atac_regions.tsv.gz"],
            r"\2")
def run_atac_de(infiles, outfiles, grouping):
    '''Use DESeq to call differentially open regions for both pan and subgroup peaks.
    Bit of a cludge because the two calling processes are quite different, and input
    files to scripts are hard coded'''

    script = os.path.join(os.path.dirname(__file__), "scripts/R_atac_%s_de.R" % grouping)
    statement = '''Rscript %(script)s'''

    job_memory="16G"
    P.run(statement)

##########################################################################################
# RNA Processing
##########################################################################################
@mkdir('geneset.dir')
@transform(PARAMS['geneset'],
           regex("(\S+).gtf.gz"),
           r"geneset.dir/\1.fa")
def buildReferenceTranscriptome(infile, outfile):
    '''
    Builds a reference transcriptome from the provided GTF geneset - generates
    a fasta file containing the sequence of each feature labelled as
    "exon" in the GTF.
    --fold-at specifies the line length in the output fasta file

    Parameters
    ----------
    infile: str
        path to the GTF file containing transcript and gene level annotations
    genome_dir: str
        :term: `PARAMS` the directory of the reference genome
    genome: str
        :term: `PARAMS` the filename of the reference genome (without .fa)
    outfile: str
        path to output file
    '''

    genome_file = os.path.abspath(
        os.path.join(PARAMS["genome_dir"], PARAMS["genome"] + ".fa"))

    statement = '''
    zcat %(infile)s |
    awk '$3=="exon"'|
    cgat gff2fasta
    --is-gtf --genome-file=%(genome_file)s --fold-at=60 -v 0
    --log=%(outfile)s.log > %(outfile)s &&
    samtools faidx %(outfile)s
    '''

    P.run(statement)


#-----------------------------------------------------------------------------------------
@transform(buildReferenceTranscriptome,
           suffix(".fa"),
           ".salmon.index")
def buildSalmonIndex(infile, outfile):
    '''
    Builds a salmon index for the reference transriptome
    Parameters
    ----------
    infile: str
       path to reference transcriptome - fasta file containing transcript
       sequences
    salmon_kmer: int
       :term: `PARAMS` kmer size for salmon.  Default is 31.
       Salmon will ignores transcripts shorter than this.
    salmon_index_options: str
       :term: `PARAMS` string to append to the salmon index command to
       provide specific options e.g. --force --threads N
    outfile: str
       path to output file
    '''

    job_memory = "16G"
    job_threads = 2
    # need to remove the index directory (if it exists) as ruffus uses
    # the directory timestamp which wont change even when re-creating
    # the index files
    statement = '''
    rm -rf %(outfile)s;
    salmon index -k 31 --threads 2 -t %(infile)s -i %(outfile)s
    '''

    P.run(statement)


#-----------------------------------------------------------------------------------------
@transform(PARAMS["geneset"],
           formatter(),
           "transcript2geneMap.tsv")
def getTranscript2GeneMap(infile, outfile):
    ''' Extract a 1:1 map of transcript_id to gene_id from the geneset '''

    pipelineAtacseq.gene_to_transcript_map(infile, outfile)
    

#------------------------------------------------------------------------------
@follows(mkdir("salmon.dir"))
@collate(("input_rna.dir/*.fastq.1.gz",
          "input_rna.dir/*.fastq.gz",
          "input_rna.dir/*.sra"),
         regex(".+/(\S+).(fastq.1.gz|fastq.gz|sra)"),
         add_inputs(buildSalmonIndex, getTranscript2GeneMap),
         [r"salmon.dir/\1/transcripts.tsv.gz",
         r"salmon.dir/\1/genes.tsv.gz"])
def runSalmon(infiles, outfiles):
    '''
    Computes read counts across transcripts and genes based on a fastq
    file and an indexed transcriptome using Salmon.

    Runs the salmon "quant" function across transcripts with the specified
    options.  Read counts across genes are counted as the total in all
    transcripts of that gene (based on the getTranscript2GeneMap table)

    Parameters
    ----------
    infiles: list
        list with three components
        0 - list of strings - paths to fastq files to merge then quantify
        across using salmon
        1 - string - path to salmon index file
        2 - string - path to table mapping transcripts to genes

    salmon_threads: int
       :term: `PARAMS` the number of threads for salmon
    salmon_memory: str
       :term: `PARAMS` the job memory for salmon
    salmon_options: str
       :term: `PARAMS` string to append to the salmon quant command to
       provide specific
       options, see https://salmon.readthedocs.io/en/latest/salmon.html#description-of-important-optionsç∂
    salmon_bootstrap: int
       :term: `PARAMS` number of bootstrap samples to run.
       Note, you need to bootstrap for differential expression with sleuth
       if there are no technical replicates. If you only need point estimates,
       set to 1.
    salmon_libtype: str
       :term: `PARAMS` salmon library type, https://salmon.readthedocs.io/en/latest/library_type.html
    outfiles: list
       paths to output files for transcripts and genes
    '''
    fastqfile = [x[0] for x in infiles]
    index = infiles[0][1]
    transcript2geneMap = infiles[0][2]

    transcript_outfile, gene_outfile = outfiles
    Quantifier = rnaseq.SalmonQuantifier(
        infile=fastqfile,
        transcript_outfile=transcript_outfile,
        gene_outfile=gene_outfile,
        annotations=index,
        job_threads=PARAMS["salmon_threads"],
        job_memory=PARAMS["salmon_memory"],
        options="--writeUnmappedNames --gcBias",
        bootstrap=0,
        libtype="A",
        kmer=31,
        transcript2geneMap=transcript2geneMap)

    Quantifier.run_all()


#-----------------------------------------------------------------------------
@collate(runSalmon,
         regex("(\S+).dir/(\S+)/transcripts.tsv.gz"),
         [r"\1.dir/transcripts.tsv.gz",
          r"\1.dir/genes.tsv.gz"])
def merge_rna_counts(infiles, outfiles):
    '''Merge counts from each sample into a single table'''

    pipelineAtacseq.merge_counts(infiles, outfiles, submit=True)


#------------------------------------------------------------------------------
@follows(mkdir("geneset.dir"))
@transform("input_rna.dir/*.bam", formatter(),
           add_inputs(PARAMS["geneset"]),
           "geneset.dir/{basename[0]}.gtf.gz")
def assembleWithStringTie(infiles, outfile):
    '''Generate a transcript assembly from each sample using stringtie'''

    infile, reference = infiles

    job_threads = PARAMS["stringtie_threads"]
    job_memory = PARAMS["stringtie_memory"]

    statement = '''stringtie %(infile)s
                           -p %(stringtie_threads)s
                           -G <(zcat %(reference)s)
                           2> %(outfile)s.log
                   | gzip > %(outfile)s '''

    P.run(statement)


# ----------------------------------------------------------------------------
@merge((assembleWithStringTie,
        PARAMS["geneset"]),
       "geneset.dir/agg-agg-agg.gtf.gz")
def mergeAllAssemblies(infiles, outfile):
    '''Merge the assemblies from each sample together into a single merged
    assembly'''

    infiles = ["<(zcat %s)" % infile for infile in infiles]
    infiles, reference = infiles[:-1], infiles[-1]

    job_threads = PARAMS["stringtie_merge_threads"]
    job_memory = "%iG" % (24.0/PARAMS["stringtie_merge_threads"])
    infiles = " ".join(infiles)

    statement = '''stringtie --merge
                             -G %(reference)s
                             -p %(stringtie_merge_threads)s
                             %(infiles)s
                            2> %(outfile)s.log
                   | cgat gtf2gtf --method=sort
                           --sort-order=gene+transcript
                            -S %(outfile)s -L %(outfile)s.log'''

    P.run(statement) 


# ----------------------------------------------------------------------------
@transform(mergeAllAssemblies,
           regex(".+"),
           "geneset.dir/denovo_promoters.bed.gz")
def get_denovo_promoters(infile, outfile):
    '''Convert the denovo generated geneset to BED12 format,
    filter out the one exon transcripts, and get a region 2kb uptream
    and 100bp downstream'''

    contigs = os.path.join(os.path.dirname(__file__), "contigs.tsv")
    statement = '''cgat gff2bed
                           --bed12-from-transcripts
                            -I %(infile)s
                            -L %(outfile)s.log 
                    | awk '$10>1' 
                    | bedtools flank -i - -l 2000 -r 0 -s -g %(contigs)s
                    | bedtools slop -r 100 -l 0 -s -g %(contigs)s
                    | sort -k1,1 -k2,2n
                    | gzip > %(outfile)s '''

    P.run(statement)


# ----------------------------------------------------------------------------
@transform(PARAMS["geneset"],
           formatter(),
           "geneset.dir/reference_promoters.bed.gz")
def get_reference_promoters(infile, outfile):
    '''Get a region 2kb upstream and 100bp downstream of each transcripts. Note
    that in contrast to the denovo transcripts processed above, we keep 1 exon
    transcripts.'''

    contigs = os.path.join(os.path.dirname(__file__), "contigs.tsv")
    statement = '''cgat gff2bed
                           --bed12-from-transcripts
                            -I %(infile)s
                            -L %(outfile)s.log 
                    | bedtools flank -i - -l 2000 -r 0 -s -g %(contigs)s
                    | bedtools slop -r 100 -l 0 -s -g %(contigs)s
                    | sort -k1,1 -k2,2n
                    | gzip > %(outfile)s '''

    P.run(statement)


# ----------------------------------------------------------------------------
@merge([get_denovo_promoters, get_reference_promoters],
        "all_promoters.bed.gz")
def merge_promoters(infiles, outfile):
    '''Merge the reference and denovo TSSes'''

    infiles = " ".join(infiles)
    statement = '''zcat %(infiles)s
                 | sort -k1,1 -k2,2n
                 | bedtools merge -i -
                 | gzip > %(outfile)s'''

    P.run(statement)


@follows(merge_promoters,
         merge_rna_counts)
def rna_preprocessing():
    pass

##########################################################################################
##                   Targets          
##########################################################################################

@follows(getFirstFilteringStats,
        getInitialMappingStats,
        getPostDuplicationStats,
        calculateLibrarycomplexity)
def get_stats():
    pass

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))