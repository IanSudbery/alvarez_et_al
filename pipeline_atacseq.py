from ruffus import follows, transform, add_inputs, mkdir, regex, formatter
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
#sys.path.insert(0, "/home/mbp15ja/dev/pipeline_chromHMM/src/")
#import PipelineChromHMM
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
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"])


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
                       submit=True,
                       job_memory="6G")
    
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
    


    statement = '''samtools view -F 524 -f 2 -u %(infile)s | 
                   samtools sort -n - -o %(first_filtering_bam_output)s 
                                  -T %(samtools_temp_file)s 2> %(log_file)s &&
                   samtools view -h %(first_filtering_bam_output)s | 
                    scripts/assign_multimappers.py -k %(allowed_multimappers)s --paired-end | 
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
    five_prime_trim = pipelineAtacseq.getSampleQCShift(sample_name, sample_details)
    
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
@transform(markDuplicates,
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
@transform("dedupped.dir/*_pos_sorted.bam",
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
@transform("dedupped.dir/*_pos_sorted.bam",
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
    gzip -c > %(temp_file)s 2> %(log_file)s;
    checkpoint;
    
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
    
    statement += '''
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
    five_prime_trim = pipelineAtacseq.getSampleQCShift(sample_name, sample_details)
    
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
    contigs = PARAMS["general_contigs"]
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
                         $2 = $2 + 4%(positive_strand_correction)s
                        } else if ($6 == "-") {
                         $3 = $3 - 5%(negative_strand_correction)s} print $0}'
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
    
    statement += ''' checkpoint;
                    mv %(temp_file)s %(outfile)s'''

    P.run(statement)
