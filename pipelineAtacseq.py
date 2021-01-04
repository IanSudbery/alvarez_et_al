import copy
import pysam
import collections
import os
import cgatcore.experiment as E
from cgatcore.pipeline import cluster_runnable
import cgatcore.iotools as IOTools
import cgatcore.pipeline as P
from cgat import GTF
import tempfile
import re
from cgat import Bed

import pandas as pd
#from pandas.core.frame import DataFrame
#import numpy as np
#import matplotlib as mpl
#mpl.use('Agg') # So that matplotlib.pyplot doesn't try to find the X-server and give an error
#import matplotlib.pyplot as plt
#import sys
#sys.path.insert(0, "/home/mbp15ja/dev/AuxiliaryPrograms/StringOperations/")
#import StringOperations

@cluster_runnable
def getMappedUnmappedReads(infile, outfile):
    
    # Determine the extension of the file (BAM or SAM), done this for SAM/BAM compatibility
    # At the moment, all calls should only pass SAM, if it stops working just leave 
    # samfile = pysam.AlignmentFile(infile, "r")
    _ , file_extension = os.path.splitext(infile)
    
    if (file_extension.lower() == '.sam'):
        samfile = pysam.AlignmentFile(infile, "r")   
    elif (file_extension.lower() == '.bam'):
        samfile = pysam.AlignmentFile(infile, "rb")
    
    # Create the name of the files for the unmapped and mapped reads
    out_unmapped_basename = P.snip(os.path.basename(outfile), '.tsv') + '.unmapped.txt.gz'
    out_unmapped = os.path.join(os.path.dirname(outfile), (out_unmapped_basename))
      
    out_mapped_basename = P.snip(os.path.basename(outfile), '.tsv') + '.mapped.txt.gz'
    out_mapped = os.path.join(os.path.dirname(outfile), (out_mapped_basename))
    
    # Unmapped readnames
    unmapped = IOTools.open_file(out_unmapped, "w")
    
    # Mapped readnames
    mapped = IOTools.open_file(out_mapped, "w")
    
    dictionary = {'read_name': 'first', 'mapped': 'no'}
    
    c = collections.Counter()
       
    # Go through the reads
    for read in samfile:
           
        # Get the read name
        read_name = read.query_name
           
        # First case, assign the read name to the dictionary
        if dictionary['read_name'] == 'first':
            dictionary['read_name'] = read_name
            # Add one read to the total
            c["total"] += 1
                 
        # If a fragment from the read has been seen before
        if dictionary['read_name'] == read_name:  
               
            # If that fragment was mapped
            if dictionary['mapped'] == 'yes':
                # Skip to the next data line
                continue
               
            # If all the prior fragments from the read haven't been mapped
            else:
                # If it is mapped, change the state to mapped (the default beginning state is unmapped)
                if not read.is_unmapped:
                    # Put the read as mapped
                    dictionary['mapped'] = 'yes'
    
        # If a fragment from the same read hasn't been seen before (new read)        
        else:
  
            # If the last read wasn't mapped, store it
            if dictionary['mapped'] == 'no':
                unmapped.write(dictionary['read_name']+"\n")   
                # Add one read to the unmapped
                c["unmapped"] += 1
                  
            else:
                  
                mapped.write(dictionary['read_name'] + "\n") 
                # Add one read to the mapped
                c["mapped"] += 1
               
               
            # Add one read to the unmapped
            c["total"] += 1
               
            # Create a new dictionary with unmapped and the current read name
            dictionary = {'read_name': read_name, 'mapped': 'no'}
               
            # If it is mapped, change the state to mapped (the default beginning state is unmapped)
            if not read.is_unmapped: 
                # Put the read as mapped
                dictionary['mapped'] = 'yes'

    # The last read is hanging
    # If it is unmapped, store it
    if dictionary['mapped'] == 'no':
        unmapped.write(read_name + "\n")  
        # Add one read to the unmapped
        c["unmapped"] += 1    
    else:          
        mapped.write(dictionary['read_name'] + "\n")

        # Add one read to the mapped
        c["mapped"] += 1
  
    # Store the counters        
     
    # Create a name for the counters
    with IOTools.open_file(outfile, "w") as output_file_write:   
        for counter in c:  
            output_file_write.write(counter + ':\t' + str(c[counter]) + '\n')
             

#------------------------------------------------------------------------------
# Looks into sample_table for the sample specified by sample_name
# (it performs a case insensitive comparison).
# If it finds it, returns the shift specified for that sample.
# If it finds an empty string, it returns "0".
# If it doesn't find it, or if it is not an integer or empty returns an exception
# The shift is defined as the previously 5' constant (performed on all reads within a read pair) 
# trimming.
# For example, if 2bp have been trimmed from the 5' end in qc, 
# five_prime_quality_trimming=2
# Negative numbers are also allowed if the reads have been extended.
# Inputs:
#     -sample_name: The name of the sample
#     -sample_table: The tab separated table specifying the following 
#        (doesn't contain a header):
#    
#        Sample_name    condition(control/treatment)    five_prime_quality_trimming(number of base pairs trimmed from both reads' 5' end in a pair)
#        nd2345    control    2
#    
#
# Outputs:
#     -returns the shift specified for that sample or an exception
def getSampleQCShift(sample_name, sample_table):
    
    shift = ""
    
    # To make sure the sample is found
    found = False
    
    with open(sample_table, 'r') as reader:
    
        for line in reader:
            
            # Remove the new line character from the line 
            # (avoids taking this into account downstream)
            line = line.rstrip('\n')
            
            fields = line.split("\t")
            
            # If the name of the sample is equal to the name specified in the 
            # file
            if fields[0].lower() == sample_name.lower():
                shift = fields[2]

                # If it is not empty test whether it is an integer
                if shift != "":
                    try:
                        int(shift)
                    except ValueError:
                        raise Exception("The sample specified: " + sample_name + 
                                        " contains a non empty non number string")
                # If it is empty return 0
                else:
                    shift = "0"
                found = True
                break 
    
    if not found:
        raise Exception("The sample specified: "+ sample_name +
                         " has not been found")
    else:
        return shift


#------------------------------------------------------------------------------
# inspired pipelineCaptCPerl.filterReadMappings
# Assuming a Sam/Bam file resulting from read pairs mapping
# and containing ONLY UNIQUE mappings per pair in each read
# and SORTED BY READNAME from any mapper output. 
# Because a pair R1 and R2 can only have one alignment each, 
# at most there would be 2 alignments, both alignments must 
# reflect either a proper pair or not a proper pair, same with
# pcr_duplicates. Therefore this information can be extracted 
# after both reads have been seen  
# Obtains statistics from the reads and read pairs.
# For individual reads:
#    -Unmapped reads
#    -Reads which are supplementary alignment (subread parts mapping to
#        multiple places.
#    -Exception if any non primary alignments found
# For read pairs:
#    -Read is not paired (will count 1 read per supposed read pair: not found)
#    -Read pair is not mapped in proper pair
#    -PCR duplicate
#    -Read pairs which are all of the below: 
#        -Read pair mapped in proper pair
#        -Doesn't contain supplementary alignments
#        -Doesn't contain secondary alignments
#        -Not PCR duplicate
#        -Contains R1 and R2
#    -Read pairs which are all of the below: 
#        -Read pair mapped in proper pair
#        -Contains supplementary alignments (maybe translocations, inversions, etc..)
#        -Doesn't contain secondary alignments
#        -Not PCR duplicate
#        -Contains R1 and R2
#
# Inputs:
#     -infile: Sam/Bam file containing only unique mappings per pair in each read
#         sorted by readname sorted by readname from any mapper output
#     -stats_outfile: Outfile with the statistics
#
# Output:
#     -returns printable output with the stats
#
# Exception:
#     -If a multimapping (read alignment not being primary is found).
@cluster_runnable
def getUniquelyMappedPairsNoMultimapping(infile, stats_outfile):
    
    c = collections.Counter()
    
    # Determine the extension of the file (BAM or SAM), done this for SAM/BAM 
    # compatibility At the moment, all calls should only pass SAM, if it stops
    # working just leave samfile = pysam.AlignmentFile(infile, "r")
    _ , file_extension = os.path.splitext(infile)
    
    if (file_extension.lower() == '.sam'):
        samfile = pysam.AlignmentFile(infile, "r")
    
    elif (file_extension.lower() == '.bam'):
        samfile = pysam.AlignmentFile(infile, "rb")
    
    
    # We need to know when the read name changes
    last_read_name = "first"
    
    # We need to know the last read when the read name changes to assign
    # the read_pair details from the last read when we have seen all the read 
    # pair alignments
    last_read = ""
    
    # Create new dictionary template to store the details
    read_pair_dict_template = {'read_pair_is_paired':"", # Read pair contains R1 and R2    
            'read_pair_mapped_in_proper_pair':"", # Read pair contains proper mapped pair
            'read_pair_is_pcr_duplicate':"", # Read pair is a PCR duplicate
            'read_number_unmapped':0, # Number of R1 and R2 unmapped
            'read_number_supplementary':0, # Number of R1 and R2 supplementary (Eg. Translocations) 
            }
    
    # Create new copy of the template dictionary to store
    read_pair_dict = copy.deepcopy(read_pair_dict_template)
        
    # Process the reads
    for read in samfile:
        
        # If it is the first read, assign the name (it will consider
        # it as read name hasn't changed)
        if last_read_name == "first":
            last_read_name = read.query_name
            

        # The read has changed excepting the first read from the file
        if (last_read_name != read.query_name):       
            # We have seen all the reads from a pair, fill the details
            # of the pair
            
            # reads contain a pair (R1 and R2)
            read_pair_dict['read_pair_is_paired'] = last_read.is_paired
            
            # Read pair contains proper mapped pair
            read_pair_dict['read_pair_mapped_in_proper_pair'] = last_read.is_proper_pair
            
            # Read pair is a PCR duplicate
            read_pair_dict['read_pair_is_pcr_duplicate'] = last_read.is_duplicate
            
            # Store everything updating the counters
            c = getUniquelyMappedPairsUpdateCounters(c, read_pair_dict)
 
            # Create new copy of the template dictionary to store
            read_pair_dict = copy.deepcopy(read_pair_dict_template)
            

        # Fill the dictionary with the details for the individual reads
        if read.is_unmapped:      
            read_pair_dict['read_number_unmapped'] += 1
  
        if read.is_supplementary: 
            read_pair_dict['read_number_supplementary'] += 1
            
        # If the read alignment is secondary, it means multimapping is enabled
        # exit with exception
        if read.is_secondary:
            raise Exception("Secondary read found: "+ read.query_name + 
                            " Multimapping was enabled")
             
        # Last command in the loop: update the last_read with
        # the read currently viewed
        last_read = read
        last_read_name = read.query_name
              
    # Process the last read
    # We have seen all the reads from a pair, fill the details
    # of the pair
    
    # reads contain a pair (R1 and R2)
    read_pair_dict['read_pair_is_paired'] = last_read.is_paired
    
    # Read pair contains proper mapped pair
    read_pair_dict['read_pair_mapped_in_proper_pair'] = last_read.is_proper_pair
    
    # Read pair is a PCR duplicate
    read_pair_dict['read_pair_is_pcr_duplicate'] = last_read.is_duplicate
    
    
    # Store everything updating the counters
    c = getUniquelyMappedPairsUpdateCounters(c, read_pair_dict)
    
    createOutputStringFromCounters(c, stats_outfile)


#------------------------------------------------------------------------------
# Updates the counters of reads taking into account the guidelines from
# getUniquelyMappedPairs  
# Inputs:
#    -counter: Counter with the different keys and the counts
#    -dictionary: Dictionary containing the counts for the different keys
#        of the last read pair
#
# Outputs:
#    -updated_counter: Counter after adding the different keys and corresponding counts
def getUniquelyMappedPairsUpdateCounters(counter, dictionary):
    
    # Get proper mapped pairs
    if dictionary["read_pair_is_paired"] and dictionary["read_pair_mapped_in_proper_pair"]:
        
        counter["total_proper_mapped_pairs"] += 1
        
        # For supplementary, for now, all we want to know is whether a read
        # has these alignments or not
        if dictionary['read_number_supplementary'] != 0:  
            counter["proper_mapped_pairs_with_supp_alignment"] += 1

        if dictionary["read_pair_is_pcr_duplicate"]:    
            counter["proper_mapped_pairs_pcr_duplicate"] += 1
            
        # -Read pairs which are all of the below (group1): 
        #        -Read pair mapped in proper pair
        #        -Doesn't contain supplementary alignments
        #        -Doesn't contain secondary alignments
        #        -Not PCR duplicate
        #        -Contains R1 and R2
        if ((dictionary['read_number_supplementary'] == 0) and \
             not dictionary["read_pair_is_pcr_duplicate"]):
            counter["proper_mapped_pairs_with_no_supp_or_pcr_duplicate"] += 1
        
        #    -Read pairs which are all of the below (group2): 
        #        -Read pair mapped in proper pair
        #        -Contains supplementary alignments (maybe translocations, inversions, etc..)
        #        -Doesn't contain secondary alignments
        #        -Not PCR duplicate
        #        -Contains R1 and R2 
        elif ((dictionary['read_number_supplementary'] != 0) and \
             not dictionary["read_pair_is_pcr_duplicate"]):
            counter["proper_mapped_pairs_with_supp_or_pcr_duplicate"] += 1
        
    return counter    


#------------------------------------------------------------------------------
# Writes the counters in the stats_outfile
# Inputs:
#     -counters: Counters with different keys and counts
#     -stats_outfile: File with the output of the stats
def createOutputStringFromCounters(counters, stats_outfile):
    
    with open(stats_outfile, 'w') as writer:
        for key in counters:  
            writer.write(key + ":\t" +str(counters[key]) + "\n")


#--------------------------------------------------------------------------------
# Updates the counters of reads taking into account the guidelines from
# getCorrectReadPairs  
# Inputs:
#    -counter: Counter with the different keys and the counts
#    -dictionary: Dictionary containing the counts for the different keys
#        of the last read pair
#
# Outputs:
#    -updated_counter: Counter after adding the different keys and corresponding counts
def getCorrectReadPairsUpdateCounters(counter, dictionary):
    
    # Variable to control the pair isn't into any "incorrect category"
    correct_pair = True
    
    if (dictionary["read_pair_mapped_in_proper_pair"] and 
        dictionary["read_pair_is_paired"]): 
       
        counter["read_pairs_mapped_in_proper_pair"] += 1
        if dictionary["read_pair_is_pcr_duplicate"]:  
            counter["proper_read_pair_is_pcr_duplicate"] += 1   
            correct_pair = False

        if dictionary['proper_pair_additional_hit_equally_good']:   
            counter["proper_read_pair_additional_hit_equally_good"] += 1     
            correct_pair = False
                
        if dictionary['proper_pair_with_chimeric_read']:    
            counter["proper_read_pair_with_chimeric_read"] += 1       
            correct_pair = False
        
        # After performing all the checks, see which are "correct" pairs    
        if correct_pair:  
            counter["proper_read_pair_no_dup_no_additional_hit_equally_good_no_chimeric"] += 1
    
    # Not proper pairs    
    else:
        counter["read_pairs_mapped_not_in_proper_pair"] += 1

    return counter 


#------------------------------------------------------------------------------
# inspired pipelineCaptCPerl.filterReadMappings
# Assuming a Sam/Bam file resulting from read pairs mapping with BWA mem -M
# and containing unique or multiple mappings (1 pair - multiple alignments)
# per pair in each read and SORTED BY READNAME from any mapper
# output. 
# Obtains the number of correctly mapped read pairs and the number of those pairs
# correctly mapped read pairs which are marked as PCR duplicates
# Note: For each read pair, it only needs one proper read pair mapping in it's alignments
#     to be considered a correctly mapped read pair.
#     As long as one correctly mapped read pair is classified as a PCR duplicate,
#     the read pair will be considered a PCR duplicate.
# 
#
# Inputs:
#     -infile: Sam/Bam file containing only unique mappings per pair in each read
#         sorted by readname sorted by readname from BWA mem -M
#     -stats_outfile: Outfile with the statistics
#
# Output:
#     -returns printable output with the stats
@cluster_runnable
def getCorrectReadPairs(infile, stats_outfile):
    
    c = collections.Counter()
    
    # Determine the extension of the file (BAM or SAM), done this for SAM/BAM compatibility
    # At the moment, all calls should only pass SAM, if it stops working just leave 
    # samfile = pysam.AlignmentFile(infile, "r")
    filename, file_extension = os.path.splitext(infile)
    
    if (file_extension.lower() == '.sam'):
        samfile = pysam.AlignmentFile(infile, "r")
    
    elif (file_extension.lower() == '.bam'):
        samfile = pysam.AlignmentFile(infile, "rb")
    
    
    # We need to know when the read name changes
    last_read_name = "first"
    
    # We need to know the last read when the read name changes to assign
    # the read_pair details from the last read when we have seen all the read pair
    # alignments
    last_read = ""
    
    # Create new dictionary template to store the details
    read_pair_dict_template = {'read_pair_is_paired':False, # Read pair contains R1 and R2    
            'read_pair_mapped_in_proper_pair':False, # Read pair contains proper mapped pair
            'read_pair_is_pcr_duplicate':False, # Read pair is a PCR duplicate
            'proper_pair_additional_hit_equally_good':False, # If additional hits with equally good score are found
            'proper_pair_with_chimeric_read':False, # Chimeric read
            }
    
    # Create new copy of the template dictionary to store
    read_pair_dict = copy.deepcopy(read_pair_dict_template)
            
    
    # Process the reads
    for read in samfile:        
        
        # If it is the first read, assign the name (it will consider
        # it as read name hasn't changed
        if last_read_name == "first":
            last_read_name = read.query_name
            
        
        # The read has changed excepting the first read from the file
        if (last_read_name != read.query_name):
            
            # We have seen all the reads from a pair, update the counters
            
            # Store everything updating the counters
            c = getCorrectReadPairsUpdateCounters(c, read_pair_dict)
            
            
            # Create new copy of the template dictionary to store
            read_pair_dict = copy.deepcopy(read_pair_dict_template)
        
        
        # Checks for the read
        # For every read pair we check whether it is a proper pair, 
        # We only need one mapping from the read pair as a proper pair to
        # consider it a proper pair towards the output
        if(read.is_paired and read.is_proper_pair):
            
            read_pair_dict['read_pair_is_paired'] = read.is_paired
            
            read_pair_dict['read_pair_mapped_in_proper_pair'] = read.is_proper_pair
            
            # We only need one proper pair marked as duplicate for a read pair
            # to consider it a duplicate. 
            if read.is_duplicate:
                read_pair_dict['read_pair_is_pcr_duplicate'] = read.is_duplicate
            
            # Check additional hits when the read pair happens
            if read.has_tag("XA"):
                
                # If there are additional hits
                if (len(read.get_tag("XA").split(","))>1):
                    
                    # If the score on the additional hit (XS) is as good or better
                    # than the score of the original hit (AS)
                    if (read.get_tag("AS") <= read.get_tag("XS")):
                    
                        read_pair_dict['proper_pair_additional_hit_equally_good'] = True
                        
            # Check chimeric hits (split reads)
            if read.has_tag("SA"):
                 
                # If there are chimeric hits
                if (len(read.get_tag("SA").split(","))>1):
                    
                    read_pair_dict['proper_pair_with_chimeric_read'] = True
        
        # Last command in the loop: update the last_read with
        # the read currently viewed
        last_read = read
        
        last_read_name = read.query_name
        
        
    # Process the last read
    # We have seen all the reads from a pair, update the counters
            
    # Store everything updating the counters
    c = getCorrectReadPairsUpdateCounters(c, read_pair_dict)
    
    createOutputStringFromCounters(c, stats_outfile)


#------------------------------------------------------------------------------
# Creates a statement to exclude from the infile the chrs
# indicated in excluded_chrs and saves the processed file to outfile
#
# Inputs:
#    -infile: peaks infile (compressed or uncompressed) to process.
#    -excluded_chrs: list of chrs to exclude (separated by ,)
#    -outfile: processed file without the chrs regions.
#
# Outputs:
#    -Writes the processed file to outfile.
def createExcludingChrFromBedStatement(infile, excluded_chrs, outfile):
    
    statement = ""
    
    # If no excluded chrs provided just copy the file
    if excluded_chrs == "":
        statement = "cp "+infile+" "+outfile
   
    else:        
        
        # See the extension of the file to determine the cat/zcat command
        _, file_extension = os.path.splitext(infile)
        
        if file_extension == ".gz":
            statement = "zcat "  
        else:   
            statement = "cat "
            
        # Create the statement
        statement += infile
        statement += " | egrep -v '("
        statement += excluded_chrs
        statement += ")' | gzip -c > "
        statement += outfile

    return statement


#------------------------------------------------------------------------------
# Given a contig_name, gets the corresponding length of the contig looking into 
# contig_file for the corresponding contig_name. If it doesn't find it, returns
# -1.
#
# contig_file format:
# chr1    2500000
# chr2    170000
# 
# Inputs:
#    -contig_name: Name of the contig
#    -contig_file: Name of the contig file
#
# Outputs:
#    -contig_length: The length of the contig. If it doesn't find it, returns -1.
def getContigLength(contig_name, contig_file):
    
    contig_length = -1
    with open(contig_file, 'r') as reader:
                    
        for line in reader:  
            # Remove the new line character from the line (avoids taking this into account downstream)
            line = line.rstrip('\n') 
            fields = line.split("\t")
            
            if( fields[0] == contig_name):
                # Get the fragment length
                contig_length = fields[1]
                break
            
    return contig_length


#------------------------------------------------------------------------------
# Performs a lot faster if the infile has been previously chromosome sorted            
# First checks that all the bed segments don't surpass the start and end of the chr.
# Regions for which its ends are before the beginning of the chr, are considered empty and disregarded.
# Regions for which its starts are after the end of the chr, are considered empty and disregarded.
# After those filters, starts which are past the beginning of the chr are put at 0
# and ends which are past the end of the chr are put at the chromosome length
# Any regions where the start and end positions mean the region is empty are disregarded.
# Inputs:
#    -infile: bed file with strand information after being processed with bedtools slop.
#             In theory, any bed file should serve since all fields are outputted, this includes bedgraph files.
#             Tested formats:
#             -bed:    chr1    0    100
#             -bedgraph:    chr1    0    100    1.2345
#             -bed6 (with strand): chr1    10213   10256   HSQ-700220:198:C5WA4ACXX:5:1304:17892:7253      7       +
#    -contigs: file with contig \t length pairs for the infile
#    -outfile: processed file eliminating the problems described above
#    -log_file: Any deleted or modified regions are logged here.
#
# Outputs:
#    processed contents to the outfile
#    Exception if a contig is not found in the list of contigs.
@cluster_runnable
def correctSlopChromosomeEdges(infile, contigs, outfile, log_file):

    with IOTools.open_file(outfile, "w") as writer, \
         IOTools.open_file(log_file, "w") as log:

        last_contig = "first"    
        contig_length = -1
    
        # Go through each region
        for bed_region in Bed.iterator(IOTools.open_file(infile, "r")):
            
            # Create and output bed copying the details
            output_bed = bed_region.copy() 
            if bed_region.contig != last_contig:
                # Get the contig length
                contig_length = int(getContigLength(bed_region.contig, contigs))
                last_contig = bed_region.contig
                 
            # If a contig length is not found in the contigs file raise an exception
            if contig_length == -1:        
                raise Exception("Correcting positions, contig "+bed_region.contig+ " not found in the list of contigs: "+contigs)
 
            # Start >= length of chromosome
            if (bed_region.start >= contig_length):
                log.write("Correcting positions, file: "+infile+"."
                          +" Eliminating bed region: "+str(bed_region)+
                          " start >= contig length "+str(contig_length)+ "\n")
                
                # Skip outputting this bed
                continue
            
            # end <= 0
            if (bed_region.end <= 0):
                log.write("Correcting positions, file: "+infile+"."
                          +" Eliminating bed region: "+str(bed_region)+
                          " end <= 0\n")
                
                # Skip outputting this bed
                continue
            
            # End surpassing contig length, correct    
            if (bed_region.end > contig_length):
                
                # The last allowable position is the end of the chromosome base
                # The end coordinate is not included in the region
                output_bed.end = contig_length
                
                log.write("Correcting positions, file: "+infile+"."
                          +" End surpassing length of chromosome. "
                          +"Before: "+str(bed_region)+
                          " Now: "+str(output_bed)+"\n")
                
            
            # End surpassing contig beginning, correct
            if (bed_region.start < 0):
                
                # The last allowable position is the end of the chromosome base
                output_bed.start = 0
                
                log.write("Correcting positions, file: "+infile+"."
                          +" Beginning surpassing beginning of chromosome. "
                          +"Before: "+str(bed_region)+
                          " Now: "+str(output_bed)+"\n")
                
            
            
            
            # At this point all the points in the region can't exceed the chromosome edges, now verify if the
            # region is not empty/is adequate (start < end)
            if output_bed.start >= output_bed.end:
                
                # Print reason to disregard in the log
                log.write("Correcting positions, file: "+infile+
                          "Start equal or larger than end (empty). "
                          +"Start: " +str(output_bed.start)+
                          " End: " +str(output_bed.end)+ "\n")
            
            # If it gets here the bed file can be outputted
            else:
                
                writer.write(str(output_bed))
                writer.write("\n")


# Creates a statement to exclude from the peaks infile the bed regions
# indicated in excluded_beds and saves the processed file to outfile
#
# Inputs:
#    -infile: bed infile (compressed or uncompressed) to process.
#    -excluded_beds: list of chrs to exclude (separated by |)
#    -outfile: processed file without the chrs.
#
# Outputs:
#    -Writes the processed file to outfile.
def createExcludingBedsFromBedStatement(infile, excluded_beds, outfile):
    
    statement = ""
    
    # Get the string into a list, separating the elements
    try:
        list_excluded_beds = excluded_beds.split(",")
    except AttributeError:
        list_excluded_beds = excluded_beds
    
    # If no excluded beds provided just copy the file
    if len(list_excluded_beds) == 0:
        
        statement = "cp "+infile+" "+outfile
    
    
    else:                    
        
        # Create the statement
        # Because sometimes inputing a compressed file to bedtools intersect
        # can cause problems (https://github.com/arq5x/bedtools2/issues/513)
        # The file is zcat first and inputted to bedtools intersect -a stdin
        if infile.endswith(".gz"):
            statement += "zcat "+infile
            
        else:
            statement += "cat "+infile
        
        statement += " | bedtools intersect -v -a stdin"
        statement += " -b "
        
        for excluded_bed in list_excluded_beds:
            statement += excluded_bed
            statement += " "
        
        statement += " | gzip -c > "
        statement += outfile
            
    return statement
                 



# Creates a statement to exclude from the peaks infile the bed regions
# indicated in excluded_beds and saves the processed file to outfile
#
# Inputs:
#    -infile: bed infile (compressed or uncompressed) to process.
#    -excluded_beds: list of chrs to exclude (separated by |)
#    -outfile: processed file without the chrs.
#
# Outputs:
#    -Writes the processed file to outfile.
@cluster_runnable
def createExcludingBedsFromBedStatement(infile, excluded_beds, outfile):
    
    statement = ""
    
    # Get the string into a list, separating the elements
    if isinstance(excluded_beds, str):
        list_excluded_beds = list(excluded_beds)
    else:
        list_excluded_beds = excluded_beds

    # If no excluded beds provided just copy the file
    if len(list_excluded_beds) == 0:
        
        statement = "cp "+infile+" "+outfile
    
    
    else:                    
        
        # Create the statement
        # Because sometimes inputing a compressed file to bedtools intersect
        # can cause problems (https://github.com/arq5x/bedtools2/issues/513)
        # The file is zcat first and inputted to bedtools intersect -a stdin
        if infile.endswith(".gz"):
            
            statement += "zcat "+infile
            
        else:
            
            statement += "cat "+infile
        
        statement += " | bedtools intersect -v -a stdin"
        statement += " -b "
        
        for excluded_bed in list_excluded_beds:
            
            statement += excluded_bed
            statement += " "
        
        
        statement += " | gzip -c > "
        statement += outfile

            
    return statement

def gene_to_transcript_map(infile, outfile):
    '''Parses infile GTF and extracts mapping of transcripts to gene, outputs
    as a tsv file'''
    
    iterator = GTF.iterator(IOTools.open_file(infile))
    transcript2gene_dict = {}

    for entry in iterator:

        # Check the same transcript_id is not mapped to multiple gene_ids!
        if entry.transcript_id in transcript2gene_dict:
            if not entry.gene_id == transcript2gene_dict[entry.transcript_id]:
                raise ValueError('''multipe gene_ids associated with
                the same transcript_id %s %s''' % (
                    entry.gene_id,
                    transcript2gene_dict[entry.transcript_id]))
        else:
            transcript2gene_dict[entry.transcript_id] = entry.gene_id

    with IOTools.open_file(outfile, "w") as outf:
        outf.write("transcript_id\tgene_id\n")
        for key, value in sorted(transcript2gene_dict.items()):
            outf.write("%s\t%s\n" % (key, value))


@cluster_runnable
def merge_counts(infiles, outfiles):
    ''' Take salmon inputs from each file and produce a matrix of counts'''
    transcript_infiles = [x[0] for x in infiles]
    gene_infiles = [x[1] for x in infiles]

    transcript_outfile, gene_outfile = outfiles

    def mergeinfiles(infiles, outfile):
        final_df = pd.DataFrame()

        for infile in infiles:
            tmp_df = pd.read_table(infile, sep="\t", index_col=0)
            final_df = final_df.merge(
                tmp_df, how="outer",  left_index=True, right_index=True)

        final_df = final_df.round()
        final_df.sort_index(inplace=True)
        final_df.to_csv(outfile, sep="\t", compression="gzip")

    mergeinfiles(transcript_infiles, transcript_outfile)
    mergeinfiles(gene_infiles, gene_outfile)