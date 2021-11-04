#!/bin/bash
########## ChIPseq analysis script ##########

#1. Dependencies 
# Python 2.7.4
# bowtie2 v2.4 -> conda install bowtie2
# picard tools-> install from https://broadinstitute.github.io/picard/
# samtools -> conda install samtools
# deeptools -> pip install deeptools
# macs2 -> pip install MACS2
# homer4.9 tools -> install from http://homer.ucsd.edu/homer/download.html
# pybedtools -> pip install pybedtools

#2. Main scripts
#enter fastq files directory
export fastq_files = $fastq_files_directory
export genome_dir = $bowtie_GRCh38_genome_directory

> Stats_bamfiles_depth.table
mkdir Analysis
mkdir ./Analysis/fastqc_reports

for file in ./$fastq_files/*.fastq.gz; do
	name=$(basename $file .fastq.gz)
	#echo $file; echo $name

#1. fastqc 
	fastqc $file -o './Analysis/fastqc_reports/'

#2 Alignment
	bowtie2 -x $genome_dir -U $file  -p 8 -S ./Analysis/$name.sam 2> ./Analysis/$name.log
    
cd ./Analysis/
##3. Run main analysis
	samtools view -bSo $name.bam $name.sam 
	samtools sort $name.bam -o $name.sorted.bam
	picard MarkDuplicates INPUT=$name.sorted.bam OUTPUT=$name.nodup.bam REMOVE_DUPLICATES=true METRICS_FILE=$namenodup_metrics_file.txt
	picard BuildBamIndex INPUT=$name.nodup.bam
	samtools sort  $name.nodup.bam -o $name.nodup.sorted.bam
 	samtools index $name.nodup.sorted.bam
	samtools view -c $name.nodup.sorted.bam >> Stats_bamfiles_depth.table
	cd ..
done 

#4 Peak calling 
export ChIPseq_table = $ChIPseq_table_samples_match_table
mkdir peak_calling

while read -r treatment control ; do
echo $treatment $control 
treatment_name=$(basename "$treatment" .nodup.sorted.bam)
control_name=$(basename "$control" .nodup.sorted.bam)

macs2 callpeak -t $treatment -c $control -g hs -n ./peak_calling/$treatment_name.0.01 -B -q 0.01 --nomodel --verbose 4 --SPMR --call-summits
macs2 callpeak -t $treatment -c $control  -g hs -n ./peak_calling/$treatment_name.0.01 --nomodel -B --broad --broad-cutoff 0.01 --verbose 4 --SPMR
bamCompare --bamfile1 $treatment --bamfile2 $control --binSize 50 --normalizeUsing RPKM --scaleFactorsMethod None -e 200 --operation ratio -o peak_calling/$treatment_name.FC_over_Input.bw

done <$ChIPseq_table


#5. Super-enhancer calling, example for H3K27ac analysis
#a. First make Tag Directories for CHIP AND INPUT files you want to use#
makeTagDirectory "MM1S_H3K27ac_ChIP/" "MM1S_H3K27ac.nodup.sorted.bam"
makeTagDirectory "MM1S_Input_ChIP/" "MM1S_Input.nodup.sorted.bam"

#b. Then perform super-enhancer analysis
findPeaks "MM1S_H3K27ac_ChIP" -i MM1S_Input_ChIP -style super -o MM1S_H3K27ac_ChIP.SuperEnhancers.txt -typical MM1S_H3K27ac_ChIP.typical_enhancers.txt -L 0

sed '/^#/ d' MM1S_H3K27ac_ChIP.SuperEnhancers.txt | awk '{printf("%s\t%s\t%s\t%s\t%s\t%s\n",$2,$3,$4,$1,".",$6)}' - > MM1S_H3K27ac_ChIP.SuperEnhancers.bed
sed '/^#/ d' MM1S_H3K27ac_ChIP.typical_enhancers.txt | awk '{printf("%s\t%s\t%s\t%s\t%s\t%s\n",$2,$3,$4,$1,".",$6)}' - > MM1S_H3K27ac_ChIP.typical_enhancers.bed


#6. Motif analysis & genomic annotation
#for a $file = file.bed/file.narrowPeaks/file.broadPeaks 
findMotifsGenome.pl $file hg38 MotifOutput/ -size given -mask
annotatePeaks.pl $file hg38 > outputfile.txt


####### Footprinting analysis #######
#Example for construction of MMSET network. 
export bedfile=$merged_peaks_bedfile
export $bamfile=$merged_reads_bamfile
wellington_footprints.py $bedfile $bamfile $outputfolder -fdr 0.05 -fdrlimit -10 -A
awk $5< "-10" $$outputfolder/*"-10fdr".bed | awk '{printf("%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4":"$1":"$2":"$3":"$5,$5)}' > MMSET_FT_reformatted.bed

#obtain TF motif occurances on footprints
export motif_database=HOCOMOCOv11_full_HUMAN_mono_homer_format_0.0001.motif
perl findMotifsGenome.pl MMSET_FT_reformatted.bed hg38 $outputfolder -size given -mask -find $motif_database > MMSET_FT_motif.txt
cut -f1,6 MMSET_FT_motif.txt | sed -e '1d' | awk -F":" '$1=$1' OFS="\t" - | awk -F"-" '$1=$1' OFS="\t" - | awk '{printf("%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,"'$textname'",$4)}' -  > MMSET_FT_motif.bed
bedtools sort MMSET_FT_motif.bed >ALL_FTs_MMSET_sorted.bed

#obtain TF motif occurances on peaks (background)
export motif_database=HOCOMOCOv11_full_HUMAN_mono_homer_format_0.0001.motif
perl findMotifsGenome.pl $bedfile hg38 $outputfolder -size given -mask -find $motif_database > MMSET_peaks_motif.txt
cut -f1,6 MMSET_peaks_motif.txt | sed -e '1d' | awk -F":" '$1=$1' OFS="\t" - | awk -F"-" '$1=$1' OFS="\t" - | awk '{printf("%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,"'$textname'",$4)}' -  > MMSET_peaks_motif.bed

#calculate scores
cut -f1  MMSET_peaks_motif.bed | sort - | uniq -c - > MMSET_TFs.txt
while read TF in MMSET_TFs.txt; do
echo $TF >>Stats.txt
echo "No. of occurances on FTs"
grep $TF MMSET_FT_motif.bed| wc -l - >>Stats.txt
echo "No. of occurances on peaks"
grep $TF MMSET_peaks_motif.bed| wc -l - >>Stats.txt
done

## Construct Gene Regulatory Networks per subgroup
# Use the GRN_script.R  

######## Differential footprinting analysis ########
#comparative analysis for HD_high vs HD_low 
wellington_bootstrap.py $HD_HIGH.bam $HD_LOW.bam $HD_high_n_low_merged.bed $HD_high.bed $HD_low.bed -fdr 0.05 -fdrlimit -8 -A 













