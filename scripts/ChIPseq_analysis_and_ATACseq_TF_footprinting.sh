
########## ChIPseq analysis script ##########

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
treatment_name=$(basename "$treatment" .noall.sorted.bam)
control_name=$(basename "$control" .nodup.sorted.bam)

macs2 callpeak -t $treatment -c $control -g hs -n ./peak_calling/$treatment_name.0.01 -B -q 0.01 --nomodel --verbose 4 --SPMR --call-summits
macs2 callpeak -t $treatment -c $control  -g hs -n ./peak_calling/$treatment_name.0.01 --nomodel -B --broad --broad-cutoff 0.01 --verbose 4 --SPMR
bamCompare --bamfile1 $treatment --bamfile2 $control --binSize 50 --normalizeUsing RPKM --scaleFactorsMethod None -e 200 --operation ratio -o peak_calling/$treatment_name.FC_over_Input.bw

done <$ChIPseq_table


#5. Super-enhancer calling, example for H3K27ac analysis
#1. First make Tag Directories for CHIP AND INPUT files you want to use#
makeTagDirectory "MM1S_H3K27ac_ChIP/" "MM1S_H3K27ac.noall.sorted.bam"
makeTagDirectory "MM1S_Input_ChIP/" "MM1S_Input.noall.sorted.bam"
findPeaks "MM1S_H3K27ac_ChIP" -i MM1S_Input_ChIP -style super -o MM1S_H3K27ac_ChIP.SuperEnhancers.txt -typical MM1S_H3K27ac_ChIP.typical_enhancers.txt -L 0

sed '/^#/ d' MM1S_H3K27ac_ChIP.SuperEnhancers.txt | awk '{printf("%s\t%s\t%s\t%s\t%s\t%s\n",$2,$3,$4,$1,".",$6)}' - > MM1S_H3K27ac_ChIP.SuperEnhancers.bed
sed '/^#/ d' MM1S_H3K27ac_ChIP.typical_enhancers.txt | awk '{printf("%s\t%s\t%s\t%s\t%s\t%s\n",$2,$3,$4,$1,".",$6)}' - > MM1S_H3K27ac_ChIP.typical_enhancers.bed


#6. Motif analysis & genomic annotation
#for a $file = file.bed/file.narrowPeaks/file.broadPeaks 
findMotifsGenome.pl $file hg38 MotifOutput/ -size given -mask
annotatePeaks.pl $file hg38 > outputfile.txt


####### Footprinting analysis #######
#for each subgroup_merged files, $bedfile and $bamfile do
wellington_footprints.py $bedfile $bamfile $outputfolder -fdr 0.05 -fdrlimit -10 -A
awk $5< -10 $$outputfolder/*"-10fdr".bed | awk '{printf("%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4":"$1":"$2":"$3":"$5,$5)}' > output_reformatted


######## Differential footprinting analysis ########
#comparative analysis for HD_high vs HD_low 
wellington_bootstrap.py $HD_HIGH.bam $HD_LOW.bam $HD_high_n_low_merged.bed $HD_high.bed $HD_low.bed -fdr 0.05 -fdrlimit -8 -A 













