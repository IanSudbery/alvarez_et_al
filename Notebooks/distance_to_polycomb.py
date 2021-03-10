from cgatcore import pipeline as P
from ruffus import transform, collate, add_inputs, follows, regex, formatter, mkdir
from ruffus.combinatorics import product

PARAMS = P.get_parameters()
SAMPLES = {x.split("\t")[0]:x.strip().split()[2]
           for x in open("chromHMM_samples.txt").readlines()[1:]}

@follows(mkdir("polycomb_regions.dir"))
@transform(
    ["../BP_chromHMM/VISUALIZATION_Blueprint_release_201608/VISUALIZATION_healthy_model/%s_12_12_Blueprint_release_201608_dense.bed" % x
     for x in SAMPLES.keys()],
     regex(".+/(.+)_12_12_Blueprint_release_201608_dense.bed"),
     r"polycomb_regions.dir/\1.bed.gz")
def get_polycomb(infile, outfile):

    statement = '''awk -v OFS="\\t" '$4==1 || $4==2' %(infile)s
                 | gzip > %(outfile)s'''

    P.run(statement)


@follows(mkdir("closest.dir"))
@product(get_polycomb,
         formatter(r"polycomb_regions.dir/(?P<chromhmm>.+).bed"),
         ["../data/pan_MM_differential_ATAC-seq_regions_no_TSS.tsv",
          "../data/MM_subgroup_differential_ATAC-seq_regions_no_TSS.tsv",
          "../data/pan_MM_differential_genes_regulated_differential_enhancers.bed"],
         formatter(r"../data/(?P<peakset>.+)_differential_(?P<target>ATAC-seq|genes)"),
         r"closest.dir/{chromhmm[0][0]}_vs_{peakset[1][0]}_{target[1][0]}.bed.gz")
def get_closest(infiles, outfile):

    polycomb, peaks = infiles
    statement = '''cut -f1-3 %(peaks)s
                 | tail -n+2
                 | sort -k1,1 -k2,2n
                 | uniq
                 | bedtools closest -a - 
                               -b <( zcat %(polycomb)s | sort -k1,1 -k2,2n ) -d
                 | cut -f1-3,13 
                 | gzip > %(outfile)s'''

    P.run(statement)

@collate(get_closest,
         regex(".+/(.+)_vs_(.+).bed.gz"),
         r"\2_closest_merge.tsv")
def merge_samples(infiles, outfile):

    infiles = " ".join(infiles)
    statement = '''cgat combine_tables
                    --regex-filename='.+/(.+)_vs'
                    --cat chromHMM_sample
                    -L %(outfile)s.log
                    -S %(outfile)s
                    --no-titles
                
                    %(infiles)s'''

    P.run(statement)

@follows(merge_samples)
def full():
    pass

P.main()
