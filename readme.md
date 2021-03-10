# Pipeline and code for Alvarez-Benayas, Katsarou and Trasanidis *et al* #

This pipeline is created in the ruffus (http://www.ruffus.org.uk/) and cgatcore (https://cgat-core.readthedocs.io/en/latest/) frameworks. See there documentation to understand how pipeline script files. 

The script `pipeline_atacseq.py` contains the pipeline logic.

The module `pipelineAtacseq.py` contains various work functions that are executed by the pipeline. 

The directory `scripts` contains several ancillary scripts (mostly in R).

The directory `Notebooks` contains Rmarkdown and jupyter notebooks containing interpretive analysis and plotting code. 

## Dependency installation ##
Dependencies are listed in `environment.yml`. This can be used to build a conda environment. We recommend using `mamba` to do this, as it is significantly faster than `conda`.

```
$ mamba env create [-p PATH_TO_STORE_ENV|-n NAME_FOR_ENV] -f environment.yml
```

Next you will need to install `cgat-flow`, which is not yet available on conda. This can be achieved by cloning the github repository and installing it into your environment. 

```
$ conda activate ENV_NAME_OR_PATH
$ git clone https://github.com/cgat-developers/cgat-flow.git
$ cd cgat-flow
$ python setup.py develop
```

## Configuring the pipeline ##

1. To set up the pipeline create a folder and run `python PATH/TO/REPO/pipeline_atac.py config`
2. Edit the configuration file. In particular set the location of the geneset, a samtools indexed copy of version 38 of the human genome, and the ENCODE 3. black list and low mappability regions bed files. 
3. Copy into the directory the sample annotation file in this rep `samples.tsv`.
4. Put all the ATAC seq mapped BAM (mapped with bowtie2) files in the root of the directory. 
5. Put both fastq files and bam files (HiSat mapped) into a folder called `input_rna.dir`

## Running the pipeline ##

You will need to configure the cgatcore frame work to work with your cluster. Do this by creating a `.cgat.yml` file in your home directory. An example is provided below

    jobs_limit_db: 1
    shared_tmpdir: /fastdata/mb1ims/tmp
    tmpdir: /scratch

    cluster:
        queue_manager: sge
        queue: NONE
        parallel_environment: smp
        memory_resource: rmem
        options: -P gen2reg

You can now run show the tasks to be run with `python PATH/TO/REPO/pipeline_atacseq.py show full` or run it with `python PATH/TO/REPO/pipeline_atacseq.py make full`.

On our cluster, a complete run takes around 2.5 days (node availability is never limiting).
