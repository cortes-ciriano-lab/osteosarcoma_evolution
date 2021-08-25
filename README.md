# hmf-pipeline
Run the HMF pipeline tools on a tumor-normal pair on the EBI cluster. It includes somatic SNV calling with SAGE, somatic SV calling with GRIDSS (and filtering with GRIPSS), somatic CNV calling with PURPLE and somatic SV interpretation with LINX. It will submit a bunch of jobs and take care of the dependencies. It does not include mapping (yet), so you need to start with BAM files. You can read more about HMFtools at https://github.com/hartwigmedical/hmftools

## Install
Create a conda environment with:
```
conda env create -f hmf.yml
```
(this step might take a while).

When the environment is installed, there is a small hack needed to fix a library dependency (circos is picky about the version but conda does not seem to realize). To get around it you need to find your conda environment and make a symlink:
```
conda activate hmf
condaBin=$(dirname $(which PURPLE))
ln -sf ${condaBin}/../lib/libwebp.so.7 ${condaBin}/../lib/libwebp.so.6
```

## Run
You will need to provide the path to the tumor and normal BAM files and an output directory:
```
sh runHMF.sh -t path/to/tumorBam -n path/to/normalBam -o path/to/outputDir
```
All appropriate sub-directories will be created in outputDir. If needed, you can change a lot of memory, threads and parameters in an ini file. Copy hmf.ini to your preferred dir, modify what you need and then provide the path with:
```
-i path/to/iniFile
```
You can load your conda environment with the hmftools prior to running the script, then add the
```
--condaLoaded
```
argument. Otherwise, you can specify the path to your conda installation and the name of the environment in the ini file and the script will load it for you.
You can also add your email address if you want to receive a message when the pipeline finishes, detailing if the steps were succesfully completed:
```
-m email@ebi.ac.uk
```

To run in the codon cluster, you need to add the path to the codon ini file:
```
-i /path/to/hmf_codon.ini
```

Full usage is:
```
Usage: run_HMF.sh [options] -t <tumor.bam> -n <normal.bam> -o <outputDir>

Required parameters:
	-t/--tumorBam: path to tumor BAM file.
	-n/--normalBam: path to normal BAM file.
	-o/--outputDir: path to output directory (will be created).

Optional parameters:
	-h/--help: show this usage help.
	-i/--iniFile: path to ini file [/hps/research1/icortes/jespejo/hmf-pipeline/hmf.ini]
	-m/--mail: Add an email to send a final report on the pipeline []
	-r/--reference: reference genome to use [/hps/research1/icortes/DATA/hg38/Homo_sapiens_assembly38.fasta]
	--germline: Perform SAGE germline calling
	--id: Specific ID to append to job names [random string]
	--condaLoaded: flag; your env is already pre-loaded, don't load another one [false]
```
Please note that if one step downstream fails, when re-running the pipeline will pick up where it failed.
Here an example of how I ran a PCAWG tumor-normal pair (I used the default ini-file):
```
conda activate hmf
sh /hps/research1/icortes/jespejo/hmf-pipeline/run_HMF.sh \
-t test/DO220842/bam/Tumor_SA557318.sorted.bam \
-n test/DO220842/bam/Normal_SA557554.sorted.bam \
-o test/DO220842/hmf_full/ \
-m jespejo@ebi.ac.uk \
--condaLoaded
```


## What does it do:
The HMF pipeline runs a bunch of specific tools and then everything comes together with PURPLE. Therefore, that's primarily where you need to go for the end files.
#### 1. Preparations
The first step is to check that all the tools needed are in path. It also assigns a random 10-character string to each individual run, to ensure dependencies don't collide. It creates the output directory with a log directory inside, and will write a version.log file with the package versions. It also gets the tumor and normal sample names from the corresponding BAM files.

#### 2. SAGE
SAGE is a somatic SNV, MNV and indel caller. Details are in: https://github.com/hartwigmedical/hmftools/blob/master/sage/README.md
The SAGE output is also annotated with SnpEff and with the HMF-PON.

#### 3. SAGE-GERMLINE
SAGE can also perform germline calling. Details are in: https://github.com/hartwigmedical/hmftools/blob/master/sage/GERMLINE.md
The SAGE output is filtered like mentioned in the link and annotated with SnpEff. 

#### 4. AMBER
Amber checks the BAF of the tumor/normal pair of likely heterozygous loci. Details in: https://github.com/hartwigmedical/hmftools/blob/master/amber/README.md

#### 5. COBALT
Cobalt checks the read depth of the tumor/normal pairs while taking into account GC content. Read more at https://github.com/hartwigmedical/hmftools/blob/master/cobalt/README.md

#### 6. GRIDSS
GRIDSS is an SV caller. It will call SVs on the tumor and the normal jointly, than will be then filtered for somatic SVs downstream. The output is also annotated with repeat regions and viral integration evidence. You can read more about GRIDSS in: https://github.com/PapenfussLab/gridss

#### 7. GRIPSS
GRIPSS applies a set of filtering and post processing steps on GRIDSS paired tumor-normal output to produce a high confidence set of somatic SV for a tumor sample. GRIPSS processes the GRIDSS output and produces a somatic vcf. You can read more at: https://github.com/hartwigmedical/hmftools/blob/master/gripss/README.md

#### 8. PURPLE
PURPLE is a purity-ploidy estimator, but also a CNA caller and integrates all the data from the tools upstream. It will generate annotated somatic SNV and SV VCF files with CNA information. It also generates sample QC (purity, ploidy, WGD, microsatellite status, contamination), detailed purity estimation files, segmented copy number estimation, copy number per gene and a driver catalog. It generates also circos plots with a lot of information, and informative model-fitting charts. Everything is well-explained here: https://github.com/hartwigmedical/hmftools/blob/master/purple/README.md

#### 9. LINX
LINX is an annotation, interpretation and visualisation tool for structural variants. The primary function of LINX is grouping together individual SV calls into distinct events and properly classify and annotating the event to understand both its mechanism and genomic impact. Read more at: https://github.com/hartwigmedical/hmftools/blob/master/sv-linx/README.md

### Clean up
The pipeline will generate Gbs of intermediate files, mostly through GRIDSS. If you are sure your pipeline succesfully finished, you can clean up the output directory removing such files with:
```
sh clean_HMF.sh -o <outputDir>
```
```
Usage: clean_HMF.sh [options] -o <outputDir>
Cleans up the HMF run directory, leaving just the input files needed to regenerate PURPLE and LINX output.

Required parameters:
	-o/--outputDir: path to output directory of the HMF run

Optional parameters:
	-h/--help: show this usage help.
	-f/--force: Ignore done files sanity check and force removal.
```

### Common problems
#### libwep dependency error in circos/purple
Some library is more advanced that circos wants. To get around it you need to find your conda environment and make a symlink:
```
conda activate hmf
condaBin=$(dirname $(which PURPLE))
ln -sf ${condaBin}/../lib/libwebp.so.7 ${condaBin}/../lib/libwebp.so.6
```

#### Something runs out of memory and jobs are trapped in PEND status due to dependency failure
You can always check the status of the dependencies for a particular job with:
```
bjdepinfo <job ID>
```
If a job upstream has run out of memory or failed for whatever other reason there are two options:

1. Remove the job from the queue with
```
bkill <jobID>
```  
and resubmit the pipeline giving it more memory or after solving the error.

2. Run the job manually (commands are in the log folder) and remove the dependencies for the stuck job using:
```
bmodify -wn <jobID>
```

Sometimes the second option is nicer if the last steps of the pipeline fail, if its more upstream it is maybe more complex.
