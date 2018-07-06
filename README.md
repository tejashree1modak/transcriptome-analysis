# transcriptome-analysis

Scripts to analyze transcriptomic data produced to study the response of Eastern oyster larvae to bacterial exposures. 

| | Season 2016 | Season 2017 |
|:-:|:-|:-|
|**Site**|Bount shellfish hatchery, RWU, RI|Lab|
|**Species**|*Crassostrea virginica*|*Crassostrea virginica*|
|**Age**|5, 12 and 16 days|8-12 days|
|**Bacterial exposure**|Probiotic Bacillus pumilus RI-0695|Probiotics *Bacillus pumilus* RI-0695, *Phaeobacter inhibens* S4 and pathogen *Vibrio coraliillyticus* RE22|
|**Time of exposure**|5, 12 and 16 days| 6, 24 hours|

### Pipeline for transcriptome data analysis follows steps below:
| Step | Pipelines |
| :- | :- |
| QC | [FASTQC](#fastqc) and [BBTools](#bbtools)| 
| Alignment| [HISAT2](#hisat2) |
| Assembly|  Stringtie |
| Differential expression analysis| DESeq |
| GO annotation | BLAST2GO|
| GO enrichment | topGO  |
| Pathway mapping | KEGG annotation |

---
  
## Pipeline details

### QC 

The QC pipeline supports to scripts - 

#### FASTQC

---

#### BBTools 
  - The [qc_bbtools](qc_bbtools.sh) script contains the pipeline for BBTools 
  - _Usage_: 
    - [SLURM](https://slurm.schedmd.com/) is required to run this script
    - To speed up processing of large jobs, this script requires the jobs to be submitted as [array jobs](https://slurm.schedmd.com/job_array.html)
  - _Input Parameters_:
    - The following paremeters should be set using the `--export=<variable>=<value>` notation of SLURM
 
| Variable | Value                                                                       | Comments                                                                                                                 |
|----------|-----------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------|
| RUN_STEP | One or more of `stat adapter_trim quality_trim force_trim quality_filter` | This specifies which pipeline step to run. Multiple steps can be run by separating them by `:`, e.g. `stat:adapter_trim` |
| SRA_DIR  | The directory containing the paired SRA files                               | These files need to be paired ones with `_1.fastq` and `_2.fastq names`                                                      |
| OUT_DIR  | The directory to generate the output files in                               | For subsequent pipeline steps, the input files would be output of prior steps and would also be in this directory        |

  - _Details of pipeline steps_: 
  
| Step | Details / comments |
| - | - |
| stat | |
| adapter_trim | | 
| quality_filter | | 
| force_trim | |
| quality_trim | |
| force_trim_right | |

  - Example: assume that there are 4 paired SRA files under `/path/to/sra` directory like this - 

```shell
$ cd /path/to
$ tree sra
sra
├── fileA_1.fastq
├── fileA_2.fastq
├── fileB_1.fastq
├── fileB_2.fastq
├── fileC_1.fastq
├── fileC_2.fastq
├── fileD_1.fastq
└── fileD_2.fastq
```

  - Run the pipeline

```shell  
# to run stat and force_trim steps
sbatch -a 1-4 --export=RUN_STEP=stat:force_trim --export=SRA_DIR=/path/to/sra \
        --export=OUT_DIR=/path/to/out qc_bbtools.sh
        
# to run just adapter trim step
sbatch -a 1-4 --export=RUN_STEP=adapter_trim --export=SRA_DIR=/path/to/sra \
        --export=OUT_DIR=/path/to/out qc_bbtools.sh
 ```
 
   - This will generate the output in `OUT_DIR` and the SLURM out files are stored in current directory
---

### Alignment 

Alignment pipeline supports the following scripts - 

#### HISAT2
  - The [align_hisat2](align_hisat2.sh) script contains the pipeline for hisat alignment
    - _Usage_: 
    - [SLURM](https://slurm.schedmd.com/) is required to run this script
    - To speed up processing of large jobs, this script requires the jobs to be submitted as [array jobs](https://slurm.schedmd.com/job_array.html)
  - _Input Parameters_:
    - The following paremeters should be set using the `--export=<variable>=<value>` notation of SLURM
    - hisat

---
### Assembly 

Assembly pipeline supports the following scripts - 

#### Stringtie
  - The [assembly_stringtie](assembly_stringtie.sh) script contains the pipeline for assembly
    - _Usage_: 
    - [SLURM](https://slurm.schedmd.com/) is required to run this script
    - To speed up processing of large jobs, this script requires the jobs to be submitted as [array jobs](https://slurm.schedmd.com/job_array.html) for some steps.
  - _Input Parameters_:
    - The following paremeters should be set using the `--export=<variable>=<value>` notation of SLURM
    - assembly, merge, reestimate
  - Running the last step in this pipeline will generate the count matrices that can be used to perform differential    expression analysis. 
   
---
### Differential expression 

This script is in R and uses DESeq to obtain DEGs from the counts file obtained from stringtie

#### DESeq

- THe [deseq_pipe](deseq_pipe.R) script contains the pipeline for 
  - PCA analysis 
  - DEGs
- It takes as input two files
  - Count matrix file
  - Pheno data: This is a csv file with metadata for the samples
  
---
### GO annotation 

#### blast2GO needs to be run locally

- Download [blast2GO](https://www.blast2go.com/) 
- Load fasta file of DEGs
- Choose using Blast first 
- Mapping to GO IDs
- Export results 

---
### GO enrichment 

GO enrichment will allow us to see significantly enriched GO terms 

#### topGO 

This piipeline runs the Classic Fisher test to find significantly enriched GO terms
- It takes as input: 
  - A tab separated file containing all the genes with the mapped GO IDs one per line 
  - Interesting genes file with p values from DESeq. 
- It will output a csv file with GO terms and the p values after running the Fishers test. 
- You can use this data to generate a figure of enriched GO terms in [REVIGO](http://revigo.irb.hr/)

---
### Pathway mapping  

The online annotation tool KAAS can be used  

#### KAAS annotation server for KO mapping
- This server can be used to first map the DEGs to KO IDs [KAAS](http://www.genome.jp/tools/kaas/)
  - Upload a fasta file of the DEGs. If uploading in nucleotide form select nucleiotide in the form.
  - Choose organism to map to from the provided list
  - Submit job. 
- Use [KEGG mapper](http://www.genome.jp/kegg/tool/map_pathway.html) to observe KO mapping directly on color coded pathways. Can show upregulated and downregulated genes in different color using this tool. 

---


