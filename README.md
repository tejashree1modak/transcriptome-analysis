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
| QC | [FASTQC](#FASTQC) and [BBTools](#BBTools)| 
| Alignment| [HISAT2](#HISAT2) |
| Assembly|  Stringtie |
| Differential expression analysis| DESeq |
| GO annotation | BLAST2GO|
| GO enrichment | topGO  |
| Pathway analysis | KEGG annotation |

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

---
