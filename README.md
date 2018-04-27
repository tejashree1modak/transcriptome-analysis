# transcriptome-analysis
##Scripts to analyze transcriptomic data produced to study the response of Eastern oyster larvae to bacterial exposures. 

### Season 2016
- **Site**: Bount shellfish hatchery, RWU, RI 
- **Species**: *Crassostrea virginica*
- **Age**: 5, 12 and 16 days 
- **Bacterial exposure**: Probiotic Bacillus pumilus RI-0695
- **Time of exposure**: 5, 12 and 16 days 

### Season 2017
- **Site**: Lab
- **Species**: *Crassostrea virginica*
- **Age**: 8-12 days
- **Bacterial exposure**: Probiotics *Bacillus pumilus* RI-0695, *Phaeobacter inhibens* S4 and pathogen *Vibrio coraliillyticus* RE22
- **Time of exposure**: 6, 24 hours 

##Pipeline for transcriptome data analysis follows steps below: 
- **QC**: FASTQC, BBTools 
- **Alignment**: HISAT2
- **Assembly**: Stringtie
- **Differential expression analysis**: DESeq
- **GO annotation**: BLAST2GO
- **GO enrichment**: topGO 
-**Pathway analysis**: KEGG annotation

  
