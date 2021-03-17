# BEIDOU 
**B**ase/Prime **E**ditor **I**nduced **D**NA **O**ff-target site identification **U**nified toolkit

Version: 1.0.0

*About the name: "BeiDou" is also the name of China's navigation satellite system*

-----------------------------------

## Schema
![image](doc/BEIDOU_workflow.001.png)
Authors: Wei Xue(xuewei@picb.ac.cn), Zhi-Can Fu(fuzhican@picb.ac.cn) and Li Yang (liyang@picb.ac.cn)
Maintainer: Zhi-Can Fu(fuzhican@picb.ac.cn)

-----------------------------------

## Installation requirements [mandatory]
* Software
    - [bwa](https://github.com/lh3/bwa) [(version 0.7.17-r1188)](https://github.com/lh3/bwa/releases/tag/v0.7.17)
    - perl (version 5.26.2)
    - samtools (version 1.9)
    - bedtools (version 2.28.0)
    - picard (version 2.21.2)
    - bamtools (version 2.5.1)
    - bcftools (version 1.9)
    - GATK (version 4.1.3.0-0)
    - lofreq (version 2.1.3.1) 
    - Strelka2 (version 2.9.10)
    - Scalpel (version 0.5.4)
    - Manta (version 1.6.0)
    - GNU Parallel (version 20200722)
    - R (version 3.5.1)

## Data requirements [mandatory]
* reference sequences
    - **[hg38.fa](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz)|[mm10.fa](https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz)**
    - **hg38_all.fai|mm10_all.fai** (Created by "samtools faidx")
    - **hg38_all.dict|mm10_all.dict** (Created by "picard CreateSequenceDictionary")
* vcfs for GATK BaseRecalibrator
    - **[NCBI_dbSNP_all_hg38.vcf.gz](https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz)|[EVA_SNP_all_mm10.vcf.gz](http://ftp.ebi.ac.uk/pub/databases/eva/rs_releases/release_1/by_species/Mouse_10090/GRCm38.p4/GCA_000001635.6_current_ids.vcf.gz)** (Very importantly, chromosome names in the annotations GTF file have to match chromosome names in the FASTA genome sequence files)
* vcfs for GATK VariantRecalibrator (Human)
    - **[hapmap_3.3.hg38.vcf.gz](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz)** 
    - **[1000G_omni2.5.hg38.vcf.gz](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz)** 
    - **[1000G_phase1.snps.high_confidence.hg38.vcf.gz](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz)** 
    - **[NCBI_dbSNP_all_hg38.vcf.gz](https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz)** 
    - **[Mills_and_1000G_gold_standard.indels.hg38.vcf.gz](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz)** 
* vcfs for GATK VariantRecalibrator (Mouse)
    - **[EVA_SNP_all_mm10.vcf.gz](http://ftp.ebi.ac.uk/pub/databases/eva/rs_releases/release_1/by_species/Mouse_10090/GRCm38.p4/GCA_000001635.6_current_ids.vcf.gz)** 
    - **MGP_SNP_indel_v5.vcf.gz** (merged from [MGP_SNP_v5.vcf.gz](ftp-mouse.sanger.ac.uk/current_snps/mgp.v5.merged.snps_all.dbSNP142.vcf.gz) and [MGP_indel_v5.vcf.gz](ftp-mouse.sanger.ac.uk/current_snps/mgp.v5.merged.indels.dbSNP142.normed.vcf.gz))
    - **[MGP_indel_v5.vcf.gz](ftp-mouse.sanger.ac.uk/current_snps/mgp.v5.merged.indels.dbSNP142.normed.vcf.gz)** 
* files used to filter out the background variants (Human)
    - **[NCBI_dbSNP_all_hg38.vcf.gz](https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz)** 
    - **[UCSC_RepeatMask_hg38.bed](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz)** 
* files used to filter out the background variants (Mouse)
    - **[EVA_SNP_all_mm10.vcf.gz](http://ftp.ebi.ac.uk/pub/databases/eva/rs_releases/release_1/by_species/Mouse_10090/GRCm38.p4/GCA_000001635.6_current_ids.vcf.gz)** 
    - **[UCSC_RepeatMask_mm10.bed](http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/rmsk.txt.gz)** 

## Installation
```bash
git clone https://github.com/YangLab/BEIDOU.git
```

## Usage
```bash
BEIDOU -f Function -1 Path_of_fastq1 -2 Path_of_fastq2 -o Output_path -n Output_name -c Path_of_config_file -t number_of_maximum_threads -g genome_build_version -d tmp_folder
       [-f Function, "SNV", "Indel" or "all_steps"(default all_steps)]
       [-1 Path of fastq1]
       [-2 Path of fastq2]
       [-o Output directory(default current directory)]
       [-n Output name]
       [-c Path of config file(default ./BEIDOU_config_GENOME_BUILD_VERSION)]
       [-t Maximum_threads]
       [-g Genome build version, "hg38" or "mm10"]
       [-d wirtable temporary folder(default [Output directory]/BEIDOU_tmp)]
```


## Output
* **[Output directory]/[Output name]_Novel_SNVs** and/or **[Output directory]/[Output name]_Novel_Indels** is the result of BEIDOU pipeline.

-----------------------------------

## Example of BEIDOU_config file
```bash
#[Software]
dir_of_bwa=~/bin
dir_of_samtools=~/bin
dir_of_gatk=~/bin
dir_of_picard=~/picard/build/libs
dir_of_bamtools=~/bin
dir_of_bcftools=~/bin
dir_of_lofreq=~/lofreq_star-2.1.3.1/bin
dir_of_Strelka2=~/strelka-2.9.10.centos6_x86_64/bin
dir_of_Scalpel=~/scalpel-0.5.4
dir_of_Manta=~/bin
dir_of_perl=~/bin
dir_of_parallel=~/bin
dir_of_intersectBed=~/bin
#[ref_genome]
ref_genome_path=path/to/hg38_all.fa
dict_of_ref_genome_path=path/to/hg38_all.dict
#[vcf_files_for_BaseRecalibrator]
dbsnp_vcf_for_BaseRecalibrator=path/to/NCBI_dbSNP_all_hg38.vcf.gz
#[vcf_files_for_VariantRecalibrator]
hapmap_vcf=path/to/hapmap_3.3.hg38.vcf.gz
file_1000G_omni_vcf=path/to/1000G_omni2.5.hg38.vcf.gz
file_1000G_phase1_vcf=path/to/1000G_phase1.snps.high_confidence.hg38.vcf.gz
dbsnp_vcf=path/to/dbsnp_146.hg38.vcf.gz
Mills_and_1000G_vcf=path/to/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
#[filtering file]
filtering_dbSNP_vcf=path/to/NCBI_dbSNP_all_hg38.vcf.gz
UCSC_RepeatMask_bed=path/to/UCSC_RepeatMask_hg38.bed
#[optional files(these files can be created automatic)]
dir_of_individual_chr_genome_range_bed=path/to/dir
dir_of_individual_chr_ref_genome_path=path/to/dir
```
-----------------------------------

## Citation
Runze Gao#, Zhi-Can Fu#, Xiangyang Li#, Ying Wang#, Jia Wei, Guangye Li, Lijie Wang, Jing Wu, Wei Xue, Xingxu Huang\*, Li Yang\*, Jia Chen\*. Background levels of genome-wide off-target mutations induced by prime editor. 2020, xxxxxx

-----------------------------------

## License
Copyright (C) 2020 YangLab. Licensed GPLv3 for open source use or contact YangLab (yanglab@picb.ac.cn) for commercial use.
