#!/usr/bin/env bash  
set -eo pipefail  
# File Name: xw_work_PGM3_01_bwa_mem_wgs_hg38_20_02_20.sh  
# Author: XueWei  
# mail: xuewei@picb.ac.cn  
# Created time:  
# Last modified: Tue May 26 22:24:19 CST 2020  
###### Fri Jul 17 13:24:44 CST 2020  
###############Usage############  
#bash xw_work_PGM3_01_bwa_mem_wgs_hg38_20_02_20.sh -1 $fq1 -2 $fq2 -o $work_path -n $name  
#ex:  
#bash xw_work_PGM3_01_bwa_mem_wgs_hg38_20_02_20.sh -1 /picb/rnomics3/reads/2019/20191016_WGS_WLJ_NEB_Xten/20191016_0919_WLJ_WGS_01_19530_R1.fastq.gz -2 /picb/rnomics3/reads/2019/20191016_WGS_WLJ_NEB_Xten/20191016_0919_WLJ_WGS_01_19530_R2.fastq.gz -o /data/rnomics6/fuzhican/project/PE_off_target/20_02_25_miniature_run -n 20191016_0919_WLJ_WGS_01_19530  
###############Usage###########  
###Prerequisites  
    ##Softwares and Packages  
    #bwa (0.7.17-r1188)  
    #samtools (1.9)  
    #gatk (gatk4-4.1.3.0-0)  
    #picard (2.21.2)  
    #lofreq (2.1.3.1)  
    #Strelka2 (2.9.10)  
    #Scalpel (v0.5.4)  
    #bamtools (2.5.1)  
    #Manta (1.6.0)  

    ##ref_genome_path  
    #hg38_all.fa  
    ##vcf files for VariantRecalibrator(ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/)  
    #/picb/rnomics3/xuew/Human/backup/GATK_bundle/hg38/hapmap_3.3.hg38.vcf.gz  
    # /picb/rnomics3/xuew/Human/backup/GATK_bundle/hg38/1000G_omni2.5.hg38.vcf.gz  
    # /picb/rnomics3/xuew/Human/backup/GATK_bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz  
    # /picb/rnomics3/xuew/Human/backup/GATK_bundle/hg38/dbsnp_146.hg38.vcf.gz  
    #/picb/rnomics3/xuew/Human/backup/GATK_bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz  

    ##bed for scalpel-discovery  
    #/picb/rnomics3/xuew/Human/backup/hg38_common/hg38_all.bed  
    #/picb/rnomics3/xuew/Human/backup/hg38_common/${chr_n}.bed  

while getopts :1:2:o:n:m:c:t:p:s:g:f:h: ARGS  
    do  
    case $ARGS in   
        1)  
            fq1=$OPTARG  
            ;;  
        2)  
            fq2=$OPTARG  
            ;;  
        o)  
            work_path=$OPTARG  
            ;;  
        n)  
            name=$OPTARG  
            ;;  
        c)  
            tmp_config_file=$OPTARG  
            ;;  
        g)  
            genome_build_version=$OPTARG  
            ;;  
        t)  
            threads=$OPTARG  
            ;;  
        p)  
            Patch=$OPTARG  
            ;;  
        s)  
            Patch_For_Scapel=$OPTARG  
            ;; 
        m)  
            mutation_type=$OPTARG  
            ;;  
        *)  
            echo "Unknown option: $ARGS"  
            ;;  
        \?)  
        echo "Invalid option: -$OPTARG"  
        ;;  
    esac  
    done  
#fq1="/picb/rnomics3/reads/2019/20191016_WGS_WLJ_NEB_Xten/20191016_0919_WLJ_WGS_01_19530_R1.fastq.gz"  
#fq2="/picb/rnomics3/reads/2019/20191016_WGS_WLJ_NEB_Xten/20191016_0919_WLJ_WGS_01_19530_R2.fastq.gz"  
#work_path="/data/rnomics6/fuzhican/project/PE_off_target/20_02_25_miniature_run"  
#name="20191016_0919_WLJ_WGS_01_19530"  
this_software_name=BEIDOU  
tmp_a=aa${Patch}aa  
if [ "${tmp_a}" == "aaTrueaa" ];then  
Patch="True"  
else  
Patch="False"  
fi  
if [ "aa${mutation_type}aa" == "aaaa" ];then  
mutation_type="all"  
fi  
if [ "aa${Patch_For_Scapel}aa" == "aaaa" ];then  
Patch_For_Scapel="False"  
fi 

echo "mutation_type:$mutation_type"  
ref_genome=$genome_build_version  
##To do  
##modification note  
#install.packages("proto", repos="https://cran.rstudio.com/")  
#packageurl <- "https://cran.r-project.org/src/contrib/Archive/ggplot2/ggplot2_1.0.1.tar.gz"  
#install.packages(packageurl, repos=NULL, type="source")  
metc(){  
    max_threads=$1  
    if [ ! "$thread_ctrl_n" ];then  
    thread_ctrl_n=0  
    fi  
    let thread_ctrl_n+=1  
    if [ "$thread_ctrl_n" -eq "$max_threads" ];then  
    wait  
    thread_ctrl_n=0  
    fi  
}  
memkdir(){  
    for path in $@  
    do  
    test -d $path || mkdir -p $path  
    done  
}  
merm(){  
    for file1 in $@  
    do  
    test -e $file1 && rm -r $file1 || echo " " 
    done  
}  
#test -z "$dir_of_individual_chr_genome_range_bed" && dir_of_individual_chr_genome_range_bed=  
step01_06_construct_bam(){  
    ###### Sun Mar 15 19:49:45 CST 2020  
    block_output=${work_path}/${name}_02_bwa_mem_PE.bam  
    test -s $block_output -a "${Patch}" == "True" ||{  
    # /data/rnomics6/xuew/Human/BE_WangLiJie/SBE_genome_20191023/mapping/01_C9_11/  
    # 1. bwa mem  
    ${dir_of_bwa}/bwa mem -t $threads -R '@RG\tID:2019\tPL:illumina\tLB:library\tSM:'"${name}" ${ref_genome_path} $fq1 $fq2 > ${work_path}/${name}_01_bwa_mem_PE.sam 2> >(tee ${work_path}/${name}_01_bwa_mem_PE.log >&2)  
    # 2. sam2bam  
    ${dir_of_samtools}/samtools view -bh -F 4 -q 30 ${work_path}/${name}_01_bwa_mem_PE.sam |samtools sort -@ $threads -o ${work_path}/${name}_02_bwa_mem_PE.bam  
    ${dir_of_samtools}/samtools flagstat -@ 3 ${work_path}/${name}_02_bwa_mem_PE.bam > ${work_path}/${name}_02_bwa_mem_PE_flagstat.log &  
    merm ${work_path}/${name}_01_bwa_mem_PE.sam  
    }  

    # 5. Duplicates Marking  

    # https://www.jianshu.com/p/938d362fc48d  
    # https://www.plob.org/article/11698.html  
    block_output=${work_path}/${name}_06_BQSR.bam  
    test -s $block_output -a "${Patch}" == "True" ||{  
    java -jar ${dir_of_picard}/picard.jar MarkDuplicates REMOVE_DUPLICATES=false I=${work_path}/${name}_02_bwa_mem_PE.bam O=${work_path}/${name}_05_2_markdup.bam M=${work_path}/${name}_05_2_markdup_metrics.txt  
    }  
    # 6. BaseRecalibrator  
    block_output=${work_path}/${name}_06_BQSR.bam  
    test -s $block_output -a "${Patch}" == "True" ||{  
    ${dir_of_gatk}/gatk BaseRecalibrator -I ${work_path}/${name}_05_2_markdup.bam -R ${ref_genome_path} -O ${work_path}/${name}_06_BQSR.table --known-sites $dbsnp_vcf_for_BaseRecalibrator  
    ${dir_of_gatk}/gatk ApplyBQSR -I ${work_path}/${name}_05_2_markdup.bam -R ${ref_genome_path} -bqsr ${work_path}/${name}_06_BQSR.table -O ${work_path}/${name}_06_BQSR.bam  
    }  
    if [ -s ${work_path}/${name}_06_BQSR.bam ] && [ -s ${work_path}/${name}_05_2_markdup.bam ];then
    merm ${work_path}/${name}_05_2_markdup.bam 
    fi

}  
cpu_info(){  
    ###### Sat Sep 19 17:10:37 CST 2020  
    total_threads=`grep 'processor' /proc/cpuinfo | sort -u | wc -l`  
    run_threads=`cat /proc/loadavg|cut -f1 -d" "`  
    idle_threads=`echo $total_threads-$run_threads|bc`  
    #echo `echo $run_threads/$total_threads|bc`  
    echo $idle_threads  
}  
threads_ctrl_for_step7_HaplotypeCaller(){  
    ###### Sat Sep 19 17:19:26 CST 2020  
    count_hap=$(echo "scale=0;($(cpu_info)-10)/2"|bc)  
    if [[ $count_hap -gt 10 ]];then  
    count_threads_ctrl_for_step7_HaplotypeCaller=10  
    else  
    count_threads_ctrl_for_step7_HaplotypeCaller=`echo "scale=0;$count_hap/1"|bc`  
    fi  
    echo $count_threads_ctrl_for_step7_HaplotypeCaller  
}  
threads_ctrl_for_step12_scalpel_indels(){  
    ###### Sat Sep 19 17:19:26 CST 2020  
    count_idles=$(echo "scale=0;($(cpu_info)-20)/3"|bc)  
    if [[ $count_idles -gt $threads ]];then  
    count_threads_ctrl_for_step12_scalpel_indels=$threads  
    else  
    count_threads_ctrl_for_step12_scalpel_indels=`echo "scale=0;$count_idles/1"|bc`  
    fi  
    echo $count_threads_ctrl_for_step12_scalpel_indels  
}  

step7_HaplotypeCaller(){  
    # 7. gatk HaplotypeCaller  
    # 7.1 for only one bam  
    block_output=${work_path}/${name}_07_HC.vcf.gz  
    test -s $block_output -a "${Patch}" == "True" ||{  
    echo "##Time `date +"%R %d %m"` hap begining"  

    if [ "$genome_build_version" == "hg38" ];then  
    {  
    chrn_list=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM)  
    }  
    elif [ "$genome_build_version" == "mm10" ];then  
    chrn_list=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY chrM)  
    fi  
    #echo ${chrn_list[@]}|awk 'BEGIN{RS=" ";ORS="\" -L \""}{print}'  
    

    hap_out=$block_output  
    input_variant_files=${work_path}/hap_split_chr/${name}_input_variant_files  
    memkdir ${work_path}/hap_split_chr  

    for chrn in ${chrn_list[@]}  
    do  
    hap_out_split_chr=${work_path}/hap_split_chr/${name}_${chrn}_07_HC.vcf.gz  
    #echo ${dir_of_gatk}/gatk --java-options \"-Xmx30g -Xms10g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true\" HaplotypeCaller -I $hap_in_bam -R ${ref_genome_path} -O $hap_out_split_chr -L $chrn  
    ${dir_of_gatk}/gatk --java-options "-Xmx30g -Xms10g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" HaplotypeCaller -I ${work_path}/${name}_06_BQSR.bam -R ${ref_genome_path} -O $hap_out_split_chr -L $chrn &  
    metc `threads_ctrl_for_step7_HaplotypeCaller`  
    done  
    wait  
    ls ${work_path}/hap_split_chr/${name}_*_07_HC.vcf.gz >$input_variant_files  
    java -jar ${dir_of_picard}/picard.jar MergeVcfs I=$input_variant_files O=$hap_out  

    chrn_sum_infile=`${dir_of_bcftools}/bcftools view -H $hap_out|cut -f1|uniq|sort|uniq|wc -l`  
    [[ $chrn_sum_infile -eq ${#chrn_list[@]} ]] && rm -r ${work_path}/hap_split_chr || {  
        mv $hap_out ${hap_out}_chrn_error  
    }  

    #${dir_of_gatk}/gatk --java-options "-Xmx30g -Xms10g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" HaplotypeCaller -I ${work_path}/${name}_06_BQSR.bam -R ${ref_genome_path} -O ${work_path}/${name}_07_HC.vcf.gz  
    echo "##Time `date +"%R %d %m"` hap finished"  
    }  


    if ([ ! -s "${work_path}/${name}_09_SNVs_VQSR.vcf.gz" ] || [ ! -s "${work_path}/${name}_09_indels_VQSR.vcf.gz" ]) || ( [ "${Patch}" != "True" ]) ;then  
    {  
    # 8. VQSR  
    # 8.1 SNP  
    # VariantRecalibrator SNP  
    if [ "$genome_build_version" == "hg38" ];then  
    {  
        echo "##Time `date +"%R %d %m"` VariantRecalibrator begining"  
        if ([ ! -s "${work_path}/${name}_09_SNVs_VQSR.vcf.gz" ]) || ([ "${Patch}" != "True" ]) ;then  
        ${dir_of_gatk}/gatk --java-options "-Xmx60g -Xms10g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" VariantRecalibrator -R ${ref_genome_path} -V ${work_path}/${name}_07_HC.vcf.gz -O ${work_path}/${name}_08_snps.recal --resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap_vcf --resource:omni,known=false,training=true,truth=false,prior=12.0 $file_1000G_omni_vcf --resource:1000G,known=false,training=true,truth=false,prior=10.0 $file_1000G_phase1_vcf --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp_vcf -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP --max-gaussians 4 --tranches-file ${work_path}/${name}_08_snps.tranches --rscript-file ${work_path}/${name}_08_snps.R  -L "chr1" -L "chr2" -L "chr3" -L "chr4" -L "chr5" -L "chr6" -L "chr7" -L "chr8" -L "chr9" -L "chr10" -L "chr11" -L "chr12" -L "chr13" -L "chr14" -L "chr15" -L "chr16" -L "chr17" -L "chr18" -L "chr19" -L "chr20" -L "chr21" -L "chr22" -L "chrX" -L "chrY"  -L "chrM"  2> >(tee ${work_path}/${name}_08_snps.log >&2)  ###### Thu Apr 23 08:56:09 CST 2020 fzc add -L ;###### Wed Sep 30 08:50:47 CST 2020 remove -L "chrM" ###### Thu Oct 8 10:54:57 CST 2020 add back -L "chrM"  
        # ApplyVQSR SNP  
        echo "##Time `date +"%R %d %m"` ApplyVQSR begining"  
        ${dir_of_gatk}/gatk --java-options "-Xmx30g -Xms10g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" ApplyVQSR -R ${ref_genome_path} -V ${work_path}/${name}_07_HC.vcf.gz --recal-file ${work_path}/${name}_08_snps.recal -O ${work_path}/${name}_08_snps_VQSR.vcf --tranches-file ${work_path}/${name}_08_snps.tranches -mode SNP -ts-filter-level 95 -L "chr1" -L "chr2" -L "chr3" -L "chr4" -L "chr5" -L "chr6" -L "chr7" -L "chr8" -L "chr9" -L "chr10" -L "chr11" -L "chr12" -L "chr13" -L "chr14" -L "chr15" -L "chr16" -L "chr17" -L "chr18" -L "chr19" -L "chr20" -L "chr21" -L "chr22" -L "chrX" -L "chrY"  -L "chrM" 
        fi  
        if [ "$mutation_type" != "SNV" ];then  
        # 8.2 INDEL  
        # VariantRecalibrator INDEL  
        echo "##Time `date +"%R %d %m"` VariantRecalibrator indel begining"  
        ${dir_of_gatk}/gatk --java-options "-Xmx30g -Xms10g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" VariantRecalibrator -R ${ref_genome_path} -V ${work_path}/${name}_08_snps_VQSR.vcf -O ${work_path}/${name}_08_indels.recal -resource:mills,known=true,training=true,truth=true,prior=12.0 $Mills_and_1000G_vcf -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP --max-gaussians 4 -mode INDEL --tranches-file ${work_path}/${name}_08_indels.tranches --rscript-file ${work_path}/${name}_08_indels.R -L "chr1" -L "chr2" -L "chr3" -L "chr4" -L "chr5" -L "chr6" -L "chr7" -L "chr8" -L "chr9" -L "chr10" -L "chr11" -L "chr12" -L "chr13" -L "chr14" -L "chr15" -L "chr16" -L "chr17" -L "chr18" -L "chr19" -L "chr20" -L "chr21" -L "chr22" -L "chrX" -L "chrY"   -L "chrM"  2> >(tee ${work_path}/${name}_08_indels.log >&2)  
        # ApplyVQSR INDEL  
        echo "##Time `date +"%R %d %m"` ApplyVQSR indel begining"  
        ${dir_of_gatk}/gatk --java-options "-Xmx30g -Xms10g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" ApplyVQSR -R ${ref_genome_path} -V ${work_path}/${name}_08_snps_VQSR.vcf --recal-file ${work_path}/${name}_08_indels.recal -O ${work_path}/${name}_08_indels_VQSR.vcf --tranches-file ${work_path}/${name}_08_indels.tranches -mode INDEL -ts-filter-level 95 -L "chr1" -L "chr2" -L "chr3" -L "chr4" -L "chr5" -L "chr6" -L "chr7" -L "chr8" -L "chr9" -L "chr10" -L "chr11" -L "chr12" -L "chr13" -L "chr14" -L "chr15" -L "chr16" -L "chr17" -L "chr18" -L "chr19" -L "chr20" -L "chr21" -L "chr22" -L "chrX" -L "chrY"  -L "chrM"
        fi  
    }  
    elif [ "$genome_build_version" == "mm10" ];then  
    {  
        if ([ ! -s "${work_path}/${name}_09_SNVs_VQSR.vcf.gz" ]) || ( [ "${Patch}" != "True" ]) ;then  
        #Raw vcf files from variant calling step for all chromosomes except chromosome Y were pooled together for variant quality score recalibration (VQSR) using GATK’s VariantRecalibrator under SNP mode. Training, known and true sets for building the positive model are the SNPs which segregate among the classical laboratory strains of the Mouse Genomes Project 11 (2011 release REL-1211) on all chromosomes except chromosome Y. #Nat Genet. 2016 PMID: 27376238  
        #All animals were jointly genotyped using GenotypeGVCFs, and variant quality score recalibration was performed separately for SNVs and indels using the VariantRecalibrator and ApplyRecalibration, using dbSNP for mouse v142 as the truth set of known variants. # Sci Rep. 2019 PMID: 31551502  
        ${dir_of_gatk}/gatk --java-options "-Xmx60g -Xms10g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" VariantRecalibrator -R ${ref_genome_path} -V ${work_path}/${name}_07_HC.vcf.gz -O ${work_path}/${name}_08_snps.recal --resource:MGP,known=false,training=true,truth=true,prior=15.0  $MGP_vcf --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp_vcf -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP --max-gaussians 4 --tranches-file ${work_path}/${name}_08_snps.tranches --rscript-file ${work_path}/${name}_08_snps.R  -L "chr1" -L "chr2" -L "chr3" -L "chr4" -L "chr5" -L "chr6" -L "chr7" -L "chr8" -L "chr9" -L "chr10" -L "chr11" -L "chr12" -L "chr13" -L "chr14" -L "chr15" -L "chr16" -L "chr17" -L "chr18" -L "chr19" -L "chrX" -L "chrY" -L "chrM" 2> >(tee ${work_path}/${name}_08_snps.log >&2)   ###### Thu Apr 23 08:56:09 CST 2020 fzc add -L ###### Tue Jul 21 14:04:29 CST 2020 deletion -L "chrM"  

        # ApplyVQSR SNP  
        ${dir_of_gatk}/gatk --java-options "-Xmx30g -Xms10g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" ApplyVQSR -R ${ref_genome_path} -V ${work_path}/${name}_07_HC.vcf.gz --recal-file ${work_path}/${name}_08_snps.recal -O ${work_path}/${name}_08_snps_VQSR.vcf --tranches-file ${work_path}/${name}_08_snps.tranches -mode SNP -ts-filter-level 95 -L "chr1" -L "chr2" -L "chr3" -L "chr4" -L "chr5" -L "chr6" -L "chr7" -L "chr8" -L "chr9" -L "chr10" -L "chr11" -L "chr12" -L "chr13" -L "chr14" -L "chr15" -L "chr16" -L "chr17" -L "chr18" -L "chr19" -L "chrX" -L "chrY" -L "chrM" ###### Tue Jul 21 14:04:29 CST 2020 deletion -L "chrM" added backed ###### Mon Aug 17 09:20:44 CST 2020  added backed -L "chrM"  
        fi  
        if [ "$mutation_type" != "SNV" ];then  
        # 8.2 INDEL  
        # VariantRecalibrator INDEL  
        ${dir_of_gatk}/gatk --java-options "-Xmx30g -Xms10g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" VariantRecalibrator -R ${ref_genome_path} -V ${work_path}/${name}_08_snps_VQSR.vcf -O ${work_path}/${name}_08_indels.recal -resource:MGP,known=true,training=true,truth=true,prior=12.0 $MGP_indel_vcf -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP --max-gaussians 4 -mode INDEL --tranches-file ${work_path}/${name}_08_indels.tranches --rscript-file ${work_path}/${name}_08_indels.R -L "chr1" -L "chr2" -L "chr3" -L "chr4" -L "chr5" -L "chr6" -L "chr7" -L "chr8" -L "chr9" -L "chr10" -L "chr11" -L "chr12" -L "chr13" -L "chr14" -L "chr15" -L "chr16" -L "chr17" -L "chr18" -L "chr19" -L "chrX" -L "chrY"  -L "chrM" 2> >(tee ${work_path}/${name}_08_indels.log >&2) ###### Tue Jul 21 14:04:29 CST 2020 deletion -L "chrM" ###### Mon Aug 17 09:20:44 CST 2020  added backed -L "chrM"  
        # ApplyVQSR INDEL  
        ${dir_of_gatk}/gatk --java-options "-Xmx30g -Xms10g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" ApplyVQSR -R ${ref_genome_path} -V ${work_path}/${name}_08_snps_VQSR.vcf --recal-file ${work_path}/${name}_08_indels.recal -O ${work_path}/${name}_08_indels_VQSR.vcf --tranches-file ${work_path}/${name}_08_indels.tranches -mode INDEL -ts-filter-level 95 -L "chr1" -L "chr2" -L "chr3" -L "chr4" -L "chr5" -L "chr6" -L "chr7" -L "chr8" -L "chr9" -L "chr10" -L "chr11" -L "chr12" -L "chr13" -L "chr14" -L "chr15" -L "chr16" -L "chr17" -L "chr18" -L "chr19" -L "chrX" -L "chrY" -L "chrM" ###### Tue Jul 21 14:04:29 CST 2020 deletion -L "chrM" ###### Mon Aug 17 09:20:44 CST 2020  added backed -L "chrM"  
        fi  
    }  
    else  
    echo "unsupport genome_build_version: $genome_build_version"  
    fi  
    # 9. SelectVariants  
    # 9.1 select SNP  
    if ([ ! -s "${work_path}/${name}_09_SNVs_VQSR.vcf.gz" ]) || ( [ "${Patch}" != "True" ]) ;then  
    ${dir_of_gatk}/gatk --java-options "-Xmx30g -Xms10g" SelectVariants -V ${work_path}/${name}_08_snps_VQSR.vcf -O ${work_path}/${name}_09_SNVs_VQSR.vcf.gz -select-type SNP  
    test -s ${work_path}/${name}_09_SNVs_VQSR.vcf.gz && merm ${work_path}/${name}_08_snps.recal ${work_path}/${name}_08_snps_VQSR.vcf  
    fi  
    if [ "$mutation_type" != "SNV" ];then  
    # 9.2 select INDEL  
    ${dir_of_gatk}/gatk --java-options "-Xmx30g -Xms10g" SelectVariants -V ${work_path}/${name}_08_indels_VQSR.vcf -O ${work_path}/${name}_09_indels_VQSR.vcf.gz -select-type INDEL  
    test -s ${work_path}/${name}_09_indels_VQSR.vcf.gz && merm ${work_path}/${name}_08_indels.recal ${work_path}/${name}_08_indels_VQSR.vcf  
    fi  

    }  
    fi  
}  
step8_SNVs_strelka2(){  
    if [ ! -s "${work_path}/${name}_SNVs_strelka2.vcf.gz" ] || [ "${Patch}" != "True" ];then  
    {  
        # Strelka2 2.9.10 [SNVs @3 ok]  
        # https://github.com/Illumina/strelka  
        # https://github.com/Illumina/strelka/blob/v2.9.x/docs/userGuide/quickStart.md  
        # /picb/rnomics3/xuew/software/WGS/strelka-2.9.10.centos6_x86_64/bin/configureStrelkaGermlineWorkflow.py  

            #20:47  
        merm ${work_path}/SNVs_strelka2
        ${dir_of_Strelka2}/configureStrelkaGermlineWorkflow.py --bam ${work_path}/${name}_06_BQSR.bam --referenceFasta ${ref_genome_path} --runDir ${work_path}/SNVs_strelka2  
        ${work_path}/SNVs_strelka2/runWorkflow.py -m local -j $threads  
        # output [./SNVs_strelka2/results/variants/variants.vcf.gz]  
        #zcat ${work_path}/SNVs_strelka2/results/variants/variants.vcf.gz > ${work_path}/${name}_SNVs_strelka2.vcf  
        mv ${work_path}/SNVs_strelka2/results/variants/variants.vcf.gz ${work_path}/${name}_SNVs_strelka2.vcf.gz  
        test -s ${work_path}/${name}_SNVs_strelka2.vcf.gz && rm -r ${work_path}/SNVs_strelka2  
        #svr3 100 21:46  
        }  
    fi  
}  
step9_SNVs_lofreq(){  
    if [ ! -s "${work_path}/${name}_SNVs_lofreq.vcf" ] || [ "${Patch}" != "True" ];then  
        # lofreq 2.1.3.1 [SNVs @2 ok]  
        # https://github.com/CSB5/lofreq  
        # https://csb5.github.io/lofreq/  
        # /picb/rnomics3/xuew/software/WGS/lofreq_star-2.1.3.1/bin/lofreq  
        # lofreq call -f ref.fa -o vars.vcf aln.bam  
        # lofreq call-parallel --pp-threads 8 -f ref.fa -o vars.vcf aln.bam  
        test -s ${work_path}/${name}_SNVs_lofreq.vcf || ${dir_of_lofreq}/lofreq call-parallel --pp-threads $threads -f ${ref_genome_path} -o ${work_path}/${name}_SNVs_lofreq.vcf ${work_path}/${name}_06_BQSR.bam  
        # output [SNVs_lofreq.vcf]  
    fi  
}  

step10_indels_Manta(){  
if [ ! -s "${work_path}/${name}_09_indels_Manta.vcf.gz" ] || [ "${Patch}" != "True" ] ;then  
{  
    # Manta 1.6.0 [indels @3.1 ok]  
    # /picb/rnomics3/wangying/tools/manta-1.6.0.centos6_x86_64/bin/configManta.py  
    # /picb/rnomics3/xuew/software/WGS/manta-1.6.0.centos6_x86_64/bin/configManta.py  
        #20-03-15 20:56:26  
    ${dir_of_Manta}/configManta.py --bam ${work_path}/${name}_06_BQSR.bam --referenceFasta ${ref_genome_path} --runDir ${work_path}/indels_Manta  
    ${work_path}/indels_Manta/runWorkflow.py -j $threads  
    #zcat ${work_path}/indels_Manta/results/variants/candidateSmallIndels.vcf.gz > ${work_path}/${name}_09_indels_Manta.vcf  
    cp ${work_path}/indels_Manta/results/variants/candidateSmallIndels.vcf.gz ${work_path}/${name}_09_indels_Manta.vcf.gz  
    #20-03-15 20:59:19 60  
}  
fi  
}  

step11_indels_strelka2(){  
    if [ ! -s "${work_path}/${name}_09_indels_strelka2.vcf.gz" ] || [ "${Patch}" != "True" ];then  
    # output [./indels_Manta/results/variants/candidateSmallIndels.vcf.gz]  
    # Strelka2 2.9.10 [indels @3.2 ok]  
    # https://github.com/Illumina/strelka  
    # https://github.com/Illumina/strelka/blob/v2.9.x/docs/userGuide/quickStart.md  
    # /picb/rnomics3/xuew/software/WGS/strelka-2.9.10.centos6_x86_64/bin/configureStrelkaGermlineWorkflow.py  
    merm ${work_path}/indels_strelka2  
    ${dir_of_Strelka2}/configureStrelkaGermlineWorkflow.py --bam ${work_path}/${name}_06_BQSR.bam --referenceFasta ${ref_genome_path} --indelCandidates ${work_path}/indels_Manta/results/variants/candidateSmallIndels.vcf.gz --runDir ${work_path}/indels_strelka2  
    
    ${work_path}/indels_strelka2/runWorkflow.py -m local -j $threads  
    #zcat ${work_path}/indels_strelka2/results/variants/variants.vcf.gz > ${work_path}/${name}_09_indels_strelka2.vcf  
    mv ${work_path}/indels_strelka2/results/variants/variants.vcf.gz ${work_path}/${name}_09_indels_strelka2.vcf.gz  
    # output [./indels_strelka2/results/variants/genome.S1.vcf.gz]  
    fi  
}  
step12_1_split_bam(){  
    if ( [ ! -s ${work_path}/${name}_scalpel_indels.vcf ] && [ ! -e "${work_path}/06_split_chr_bam/split_bam_ok" ] ) || [ "${Patch}" != "True" ];then  
    #if [ ! -s "${work_path}/${name}_scalpel_indels.vcf" ] || [ "${Patch}" != "True" ];then  
    memkdir ${work_path}/06_split_chr_bam  
    ${dir_of_bamtools}/bamtools split -in ${work_path}/${name}_06_BQSR.bam -stub ${work_path}/06_split_chr_bam/${name}_06_BQSR -reference  
    touch ${work_path}/06_split_chr_bam/split_bam_ok  
    fi  
}  
step12_scalpel_indels(){  
    if [ ! -s ${work_path}/${name}_scalpel_indels.vcf ] || [ "${Patch}" != "True" ];then  
    {  
        #wait  
        # Scalpel v0.5.4 [indels @2]  
        # http://scalpel.sourceforge.net/index.html  
        # /picb/rnomics3/xuew/software/WGS/scalpel-0.5.4/scalpel-discovery  
        # scalpel-discovery --single --bam file.bam --ref genome.fa --bed 22:1-51304566 --window 600 --numprocs 10 [OPTIONS]  
            #20-03-15 21:00:51  
        #${dir_of_Scalpel}/scalpel-discovery --single --bam ${work_path}/${name}_06_BQSR.bam --ref ${ref_genome_path} --bed $genome_range_bed --window 600 --numprocs 10 --dir ${work_path}/indels_scalpel ###### Sun Mar 15 21:47:29 CST 2020 annotated by fzc confirmed by xw  
        #  
        #60  
        # output [./indels_scalpel/variants.indel.vcf]  

        #memkdir ${work_path}/06_split_chr_bam  
        #${dir_of_bamtools}/bamtools split -in ${work_path}/${name}_06_BQSR.bam -stub ${work_path}/06_split_chr_bam/${name}_06_BQSR -reference  


        # sh ../script_shell/indel_Scalpel_01_split_bam.sh  
        #sh ../../../script_shell/indel_Scalpel_02_del_chrUn_random_alt_hg38.sh ###### Tue Feb 25 22:11:02 CST 2020 fzc;  
        #sh ../../../script_shell/indel_Scalpel_03_samtools_quickcheck.sh  ###### Tue Feb 25 22:27:41 CST 2020 fzc;  
        #sh ../../../script_shell/indel_Scalpel_04_samtools_index_hg38.sh ###### Tue Feb 25 22:27:48 CST 2020 fzc;  
        #sh ../../../script_shell/indel_Scalpel_05_discovery_split_chr_pu32_hg38.sh ###### Tue Feb 25 22:27:54 CST 2020 fzc;  


        while [ ! -e "${work_path}/06_split_chr_bam/split_bam_ok" ];do  
        sleep 300  
        done  


        for file1 in `ls ${work_path}/06_split_chr_bam/*.bam|grep -v "chr[^_]\{1,2\}\.bam"`  
        do  
        merm $file1  
        done ###### Tue Feb 25 22:11:10 CST 2020 fzc;  
        scalpel_indels_list=${work_path}/${name}_scalpel_indels_list  
        if [ -e ${scalpel_indels_list} ];then  
        ok_chrn_num=`cat ${scalpel_indels_list}|sort|uniq|wc -l`  
        else  
        ok_chrn_num=0  
        fi  

        merm $scalpel_indels_list  
        #20-03-15 22:04:25 100  
        chrn_sum=`ls -Sr ${work_path}/06_split_chr_bam/*.bam|wc -l`  
        need_do_chrn_num=`echo $chrn_sum - $ok_chrn_num|bc`  
        one_shot_chrn_num=`awk 'BEGIN{print int('${need_do_chrn_num}'/2+1)}'`  
        if [ "$need_do_chrn_num" -eq 1 ];then  
        threads_1=30 #`threads_ctrl_for_step12_scalpel_indels`  
        else  
        threads_1=3  
        fi  
        
        
        for bam_file in `ls -Sr ${work_path}/06_split_chr_bam/*.bam` ;do ###### Tue Feb 25 22:26:59 CST 2020 added by fzc;  
        #for bam_file in `ls -S ${work_path}/06_split_chr_bam/*.bam|head -n10` ;do ###### Tue Feb 25 22:26:59 CST 2020 added by fzc;  
            #06_BQSR.REF_chr10.bam  
            chr_n=$(echo $bam_file|awk -F ".REF_" '{print $2}'|cut -d. -f1)  
            output_file=${work_path}/06_split_chr_bam/scalpel_${chr_n}/variants.indel.vcf  
            test -s $output_file  &&{  
                echo $output_file >>$scalpel_indels_list  
                continue  
            }|| echo "BEGIN $output_file"  
            test -d ${work_path}/06_split_chr_bam/scalpel_${chr_n} && rm -rf ${work_path}/06_split_chr_bam/scalpel_${chr_n}  
            ${dir_of_samtools}/samtools quickcheck $bam_file || {  
                echo $bam_file error;continue  
                }  
        {  
            rm -f ${bam_file}.bai  
            ${dir_of_samtools}/samtools index $bam_file  
            ${dir_of_Scalpel}/scalpel-discovery --single --bam $bam_file --ref ${dir_of_individual_chr_ref_genome_path}/${chr_n}.fa --bed ${dir_of_individual_chr_genome_range_bed}/${chr_n}.bed --window 600 --numprocs ${threads_1} --dir ${work_path}/06_split_chr_bam/scalpel_${chr_n}  
            echo $output_file >>$scalpel_indels_list  
        }&  
        metc `threads_ctrl_for_step12_scalpel_indels`  
        done  
        wait  
        ##awk '{print "06_split_chr_bam/scalpel_"$0"/variants.indel.vcf"}' /picb/rnomics3/xuew/Human/backup/hg38_common/hg38_all.txt > scalpel_indels.list  
        java -jar ${dir_of_picard}/picard.jar MergeVcfs I=$scalpel_indels_list O=${work_path}/${name}_scalpel_indels.vcf SEQUENCE_DICTIONARY=${dict_of_ref_genome_path}  
        ## /picb/rnomics3/xuew/software/WGS/scalpel-0.5.4/scalpel-discovery --single --bam ${work_path}/${name}_06_BQSR.REF_chr22.bam --ref /picb/rnomics3/xuew/Human/backup/hg38_common/chr22.fa --bed /picb/rnomics3/xuew/Human/backup/hg38_common/chr22.bed --window 600 --numprocs 10 --dir scalpel_chr22  
        ###### Tue Jul 21 14:04:29 CST 2020 deletion -L "chrM"  
        chrn_sum_infile=`${dir_of_bcftools}/bcftools view -H ${work_path}/${name}_scalpel_indels.vcf|cut -f1|sort|uniq|wc -l`  
        chrn_sum_inlist=`cat $scalpel_indels_list|wc -l`  
        if [ "$chrn_sum_inlist" != "$chrn_sum" ];then  
        echo -e 'ERROR! ! ! ! ! #'" \n $scalpel_indels_list do not have correct chromosome number; please check it! "  
        mv ${work_path}/${name}_scalpel_indels.vcf ${work_path}/${name}_scalpel_indels.vcf_chrn_error  
        fi  

        if [ "$chrn_sum" != "$chrn_sum_infile" ] ;then  
        echo -e 'ERROR! ! ! ! ! #'" \n ${work_path}/${name}_scalpel_indels.vcf do not have correct chromosome number; please check it! \n ${dir_of_bcftools}/bcftools view -H ${work_path}/${name}_scalpel_indels.vcf|cut -f1|sort|uniq|wc -l"  
        mv ${work_path}/${name}_scalpel_indels.vcf ${work_path}/${name}_scalpel_indels.vcf_chrn_error  
        fi  
        #_scalpel_indels  
        test -s ${work_path}/${name}_scalpel_indels.vcf && rm -r ${work_path}/06_split_chr_bam  
        }  

    fi  
}  
skipping_steps(){  
    ###### Fri Feb 21 10:13:22 CST 2020 xw recommend skipping next two steps, added by fzc  
    NOT_RUN=True  
    if [ "$NOT_RUN" != "True" ];then  
    {  
        
        # /data/rnomics6/xuew/Human/BE_WangLiJie/SBE_genome_20191023/mapping/01_C9_11/  
        bcftools mpileup -d 8000 --threads 12 -f ${ref_genome_path} ${work_path}/${name}_06_BQSR.bam > ${work_path}/${name}_07_bcftools.vcf  
        bcftools call ${work_path}/${name}_07_bcftools.vcf -o ${work_path}/${name}_08_bcftools_indel.vcf -O v -V snps -v -c  
    }  
    fi  
    NOT_RUN=True  
    if [ "$NOT_RUN" != "True" ];then  
    {  
        # annovar Version 2018-04-16  
        # http://annovar.openbioinformatics.org/en/latest/user-guide/startup/  
        perl /picb/rnomics3/xuew/software/SNV_annovar/annotate_variation.pl  

        perl /picb/rnomics3/xuew/software/SNV_annovar/convert2annovar.pl -format vcf4 10_SNVs_gatk_d10_Q20.vcf > 11_SNVs_gatk_d10_Q20.avinput  
        perl /picb/rnomics3/xuew/software/SNV_annovar/annotate_variation.pl -filter -buildver hg19 -dbtype 1000g2015aug_all -maf 0.001 11_SNVs_gatk_d10_Q20.avinput /picb/rnomics3/xuew/software/SNV_annovar/hg19/ --outfile 12_SNVs_gatk_d10_Q20_snps  
        perl /picb/rnomics3/xuew/software/SNV_annovar/table_annovar.pl 10_SNVs_gatk_d10_Q20.vcf /picb/rnomics3/xuew/software/SNV_annovar/hg19/ -buildver hg19 -out 13_SNVs_gatk_d10_Q20 -remove -protocol refGene,genomicSuperDups,cytoBand,exac03,avsnp150,dbnsfp30a,clinvar_20180603,esp6500siv2_all -operation g,r,r,f,f,f,f,f -nastring . -vcfinput  
        perl /picb/rnomics3/xuew/software/SNV_annovar/convert2annovar.pl -format vcf4 10_indels_gatk_d10_Q20.vcf > 11_indels_gatk_d10_Q20.avinput  
        perl /picb/rnomics3/xuew/software/SNV_annovar/annotate_variation.pl -filter -buildver hg19 -dbtype 1000g2015aug_all -maf 0.001 11_indels_gatk_d10_Q20.avinput /picb/rnomics3/xuew/software/SNV_annovar/hg19/ --outfile 12_indels_gatk_d10_Q20_snps  
        perl /picb/rnomics3/xuew/software/SNV_annovar/table_annovar.pl 10_indels_gatk_d10_Q20.vcf /picb/rnomics3/xuew/software/SNV_annovar/hg19/ -buildver hg19 -out 13_indels_gatk_d10_Q20 -remove -protocol refGene,genomicSuperDups,cytoBand,exac03,avsnp150,dbnsfp30a,clinvar_20180603,esp6500siv2_all -operation g,r,r,f,f,f,f,f -nastring . -vcfinput  
        # vcf  
        # /data/rnomics6/xuew/Human/BE_WangLiJie/SBE_genome_20191023/backup/vcf151/  
        ln -s /picb/rnomics4/wangying/DNA_Cas9_ChenJia/Cpf1_paper/ClinVar_Database/Download/dbSNP/GRCh37_b151_p13/common_all_20180423_chr.vcf ${work_path}/${name}_01_GRCh37_b151_p13_common.vcf  
        ln -s /picb/rnomics4/wangying/DNA_Cas9_ChenJia/Cpf1_paper/ClinVar_Database/Download/dbSNP/GRCh37_b151_p13/All/All_20180423.vcf ${work_path}/${name}_01_GRCh37_b151_p13_all.vcf  
        sed '1,58d' ${work_path}/${name}_01_GRCh37_b151_p13_common.vcf |awk '{print $1":"$2"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5}' > ${work_path}/${name}_01_GRCh37_b151_p13_common.txt  
        sed '1,58d' ${work_path}/${name}_01_GRCh37_b151_p13_all.vcf |awk '{print $1":"$2"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5}' > ${work_path}/${name}_01_GRCh37_b151_p13_all.txt  
    }  
    fi  
}  
memkdir $work_path  
source $tmp_config_file  
echo "[`date`] Main flow begining..."  
echo "Patch: " ${Patch}  
echo "Patch_For_Scapel: " ${Patch_For_Scapel}  
if [ "$Patch_For_Scapel" == True ];then
if [[ $(uname -n) == "liyang-svr9" ]] || [[ $(uname -n) == "liyang-svr8" ]] ||[[ $(uname -n) == "liyang-svr5.icb.ac.cn" ]]||([[ $(uname -n) == "liyang-svr6.icb.ac.cn" ]] && [[ $(id -u) == "4608" ]]) ;then  
if [ -e ${work_path}/06_split_chr_bam/split_bam_ok  ];then  
step12_scalpel_indels 
wait  
exit  
else  
echo ${work_path}/06_split_chr_bam/split_bam_ok " not exists" 
exit  
fi  
fi  
fi

step01_06_construct_bam  

if [ "$mutation_type" != "SNV" ];then  
step12_1_split_bam &  

fi  

step7_HaplotypeCaller &  


if [ "$mutation_type" != "Indel" ];then  
step8_SNVs_strelka2  
step9_SNVs_lofreq  
fi  
if [ "$mutation_type" != "SNV" ];then  
step10_indels_Manta  
step11_indels_strelka2  
if [ $(id -u) != "5158" ]||([ $(uname -n) != "liyang-svr6.icb.ac.cn" ]&&[ $(uname -n) != "liyang-svr3.icb.ac.cn" ]);then  

step12_scalpel_indels  
fi  
fi  

wait  
test -s ${work_path}/${name}_scalpel_indels.vcf  && touch ${work_path}/${name}_Main_stream_ok || rm -f ${work_path}/${name}_Main_stream_ok  
###### Fri Feb 21 10:02:08 CST 2020 xw recommended end, added by fzc  


