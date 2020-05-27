MYDIR=`dirname $0`
# File Name: xw_work_PGM3_02_SNVs_d10_af10.sh
# Author: XueWei
# mail: xuewei@picb.ac.cn
# Created time: 
# Last modified: Tue May 26 22:24:02 CST 2020
while getopts :o:n:c: ARGS  
    do  
    case $ARGS in   
        o)  
            work_path=$OPTARG
            ;;  
        n)  
            name=$OPTARG
            ;; 
        c)
            tmp_config_file=$OPTARG
            ;;
        *)  
            echo "Unknown option: $ARGS"
            ;;
        \?)
        echo "Invalid option: -$OPTARG" 
        ;;
    esac
    done
memkdir(){
    for path in $@
    do
    test -d $path ||mkdir -p $path
    done
}
merm(){
    for file1 in $@
    do
    test -e $file1 && rm $file1 || echo "$file1" not exist!
    done
}
check_status(){
    ###### Sun Mar 1 21:06:32 CST 2020 fzc
    bin_file=$1
    status=0
    test -e $bin_file && { 
        let status+=1
        test -x $bin_file && let status+=1
    }
    test "$status" == "0" && {
        echo "$bin_file is non-exist"
        check_token="F"
    }
    test "$status" == "1" && {
        echo "$bin_file is not executable!"
        check_token="F" 
    }
		a=a
}
parse_parameters(){
    ###### Sun Mar 1 21:06:21 CST 2020 fzc
    tmp_a="aa${tmp_config_file}aa"
    if [ "$tmp_a" != "aaaa" ];then
    config_file=$tmp_config_file
    test -e $config_file || { 
        echo "The config file($config_file) is non-exist!" 
        exit 1
    } 
    else
    config_file=$(dirname $0)/${this_software_name}_config
    test -e $config_file || { 
        echo "Config_file not specified in CMD, and the default config file not find in $(dirname $0)" 
        exit 1
    }
    fi
    test -e $config_file && echo "use $config_file as config_file" 
    source $config_file
    para_list=(dir_of_gatk dir_of_perl dir_of_intersectBed)
    non_specified_para_list=""
    non_exists_para_list=""
    for para1 in ${para_list[@]}
    do
    {
        para1_value=`eval echo "$"$para1`
        tmp_a="aa${para1_value}aa"
        #echo $para1 "a${para1_value}a"
        test "${tmp_a}" == "aaaa" && {
            non_specified_para_list=$non_specified_para_list"\n"${para1}
            continue
        }
        test -e ${para1_value} || {
            non_exists_para_list=$non_exists_para_list"\n"${para1_value}
            continue
        }
        if [[ $(echo $para1|awk -F "_" '{print $1}' ) == "dir" ]];then
        {
            software_name=$(echo $para1|awk -F "_" '{print $NF}')
            check_status ${para1_value}/$software_name
        }
        fi
    }
    done
    test "$non_specified_para_list" != "" &&{
        echo "The following para(s) are not specified but necessary:$non_specified_para_list"
        check_token="F"
    }   
    test "$non_exists_para_list" != "" &&{
        echo -e "The following specified para are not exist:$non_exists_para_list"
        check_token="F"
    }   
    if [ "$check_token" == "F" ];then exit 1; fi
echo "[`date`]Finish parse_parameters"
}
join_ID_pl=${MYDIR}/join_ID.pl
select_v_ID_pl=${MYDIR}/select_v_ID.pl

# 10. vcfutils.pl varFilter
ln -sf ${work_path}/${name}_09_SNVs_VQSR.vcf ${work_path}/${name}_10_SNVs_gatk_total.vcf
ln -sf ${work_path}/${name}_SNVs_lofreq.vcf ${work_path}/${name}_10_SNVs_lofreq_total.vcf
ln -sf ${work_path}/${name}_SNVs_strelka2.vcf ${work_path}/${name}_10_SNVs_strelka2_total.vcf

memkdir ${work_path}/Novel_SNVs_d10_af10 #Novel_SNVs_d10_af10
# 11. PASS
# 11.1 gatk
${dir_of_gatk}/gatk --java-options "-Xmx30g -Xms10g" SelectVariants -V ${work_path}/${name}_10_SNVs_gatk_total.vcf -O ${work_path}/${name}_11_SNVs_gatk_PASS.vcf --exclude-filtered true -select-type SNP
awk '{if($7=="PASS"&&$5!~",")print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$NF}' ${work_path}/${name}_11_SNVs_gatk_PASS.vcf |awk '{split($6,x,":"); print $0"\t"x[2]"\t"x[3]}' |cut -f 1-2,4-5,7-8 |sed 's/,/\t/g' |awk '{if($7>0)print$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$6/$7}' > ${work_path}/${name}_11_SNVs_gatk_PASS.txt
# 11.2 lofreq
${dir_of_gatk}/gatk --java-options "-Xmx30g -Xms10g" SelectVariants -V ${work_path}/${name}_10_SNVs_lofreq_total.vcf -O ${work_path}/${name}_11_SNVs_lofreq_PASS.vcf --exclude-filtered true -select-type SNP
awk '{if($7=="PASS")print $0}' ${work_path}/${name}_11_SNVs_lofreq_PASS.vcf |awk '{split($8,x,";"); print $1"\t"$2"\t"$4"\t"$5"\t"x[1]"\t"x[2]"\t"x[3]}' |sed 's/AF=//g' |sed 's/DP=//g' |sed 's/DP4=//g' |awk '{split($7,x,","); print $1"\t"$2"\t"$3"\t"$4"\t"x[1]+x[2]"\t"x[3]+x[4]"\t"$6"\t"$5}' > ${work_path}/${name}_11_SNVs_lofreq_PASS.txt
# 11.3 Strelka2
${dir_of_gatk}/gatk --java-options "-Xmx30g -Xms10g" SelectVariants -V ${work_path}/${name}_10_SNVs_strelka2_total.vcf -O ${work_path}/${name}_11_SNVs_strelka2_PASS.vcf --exclude-filtered true -select-type SNP
awk '{if($7=="PASS"&&$5!~",")print $1"\t"$2"\t"$4"\t"$5"\t"$NF}' ${work_path}/${name}_11_SNVs_strelka2_PASS.vcf |awk '{split($5,x,":"); print $1"\t"$2"\t"$3"\t"$4"\t"x[2]"\t"x[5]}' |sed 's/,/\t/g' |awk '{print $0"\t"$6/$7}' > ${work_path}/${name}_11_SNVs_strelka2_PASS.txt

# 12. depth ≥ 10/20 Allele Frequency ≥ 0.1
# 12.1 gatk d10_af10/d10_af10
awk '{if($7>=10 && $8>=0.1)print $1"\t"$2-1"\t"$2"\t"$1":"$2"\t"$0}' ${work_path}/${name}_11_SNVs_gatk_PASS.txt > ${work_path}/Novel_SNVs_d10_af10/${name}_12_SNVs_gatk_d10_af10.txt
# 12.2 lofreq d10_af10/d10_af10
awk '{if($7>=10 && $8>=0.1)print $1"\t"$2-1"\t"$2"\t"$1":"$2"\t"$0}' ${work_path}/${name}_11_SNVs_lofreq_PASS.txt > ${work_path}/Novel_SNVs_d10_af10/${name}_12_SNVs_lofreq_d10_af10.txt
# 12.3 Strelka2 d10_af10/d10_af10
awk '{if($7>=10 && $8>=0.1)print $1"\t"$2-1"\t"$2"\t"$1":"$2"\t"$0}' ${work_path}/${name}_11_SNVs_strelka2_PASS.txt > ${work_path}/Novel_SNVs_d10_af10/${name}_12_SNVs_strelka2_d10_af10.txt

# 13. remove dbSNP
# dbSNP151 all
{
# gatk_d10_af10
${dir_of_perl}/perl ${join_ID_pl} ${db_path}/backup/vcf151/01_GRCh38_b151_p7_all.txt ${work_path}/Novel_SNVs_d10_af10/${name}_12_SNVs_gatk_d10_af10.txt 1 4 |cut -f 1,4,7-9,13-18 > ${work_path}/Novel_SNVs_d10_af10/${name}_13_all_gatk_d10_af10.txt
cut -f 1 ${work_path}/Novel_SNVs_d10_af10/${name}_13_all_gatk_d10_af10.txt |sort -u > ${work_path}/Novel_SNVs_d10_af10/${name}_13_all_gatk_d10_af10.sites
${dir_of_perl}/perl ${select_v_ID_pl} ${work_path}/Novel_SNVs_d10_af10/${name}_12_SNVs_gatk_d10_af10.txt ${work_path}/Novel_SNVs_d10_af10/${name}_13_all_gatk_d10_af10.sites 4 > ${work_path}/Novel_SNVs_d10_af10/${name}_13_all_gatk_d10_af10_rmdbSNP151.txt
}&
{
# lofreq_d10_af10
${dir_of_perl}/perl ${join_ID_pl} ${db_path}/backup/vcf151/01_GRCh38_b151_p7_all.txt ${work_path}/Novel_SNVs_d10_af10/${name}_12_SNVs_lofreq_d10_af10.txt 1 4 |cut -f 1,4,7-9,13-18 > ${work_path}/Novel_SNVs_d10_af10/${name}_13_all_lofreq_d10_af10.txt
cut -f 1 ${work_path}/Novel_SNVs_d10_af10/${name}_13_all_lofreq_d10_af10.txt |sort -u > ${work_path}/Novel_SNVs_d10_af10/${name}_13_all_lofreq_d10_af10.sites
${dir_of_perl}/perl ${select_v_ID_pl} ${work_path}/Novel_SNVs_d10_af10/${name}_12_SNVs_lofreq_d10_af10.txt ${work_path}/Novel_SNVs_d10_af10/${name}_13_all_lofreq_d10_af10.sites 4 > ${work_path}/Novel_SNVs_d10_af10/${name}_13_all_lofreq_d10_af10_rmdbSNP151.txt
}&
{
# strelka2_d10_af10
${dir_of_perl}/perl ${join_ID_pl} ${db_path}/backup/vcf151/01_GRCh38_b151_p7_all.txt ${work_path}/Novel_SNVs_d10_af10/${name}_12_SNVs_strelka2_d10_af10.txt 1 4 |cut -f 1,4,7-9,13-18 > ${work_path}/Novel_SNVs_d10_af10/${name}_13_all_strelka2_d10_af10.txt
cut -f 1 ${work_path}/Novel_SNVs_d10_af10/${name}_13_all_strelka2_d10_af10.txt |sort -u > ${work_path}/Novel_SNVs_d10_af10/${name}_13_all_strelka2_d10_af10.sites
${dir_of_perl}/perl ${select_v_ID_pl} ${work_path}/Novel_SNVs_d10_af10/${name}_12_SNVs_strelka2_d10_af10.txt ${work_path}/Novel_SNVs_d10_af10/${name}_13_all_strelka2_d10_af10.sites 4 > ${work_path}/Novel_SNVs_d10_af10/${name}_13_all_strelka2_d10_af10_rmdbSNP151.txt
}&
wait
#cat ${work_path}/Novel_SNVs_d10_af10/${name}_13*_d10_af10_rmdbSNP151.txt |cut -f 4 |sort |uniq -c |awk '{print $2"\t"$1}' |awk '{if($2==3)print $1}' |wc -l

# 14. remove repeat region
# gatk_d10_af10
${dir_of_intersectBed}/intersectBed -a ${work_path}/Novel_SNVs_d10_af10/${name}_13_all_gatk_d10_af10_rmdbSNP151.txt -b ${db_path}/backup/vcf151/UCSC_RepeatMask_hg38.bed -wa -wb |cut -f 4 |sort -u > ${work_path}/Novel_SNVs_d10_af10/${name}_14_all_gatk_d10_af10_RepeatMask.sites
${dir_of_perl}/perl ${select_v_ID_pl} ${work_path}/Novel_SNVs_d10_af10/${name}_13_all_gatk_d10_af10_rmdbSNP151.txt ${work_path}/Novel_SNVs_d10_af10/${name}_14_all_gatk_d10_af10_RepeatMask.sites 4 |cut -f 1-4,7-12 > ${work_path}/Novel_SNVs_d10_af10/${name}_14_all_gatk_d10_af10_novel.txt
cut -f 4 ${work_path}/Novel_SNVs_d10_af10/${name}_14_all_gatk_d10_af10_novel.txt |sort -u > ${work_path}/Novel_SNVs_d10_af10/${name}_14_all_gatk_d10_af10_novel.sites
# lofreq_d10_af10
${dir_of_intersectBed}/intersectBed -a ${work_path}/Novel_SNVs_d10_af10/${name}_13_all_lofreq_d10_af10_rmdbSNP151.txt -b ${db_path}/backup/vcf151/UCSC_RepeatMask_hg38.bed -wa -wb |cut -f 4 |sort -u > ${work_path}/Novel_SNVs_d10_af10/${name}_14_all_lofreq_d10_af10_RepeatMask.sites
${dir_of_perl}/perl ${select_v_ID_pl} ${work_path}/Novel_SNVs_d10_af10/${name}_13_all_lofreq_d10_af10_rmdbSNP151.txt ${work_path}/Novel_SNVs_d10_af10/${name}_14_all_lofreq_d10_af10_RepeatMask.sites 4 |cut -f 1-4,7-12 > ${work_path}/Novel_SNVs_d10_af10/${name}_14_all_lofreq_d10_af10_novel.txt
cut -f 4 ${work_path}/Novel_SNVs_d10_af10/${name}_14_all_lofreq_d10_af10_novel.txt |sort -u > ${work_path}/Novel_SNVs_d10_af10/${name}_14_all_lofreq_d10_af10_novel.sites
# strelka2_d10_af10
${dir_of_intersectBed}/intersectBed -a ${work_path}/Novel_SNVs_d10_af10/${name}_13_all_strelka2_d10_af10_rmdbSNP151.txt -b ${db_path}/backup/vcf151/UCSC_RepeatMask_hg38.bed -wa -wb |cut -f 4 |sort -u > ${work_path}/Novel_SNVs_d10_af10/${name}_14_all_strelka2_d10_af10_RepeatMask.sites
${dir_of_perl}/perl ${select_v_ID_pl} ${work_path}/Novel_SNVs_d10_af10/${name}_13_all_strelka2_d10_af10_rmdbSNP151.txt ${work_path}/Novel_SNVs_d10_af10/${name}_14_all_strelka2_d10_af10_RepeatMask.sites 4 |cut -f 1-4,7-12 > ${work_path}/Novel_SNVs_d10_af10/${name}_14_all_strelka2_d10_af10_novel.txt
cut -f 4 ${work_path}/Novel_SNVs_d10_af10/${name}_14_all_strelka2_d10_af10_novel.txt |sort -u > ${work_path}/Novel_SNVs_d10_af10/${name}_14_all_strelka2_d10_af10_novel.sites

cat ${work_path}/Novel_SNVs_d10_af10/${name}_14_all_gatk_d10_af10_novel.sites ${work_path}/Novel_SNVs_d10_af10/${name}_14_all_lofreq_d10_af10_novel.sites ${work_path}/Novel_SNVs_d10_af10/${name}_14_all_strelka2_d10_af10_novel.sites |sort |uniq -c |awk -F" " '$1==3{print $2}' > ${work_path}/Novel_SNVs_d10_af10/${name}_14_overlap_d10_af10_novel.sites

# 15. variants
# C-to-T
awk '{if(($5~"[Cc]")&&($6~"[Tc]")) print $0}' ${work_path}/Novel_SNVs_d10_af10/${name}_14_all_gatk_d10_af10_novel.txt > ${work_path}/Novel_SNVs_d10_af10/${name}_15_all_gatk_d10_af10_C2T.txt
awk '{if(($5~"[Cc]")&&($6~"[Tc]")) print $0}' ${work_path}/Novel_SNVs_d10_af10/${name}_14_all_lofreq_d10_af10_novel.txt > ${work_path}/Novel_SNVs_d10_af10/${name}_15_all_lofreq_d10_af10_C2T.txt
awk '{if(($5~"[Cc]")&&($6~"[Tc]")) print $0}' ${work_path}/Novel_SNVs_d10_af10/${name}_14_all_strelka2_d10_af10_novel.txt > ${work_path}/Novel_SNVs_d10_af10/${name}_15_all_strelka2_d10_af10_C2T.txt
cat <(cut -f4 ${work_path}/Novel_SNVs_d10_af10/${name}_15_all_gatk_d10_af10_C2T.txt|sort -u) <(cut -f4 ${work_path}/Novel_SNVs_d10_af10/${name}_15_all_lofreq_d10_af10_C2T.txt|sort -u) <(cut -f4  ${work_path}/Novel_SNVs_d10_af10/${name}_15_all_strelka2_d10_af10_C2T.txt|sort -u)|sort |uniq -c  |awk -F" " '$1==3{print $2}' > ${work_path}/Novel_SNVs_d10_af10/${name}_15_overlap_d10_af10_C2T.sites
# G-to-A
awk '{if(($5~"[Gg]")&&($6~"[Aa]")) print $0}' ${work_path}/Novel_SNVs_d10_af10/${name}_14_all_gatk_d10_af10_novel.txt > ${work_path}/Novel_SNVs_d10_af10/${name}_15_all_gatk_d10_af10_G2A.txt
awk '{if(($5~"[Gg]")&&($6~"[Aa]")) print $0}' ${work_path}/Novel_SNVs_d10_af10/${name}_14_all_lofreq_d10_af10_novel.txt > ${work_path}/Novel_SNVs_d10_af10/${name}_15_all_lofreq_d10_af10_G2A.txt
awk '{if(($5~"[Gg]")&&($6~"[Aa]")) print $0}' ${work_path}/Novel_SNVs_d10_af10/${name}_14_all_strelka2_d10_af10_novel.txt > ${work_path}/Novel_SNVs_d10_af10/${name}_15_all_strelka2_d10_af10_G2A.txt
cat <(cut -f4  ${work_path}/Novel_SNVs_d10_af10/${name}_15_all_gatk_d10_af10_G2A.txt|sort -u) <(cut -f4  ${work_path}/Novel_SNVs_d10_af10/${name}_15_all_lofreq_d10_af10_G2A.txt|sort -u) <(cut -f4  ${work_path}/Novel_SNVs_d10_af10/${name}_15_all_strelka2_d10_af10_G2A.txt|sort -u) |sort |uniq -c |awk -F" " '$1==3{print $2}' > ${work_path}/Novel_SNVs_d10_af10/${name}_15_overlap_d10_af10_G2A.sites


#mv ${work_path}/${name}_11* ${work_path}/Novel_SNVs_d10_af10/${name}_12* ${work_path}/Novel_SNVs_d10_af10/${name}_13* ${work_path}/Novel_SNVs_d10_af10/${name}_14* ${work_path}/Novel_SNVs_d10_af10/${name}_15* Novel_SNVs_d10_af10
NOT_RUN=True
if [ "$NOT_RUN" != "True" ];then
{
# 16. PGM3
# PGM3 DNA region (hg38) chr6:83188483-83188924
# PGM3 sgRNA region (hg38) chr6:83188738-83188760
# Ctrl, Cas9, BE3, Y130F and iBE
awk '{if($1=="chr6" && $2>=83188483 && $3<=83188924)print $0}' 12_SNVs_gatk_d10_af10.txt
awk '{if($1=="chr6" && $2>=83188483 && $3<=83188924)print $0}' 12_SNVs_lofreq_d10_af10.txt
awk '{if($1=="chr6" && $2>=83188483 && $3<=83188924)print $0}' 12_SNVs_strelka2_d10_af10.txt
# Ctrl
samtools view 06_BQSR.bam chr6:83188770-83188771 |wc -l
samtools view 06_BQSR.bam chr6:83188769-83188770 |wc -l
samtools view 06_BQSR.bam chr6:83188754-83188755 |wc -l
samtools view 06_BQSR.bam chr6:83188751-83188752 |wc -l
samtools view 06_BQSR.bam chr6:83188749-83188750 |wc -l
samtools view 06_BQSR.bam chr6:83188748-83188749 |wc -l
samtools view 06_BQSR.bam chr6:83188745-83188746 |wc -l

}
fi
Patch_1(){
    ###### Thu Apr 2 10:01:21 CST 2020
    file1="${work_path}/Novel_SNVs_d10_af10/${name}_15_overlap_d10_af10_G2A.sites"
    file2="${work_path}/Novel_SNVs_d10_af10/${name}_15_overlap_d10_af10_C2T.sites"
    output_file="${work_path}/Novel_SNVs_d10_af10/${name}_15_overlap_d10_af10_G2A-C2T.sites"
    #check_file="${work_path}/Novel_SNVs_d10_af10/${name}_15_overlap_d10_af10_G2A.sites"
    #check_file="${work_path}/${name}$file_suffix"
    cat $file1 $file2 >$output_file
}
Patch_1
