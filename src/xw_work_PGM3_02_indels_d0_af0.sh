MYDIR=`dirname $0`
# File Name: xw_work_PGM3_02_indels_d0_af0.sh
# Author: XueWei
# mail: xuewei@picb.ac.cn
# Created time: 
# Last modified: Tue May 26 22:24:09 CST 2020
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
source $tmp_config_file
echo "[`date`] Indel branch begining..."
join_ID_pl=${MYDIR}/join_ID.pl
select_v_ID_pl=${MYDIR}/select_v_ID.pl
memkdir ${work_path}/Novel_Indels_d0_af0
# 10. vcfutils.pl varFilter
ln -sf ${work_path}/${name}_09_indels_VQSR.vcf.gz ${work_path}/${name}_10_indels_gatk_total.vcf.gz
test -s ${work_path}/${name}_10_indels_scalpel_total.vcf ||ln -sf ${work_path}/${name}_scalpel_indels.vcf ${work_path}/${name}_10_indels_scalpel_total.vcf
ln -sf ${work_path}/${name}_09_indels_strelka2.vcf.gz  ${work_path}/${name}_10_indels_strelka2_total.vcf.gz

#%s/_scalpel_indels.vcf/_10_indels_scalpel_total.vcf/g
tabix -fp vcf ${work_path}/${name}_10_indels_gatk_total.vcf.gz
tabix -fp vcf ${work_path}/${name}_10_indels_scalpel_total.vcf 
tabix -fp vcf ${work_path}/${name}_10_indels_strelka2_total.vcf.gz 

# 11. PASS
# 11.1 gatk
${dir_of_gatk}/gatk --java-options "-Xmx30g -Xms10g" SelectVariants -V ${work_path}/${name}_10_indels_gatk_total.vcf.gz -O ${work_path}/Novel_Indels_d0_af0/${name}_11_indels_gatk_PASS.vcf.gz --exclude-filtered true -select-type INDEL
awk '{if($7=="PASS"&&$5!~",")print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$NF}' <(bcftools view -H ${work_path}/Novel_Indels_d0_af0/${name}_11_indels_gatk_PASS.vcf.gz) |awk '{split($6,x,":"); print $0"\t"x[2]"\t"x[3]}' |cut -f 1-2,4-5,7-8 |sed 's/,/\t/g' |awk '{if($7>0)print$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$6/$7}' > ${work_path}/Novel_Indels_d0_af0/${name}_11_indels_gatk_PASS.txt
# 11.2 scalpel
${dir_of_gatk}/gatk --java-options "-Xmx30g -Xms10g" SelectVariants -V ${work_path}/${name}_10_indels_scalpel_total.vcf -O ${work_path}/Novel_Indels_d0_af0/${name}_11_indels_scalpel_PASS.vcf --exclude-filtered true -select-type INDEL
awk '{if($7=="PASS"&&$5!~",")print $1"\t"$2"\t"$4"\t"$5"\t"$NF}' ${work_path}/Novel_Indels_d0_af0/${name}_11_indels_scalpel_PASS.vcf |awk '{split($5,x,":"); print $1"\t"$2"\t"$3"\t"$4"\t"x[2]}' |sed 's/,/\t/g' |awk '{print $0"\t"$5+$6"\t"$6/($5+$6)}' > ${work_path}/Novel_Indels_d0_af0/${name}_11_indels_scalpel_PASS.txt
# 11.3 Strelka2
${dir_of_gatk}/gatk --java-options "-Xmx30g -Xms10g" SelectVariants -V ${work_path}/${name}_10_indels_strelka2_total.vcf.gz -O ${work_path}/Novel_Indels_d0_af0/${name}_11_indels_strelka2_PASS.vcf.gz --exclude-filtered true -select-type INDEL
awk '{if($7=="PASS"&&$5!~",")print $1"\t"$2"\t"$4"\t"$5"\t"$NF}' <(bcftools view -H ${work_path}/Novel_Indels_d0_af0/${name}_11_indels_strelka2_PASS.vcf.gz) |awk '{split($5,x,":"); print $1"\t"$2"\t"$3"\t"$4"\t"x[2]}' |sed 's/,/\t/g' |awk '{print $0"\t"$5+$6"\t"$6/($5+$6)}' > ${work_path}/Novel_Indels_d0_af0/${name}_11_indels_strelka2_PASS.txt

# 12. depth ≥ 10/20 Allele Frequency ≥ 0.1
# 12.1 gatk d0_af0/d0_af0
#%s/d10_af10/d0_af0/g
awk '{if($7>=0 && $8>=0 )print $1"\t"$2-1"\t"$2"\t"$1":"$2"\t"$0}' ${work_path}/Novel_Indels_d0_af0/${name}_11_indels_gatk_PASS.txt > ${work_path}/Novel_Indels_d0_af0/${name}_12_indels_gatk_d0_af0.txt
# 12.2 scalpel d0_af0/d0_af0
awk '{if($7>=0 && $8>=0 )print $1"\t"$2-1"\t"$2"\t"$1":"$2"\t"$0}' ${work_path}/Novel_Indels_d0_af0/${name}_11_indels_scalpel_PASS.txt > ${work_path}/Novel_Indels_d0_af0/${name}_12_indels_scalpel_d0_af0.txt
# 12.3 Strelka2 d0_af0/d0_af0
awk '{if($7>=0 && $8>=0 )print $1"\t"$2-1"\t"$2"\t"$1":"$2"\t"$0}' ${work_path}/Novel_Indels_d0_af0/${name}_11_indels_strelka2_PASS.txt > ${work_path}/Novel_Indels_d0_af0/${name}_12_indels_strelka2_d0_af0.txt
# 13. remove dbSNP
# dbSNP151 all
{
# gatk_d0_af0
test -s ${work_path}/Novel_Indels_d0_af0/${name}_12_indels_gatk_d0_af0.txt && {
    ${dir_of_perl}/perl ${join_ID_pl} <(bcftools view $filtering_dbSNP_vcf -H|awk 'BEGIN{OFS="\t"}{print $1":"$2,$1,$2,$3,$4}') ${work_path}/Novel_Indels_d0_af0/${name}_12_indels_gatk_d0_af0.txt 1 4 |cut -f 1,4,7-9,13-18 > ${work_path}/Novel_Indels_d0_af0/${name}_13_all_gatk_d0_af0.txt
    cut -f 1 ${work_path}/Novel_Indels_d0_af0/${name}_13_all_gatk_d0_af0.txt |sort -u > ${work_path}/Novel_Indels_d0_af0/${name}_13_all_gatk_d0_af0.sites
    ${dir_of_perl}/perl ${select_v_ID_pl} ${work_path}/Novel_Indels_d0_af0/${name}_12_indels_gatk_d0_af0.txt ${work_path}/Novel_Indels_d0_af0/${name}_13_all_gatk_d0_af0.sites 4 > ${work_path}/Novel_Indels_d0_af0/${name}_13_all_gatk_d0_af0_rmdbSNP151.txt
}
}&
{
# scalpel_d0_af0
test -s ${work_path}/Novel_Indels_d0_af0/${name}_12_indels_scalpel_d0_af0.txt && {
${dir_of_perl}/perl ${join_ID_pl} <(bcftools view $filtering_dbSNP_vcf -H|awk 'BEGIN{OFS="\t"}{print $1":"$2,$1,$2,$3,$4}') ${work_path}/Novel_Indels_d0_af0/${name}_12_indels_scalpel_d0_af0.txt 1 4 |cut -f 1,4,7-9,13-18 > ${work_path}/Novel_Indels_d0_af0/${name}_13_all_scalpel_d0_af0.txt
cut -f 1 ${work_path}/Novel_Indels_d0_af0/${name}_13_all_scalpel_d0_af0.txt |sort -u > ${work_path}/Novel_Indels_d0_af0/${name}_13_all_scalpel_d0_af0.sites
${dir_of_perl}/perl ${select_v_ID_pl} ${work_path}/Novel_Indels_d0_af0/${name}_12_indels_scalpel_d0_af0.txt ${work_path}/Novel_Indels_d0_af0/${name}_13_all_scalpel_d0_af0.sites 4 > ${work_path}/Novel_Indels_d0_af0/${name}_13_all_scalpel_d0_af0_rmdbSNP151.txt
}
}&
{
# strelka2_d0_af0
test -s ${work_path}/Novel_Indels_d0_af0/${name}_12_indels_strelka2_d0_af0.txt && {
${dir_of_perl}/perl ${join_ID_pl} <(bcftools view $filtering_dbSNP_vcf -H|awk 'BEGIN{OFS="\t"}{print $1":"$2,$1,$2,$3,$4}') ${work_path}/Novel_Indels_d0_af0/${name}_12_indels_strelka2_d0_af0.txt 1 4 |cut -f 1,4,7-9,13-18 > ${work_path}/Novel_Indels_d0_af0/${name}_13_all_strelka2_d0_af0.txt
cut -f 1 ${work_path}/Novel_Indels_d0_af0/${name}_13_all_strelka2_d0_af0.txt |sort -u > ${work_path}/Novel_Indels_d0_af0/${name}_13_all_strelka2_d0_af0.sites
${dir_of_perl}/perl ${select_v_ID_pl} ${work_path}/Novel_Indels_d0_af0/${name}_12_indels_strelka2_d0_af0.txt ${work_path}/Novel_Indels_d0_af0/${name}_13_all_strelka2_d0_af0.sites 4 > ${work_path}/Novel_Indels_d0_af0/${name}_13_all_strelka2_d0_af0_rmdbSNP151.txt
}
}&
wait
#cat 13*_d0_af0_rmdbSNP151.txt |cut -f 4 |sort |uniq -c |awk '{print $2"\t"$1}' |awk '{if($2==3)print $1}' |wc -l

# 14. remove repeat region
# gatk_d0_af0
${dir_of_intersectBed}/intersectBed -a ${work_path}/Novel_Indels_d0_af0/${name}_13_all_gatk_d0_af0_rmdbSNP151.txt -b ${UCSC_RepeatMask_bed} -wa -wb |cut -f 4 |sort -u > ${work_path}/Novel_Indels_d0_af0/${name}_14_all_gatk_d0_af0_RepeatMask.sites
${dir_of_perl}/perl ${select_v_ID_pl} ${work_path}/Novel_Indels_d0_af0/${name}_13_all_gatk_d0_af0_rmdbSNP151.txt ${work_path}/Novel_Indels_d0_af0/${name}_14_all_gatk_d0_af0_RepeatMask.sites 4 |cut -f 1-4,7-12 > ${work_path}/Novel_Indels_d0_af0/${name}_14_all_gatk_d0_af0_novel.txt
cut -f 4 ${work_path}/Novel_Indels_d0_af0/${name}_14_all_gatk_d0_af0_novel.txt |sort -u > ${work_path}/Novel_Indels_d0_af0/${name}_14_all_gatk_d0_af0_novel.sites
# scalpel_d0_af0
${dir_of_intersectBed}/intersectBed -a ${work_path}/Novel_Indels_d0_af0/${name}_13_all_scalpel_d0_af0_rmdbSNP151.txt -b ${UCSC_RepeatMask_bed} -wa -wb |cut -f 4 |sort -u > ${work_path}/Novel_Indels_d0_af0/${name}_14_all_scalpel_d0_af0_RepeatMask.sites
${dir_of_perl}/perl ${select_v_ID_pl} ${work_path}/Novel_Indels_d0_af0/${name}_13_all_scalpel_d0_af0_rmdbSNP151.txt ${work_path}/Novel_Indels_d0_af0/${name}_14_all_scalpel_d0_af0_RepeatMask.sites 4 |cut -f 1-4,7-12 > ${work_path}/Novel_Indels_d0_af0/${name}_14_all_scalpel_d0_af0_novel.txt
cut -f 4 ${work_path}/Novel_Indels_d0_af0/${name}_14_all_scalpel_d0_af0_novel.txt |sort -u > ${work_path}/Novel_Indels_d0_af0/${name}_14_all_scalpel_d0_af0_novel.sites
# strelka2_d0_af0
${dir_of_intersectBed}/intersectBed -a ${work_path}/Novel_Indels_d0_af0/${name}_13_all_strelka2_d0_af0_rmdbSNP151.txt -b ${UCSC_RepeatMask_bed} -wa -wb |cut -f 4 |sort -u > ${work_path}/Novel_Indels_d0_af0/${name}_14_all_strelka2_d0_af0_RepeatMask.sites
${dir_of_perl}/perl ${select_v_ID_pl} ${work_path}/Novel_Indels_d0_af0/${name}_13_all_strelka2_d0_af0_rmdbSNP151.txt ${work_path}/Novel_Indels_d0_af0/${name}_14_all_strelka2_d0_af0_RepeatMask.sites 4 |cut -f 1-4,7-12 > ${work_path}/Novel_Indels_d0_af0/${name}_14_all_strelka2_d0_af0_novel.txt
cut -f 4 ${work_path}/Novel_Indels_d0_af0/${name}_14_all_strelka2_d0_af0_novel.txt |sort -u > ${work_path}/Novel_Indels_d0_af0/${name}_14_all_strelka2_d0_af0_novel.sites

cat ${work_path}/Novel_Indels_d0_af0/${name}_14_all_gatk_d0_af0_novel.sites ${work_path}/Novel_Indels_d0_af0/${name}_14_all_scalpel_d0_af0_novel.sites ${work_path}/Novel_Indels_d0_af0/${name}_14_all_strelka2_d0_af0_novel.sites |sort |uniq -c  |awk -F" " '$1==3{print $2}' > ${work_path}/Novel_Indels_d0_af0/${name}_14_overlap_d0_af0_novel.sites

# 15. variants
# C-to-T
awk '{if(($5~"[Cc]")&&($6~"[Tc]")) print $0}' ${work_path}/Novel_Indels_d0_af0/${name}_14_all_gatk_d0_af0_novel.txt > ${work_path}/Novel_Indels_d0_af0/${name}_15_all_gatk_d0_af0_C2T.txt
awk '{if(($5~"[Cc]")&&($6~"[Tc]")) print $0}' ${work_path}/Novel_Indels_d0_af0/${name}_14_all_scalpel_d0_af0_novel.txt > ${work_path}/Novel_Indels_d0_af0/${name}_15_all_scalpel_d0_af0_C2T.txt
awk '{if(($5~"[Cc]")&&($6~"[Tc]")) print $0}' ${work_path}/Novel_Indels_d0_af0/${name}_14_all_strelka2_d0_af0_novel.txt > ${work_path}/Novel_Indels_d0_af0/${name}_15_all_strelka2_d0_af0_C2T.txt
cat <(cut -f4  ${work_path}/Novel_Indels_d0_af0/${name}_15_all_gatk_d0_af0_C2T.txt|sort -u) <(cut -f4  ${work_path}/Novel_Indels_d0_af0/${name}_15_all_scalpel_d0_af0_C2T.txt|sort -u) <(cut -f4  ${work_path}/Novel_Indels_d0_af0/${name}_15_all_strelka2_d0_af0_C2T.txt|sort -u)  |sort |uniq -c  |awk -F" " '$1==3{print $2}' > ${work_path}/Novel_Indels_d0_af0/${name}_15_overlap_d0_af0_C2T.sites
# G-to-A
awk '{if(($5~"[Gg]")&&($6~"[Aa]")) print $0}' ${work_path}/Novel_Indels_d0_af0/${name}_14_all_gatk_d0_af0_novel.txt > ${work_path}/Novel_Indels_d0_af0/${name}_15_all_gatk_d0_af0_G2A.txt
awk '{if(($5~"[Gg]")&&($6~"[Aa]")) print $0}' ${work_path}/Novel_Indels_d0_af0/${name}_14_all_scalpel_d0_af0_novel.txt > ${work_path}/Novel_Indels_d0_af0/${name}_15_all_scalpel_d0_af0_G2A.txt
awk '{if(($5~"[Gg]")&&($6~"[Aa]")) print $0}' ${work_path}/Novel_Indels_d0_af0/${name}_14_all_strelka2_d0_af0_novel.txt > ${work_path}/Novel_Indels_d0_af0/${name}_15_all_strelka2_d0_af0_G2A.txt
cat <(cut -f4 ${work_path}/Novel_Indels_d0_af0/${name}_15_all_gatk_d0_af0_G2A.txt|sort -u) <(cut -f4 ${work_path}/Novel_Indels_d0_af0/${name}_15_all_scalpel_d0_af0_G2A.txt|sort -u) <(cut -f4  ${work_path}/Novel_Indels_d0_af0/${name}_15_all_strelka2_d0_af0_G2A.txt|sort -u) |sort |uniq -c |awk -F" " '$1==3{print $2}' > ${work_path}/Novel_Indels_d0_af0/${name}_15_overlap_d0_af0_G2A.sites

#mv ${work_path}/Novel_Indels_d0_af0/${name}_12* ${work_path}/Novel_Indels_d0_af0/${name}_13* ${work_path}/Novel_Indels_d0_af0/${name}_14* ${work_path}/Novel_Indels_d0_af0/${name}_15* ${work_path}/Novel_Indels_d0_af0
#%s#\${work_path}/${name}_1[2345]#\${work_path}/Novel_Indels_d0_af0/${name}_


