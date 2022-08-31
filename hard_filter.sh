#########  HARD FILTER   ###################
for sample in `ls ./no_filter/*.gz`
do
ss=${sample#*er/}
bcftools view -Oz -o ./hard2/${ss%_*}.hard2.vcf.gz --threads 8 \
    -i 'FILTER="PASS" & FORMAT/AF[0:0]<0.5 & %TYPE="snp" & FORMAT/DP>=15 & FORMAT/AD[0:1]>3 & FORMAT/DP<=120' \
    -P \
    -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX \
    -T ^/home/owht/Data/KIZ/data/Basic/RefNome/Human/bed/genomicSuperDups_hg38.reg \
    ${sample}
done

for sample in `ls ./hard2/*.gz`
do
bcftools index -t ${sample}
done
#######   END #################

###
mkdir gnomad30_filter2
for sample in `ls ./hard2/*.gz`
do
ss=${sample#*rd2/}
sn=${ss%%.*}
convert2annovar.pl -format vcf4 ${sample} > tmp.${sn}.avinput
annotate_variation.pl -filter -dbtype gnomad30_genome -buildver hg38 --otherinfo tmp.${sn}.avinput ~/Data/KIZ/data/Basic/RefNome/Human/annovar/
cat tmp.${sn}.avinput.hg38_gnomad30_genome_dropped | awk -F ',' '$9>0.01{print $0}' | cut -f '3,4' > tmp.reg
bcftools view -Oz -o ./gnomad30_filter2/${sn%_*}.gnomad30.vcf.gz \
    -T ^./tmp.reg \
    ${sample}
rm tmp*
done

for i in `ls ./gnomad30_filter2/*.gz`
do
bcftools index -t ${i}
done

ls ./gnomad30_filter2/*.gz > tm
bcftools merge -0 -Oz -o ./gnomad30_filter2/Allsample184_gnomad30.vcf.gz -l ./tm
bcftools view -S ../sampleinfo124 -Oz -o ./gnomad30_filter2/Allsample124_gnomad30.vcf.gz ./gnomad30_filter2/Allsample184_gnomad30.vcf.gz
bcftools index -t ./gnomad30_filter2/Allsample124_gnomad30.vcf.gz
bcftools view -C 2 -c 1 -Ov ./gnomad30_filter2/Allsample124_gnomad30.vcf.gz | grep -v '#' | cut -f '1,2' > max2.list

for sample in `cat ../sampleinfo124`
do
sn=${sample}
bcftools view -Oz -o ./gnomad30_filter2/nomore2/${sn}.gnomad30.nomore2.vcf.gz -R ./max2.list ./gnomad30_filter2/${sample}.gnomad30.vcf.gz
done

for sample in `cat ../sampleinfo124`
do
sn=${sample}
bcftools view -Oz -o ./gnomad30_filter2/raw/${sn}.gnomad30.raw.vcf.gz ./gnomad30_filter2/${sample}.gnomad30.vcf.gz
done
