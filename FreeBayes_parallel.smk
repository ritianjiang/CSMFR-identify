## A simple example snakemake .smk file for parallelising freebayes
## Uses a fasta_generate_regions to split the genome into regions of equal size based on the .fai index
## As snakemake automatically moves each cpu core to the next genome chunk, this works out faster
## than the freebayes-parallel wrapper.
## This .smk file assumes we have a list of the bam files called bam.list
## This .smk file splits the genome by chromosome, which of course, is not necessary.
## One will want to edit the paths (for example, the path to bam files)


# these parameters should usually be stored in the snakemake configuration file (config.yaml) and accessed e.g. config['ref']
#samples = ['sample101']
reference = "/home/wanght/refNome/BasicResource/Ref/Homo_sapiens_assembly38_gatk.fasta"
chroms = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"]
nchunks = 9
#samples = ["sample101"]

chunks = list(range(1,nchunks+1))

rule all:
    input: "results/variants/{sample}.raw.fb.sorted.vcf.gz"

rule GenerateFreebayesRegions:
    input:
        ref_idx = reference,
        index = reference + ".fai",
       # bams = "resources/alignments/{sample}.dedup.BQSR.sorted.bam"
    output:
        regions = expand("resources/regions/genome.{chrom}.region.{i}.bed", chrom=chroms, i = chunks)
    log:
        "logs/GenerateFreebayesRegions.log"
    params:
        chroms = chroms,
        nchunks = nchunks
    shell:
        # "../scripts/GenerateFreebayesRegions.R" # This is located in the scripts/ directory of freebayes
        "fasta_generate_regions.py --chunks --chromosomes {params.chroms} --bed resources/regions/genome {input.index} {params.nchunks}"


rule VariantCallingFreebayes:
    input:
        bams = "resources/alignments/{sample}.dedup.BQSR.sorted.bam",
        index = "resources/alignments/{sample}.dedup.BQSR.sorted.bam.bai", 
        ref = reference,
        regions = "resources/regions/genome.{chrom}.region.{i}.bed"
    output:
        "results/variants/vcfs/{sample}/{chrom}/variants.{i}.vcf"
    log:
        "logs/VariantCallingFreebayes/{chrom}.{i}.log"
    threads:1
    shell:  "freebayes -C 3 --min-alternate-qsum 40 --pvar 0.0001 --use-mapping-quality --posterior-integration-limits 1,3 --genotype-variant-threshold 4 --site-selection-max-iterations 3 --genotyping-max-iterations 25 --max-complex-gap 3 -f {input.ref} -t {input.regions} {input.bams} > {output} 2> {log}"

#rule GenerateVcfList:
#    input: expand("results/variants/vcfs/tmp/{chrom}/variants.{i}.vcf", chrom=chroms, i=chunks)
#    output: expand("results/variants/vcfs/tmp/{sample}.vcflist",sample=samples)
#    threads:1
#    shell: "echo {input} > {output}"

#rule ModifyVcfList:
#    input: expand("results/variants/vcfs/tmp/{sample}.vcflist",sample=samples)
#    output: expand("results/variants/vcfs/tmp/{sample}.vcflist2",sample=samples)
#    threads:1
#    shell: "sed 's/\s/\\n/g' {input} > {output}"

rule ConcatVCFs:
    input:
        calls = expand("results/variants/vcfs/{wildcards.sample}/{chrom}/variants.{i}.vcf", chrom=chroms, i=chunks)
    output:
        "results/variants/vcfs/{sample}/{sample}.vcf"
    log:
        "logs/ConcatVCFs/{sample}_Concat.log"
    threads:1
    shell:  
        "bcftools concat -f {input.calls} | vcfuniq > {output} 2> {log}"

rule ConverVcfgzs:
    input:
        "results/variants/vcfs/tmp/{sample}.vcf"
    output:
        temp("results/variants/vcfs/tmp/{sample}.vcf.gz")
    log:
        "logs/ConcatVCFs/{sample}_convert.log"
    threads:1
    shell:
        "bcftools view -Oz -o {output} {input} 2> {log}"

rule SortVcfgzs:
    input:
        "results/variants/vcfs/tmp/{sample}.vcf.gz"
    output:
        "results/variants/{sample}.raw.fb.sorted.vcf.gz"
    log:
        "ConcatVCFs/{sample}_sortvcfgz.log"
    threads:1
    shell:
        "bcftools sort -Oz -o {output} {input} 2> {log}"
