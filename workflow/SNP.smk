"""TODO:
add copyright
https://gencore.bio.nyu.edu/variant-calling-pipeline-gatk4/
"""
configfile: "config.yaml"


""" Main pipeline """
rule snp_fastp:
  input:
    r1 = "/datf/enzedeng/data/Lung/data/SNP/{sample}.R1.fastq.gz",
    r2 = "/datf/enzedeng/data/Lung/data/SNP/{sample}.R2.fastq.gz",
  output:
    r1 = "SNP/01.fastp/{sample}.R1.fq.gz",
    r2 = "SNP/01.fastp/{sample}.R2.fq.gz",
  threads: 20
  shell:
    "{config[fp]} --thread={threads} "
    "--in1={input.r1} --in2={input.r2} "
    "--out1={output.r1} --out2={output.r2}"

rule snp_bwa:
  input:
    r1 = "SNP/01.fastp/{sample}.R1.fq.gz",
    r2 = "SNP/01.fastp/{sample}.R2.fq.gz",
  output:
    "SNP/02.bwa/{sample}.bam"
  params:
    "@RG\\\\tID:{sample}\\\\tSM:{sample}\\\\tPL:ILM\\\\tLB:{sample}"
  threads: 40
  shell:
    "{config[bwa]} mem -t {threads} "
    "-R {params} "
    "{config[bwa_ref]} "
    "{input.r1} {input.r2} "
    "| {config[st]} view -bS > {output}"

rule snp_mark_duplicates:
  input:
    "SNP/02.bwa/{sample}.bam"
  output:
    "SNP/02.bwa/{sample}.dedup.bam"
  params:
    metrics = "SNP/02.bwa/{sample}.marded_dup_metrics.txt",
  threads: 20
  shell:
    "export PATH=/datf/enzedeng/mambaforge/envs/sc/bin/:$PATH; "
    "{config[gatk]} MarkDuplicatesSpark "
    "-I {input} -O {output} -M {params.metrics} "
    "--spark-master local[{threads}]"

rule snp_haplotype_caller:
  input:
    "SNP/02.bwa/{sample}.dedup.bam"
  output:
    "SNP/03.snp/{sample}.raw_variants.vcf"
  shell:
    "export PATH=/datf/enzedeng/mambaforge/envs/sc/bin/:$PATH; "
    "{config[gatk]} HaplotypeCaller "
    "-R {config[bwa_ref]} -I {input} -O {output}"

rule snp_extract_SNPs_and_Indels:
  input:
    "SNP/03.snp/{sample}.raw_variants.vcf"
  output:
    snps = "SNP/03.snp/{sample}.raw_snps.vcf",
    indels = "SNP/03.snp/{sample}.raw_indels.vcf",
  shell:
    "export PATH=/datf/enzedeng/mambaforge/envs/sc/bin/:$PATH; "
    "{config[gatk]} SelectVariants -select-type SNP "
    "-R {config[bwa_ref]} -V {input} -O {output.snps} "
    "&& {config[gatk]} SelectVariants -select-type INDEL "
    "-R {config[bwa_ref]} -V {input} -O {output.indels}"

rule snp_filter_snps:
  input:
    "SNP/03.snp/{sample}.raw_snps.vcf"
  output:
    "SNP/03.snp/{sample}.filtered_snps.vcf"
  shell:
    'export PATH=/datf/enzedeng/mambaforge/envs/sc/bin/:$PATH; '
    '{config[gatk]} VariantFiltration '
    '-R {config[bwa_ref]} -V {input} -O {output} '
    '-filter-name "QD_filter" -filter "QD < 2.0" '
    '-filter-name "FS_filter" -filter "FS > 60.0" '
    '-filter-name "MQ_filter" -filter "MQ < 40.0" '
    '-filter-name "SOR_filter" -filter "SOR > 4.0" '
    '-filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" '
    '-filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"'

rule snp_filter_indels:
  input:
    "SNP/03.snp/{sample}.raw_indels.vcf"
  output:
    "SNP/03.snp/{sample}.filtered_indels.vcf"
  shell:
    'export PATH=/datf/enzedeng/mambaforge/envs/sc/bin/:$PATH; '
    '{config[gatk]} VariantFiltration '
    '-R {config[bwa_ref]} -V {input} -O {output} '
    '-filter-name "QD_filter" -filter "QD < 2.0" '
    '-filter-name "FS_filter" -filter "FS > 200.0" '
    '-filter-name "SOR_filter" -filter "SOR > 10.0"'

rule snp_exclude_filtered_variants:
  input:
    snps = "SNP/03.snp/{sample}.filtered_snps.vcf",
    indels = "SNP/03.snp/{sample}.filtered_indels.vcf",
  output:
    snps = "SNP/03.snp/{sample}.bqsr_snps.vcf",
    indels = "SNP/03.snp/{sample}.bqsr_indels.vcf",
  shell:
    'export PATH=/datf/enzedeng/mambaforge/envs/sc/bin/:$PATH; '
    '{config[gatk]} SelectVariants '
    '--exclude-filtered -V {input.snps} -O {output.snps} '
    '&& {config[gatk]} SelectVariants '
    '--exclude-filtered -V {input.indels} -O {output.indels}'

rule snp_base_quality_score_recalibration:
  input:
    dedup = "SNP/02.bwa/{sample}.dedup.bam",
    snps = "SNP/03.snp/{sample}.bqsr_snps.vcf",
    indels = "SNP/03.snp/{sample}.bqsr_indels.vcf",
  output:
    "SNP/03.snp/{sample}.recal_data.table"
  shell:
    'export PATH=/datf/enzedeng/mambaforge/envs/sc/bin/:$PATH; '
    '{config[gatk]} BaseRecalibrator '
    '--known-sites {input.snps} --known-sites {input.indels} '
    '-R {config[bwa_ref]} -I {input.dedup} -O {output} '

rule snp_apply_BQSR:
  input:
    dedup = "SNP/02.bwa/{sample}.dedup.bam",
    table = "SNP/03.snp/{sample}.recal_data.table",
  output:
    "SNP/03.snp/{sample}.recal_reads.bam"
  shell:
    'export PATH=/datf/enzedeng/mambaforge/envs/sc/bin/:$PATH; '
    '{config[gatk]} ApplyBQSR '
    '-R {config[bwa_ref]} -I {input.dedup} -bqsr {input.table} -O {output} '

rule snp_call_variants:
  input:
    "SNP/03.snp/{sample}.recal_reads.bam"
  output:
    "SNP/03.snp/{sample}.raw_variants_recal.vcf"
  shell:
    'export PATH=/datf/enzedeng/mambaforge/envs/sc/bin/:$PATH; '
    '{config[gatk]} HaplotypeCaller '
    '-R {config[bwa_ref]} -I {input} -O {output} '
    
rule snp_extract_SNPs_and_indels_again:
  input:
    "SNP/03.snp/{sample}.raw_variants_recal.vcf"
  output:
    snps = "SNP/03.snp/{sample}.raw_snps_recal.vcf",
    indels = "SNP/03.snp/{sample}.raw_indels_recal.vcf",
  shell:
    "export PATH=/datf/enzedeng/mambaforge/envs/sc/bin/:$PATH; "
    "{config[gatk]} SelectVariants -select-type SNP "
    "-R {config[bwa_ref]} -V {input} -O {output.snps} "
    "&& {config[gatk]} SelectVariants -select-type INDEL "
    "-R {config[bwa_ref]} -V {input} -O {output.indels}"

rule snp_filter_snps_again:
  input:
    "SNP/03.snp/{sample}.raw_snps_recal.vcf"
  output:
    "SNP/03.snp/{sample}.filtered_snps_final.vcf"
  shell:
    'export PATH=/datf/enzedeng/mambaforge/envs/sc/bin/:$PATH; '
    '{config[gatk]} VariantFiltration '
    '-R {config[bwa_ref]} -V {input} -O {output} '
    '-filter-name "QD_filter" -filter "QD < 2.0" '
    '-filter-name "FS_filter" -filter "FS > 60.0" '
    '-filter-name "MQ_filter" -filter "MQ < 40.0" '
    '-filter-name "SOR_filter" -filter "SOR > 4.0" '
    '-filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" '
    '-filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"'

rule snp_filter_indels_again:
  input:
    "SNP/03.snp/{sample}.raw_indels_recal.vcf"
  output:
    "SNP/03.snp/{sample}.filtered_indels_final.vcf"
  shell:
    'export PATH=/datf/enzedeng/mambaforge/envs/sc/bin/:$PATH; '
    '{config[gatk]} VariantFiltration '
    '-R {config[bwa_ref]} -V {input} -O {output} '
    '-filter-name "QD_filter" -filter "QD < 2.0" '
    '-filter-name "FS_filter" -filter "FS > 200.0" '
    '-filter-name "SOR_filter" -filter "SOR > 10.0"'

rule snp_filter_homozygous_snps:
  input:
    "SNP/03.snp/{sample}.filtered_snps_final.vcf"
  output:
    "SNP/04.homozygous_snps/{sample}.homozygous_snps.vcf"
  shell:
    '''{config[bt]} filter -i 'FORMAT/GT=="1/1"&&FORMAT/DP>=10' {input} > {output} && '''
    '''{config[bt]} stats {output} > {output}.stats'''

rule snp_filter_gene_snps:
  input:
    "SNP/04.homozygous_snps/{sample}.homozygous_snps.vcf"
  output:
    "SNP/05.gene_snps/{sample}.gene_snps.vcf"
  params:
    "SNP/05.gene_snps/{sample}"
  shell:
    '{config[py]} {config[extract_snps_in_genes]} '
    '--gdb {config[gdb]} --gtf {config[cr_gtf]} --vcf {input} --outpre {params}'
