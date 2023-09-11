"""TODO:
add copyright
"""
configfile: "config.yaml"

rule ngs_reads_qc:
  input:
    "/datf/enzedeng/data/Lung/data/NGS/{sample}"
  output:
    "Statistics/00.quality_control/{sample}.raw.record_info.txt"
  shell:
    "{config[py]} {config[raw_length_quals]} {input} {output}"

rule ont_reads_qc:
  input:
    "/datf/enzedeng/data/Lung/data/ONT/{sample}.fq.gz"
  output:
    "Statistics/00.quality_control/{sample}.raw.record_info.txt"
  shell:
    "{config[py]} {config[raw_length_quals]} {input} {output}"

rule pac_reads_qc:
  input:
    "/datf/enzedeng/data/Lung/data/PAC/{sample}.bam"
  output:
    "Statistics/00.quality_control/{sample}.raw.record_info.txt"
  shell:
    "{config[py]} {config[raw_length_quals]} {input} {output}"

rule ont_cdna_qc:
  input:
    "ONT/01.extract/{sample}.full_length.fq"
  output:
    "Statistics/01.cDNA/{sample}.cDNA.record_info.txt"
  shell:
    "{config[py]} {config[raw_length_quals]} {input} {output}"

rule pac_cdna_qc:
  input:
    "PAC/01.longbow/{sample}.fastq"
  output:
    "Statistics/01.cDNA/{sample}.cDNA.record_info.txt"
  shell:
    "{config[py]} {config[raw_length_quals]} {input} {output}"

rule ngs_bam_flagstat:
  input:
    "NGS/01.cellranger/{sample}/outs/possorted_genome_bam.bam"
  output:
    "Statistics/02.bam_flagstat/{sample}.flagstat"
  shell:
    "{config[st]} flagstat -O tsv {input} > {output}"

rule ont_bam_flagstat:
  input:
    "ONT/02.map_reads/{sample}.sorted.bam"
  output:
    "Statistics/02.bam_flagstat/{sample}.flagstat"
  shell:
    "{config[st]} flagstat -O tsv {input} > {output}"

rule pac_bam_flagstat:
  input:
    "PAC/02.map_reads/{sample}.sorted.bam"
  output:
    "Statistics/02.bam_flagstat/{sample}.flagstat"
  shell:
    "{config[st]} flagstat -O tsv {input} > {output}"



#rule ngs_barcode_dist:
#  input:
#    "NGS/01.cellranger/{sample}/outs/possorted_genome_bam.bam"

rule tgs_barcode_dist:
  input:
    lambda wildcards: f'{wildcards.sample[:3]}/03.correct_barcode/{wildcards.sample}.barcode_counts.tsv'
  output:
    "Statistics/05.barcode_dist/{sample}.barcode_min_dist_num.tsv"
  params:
    "Statistics/05.barcode_dist/{sample}"
  shell:
    "{config[py]} {config[barcode_med]} "
    "--counts {input} --outpre {params}"
