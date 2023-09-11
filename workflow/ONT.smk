"""TODO:
add copyright
"""
configfile: "config.yaml"

""" Main pipeline """
rule pychopper:
  input:
    "/datf/enzedeng/data/Lung/data/ONT/{sample}.fq.gz"
  output:
    "ONT/01.pychopper/{sample}.full_length_with_primer.fq"
  params:
    primers = "/datf/enzedeng/data/Lung/config/pychopper/primer.fa",
    config = "/datf/enzedeng/data/Lung/config/pychopper/primer.txt",
    cutoff = 0.5,
    outpre = "ONT/01.pychopper/{sample}",
  threads: 20
  shell:
    "{config[pc]} -m edlib -b {params.primers} -c {params.config} "
    "-q {params.cutoff} "
    "-r {params.outpre}.report.pdf "
    "-S {params.outpre}.statistics.tsv "
    "-K {params.outpre}.qc_fail.fq "
    "-u {params.outpre}.unclass.fq "
    "-l {params.outpre}.len_fail.fq "
    "-w {params.outpre}.rescue.fq "
    "-A {params.outpre}.aln_hits.bed "
    "-D {params.outpre}.read_stats.txt "
    "-t {threads} -p "
    "{input} {output}"

rule extract_barcode:
  input:
    "ONT/01.pychopper/{sample}.full_length_with_primer.fq"
  output:
    "ONT/01.extract/{sample}.full_length.fq"
  params:
    rescue = "ONT/01.pychopper/{sample}.rescue.fq",
    pre = "ONT/01.extract/{sample}",
  shell:
    "cat {input} {params.rescue} | "
    "{config[py]} {config[find_polyT]} --outpre {params.pre} - {output}"

rule ont_minimap2:
  input:
    "ONT/01.extract/{sample}.full_length.fq"
  output:
    "ONT/02.map_reads/{sample}.sorted.bam"
  params:
    "ONT/02.map_reads/{sample}.bam"
  threads: 40
  shell:
    "{config[mm2]} -ax splice "
    "--MD --eqx --secondary=no -t {threads} -y -Y -L "
    "{config[mm2_ref]} {input} | {config[st]} view -bS > {params} "
    "&& {config[sb]} sort {params}"

rule ont_correct_barcode:
  input:
    "ONT/02.map_reads/{sample}.sorted.bam"
  output:
    "ONT/03.correct_barcode/{sample}.barcode_counts.tsv",
    "ONT/03.correct_barcode/{sample}.barcode_corrected.bam",
  params:
    edit_dist=2,
    outpre="ONT/03.correct_barcode/{sample}"
  shell:
    "{config[py]} {config[correct_barcodes]} --inbam {input} --outpre {params.outpre} "
    "--starcode_path {config[sc]} -w {config[wbcl]} -d {params.edit_dist}"

rule ont_spliced_bam2gff:
  input:
    "ONT/02.map_reads/{sample}.sorted.bam"
  output:
    "ONT/04.gffcompare/{sample}.gff"
  threads: 10
  shell:
    "{config[bg]} -M -S -g -t {threads} {input} > {output}"
# -M    Input is from minimap2.
# -S    Do NOT discard secondary and supplementary alignments.
# -g    Use strand tag as feature orientation then read strand if not available.

    
rule ont_gffcompare:
  input:
    "ONT/04.gffcompare/{sample}.gff"
  output:
    "ONT/04.gffcompare/{sample}.{sample}.gff.tmap"
  params:
    "ONT/04.gffcompare/{sample}"
  shell:
    "{config[gc]} -r {config[cr_gtf]} {input} -o {params}"

rule ont_filter_tmap:
  input:
    "ONT/04.gffcompare/{sample}.{sample}.gff.tmap"
  output:
    "ONT/04.gffcompare/{sample}.filter.tmap"
  shell:
    '''awk -F"\\t" '{{cnt[$4]++;l[$4]=$0;}}END{{for(k in cnt){{if(cnt[k]==1){{print l[k]}}}}}}' {input} > {output}'''

rule ont_isoquant:
  input:
    "ONT/02.map_reads/{sample}.sorted.bam"
  output:
    "ONT/04.IsoQuant/{sample}/OUT/OUT.read_assignments.tsv"
    "ONT/04.IsoQuant/{sample}/OUT/OUT.transcript_model_reads.tsv"
  params:
    outdir="ONT/04.IsoQuant/{sample}"
  threads: 40
  shell:
    "{config[py]} {config[iq]} "
    "--reference {config[cr_gfa]} --genedb {config[cr_gtf]} --complete_genedb "
    "--sqanti_output --count_exons --output {params.outdir} "
    "--threads {threads} --data_type ont --bam {input}"

rule ont_merge_isoquant_to_bam:
  input:
    bam="ONT/03.correct_barcode/{sample}.barcode_corrected.bam",
    tsv="ONT/04.IsoQuant/{sample}/OUT/OUT.read_assignments.tsv"
  output:
    "ONT/05.merge_isoquant/{sample}.isoquant.bam"
  shell:
    "{config[py]} {config[merge_isoquant_to_bam]} "
    "--gtf {config[cr_gtf]} --inbam {input.bam} --tsv {input.tsv} --outbam {output}"

rule ont_correct_umi_for_isoquant:
  input:
    "ONT/05.merge_isoquant/{sample}.isoquant.bam"
  output:
    "ONT/06.correct_umi/{sample}.UMI_corrected.bam"
  params:
    "ONT/06.correct_umi/{sample}"
  shell:
    "{config[py]} {config[correct_UMI]} "
    "--inbam {input} --model ont --outpre {params}"

rule ont_matrix_for_isoquant:
  input:
    "ONT/06.correct_umi/{sample}.UMI_corrected.bam"
  output:
    "ONT/07.matrix/{sample}.gene.filtered.h5ad",
    "ONT/07.matrix/{sample}.transcript.filtered.h5ad",
  params:
    outpre="ONT/07.matrix/{sample}",
    cells=lambda wildcards: cell_num[wildcards.sample],
  shell:
    "{config[py]} {config[create_count_matrix]} "
    "--inbam {input} --outpre {params.outpre} --expect_cells {params.cells}"

rule ont_deduplicate:
  input:
    "ONT/06.correct_umi/{sample}.UMI_corrected.bam"
  output:
    "ONT/06.correct_umi/{sample}.dedup.bam"
  shell:
    "{config[py]} {config[deduplicate]} --inbam {input} --outbam {output}"

rule ont_deepvariant:
  input:
    "ONT/06.correct_umi/{sample}.dedup.bam"
  output:
    "ONT/08.snp/{sample}.vcf.gz"
  threads: 40
  shell:
    "{config[sg]} exec "
    "--bind .:/res --bind {config[mm2_ref_dir]}:/ref "
    "{config[dv]} "
    "/opt/deepvariant/bin/run_deepvariant "
    "--model_type ONT_R104 "
    "--ref /ref/{config[mm2_ref_file]} --num_shards {threads} "
    "--reads /res/{input} --output_vcf /res/{output}"

rule ont_filter_heterozygous_snps:
  input:
    "ONT/08.snp/{sample}.vcf.gz"
  output:
    "ONT/08.snp/{sample}.heterozygous_snps.vcf"
  params:
    "ONT/08.snp/{sample}.filtered_snp.vcf"
  shell:
    '''{config[bt]} view -v snps {input} | '''
    '''{config[bt]} filter -i 'FILTER=="PASS"' | '''
    '''{config[bt]} view -Oz -o {params} && '''
    '''{config[bt]} filter -i 'FORMAT/GT=="0/1"&&FORMAT/DP>=10' {params} | '''
    '''{config[bt]} view -Oz -o {output} && '''
    '''{config[bt]} stats {output} > {output}.stats'''

rule ont_filter_gene_snps:
  input:
    "ONT/08.snp/{sample}.heterozygous_snps.vcf"
  output:
    "ONT/08.snp/{sample}.gene_snps.vcf"
  params:
    "ONT/08.snp/{sample}"
  shell:
    '{config[py]} {config[extract_snps_in_genes]} '
    '--gdb {config[gdb]} --gtf {config[cr_gtf]} --vcf {input} --outpre {params}'

rule ont_classify_reads:
  input:
    bam="ONT/06.correct_umi/{sample}.dedup.bam",
    snp="ONT/08.snp/{sample}.gene_snps.vcf",
  output:
    "ONT/09.snp_classify/{sample}.classified_reads.bam"
  params:
    outpre="ONT/09.snp_classify/{sample}",
    gene_snp="SNP/05.gene_snps/DBA2.gene_snps.tsv",
  shell:
    '{config[py]} {config[classify_reads]} '
    '--bam {input.bam} --snp {input.snp} '
    '--outpre {params.outpre} --gene_snp {params.gene_snp}'

rule ont_dedup_IsoQuant:
  input:
    bam="ONT/06.correct_umi/{sample}.dedup.bam",
    tsv="ONT/04.IsoQuant/{sample}/OUT/OUT.read_assignments.tsv"
  output:
    "ONT/10.read_assignments/{sample}.read_assignments.tsv"
  params:
    "ONT/10.read_assignments/{sample}"
  shell:
    '{config[py]} {config[dedup_IsoQuant]} '
    '--bam {input.bam} --read_assignments {input.tsv} --outpre {params}'

rule ont_IsoQuant_transcript_num:
  input:
    bam="ONT/06.correct_umi/{sample}.dedup.bam",
    tsv="ONT/04.IsoQuant/{sample}/OUT/OUT.transcript_model_reads.tsv",
    h5ad="ONT/07.matrix/{sample}.gene.filtered.h5ad",
  output:
    "ONT/10.read_assignments/{sample}.count.txt"
  params:
    "ONT/10.read_assignments/{sample}"
  shell:
    '{config[py]} {config[novel_transcript]} '
    '--bam {input.bam} --model_reads {input.tsv} --h5ad {input.h5ad} --outpre {params}'
