"""TODO:
add copyright
"""
configfile: "config.yaml"

""" Main pipeline """
rule pac_longbow_annotate:
  input:
    "/datf/enzedeng/data/Lung/data/PAC/{sample}.bam"
  output:
    "PAC/01.longbow/{sample}.annotate.bam"
  params:
    model = "/datf/enzedeng/data/Lung/config/longbow/model.json"
  threads: 20
  shell:
    "{config[lb]} annotate "
    "-m {params.model} -t {threads} -o {output} "
    "{input}"

rule pac_longbow_filter:
  input:
    "PAC/01.longbow/{sample}.annotate.bam"
  output:
    "PAC/01.longbow/{sample}.filter.bam"
  params:
    model = "/datf/enzedeng/data/Lung/config/longbow/model.json"
  shell:
    "{config[lb]} filter "
    "-m {params.model} -o {output} "
    "{input}"

rule pac_longbow_segment:
  input:
    "PAC/01.longbow/{sample}.filter.bam"
  output:
    "PAC/01.longbow/{sample}.segment.bam"
  params:
    model = "/datf/enzedeng/data/Lung/config/longbow/model.json"
  threads: 20
  shell:
    "{config[lb]} segment "
    "-m {params.model} -t {threads} -o {output} "
    "{input}"

rule pac_longbow_refine:
  input:
    "PAC/01.longbow/{sample}.segment.bam"
  output:
    "PAC/01.longbow/{sample}.segment.refine_P5.bam"
  shell:
    "{config[py]} {config[refine_segment_P5]} {input}"
  
rule pac_longbow_extract:
  input:
    "PAC/01.longbow/{sample}.segment.refine_P5.bam"
  output:
    "PAC/01.longbow/{sample}.filter_passed.bam"
  params:
    model = "/datf/enzedeng/data/Lung/config/longbow/model.json"
  shell:
    "{config[lb]} extract "
    "--force -m {params.model} -o {output} "
    "{input}"

rule pac_bam_to_fastq:
  input:
    "PAC/01.longbow/{sample}.filter_passed.bam"
  output:
    "PAC/01.longbow/{sample}.fastq"
  shell:
    "{config[st]} fastq "
    "-T CR,CY,ZU,UY,PR,PY {input} > {output}"

rule pac_minimap2:
  input:
    "PAC/01.longbow/{sample}.fastq"
  output:
    "PAC/02.map_reads/{sample}.sorted.bam"
  params:
    "PAC/02.map_reads/{sample}.bam"
  threads: 40
  shell:
    "{config[mm2]} -ax splice:hq "
    "--MD --eqx --secondary=no -t {threads} -y -Y -L "
    "{config[mm2_ref]} {input} | {config[st]} view -bS > {params} "
    "&& {config[sb]} sort {params}"

rule pac_correct_barcode:
  input:
    "PAC/02.map_reads/{sample}.sorted.bam"
  output:
    "PAC/03.correct_barcode/{sample}.barcode_counts.tsv",
    "PAC/03.correct_barcode/{sample}.barcode_corrected.bam",
  params:
    edit_dist=1,
    outpre="PAC/03.correct_barcode/{sample}"
  shell:
    "{config[py]} {config[correct_barcodes]} --inbam {input} --outpre {params.outpre} "
    "--starcode_path {config[sc]} -w {config[wbcl]} -d {params.edit_dist}"

rule pac_spliced_bam2gff:
  input:
    "PAC/02.map_reads/{sample}.sorted.bam"
  output:
    "PAC/04.gffcompare/{sample}.gff"
  threads: 10
  shell:
    "{config[bg]} -M -S -g -t {threads} {input} > {output}"
# -M    Input is from minimap2.
# -S    Do NOT discard secondary and supplementary alignments.
# -g    Use strand tag as feature orientation then read strand if not available.
    
rule pac_gffcompare:
  input:
    "PAC/04.gffcompare/{sample}.gff"
  output:
    "PAC/04.gffcompare/{sample}.{sample}.gff.tmap"
  params:
    "PAC/04.gffcompare/{sample}"
  shell:
    "{config[gc]} -r {config[cr_gtf]} {input} -o {params}"

rule pac_filter_tmap:
  input:
    "PAC/04.gffcompare/{sample}.{sample}.gff.tmap"
  output:
    "PAC/04.gffcompare/{sample}.filter.tmap"
  shell:
    '''awk -F"\\t" '{{cnt[$4]++;l[$4]=$0;}}END{{for(k in cnt){{if(cnt[k]==1){{print l[k]}}}}}}' {input} > {output}'''

# rule pac_correct_gffcompare:
#   input:
#     bam="PAC/02.map_reads/{sample}.sorted.bam",
#     tmap="PAC/04.gffcompare/{sample}.filter.tmap",
#   output:
#     "PAC/04.correct_gffcompare/{sample}.correct.bam"
#   params:
#     outpre="PAC/04.correct_gffcompare/{sample}"
#   shell:
#     "{config[py]} {config[pac_corrector]} "
#     "-w {config[wbcl]} --gtf {config[cr_gtf]} --outpre {params.outpre} "
#     "--inbam {input.bam} --tmap {input.tmap} "

rule pac_isoquant:
  input:
    "PAC/02.map_reads/{sample}.sorted.bam"
  output:
    "PAC/04.IsoQuant/{sample}/OUT/OUT.read_assignments.tsv",
    "PAC/04.IsoQuant/{sample}/OUT/OUT.transcript_model_reads.tsv",
  params:
    outdir="PAC/04.IsoQuant/{sample}"
  threads: 40
  shell:
    "{config[py]} {config[iq]} "
    "--reference {config[cr_gfa]} --genedb {config[cr_gtf]} --complete_genedb "
    "--sqanti_output --count_exons --output {params.outdir} "
    "--threads {threads} --data_type pacbio --bam {input}"

rule pac_merge_isoquant_to_bam:
  input:
    bam="PAC/03.correct_barcode/{sample}.barcode_corrected.bam",
    tsv="PAC/04.IsoQuant/{sample}/OUT/OUT.read_assignments.tsv"
  output:
    "PAC/05.merge_isoquant/{sample}.isoquant.bam"
  shell:
    "{config[py]} {config[merge_isoquant_to_bam]} "
    "--gtf {config[cr_gtf]} --inbam {input.bam} --tsv {input.tsv} --outbam {output}"

rule pac_correct_umi_for_isoquant:
  input:
    "PAC/05.merge_isoquant/{sample}.isoquant.bam"
  output:
    "PAC/06.correct_umi/{sample}.UMI_corrected.bam"
  params:
    "PAC/06.correct_umi/{sample}"
  shell:
    "{config[py]} {config[correct_UMI]} "
    "--inbam {input} --model pacbio --outpre {params}"

rule pac_matrix_for_isoquant:
  input:
    "PAC/06.correct_umi/{sample}.UMI_corrected.bam"
  output:
    "PAC/07.matrix/{sample}.gene.filtered.h5ad",
    "PAC/07.matrix/{sample}.transcript.filtered.h5ad",
  params:
    outpre="PAC/07.matrix/{sample}",
    cells=lambda wildcards: cell_num[wildcards.sample],
  shell:
    "{config[py]} {config[create_count_matrix]} "
    "--inbam {input} --outpre {params.outpre} --expect_cells {params.cells}"

rule pac_deduplicate:
  input:
    "PAC/06.correct_umi/{sample}.UMI_corrected.bam"
  output:
    "PAC/06.correct_umi/{sample}.dedup.bam"
  shell:
    "{config[py]} {config[deduplicate]} --inbam {input} --outbam {output}"

rule pac_deepvariant:
  input:
    "PAC/06.correct_umi/{sample}.dedup.bam"
  output:
    "PAC/08.snp/{sample}.vcf.gz"
  threads: 40
  shell:
    "{config[sg]} exec "
    "--bind .:/res --bind {config[mm2_ref_dir]}:/ref "
    "{config[dv]} "
    "/opt/deepvariant/bin/run_deepvariant "
    "--model_type PACBIO "
    "--ref /ref/{config[mm2_ref_file]} --num_shards {threads} "
    "--reads /res/{input} --output_vcf /res/{output}"

rule pac_filter_heterozygous_snps:
  input:
    "PAC/08.snp/{sample}.vcf.gz"
  output:
    "PAC/08.snp/{sample}.heterozygous_snps.vcf"
  params:
    "PAC/08.snp/{sample}.filtered_snp.vcf"
  shell:
    '''{config[bt]} view -v snps {input} | '''
    '''{config[bt]} filter -i 'FILTER=="PASS"' | '''
    '''{config[bt]} view -Oz -o {params} && '''
    '''{config[bt]} filter -i 'FORMAT/GT=="0/1"&&FORMAT/DP>=10' {params} | '''
    '''{config[bt]} view -Oz -o {output} && '''
    '''{config[bt]} stats {output} > {output}.stats'''

rule pac_filter_gene_snps:
  input:
    "PAC/08.snp/{sample}.heterozygous_snps.vcf"
  output:
    "PAC/08.snp/{sample}.gene_snps.vcf"
  params:
    "PAC/08.snp/{sample}"
  shell:
    '{config[py]} {config[extract_snps_in_genes]} '
    '--gdb {config[gdb]} --gtf {config[cr_gtf]} --vcf {input} --outpre {params}'

rule pac_classify_reads:
  input:
    bam="PAC/06.correct_umi/{sample}.dedup.bam",
    snp="PAC/08.snp/{sample}.gene_snps.vcf",
  output:
    "PAC/09.snp_classify/{sample}.classified_reads.bam"
  params:
    outpre="PAC/09.snp_classify/{sample}",
    gene_snp="SNP/05.gene_snps/DBA2.gene_snps.tsv",
  shell:
    '{config[py]} {config[classify_reads]} '
    '--bam {input.bam} --snp {input.snp} '
    '--outpre {params.outpre} --gene_snp {params.gene_snp}'

rule pac_dedup_IsoQuant:
  input:
    bam="PAC/06.correct_umi/{sample}.dedup.bam",
    tsv="PAC/04.IsoQuant/{sample}/OUT/OUT.read_assignments.tsv"
  output:
    "PAC/10.read_assignments/{sample}.read_assignments.tsv"
  params:
    "PAC/10.read_assignments/{sample}"
  shell:
    '{config[py]} {config[dedup_IsoQuant]} '
    '--bam {input.bam} --read_assignments {input.tsv} --outpre {params}'

rule pac_IsoQuant_transcript_num:
  input:
    bam="PAC/06.correct_umi/{sample}.dedup.bam",
    tsv="PAC/04.IsoQuant/{sample}/OUT/OUT.transcript_model_reads.tsv",
    h5ad="PAC/07.matrix/{sample}.gene.filtered.h5ad",
  output:
    "PAC/10.read_assignments/{sample}.count.txt"
  params:
    "PAC/10.read_assignments/{sample}"
  shell:
    '{config[py]} {config[novel_transcript]} '
    '--bam {input.bam} --model_reads {input.tsv} --h5ad {input.h5ad} --outpre {params}'
