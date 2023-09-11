
rule cellranger_count:
  input:
    "/datf/enzedeng/data/Lung/data/NGS/{sample}"
  output:
    "NGS/01.cellranger/{sample}/outs/possorted_genome_bam.bam"
  params:
    mem = 64,
    cells = lambda wildcards: cell_num[wildcards.sample],
    outdir = "NGS/01.cellranger/{sample}",
  threads: 40
  shell:
    "{config[cr]} count "
    "--transcriptome={config[cr_ref]} "
    "--localcores={threads} "
    "--localmem={params.mem} "
    "--id={wildcards.sample} "
    "--sample={wildcards.sample} "
    "--expect-cells={params.cells} "
    "--fastqs={input} "
    "&& mv {wildcards.sample}/* {params.outdir} "
    "&& rmdir {wildcards.sample}"

# rule ngs_correct_barcode:
#   input:
#     "NGS/01.cellranger/{sample}/outs/possorted_genome_bam.bam"
#   output:
#     "NGS/03.correct_barcode/{sample}.barcode_corrected.bam"
#   params:
#     edit_dist=1,
#     outpre="NGS/03.correct_barcode/{sample}"
#   shell:
#     "{config[py]} {config[corrector_barcodes]} --inbam {input} --outpre {params.outpre} "
#     "--starcode_path {config[sc]} -w {config[wbcl]} -d {params.edit_dist}"
