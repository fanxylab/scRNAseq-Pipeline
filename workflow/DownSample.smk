"""TODO:
add copyright
https://gencore.bio.nyu.edu/variant-calling-pipeline-gatk4/
"""
configfile: "config.yaml"

sample_dir = {'PAC_5h': 'PAC', 'ONT_5h': 'ONT', }

rule downsample_5h:
  input:
    lambda wildcards: f'{sample_dir[wildcards.sample]}/06.correct_umi/{wildcards.sample}.dedup.bam'
  output:
    'Downsample/01.downsample/{sample}.sample_0.1.bam'
  shell:
    '{config[sb]} view -h -t 10 -s 0.1 -f bam --subsampling-seed=42 '
    '{input} -o {output}'

rule downsample_create_matrix_5h:
  input:
    "Downsample/01.downsample/{sample}.sample_0.1.bam"
  output:
    "Downsample/02.matrix/{sample}.gene.filtered.h5ad",
    "Downsample/02.matrix/{sample}.transcript.filtered.h5ad",
  params:
    outpre="Downsample/02.matrix/{sample}",
    cells=500,
  shell:
    "{config[py]} {config[create_count_matrix]} "
    "--inbam {input} --outpre {params.outpre} --expect_cells {params.cells}"
