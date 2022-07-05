rule collect_allelic_counts:
  input:
    "results/alignment/recal/{sample}.bqsr.bam",
  output:
    "results/zygosity/counts/{sample}.allelicCounts.tsv",
  conda:
    "../envs/gatk.yaml",
  log:
    "logs/gatk/collectalleliccounts/{sample}.log"
  params:
    conda=config['env']['conda_shell'],
    env=directory(config['env']['preprocess']),
    ref=config['common']['genome'],
    target=config['params']['gatk']['collectalleliccounts']['target'],
  shell:
    """
    source {params.conda} && conda activate {params.env};
    gatk --java-options '-Xmx20G  -XX:ParallelGCThreads=4' CollectAllelicCounts "
    -I {input}
    -R {params.ref}
    -L {params.target}
    -O {output} > {log} 2>&1
    """

rule categorizeAD_gatk:
  input:
    "results/zygosity/counts/{sample}.allelicCounts.tsv",
  output:
    intermediate=temp("results/zygosity/counts/{sample}_out.tmp"),
    simple=temp("results/zygosity/counts/{sample}_out.tsv"),
  params:
    conda=config['env']['conda_shell'],
    env=directory(config['env']['preprocess']),
    ref=2,
    alt=3,
  shell:
    """
    source {params.conda} && conda activate {params.env};
    grep -v '^@' {input} | tail -n +2 > {output.intermediate};
    perl workflow/scripts/allelic_count_helper.pl categorize
    {output.intermediate} {params.ref} {params.alt} > {output.simple}
    """

rule aggregate_AD:
  input:
    expand("results/zygosity/counts/{sample}_out.tsv", sample=samples.index),
  output:
    "results/zygosity/AD/aggregate.csv",
  params:
    conda=config['env']['conda_shell'],
    env=directory(config['env']['preprocess']),
  shell:
    """
    source {params.conda} && conda activate {params.env};
    paste -d',' {input} > {output};
    """
rule get_AD_lines:
  input:
    "results/zygosity/AD/aggregate.csv",
  output:
    "results/zygosity/AD/aggregate_lines.txt",
  params:
    conda=config['env']['conda_shell'],
    env=directory(config['env']['preprocess']),
    min_n=2,
  shell:
    """
    source {params.conda} && conda activate {params.env};
    perl workflow/scripts/allelic_count_helper.pl setlines
    {input} {params.min_n} > {output};
    """
rule filt_AD_raw:
  input:
    filt_lines="results/zygosity/AD/aggregate_lines.txt",
    aggregate_raw="results/zygosity/AD/aggregate.csv",
  output:
    "results/zygosity/AD/aggregate_filt.csv",
  params:
    conda=config['env']['conda_shell'],
    env=directory(config['env']['preprocess']),
    samples=expand("{sample}", sample=samples.index),
  shell:
    "source {params.conda} && conda activate {params.env};"
    "echo {params.samples} | sed 's/\s/,/g' > {output}; "
    "perl workflow/scripts/allelic_count_helper.pl getlines "
    "{input.filt_lines} {input.aggregate_raw} >> {output}; "

rule filt_AD_dbsnp:
  input:
    filt_lines="results/zygosity/AD/aggregate_lines.txt",
  output:
    "results/zygosity/AD/aggregate_pos.txt",
  params:
    conda=config['env']['conda_shell'],
    env=directory(config['env']['preprocess']),
    target=config['params']['gatk']['collectalleliccounts']['target'],
  shell:
    "source {params.conda} && conda activate {params.env};"
    "perl workflow/scripts/allelic_count_helper.pl getlines "
    "{input.filt_lines} {params.target} > {output}; "

rule run_wp_zygosity:
  input:
    het_pos='results/zygosity/AD/aggregate_pos.txt',
    het_cnt='results/zygosity/AD/aggregate_filt.csv',
  output:
    pdfs=report(expand("results/zygosity/wadingpool/hmmfit_{sample}.pdf", sample=samples.index),
                caption="../report/wp_zygosity.rst", category="CNV"),
    tbls=expand("results/zygosity/wadingpool/hmmfit_{sample}.tsv", sample=samples.index),
    rdas=expand("results/zygosity/wadingpool/hmmfit_{sample}.rda", sample=samples.index),
  params:
    outdir='results/zygosity/wadingpool',
    cndir='results/cnv/ichorcna',
    genome=config['common']['build'],
    maxstate=300,
    conda=config['env']['conda_shell'],
    env=directory(config['env']['preprocess']),

  shell:
    "source {params.conda} && conda activate {params.env};"
    "wp_zygosity.R "
    "--hetposf {input.het_pos} "
    "--hetf {input.het_cnt} "
    "--cnpath {params.cndir} "
    "--outdir {params.outdir} "
    "--genome {params.genome} "
    "--maxstate {params.maxstate} "

rule run_wp_identity_autosome:
  input:
    het_cnt='results/zygosity/AD/aggregate_filt.csv',
  output:
    pdf=report("results/sampleid/autosome/similarity_autosome.pdf",
                caption="../report/autoid.rst", category="genotypeID"),
    tbl=report("results/sampleid/autosome/similarity_autosome.tsv",
                caption="../report/autoid.rst", category="genotypeID"),
  params:
    outdir='results/sampleid/autosome',
    midpoint=0.03,
    conda=config['env']['conda_shell'],
    env=directory(config['env']['preprocess']),

  shell:
    "source {params.conda} && conda activate {params.env};"
    "wp_identity.R "
    "--sample {input.het_cnt} "
    "--midpoint {params.midpoint} "
    "--mode autosome "
    "--outdir {params.outdir} "
