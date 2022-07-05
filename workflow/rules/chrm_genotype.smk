rule chrM_fa:
    input:
        fa=config['common']['genome'],
    output:
        fa="resources/chrM.fa",
        idx="resources/chrM.fa.fai",
    params:
        conda=config['env']['conda_shell'],
        env=directory(config['env']['preprocess']),
    shell:
        "source {params.conda} && conda activate {params.env};"
        "samtools faidx {input.fa} chrM > {output.fa}; "
        "samtools faidx {output.fa}; "

rule chrM_dict:
    input:
        "resources/chrM.fa",
    output:
        dict="resources/chrM.dict",
        bed="resources/chrM.bed",
    params:
        conda=config['env']['conda_shell'],
        env=directory(config['env']['preprocess']),
    shell:
        "source {params.conda} && conda activate {params.env};"
        "gatk CreateSequenceDictionary -R {input} ; "
        "echo -e 'chrM\t1\t'$(grep 'chrM' resources/chrM.dict  | cut -f3 | sed 's/LN://') > {output.bed}"

rule subset_chrM:
    input:
        "results/alignment/recal/{sample}.bqsr.bam",
    output:
        bam="results/sampleid/chrM/{sample}.chrM.bam",
    params:
        conda=config['env']['conda_shell'],
        env=directory(config['env']['preprocess']),
        param='-bh',
        chr='chrM',
    shell:
        "source {params.conda} && conda activate {params.env};"
        "samtools view {params.param} {input} {params.chr} > {output.bam} ; "

rule chrM_index:
    input:
        "results/sampleid/chrM/{sample}.chrM.bam",
    output:
        "results/sampleid/chrM/{sample}.chrM.bam.bai",
    params:
        conda=config['env']['conda_shell'],
        env=directory(config['env']['preprocess']),
    wrapper:
        "0.73.0/bio/samtools/index"

rule chrM_haplotype_caller:
    input:
        bam="results/sampleid/chrM/{sample}.chrM.bam",
        bai="results/sampleid/chrM/{sample}.chrM.bam.bai",
        ref=rules.chrM_fa.output.fa,
        refidx=rules.chrM_fa.output.fa,
        refdict=rules.chrM_dict.output.dict,
        bed=rules.chrM_dict.output.bed,
    output:
        gvcf="results/sampleid/chrM/{sample}.chrM.gvcf",
    log:
        "logs/gatk/haplotypecaller/{sample}.gvcf.chrM.log"
    params:
        extra="-L resources/chrM.bed",  # optional
        java_opts="", # optional
    resources:
        mem_mb=2048
    wrapper:
        "0.73.0/bio/gatk/haplotypecaller"

rule genotype_gvcfs:
    input:
        gvcf="results/sampleid/chrM/{sample}.chrM.gvcf",
        ref="resources/chrM.fa",
    output:
        vcf="results/sampleid/chrM/{sample}.chrM.vcf",
    log:
        "logs/gatk/haplotypecaller/{sample}.chrM.log",
    params:
        extra="--allow-old-rms-mapping-quality-annotation-data",  # optional
        java_opts="", # optional
    resources:
        mem_mb=2048
    wrapper:
        "0.73.0/bio/gatk/genotypegvcfs"

rule tabix_vcf:
    input:
        vcf="results/sampleid/chrM/{sample}.chrM.vcf"
    output:
        gz="results/sampleid/chrM/{sample}.chrM.vcf.gz",
        tbi="results/sampleid/chrM/{sample}.chrM.vcf.gz.tbi",
    log:
        "logs/gatk/genotype_checker/tabix_{sample}.log",
    params:
        conda=config['env']['conda_shell'],
        env=directory(config['env']['preprocess']),
    resources:
        mem_mb=1024
    shell:
        "source {params.conda} && conda activate {params.env};"
        "bgzip {input.vcf} ; "
        "tabix {output.gz}"

rule chrM_vcf_merge:
    input:
        vcfs=expand("results/sampleid/chrM/{sample}.chrM.vcf.gz", sample=samples.index),
    output:
        vcf="results/sampleid/chrM/merge.vcf",
    log:
        "logs/gatk/genotype_checker/merge.log",
    params:
        conda=config['env']['conda_shell'],
        env=directory(config['env']['preprocess']),
    resources:
        mem_mb=8192
    shell:
        "source {params.conda} && conda activate {params.env};"
        "vcf-merge {input.vcfs} > {output.vcf}"

rule genotype_checker:
    input:
        vcf="results/sampleid/chrM/merge.vcf",
    output:
        ntbl=report("results/sampleid/chrM/chrM_sampleid_n.tsv",
                    caption="../report/chrMid.rst", category="genotypeID"),
        jacctbl=report("results/sampleid/chrM/chrM_sampleid_jacc.tsv",
                    caption="../report/chrMid.rst", category="genotypeID"),
        plot=report("results/sampleid/chrM/chrM_sampleid.pdf",
                    caption="../report/chrMid.rst", category="genotypeID"),
    log:
        "logs/gatk/genotype_checker/merge.log",
    params:
        samples=",".join(expand("{sample}", sample=samples.index)),
        conda=config['env']['conda_shell'],
        env=directory(config['env']['preprocess']),
    resources:
        mem_mb=8192
    shell:
        "source {params.conda} && conda activate {params.env};"
        "Rscript workflow/scripts/genotypeChecker.R "
        "{input.vcf} "
        "'{params.samples}' "
        "{output.jacctbl} "
        "{output.ntbl} "
        "{output.plot} "

rule relocate_chrm_files:
    input:
        ntbl="results/sampleid/chrM/chrM_sampleid_n.tsv",
        jacctbl="results/sampleid/chrM/chrM_sampleid_jacc.tsv",
        plot="results/sampleid/chrM/chrM_sampleid.pdf",
    output:
        ntbl="results/tables/genotypeID/chrM_sampleid_n.tsv",
        jacctbl="results/tables/genotypeID/chrM_sampleid_jacc.tsv",
        plot="results/plots/genotypeID/chrM_sampleid.pdf",
    shell:
        "cp {input.plot} {output.plot}; "
        "cp {input.ntbl} {output.ntbl}; "
        "cp {input.jacctbl} {output.jacctbl}; "
