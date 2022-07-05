rule chrM_make_interval:
    input:
        bed=rules.chrM_dict.output.bed,
        dict=rules.chrM_dict.output.dict,
    output:
        bed="resources/chrM.headered.bed",
    conda:
        "../envs/picard.yaml"
    shell:
        "picard BedToIntervalList "
        "I={input.bed} "
        "O={output.bed} "
        "SD={input.dict} "

rule collect_chr_metrics:
    input:
        bam=rules.subset_chrM.output.bam,
        ref=rules.chrM_fa.output.fa,
        bed=rules.chrM_make_interval.output.bed,
    output:
        "results/chrM/copies/{sample}.wgs_metrics"
    log:
        "logs/picard/wgs_metrics/{sample}_chrM.log"
    params:
        conda=config['env']['conda_shell'],
        env=directory(config['env']['preprocess']),

    shell:
        "source {params.conda} && conda activate {params.env};"
        "picard CollectWgsMetrics "
        "I={input.bam} "
        "O={output} "
        "R={input.ref} "
        "INTERVALS={input.bed} 2> {log}"

rule get_coverages:
    input:
        chrm=rules.collect_chr_metrics.output,
        genome=rules.collect_wgs_metrics.output,
    output:
        temp("results/chrM/copies/{sample}.tmp")
    params:
        conda=config['env']['conda_shell'],
        env=directory(config['env']['preprocess']),
    shell:
        "source {params.conda} && conda activate {params.env};"
        "echo -e {wildcards.sample}'\t'$(head -8 {input.chrm} | tail -1 | cut -f2)'\t'$(head -8 {input.genome} | tail -1 | cut -f2) > {output}"

rule estimate_mt_copies:
    input:
        covs=expand("results/chrM/copies/{sample}.tmp", sample=samples.index),
    output:
        cov="results/chrM/copies/all_cov.wgs_metrics"
    params:
        conda=config['env']['conda_shell'],
        env=directory(config['env']['preprocess']),
    shell:
        "source {params.conda} && conda activate {params.env};"
        "cat {input.covs} > {output.cov};"

rule plot_chrM_copies:
    input:
        rules.estimate_mt_copies.output
    output:
        plot=report("results/chrM/copies/chrM_copies.pdf",
                    caption="../report/chrM_copies.rst", category='chrM'),
        tbl=report("results/chrM/copies/chrM_copies.tsv",
                    caption="../report/chrM_copies.rst", category='chrM')
    params:
        conda=config['env']['conda_shell'],
        env=directory(config['env']['preprocess']),
    shell:
        "source {params.conda} && conda activate {params.env};"
        "Rscript workflow/scripts/chrM_copies.R "
        "{input} {output.tbl} {output.plot}"

rule relocate_chrm_copies:
    input:
        tbl=rules.plot_chrM_copies.output.tbl,
        plot=rules.plot_chrM_copies.output.plot,
    output:
        tbl="results/tables/chrM/chrM_copies.tsv",
        plot="results/plots/chrM/chrM_copies.pdf",
    shell:
        "cp {input.plot} {output.plot}; "
        "cp {input.tbl} {output.tbl}; "
