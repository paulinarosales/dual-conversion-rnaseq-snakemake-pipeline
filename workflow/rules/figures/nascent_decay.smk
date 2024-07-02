
rule nascent_tmp_timeline_plot:
    input:
        tpmTSV = 'results/counts/all_samples/{counts}.matrix.tpm.tsv',
        goisTSV = 'resources/ref_sites/decay_gois.tsv',
        sample_manifestTSV = config['SAMPLE_MANIFEST']
    output:
       'results/figures/all_samples/GOIs_{counts}_tpm_decay_lineplot.pdf'
    params:
        filter_batch = 2023
    log:
        'logs/figures/all_samples/GOIs_{counts}_tpm_decay_lineplot.log'
    conda:
        '../../envs/downstream/r-plotter.yaml'
    threads: 12
    resources:
        mem = '12G'
    script:
        '../../scripts/figures/gois_expr.R'


rule nascent_fraction_timeline_plot:
    input:
        tpmTSV = 'results/counts/all_samples/nascent_fraction.matrix.tsv',
        goisTSV = 'resources/ref_sites/decay_gois.tsv',
        sample_manifestTSV = config['SAMPLE_MANIFEST']
    output:
       'results/figures/all_samples/GOIs_nascent_fraction_decay_lineplot.pdf'
    params:
        filter_batch = 2023
    log:
        'logs/figures/all_samples/GOIs_nascent_fraction_decay_lineplot.log'
    conda:
        '../../envs/downstream/r-plotter.yaml'
    threads: 12
    resources:
        mem = '12G'
    script:
        '../../scripts/figures/gois_expr.R'


rule nascent_tmp_timeline__old_plot:
    input:
        tpmTSV = 'results/counts/all_samples/{counts}.matrix.tpm.tsv',
        goisTSV = 'resources/ref_sites/decay_gois.tsv',
        sample_manifestTSV = config['SAMPLE_MANIFEST']
    output:
       'results/figures/all_samples/GOIs_{counts}_tpm_decay_lineplot_old.pdf'
    params:
        filter_batch = 2022
    log:
        'logs/figures/all_samples/GOIs_{counts}_tpm_decay_lineplot_old.log'
    conda:
        '../../envs/downstream/r-plotter.yaml'
    threads: 12
    resources:
        mem = '12G'
    script:
        '../../scripts/figures/gois_expr.R'


rule nascent_fraction_timeline__old_plot:
    input:
        tpmTSV = 'results/counts/all_samples/nascent_fraction.matrix.tsv',
        goisTSV = 'resources/ref_sites/decay_gois.tsv',
        sample_manifestTSV = config['SAMPLE_MANIFEST']
    output:
       'results/figures/all_samples/GOIs_nascent_fraction_decay_lineplot_old.pdf'
    params:
        filter_batch = 2022
    log:
        'logs/figures/all_samples/GOIs_nascent_fraction_decay_lineplot_old.log'
    conda:
        '../../envs/downstream/r-plotter.yaml'
    threads: 12
    resources:
        mem = '12G'
    script:
        '../../scripts/figures/gois_expr.R'