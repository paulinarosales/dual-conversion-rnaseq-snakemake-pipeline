def _input_geneset(wildcards):
    return expand('resources/external/gencode_{release}/{genome}.geneset.tsv', release=GENCODE_RELEASE, genome=GENOME)


rule labeled_fract_gene:
    input:
        sample_manifestTSV = config['SAMPLE_MANIFEST'],
        nascentFractTSV = 'results/counts/all_samples/nascent_fraction.matrix.tsv'
    output:
        'results/figures/conversions/all_samples/nascent_fraction_per_gene_density.pdf'
    log:
        'logs/figures/all_samples/fraction_boxplot.log'
    conda:
        '../../envs/downstream/r-basic.yaml'
    threads: 24
    resources:
        mem = '64G'
    script:
        '../../scripts/figures/label_fraction_dist.R'


rule labeled_fract_gene_vs_conv:
    input:
        sample_manifestTSV = config['SAMPLE_MANIFEST'],
        fraction_ctsTSV = 'results/counts/all_samples/nascent_fraction.matrix.tsv',
        gene_mutTSV = 'results/conversion_tables/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_conversionCounts.genes.tsv',
        genesetTSV = _input_geneset
    output:
        'results/figures/conversions/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_nascent_fraction_vs_conv_hist.pdf'
    params:
        target_mut = 'T_C'
    log:
        'logs/figures/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_fraction_vs_conv.log'
    conda:
        '../../envs/downstream/r-plotter.yaml'
    threads: 24
    resources:
        mem = '64G'
    script:
        '../../scripts/figures/mut_distribution_vs_fraction.R'



rule labeled_fract_gene_vs_halflife:
    input:
        sample_manifestTSV = config['SAMPLE_MANIFEST'],
        fraction_ctsTSV = 'results/counts/all_samples/nascent_fraction.matrix.diff.tsv',
        ref_halflifeCSV = 'resources/ref_sites/haflives_tab.csv',
        goisTSV = 'resources/ref_sites/decay_gois.tsv'
    output:
        'results/figures/decay/{sample_type}_{treatment}/{sample_type}_{treatment}_decay_vs_halflife_scatter.pdf'
    log:
        'logs/figures/{sample_type}_{treatment}/{sample_type}_{treatment}_decay_vs_halflife.log'
    conda:
        '../../envs/downstream/r-plotter.yaml'
    threads: 24
    resources:
        mem = '64G'
    script:
        '../../scripts/figures/halflifes_scatter.R'