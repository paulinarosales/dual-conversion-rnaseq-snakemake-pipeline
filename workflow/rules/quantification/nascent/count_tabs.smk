def _input_geneset(wildcards):
    return expand('resources/external/gencode_{release}/{genome}.geneset.tsv', release=GENCODE_RELEASE, genome=GENOME)

rule nascent_genecounts:
    input:
        bakr_metaTSV = 'results/conversion_tables/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_conversionCounts.metadata.tsv',
        feature_ctsTSV = 'results/counts/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_counts.genes.tsv',
        genesetTSV = _input_geneset
    output:
        nascent_gene_ctsTSV = 'results/counts/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_counts.genes.nascent.tsv'
    log:
        'logs/count_tabs/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_nascent_genecounts.log'
    params:
        # min_reads = config['NASCENT_GENECOUNTS']['MIN_TOTAL_READS'],
        min_conv = config['NASCENT_GENECOUNTS']['MIN_CONVERSIONS']
    threads: 24
    resources:
        mem = '32G'
    conda:
        '../../../envs/downstream/r-basic.yaml'
    script:
        '../../../scripts/quantification/nascent_genecounts.R'


rule merge_cts:
    input:
        cts_files = TARGETS['nascent_cts'],
        genesetTSV = _input_geneset
    output:
        'results/counts/all_samples/{counts}.matrix.tsv'
    log:
        'logs/count_tabs/all_samples/merge_{counts}.log'
    params:
        min_reads_th = config['NASCENT_GENECOUNTS']['MIN_TOTAL_READS']
    threads: 12
    conda:
        '../../../envs/downstream/r-basic.yaml'
    script:
        '../../../scripts/quantification/merge_nascent_cts.R'



rule tpm_cts:
    input:
        'results/counts/all_samples/{counts}.matrix.tsv'
    output:
        'results/counts/all_samples/{counts}.matrix.tpm.tsv'
    log:
        'logs/count_tabs/all_samples/tpm_{counts}.log'
    threads: 12
    conda:
        '../../../envs/downstream/r-basic.yaml'
    script:
        '../../../scripts/quantification/tpm_tabs.R'



rule nascent_fract_diff_cts:
    input:
        sample_manifestTSV = config['SAMPLE_MANIFEST'],
        fraction_ctsTSV = 'results/counts/all_samples/nascent_fraction.matrix.tsv'
    output:
        'results/counts/all_samples/nascent_fraction.matrix.diff.tsv'
    log:
        'logs/count_tabs/all_samples/diff_nascent_fraction.log'
    threads: 12
    conda:
        '../../../envs/downstream/r-basic.yaml'
    script:
        '../../../scripts/quantification/decay_rate_tpm_tabs.R'

rule nascent_content:
    input:
        cts_files = TARGETS['nascent_cts']
    output:
        'results/counts/all_samples/nascent_reads_content.tsv'
    log:
        'logs/count_tabs/all_samples/nascent_reads_content.log'
    threads: 12
    conda:
        '../../../envs/downstream/r-basic.yaml'
    script:
        '../../../scripts/quantification/nascent_percentage.R'


# rule m6a_content:
#     input:
#         cts_files = TARGETS['nascent_cts']
#     output:
#         'results/counts/all_samples/nascent_reads_content.tsv'
#     log:
#         'logs/count_tabs/all_samples/nascent_reads_content.log'
#     threads: 12
#     conda:
#         '../../../envs/downstream/r-basic.yaml'
#     script:
#         '../../../scripts/quantification/nascent_percentage.R'