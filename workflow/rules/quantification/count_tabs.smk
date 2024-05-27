def _input_transciptsInfo(wildcards):
    return expand('resources/external/gencode_{release}/{genome}.transcripts_info.tsv', release=GENCODE_RELEASE, genome=GENOME)

def _input_geneset(wildcards):
    return expand('resources/external/gencode_{release}/{genome}.geneset.tsv', release=GENCODE_RELEASE, genome=GENOME)

rule nascent_genecounts:
    input:
        bakr_mergedTSV = 'results/conversion_tables/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_conversionCounts.mergedFeatures.tsv',
        tx_infoTSV = _input_transciptsInfo,
        genesetTSV = _input_geneset
    output:
        nascent_ctsTSV = 'results/counts/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_genecounts.nascent.tsv'
    log:
        'logs/count_tabs/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_nascent_genecounts.log'
    params:
        min_reads = config['NASCENT_GENECOUNTS']['MIN_TOTAL_READS'],
        min_conv = config['NASCENT_GENECOUNTS']['MIN_CONVERSIONS']
    threads: 12
    conda:
        '../../envs/downstream/r-basic.yaml'
    script:
        '../../scripts/quantification/nascent_genecounts.R'
