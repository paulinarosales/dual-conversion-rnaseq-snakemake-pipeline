def _input_tx_infoTSV(wildcards):
    return expand('resources/external/gencode_{release}/{genome}.transcripts_info.tsv', release=GENCODE_RELEASE, genome=GENOME)

rule conv_rates:
    input:
        collapsedTSV = 'results/conversion_tables/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_conversionCounts.collpasedFeatures.tsv',
        transcriptBED = _input_tx_infoTSV
    output:
        'results/conversion_tables/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_conversionRates.tsv',
    log:
        'logs/rates/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_conversionRates.log'
    conda:
        '../../envs/downstream/r-basic.yaml'
    script:
        '../../scripts/quantification/conversion_rates.R'