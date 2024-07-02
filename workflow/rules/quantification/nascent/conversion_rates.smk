def _input_geneset(wildcards):
    return expand('resources/external/gencode_{release}/{genome}.geneset.tsv', release=GENCODE_RELEASE, genome=GENOME)

rule conv_rates:
    input:
        collapsedTSV = 'results/conversion_tables/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_conversionCounts.genes.tsv',
        genesetTSV = _input_geneset
    output:
        globalTSV = 'results/conversion_tables/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_conversionRates.global.tsv',
        geneTSV = 'results/conversion_tables/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_conversionRates.genes.tsv'
    log:
        'logs/rates/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_conversionRates.log'
    conda:
        '../../../envs/downstream/r-basic.yaml'
    script:
        '../../../scripts/quantification/conversion_rates.R'