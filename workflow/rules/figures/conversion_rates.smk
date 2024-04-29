rule conv_rates_plot:
    input:
        'results/conversion_tables/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_conversionRates.tsv'
    output:
        'results/figures/conversions/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_conversionRates.pdf'
    params:
        base_change = config['HISAT3N']['BASE_CHANGE']
    log:
        'logs/figures/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_conversionRates.log'
    conda:
        '../../envs/downstream/r-basic.yaml'
    script:
        '../../scripts/figures/conversion_rates_boxplot.R'