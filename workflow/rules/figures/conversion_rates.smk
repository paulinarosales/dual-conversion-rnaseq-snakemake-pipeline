rule conv_rates_plot:
    input:
        'results/conversion_tables/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_conversionRates.genes.tsv'
    output:
        'results/figures/conversions/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_conversionRates.pdf'
    params:
        base_change = config['HISAT3N']['BASE_CHANGE']
    log:
        'logs/figures/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_conversionRates.log'
    conda:
        '../../envs/downstream/r-basic.yaml'
    threads: 24
    resources:
        mem = '24G'
    script:
        '../../scripts/figures/conversion_rates_boxplot.R'

rule global_conv_rates_plot:
    input:
        sample_manifestTSV = config['SAMPLE_MANIFEST'],
        convFiles = TARGETS['global_conv_rates'],
        countFiles = TARGETS['featureCounts']
    output:
        all_globalRatesTSV = 'results/conversion_tables/all_samples/global_conversionRates.tsv',
        barplot_newPDF = 'results/figures/conversions/all_samples/global_conversionRates_new.pdf',
        barplot_oldPDF = 'results/figures/conversions/all_samples/global_conversionRates_old.pdf'
        # barplot_norm_newPDF = 'results/figures/conversions/all_samples/global_conversionRates_new.lib_norm.pdf',
        # barplot_norm_oldPDF = 'results/figures/conversions/all_samples/global_conversionRates_old.lib_norm.pdf'
    log:
        'logs/figures/all_samples/global_conversionRates.log'
    conda:
        '../../envs/downstream/r-basic.yaml'
    script:
        '../../scripts/figures/global_conversion_rates_barplot.R'


rule mut_dist_per_read:
    input:
        sample_manifestTSV = config['SAMPLE_MANIFEST'],
        bakR_mutTSV ='results/conversion_tables/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_conversionCounts.tsv'
    output:
        raw_distPDF = 'results/figures/conversions/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_T_C_conv_per_read_histogram.pdf',
        filter_distPDF = 'results/figures/conversions/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_T_C_conv_per_read_histogram.filter.pdf'
    params:
        target_mut = 'T_C'
    log:
        'logs/figures/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_mut_dist_histogram.log'
    conda:
        '../../envs/downstream/r-basic.yaml'
    threads: 24
    resources:
        mem = '64G'
    script:
        '../../scripts/figures/mut_distribution.R'