def _input_refGTF(wildcards):
    return expand('resources/external/gencode_{release}/{genome}.annotation.gtf', release=GENCODE_RELEASE, genome=GENOME) # de-compressed


rule featurecounts_genes:
    input:
        filterBAM = 'results/sam_files/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}.hisat3n_align.sortn.filter.bam',
        genome_gtf = _input_refGTF
    output:
        report = 'results/counts/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}.hisat3n_align.sortn.filter.bam.featureCounts',
        matrix = 'results/counts/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_counts.genes.tsv',
        summary = 'results/counts/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_counts.genes.tsv.summary'
    log: 
        'logs/featureCounts/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_counts.genes.log'
    conda:
        '../../../envs/alignment/hisat3n.yaml'
    threads: 24
    resources:
        mem = '64G'
    params:
        # count_mode = config['FEATURE_COUNTS']['COUNT_MODE'],
        ftr_type = config['FEATURE_COUNTS']['FEATURE_TYPE'],
        strandness = config['FEATURE_COUNTS']['STRANDNESS'],
        extra =config['FEATURE_COUNTS']['EXTRA']
    shell:
        """
            featureCounts -T {threads} -p --countReadPairs -R CORE\
            -t {params.ftr_type} -g gene_id -s {params.strandness} {params.extra}\
            -a {input.genome_gtf} -o {output.matrix} {input.filterBAM} 2> {log} 
        """