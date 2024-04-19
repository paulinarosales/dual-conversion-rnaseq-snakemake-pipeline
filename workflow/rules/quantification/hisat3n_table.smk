# --------- Parsing functions ----------
def _input_refGenome(wildcards):
    return expand('resources/external/gencode_{release}/{genome}.genome.fa', release=GENCODE_RELEASE, genome=GENOME) # de-compressed


# --------- Rules ----------
rule hisat3n_table:
    input:
        refGenome = _input_refGenome,
        sortedBAM = 'results/sam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.hisat3n_align.sort.bam'
    output:
        'results/conversion_tables/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.hisat3n_tab.tsv'
    log:
        'logs/hisat3n/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_table.log'
    conda:
        '../../envs/alignment/hisat3n.yaml'
    threads: 32
    resources:
        mem = '64G'
    params:
        command = f'{HISAT3N_BIN_PATH}/hisat-3n-table', # run HISAT-3N script from git clone folder
        base_change = config['HISAT3N']['BASE_CHANGE'],
        extra = config['HISAT3N']['TABLE_EXTRA']
    shell:
        """
            samtools view -h {input.sortedBAM} |\
            {params.command} -p {threads}\
            --alignments - --ref {input.refGenome}\
            --base-change {params.base_change}\
            --output-name {output} {params.extra} 2> {log}
        """
