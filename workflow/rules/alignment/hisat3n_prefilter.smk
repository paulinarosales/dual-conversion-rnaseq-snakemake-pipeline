# --------- Parsing functions ----------
def _input_for_hisat3n_align(wildcards):
    # make sure that hisat3n_build is executed first
    # INDEX_DIR = os.path.dirname(checkpoints.hisat3n_build.get(release=GENCODE_RELEASE, genome=GENOME).output[0])
    INDEX_DIR = os.path.dirname(f'resources/external/index/hisat3n/mouse_rDNA_BK000964-3/BK000964-3.rDNA')
    INDEX_FILES = []
    # print(INDEX_DIR)
    for root, directories, files in os.walk(INDEX_DIR):
        for file in files:
            # list files under hisat3n_build output directory
            INDEX_FILES.append(file)
    return expand('{index_dir}/{index_file}', index_dir=INDEX_DIR, index_file=INDEX_FILES)


def _index_basename_for_hisat3n(wildcards):
    return  'resources/external/index/hisat3n/mouse_rDNA_BK000964-3/BK000964-3.rDNA'

def _params_for_hisat3n_align(wildcards):
    # add --repeat if requested
    _repeat_index = config['HISAT3N']['REPEAT_INDEX']
    if _repeat_index == 'Yes':
        repeat_index = '--repeat'

    elif _repeat_index == 'No':
        repeat_index = ''
    return repeat_index

# --------- Rules ----------
# evaluation after run does not work, outdir doesn't match
rule get_filter_index:
    output:
        # directory('resources/external/index/hisat3n/gencode_{release}/{genome}.genome')
        directory('resources/external/index/hisat3n/mouse_rDNA_BK000964-3')
    log:
        'logs/hisat3n/build/mouse_rDNA_BK000964-3_copy_index.log'
    params:
        source_dir = config['HISAT3N']['HISAT3N_FILTER_INDEX_PATH'],
        copy_to_dir = INDEX_DIR
    shell:
        """
            rsync -Pav {params.source_dir} {params.copy_to_dir} 2> {log}
        """


rule hisat3n_prefilter_align:
    input:
        index = _input_for_hisat3n_align,
        fq1 = 'results/fastq/trimmed/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_val_1.fq.gz',
        fq2 = 'results/fastq/trimmed/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_val_2.fq.gz'
    output:
        alignSAM = 'results/sam_files/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}.hisat3n_align.prefilter.sam',
        unalignSAM_1 = 'results/sam_files/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}.hisat3n_unalign.prefilter.1.sam',
        unalignSAM_2 = 'results/sam_files/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}.hisat3n_unalign.prefilter.2.sam'
    log:
        'logs/hisat3n/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_prefilter_align.log'
    conda:
        '../../envs/alignment/hisat3n.yaml'
    threads: 32
    resources:
        mem = '300G',
        time = '1-00:00:00'
    params:
        command = f'{HISAT3N_BIN_PATH}/hisat-3n', # run HISAT-3N script from git clone folder
        index_basename = _index_basename_for_hisat3n,
        repeat = _params_for_hisat3n_align,
        base_change = config['HISAT3N']['BASE_CHANGE'],
        unalign_basename = 'results/sam_files/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}.hisat3n_unalign.prefilter.%.sam',
        extra = config['HISAT3N']['ALIGN_EXTRA']
    shell:
        """
            {params.command} -p {threads} -x {params.index_basename}\
            -q -1 {input.fq1} -2 {input.fq2}\
            -S {output.alignSAM} --base-change {params.base_change}\
            {params.repeat} {params.extra}\
            --un-conc {params.unalign_basename}\
            --new-summary --summary-file {log} 
        """


rule samtools_unaligned_sam2fq:
    input:
        'results/sam_files/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}.hisat3n_unalign.prefilter.{read}.sam'
    output:
        'results/fastq/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_R{read}.prefilter.fq'
    log:
        'logs/samtools/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_R{read}_unaligned_prefilter_sam2fq.log'
    conda:
        '../../envs/alignment/hisat3n.yaml'
    resources:
        mem= '24G'
    threads: 12
    shell:
        """
            samtools fastq -@ {threads} {input} > {output}
        """
