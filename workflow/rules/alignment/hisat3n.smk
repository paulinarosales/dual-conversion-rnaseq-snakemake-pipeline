# --------- Parsing functions ----------
# -- Build --
def _input_refGenome(wildcards):
    return expand('resources/external/gencode_{release}/{genome}.genome.fa', release=GENCODE_RELEASE, genome=GENOME) # de-compressed


def _params_for_hisat3n(wildcards):
    # add --repeat-index if requested
    _repeat_index = config['HISAT3N']['REPEAT_INDEX']
    if _repeat_index == 'Yes':
        repeat_index = '--repeat-index'
    elif _repeat_index == 'No':
        repeat_index = ''

    # add options for low memory consumption if requested
    _resources_mode = config['HISAT3N']['RESOURCES_MODE']
    if _resources_mode == 'Default':
        resources = ''

    elif _resources_mode == 'Low':
        resources = '--noauto --bmaxdivn 8 --dcv 512'

    return f'{repeat_index} {resources}'

# -- Align --

def _input_for_hisat3n_align(wildcards):
    # make sure that hisat3n_build is executed first
    # INDEX_DIR = os.path.dirname(checkpoints.hisat3n_build.get(release=GENCODE_RELEASE, genome=GENOME).output[0])
    INDEX_DIR = os.path.dirname(f'resources/external/index/hisat3n/gencode_{GENCODE_RELEASE}/{GENOME}.genome')
    INDEX_FILES = []
    # print(INDEX_DIR)
    for root, directories, files in os.walk(INDEX_DIR):
        for file in files:
            # list files under hisat3n_build output directory
            INDEX_FILES.append(file)
    return expand('{index_dir}/{index_file}', index_dir=INDEX_DIR, index_file=INDEX_FILES)


def _index_basename_for_hisat3n(wildcards):
    return expand('resources/external/index/hisat3n/gencode_{release}/{genome}.genome', release=GENCODE_RELEASE, genome=GENOME)

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
checkpoint hisat3n_build:
    input:
       'resources/external/gencode_{release}/{genome}.genome.fa'
    output:
        # directory('resources/external/index/hisat3n/gencode_{release}/{genome}.genome')
        directory('resources/external/index/hisat3n/gencode_{release}_{genome}/genome')
    log:
        'logs/hisat3n/build/gencode_{release}_{genome}_index.log'
    params:
        command = f'{HISAT3N_BIN_PATH}/hisat-3n-build', # run HISAT-3N script from git clone folder
        # index_basename = 'resources/external/index/hisat3n/gencode_{release}/{genome}.genome',
        base_change = config['HISAT3N']['BASE_CHANGE'],
        extra = _params_for_hisat3n
    conda:
        '/home/ife/paulina.rosales/miniconda3/envs/hisat-3n'
    threads: 32
    resources:
        mem = '169G',
        time = '3-00:00:00'
    shell:
        """
            {params.command} --base-change {params.base_change}\
            {params.extra} -p {threads}\
            {input} {output} > {log} 2>&1 
        """


rule hisat3n_align:
    input:
        index = _input_for_hisat3n_align,
        fq1 = 'results/fastq/trimmed/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_val_1.fq.gz',
        fq2 = 'results/fastq/trimmed/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_val_2.fq.gz'
    output:
        'results/sam_files/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}.hisat3n_align.sam'
    log:
        'logs/hisat3n/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_align.log'
    conda:
        '../../envs/alignment/hisat3n.yaml'
    threads: 32
    resources:
        mem = '140G',
        time = '1-00:00:00'
    params:
        command = f'{HISAT3N_BIN_PATH}/hisat-3n', # run HISAT-3N script from git clone folder
        index_basename = _index_basename_for_hisat3n,
        repeat = _params_for_hisat3n_align,
        base_change = config['HISAT3N']['BASE_CHANGE'],
        extra = config['HISAT3N']['ALIGN_EXTRA']
    shell:
        """
            {params.command} -p {threads} -x {params.index_basename}\
            -q -1 {input.fq1} -2 {input.fq2}\
            -S {output} --base-change {params.base_change}\
            {params.repeat} {params.extra}\
            --new-summary --summary-file {log} 
        """