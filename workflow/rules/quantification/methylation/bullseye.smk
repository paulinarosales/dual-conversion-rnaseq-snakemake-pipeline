# --------- Parsing functions ----------
def _input_refGenome(wildcards):
    return expand('resources/external/gencode_{release}/{genome}.genome.fa', release=GENCODE_RELEASE, genome=GENOME) # de-compressed

def _input_refFlat(wildcards):
    return expand('resources/external/gencode_{release}/{genome}.annotation.refFlat', release=GENCODE_RELEASE, genome=GENOME) # de-compressed

def _input_geneBED(wildcards):
    return expand('resources/external/gencode_{release}/{genome}.gene.bed', release=GENCODE_RELEASE, genome=GENOME)


rule bullseye_parseBAM:
    input:
        'results/sam_files/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}.hisat3n_align.sort.bam'
    output:
        'results/conversion_tables/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_nucContent_byPos.matrix.gz'
    log:
        'logs/bullseye/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_parseBAM.log'
    params:
        output_basename = 'results/conversion_tables/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_nucContent_byPos.matrix',
        minCov = config['BULLSEYE']['MIN_COVERAGE']
    threads: 32
    conda:
        'bullseye'
    resources:
        mem = '124G'
    shell:
        """
             perl workflow/scripts/bullseye/parseBAM.pl --input {input} --output {params.output_basename} --minCoverage {params.minCov} --removeDuplicates --verbose &> {log}
        """
        

rule bullseye_gtf2genepred:
    input:
        'resources/external/gencode_{release}/{genome}.annotation.gtf.gz'
    output:
        'resources/external/gencode_{release}/{genome}.annotation.refFlat'
    log:
        'logs/bullseye/gencode_{release}_{genome}_gtf2genepred.log'
    conda:
        'bullseye'
    shell:
        """
             perl workflow/scripts/bullseye/gtf2genepred.pl --gtf {input} --out {output}
        """
        
                
rule bullseye_findSite:
# PENDING CHANGES: make the control mtx conditional
    input:
        editMTX = 'results/conversion_tables/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_nucContent_byPos.matrix.gz',
        refGenome = _input_refGenome,
        refFlatGTF = _input_refFlat
    output:
        'results/site_calling/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_m6A_sites.bed'
    log:
        'logs/bullseye/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_find_edit_site.log'
    params:
 	    minEdit =  config['BULLSEYE']['MIN_EDIT_RATE'],  # minimal editing rates
 	    maxEdit  = config['BULLSEYE']['MAX_EDIT_RATE'],  # maximal editing rates
 	    editFoldThreshold =  config['BULLSEYE']['MIN_EDIT_FOLD_RATIO'],  # minimal editing ration over control sample
 	    minEditSites =  config['BULLSEYE']['MIN_EDIT_SITES']  # minimal number of mutations for detection of site
    threads: 32
    conda:
        'bullseye'
    resources:
        mem = '124G'
    shell:
        """
            perl workflow/scripts/bullseye/Find_edit_site.pl\
                --annotationFile {input.refFlatGTF} \
 	            --EditedMatrix {input.editMTX} \
 	            -g {input.refGenome} \
 	            --minEdit {params.minEdit} \
 	            --maxEdit {params.maxEdit} \
 	            --editFoldThreshold {params.editFoldThreshold} \
 	            --MinEditSites {params.minEditSites} \
 	            --cpu {threads} \
 	            --outfile {output} \
 	            --verbose &> {log}
        """

rule bullseye_RACfilter:
    input:
        sitesBED = 'results/site_calling/{sample_type}_wt-vs-mut_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_wt-vs-mut_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_m6A_sites.bed',
        refGenome = _input_refGenome
    output:
        'results/site_calling/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_m6A_sites.RAC.bed'
    log:
        'logs/bullseye/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_RAC_filter.log'
    params:
        window_size = config['RAC_FILTER']['WINDOW_SIZE'], # default= 0, help="How big the window should be. Use 0 for regular strict RAC filter.")
        adapt_coords = config['RAC_FILTER']['ADAPT_COORDS'] # default=True, help="Whether we should switch coordinates to A in RAC motif")
    threads: 24
    conda:
        'bullseye'
    resources:
        mem = '24G'
    shell:
        """
            python workflow/scripts/bullseye/window_rac_filter.py \
                --bed_file {input.sitesBED} --fasta_file {input.refGenome} \
                --window_size {params.window_size} --adapt_coordinates {params.adapt_coords} > {output}
        """

rule bullseye_parse_cts:
    input:
        'results/site_calling/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_m6A_sites.RAC.bed'
    output:
        parsedTSV = 'results/site_calling/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_m6A_sites.counts.tsv',
        geneListTSV = 'results/site_calling/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_m6A_sites.geneList.tsv'
    log:
        'logs/bullseye/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_parse_sites_cts.log'
    threads: 12
    conda:
        '../../../envs/downstream/r-basic.yaml'
    script:
        '../../../scripts/bullseye/parse_sites_cts.R'


rule bullseye_quantifySites:
# PENDING CHANGES: make the control mtx conditional
    input:
        editMTX = 'results/conversion_tables/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_nucContent_byPos.matrix.gz',
        genesBED = _input_geneBED
    output:
        'results/site_calling/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_m6A_sites.gene_counts.bed'
    log:
        'logs/bullseye/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_quant_edit_site.log'
    params:
 	    minEdit =  config['BULLSEYE']['MIN_EDIT_RATE'],  # minimal editing rates
 	    maxEdit  = config['BULLSEYE']['MAX_EDIT_RATE'],  # maximal editing rates
 	    editedMinCov =  config['BULLSEYE']['MIN_COVERAGE'],  # 
 	    minEditSites =  config['BULLSEYE']['MIN_EDIT_SITES']  # minimal number of mutations for detection of site
    threads: 32
    conda:
        'bullseye'
    resources:
        mem = '124G'
    shell:
        """
            perl workflow/scripts/bullseye/quantify_sites.pl\
                --bed {input.genesBED} \
 	            --EditedMatrix {input.editMTX} \
 	            --editType C2T \
 	            --minEdit {params.minEdit} \
 	            --maxEdit {params.maxEdit} \
 	            --EditedMinCoverage {params.editedMinCov} \
 	            --MinEditSites {params.minEditSites} \
 	            --cpu {threads} \
 	            --outfile {output} \
 	            --verbose &> {log}
        """