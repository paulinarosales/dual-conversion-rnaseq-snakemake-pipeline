# --------- Parsing functions ----------
def _input_refGenome(wildcards):
    return expand('resources/external/gencode_{release}/{genome}.genome.fa', release=GENCODE_RELEASE, genome=GENOME) # de-compressed

def _input_refFlat(wildcards):
    return expand('resources/external/gencode_{release}/{genome}.annotation.refFlat', release=GENCODE_RELEASE, genome=GENOME) # de-compressed



rule bullseye_parseBAM:
    input:
        'results/sam_files/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}.hisat3n_align.sort.bam'
    output:
        'results/site_calling/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_nucContent_byPos.matrix.gz'
    log:
        'logs/bullseye/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_parseBAM.log'
    params:
        output_basename = 'results/site_calling/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_nucContent_byPos.matrix',
        minCov = config['BULLSEYE']['MIN_COVERAGE']
    threads: 32
    conda:
        'bullseye'
    resources:
        mem = '124G'
    shell:
        """
             perl workflow/scripts/bullseye/parseBAM.pl --input {input} --output {params.output_basename} --minCoverage {params.minCov} --removeDuplicates --verbose
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
    input:
        editMTX = 'results/site_calling/{sample_type}_APOBEC1-YTH-wt_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_APOBEC1-YTH-wt_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_nucContent_byPos.matrix.gz',
        controlMTX = 'results/site_calling/{sample_type}_APOBEC1-YTH-mut_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_APOBEC1-YTH-mut_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_nucContent_byPos.matrix.gz',
        refGenome = _input_refGenome,
        refFlatGTF = _input_refFlat
    output:
        'results/site_calling/{sample_type}_wt-vs-mut_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_wt-vs-mut_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_m6A_sites.bed'
    log:
        'logs/bullseye/{sample_type}_wt-vs-mut_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_wt-vs-mut_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_find_edit_site.log'
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
 	            --controlMatrix  {input.controlMTX} \
 	            --minEdit {params.minEdit} \
 	            --maxEdit {params.maxEdit} \
 	            --editFoldThreshold {params.editFoldThreshold} \
 	            --MinEditSites {params.minEditSites} \
 	            --cpu {threads} \
 	            --outfile {output} \
 	            --fallback {input.refGenome} \
 	            --verbose 2>&1 {log}
        """
        
