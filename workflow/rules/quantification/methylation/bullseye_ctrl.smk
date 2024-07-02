# --------- Parsing functions ----------
def _input_refGenome(wildcards):
    return expand('resources/external/gencode_{release}/{genome}.genome.fa', release=GENCODE_RELEASE, genome=GENOME) # de-compressed

def _input_refFlat(wildcards):
    return expand('resources/external/gencode_{release}/{genome}.annotation.refFlat', release=GENCODE_RELEASE, genome=GENOME) # de-compressed

     


rule bullseye_findSite_ctrl:
# PENDING CHANGES: make the control mtx conditional
    input:
        editMTX = 'results/site_calling/{sample_type}_APOBEC1-YTH-wt_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_APOBEC1-YTH-wt_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_nucContent_byPos.matrix.gz',
        controlMTX = 'results/site_calling/{sample_type}_APOBEC1-YTH-mut_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_APOBEC1-YTH-mut_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_nucContent_byPos.matrix.gz',
        edit_c2t_snpsBED = 'results/snps/{sample_type}_APOBEC1-YTH-wt_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_APOBEC1-YTH-wt_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}.snp.c2t.bed',
        control_c2t_snpsBED = 'results/snps/{sample_type}_APOBEC1-YTH-mut_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_APOBEC1-YTH-mut_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}.snp.c2t.bed',
        refGenome = _input_refGenome,
        refFlatGTF = _input_refFlat
    output:
        'results/site_calling/against_control/{sample_type}_wt-vs-mut_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_wt-vs-mut_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_m6A_sites.bed'
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
                --filterBed {input.edit_c2t_snpsBED} \
                --filterBed {input.control_c2t_snpsBED} \
 	            --minEdit {params.minEdit} \
 	            --maxEdit {params.maxEdit} \
 	            --editFoldThreshold {params.editFoldThreshold} \
 	            --MinEditSites {params.minEditSites} \
 	            --cpu {threads} \
 	            --outfile {output} \
 	            --fallback {input.refGenome} \
 	            --verbose &> {log}
        """

rule bullseye_RACfilter_ctrl:
    input:
        sitesBED = 'results/site_calling/against_control/{sample_type}_wt-vs-mut_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_wt-vs-mut_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_m6A_sites.bed',
        refGenome = _input_refGenome
    output:
        'results/site_calling/against_control/{sample_type}_wt-vs-mut_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_wt-vs-mut_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_m6A_sites.RAC.bed'
    log:
        'logs/bullseye/{sample_type}_wt-vs-mut_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_wt-vs-mut_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_RAC_filter.log'
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

rule bullseye_parse_cts_ctrl:
    input:
        'results/site_calling/against_control/{sample_type}_wt-vs-mut_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_wt-vs-mut_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_m6A_sites.RAC.bed'
    output:
        parsedTSV = 'results/site_calling/against_control/{sample_type}_wt-vs-mut_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_wt-vs-mut_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_m6A_sites.counts.tsv',
        geneListTSV = 'results/site_calling/against_control/{sample_type}_wt-vs-mut_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_wt-vs-mut_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_m6A_sites.geneList.tsv',
        summaryTSV = 'results/site_calling/against_control/{sample_type}_wt-vs-mut_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_wt-vs-mut_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_m6A_sites.summary.tsv'
    log:
        'logs/bullseye/{sample_type}_wt-vs-mut_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_wt-vs-mut_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_parse_sites_cts.log'
    threads: 12
    conda:
        '../../../envs/downstream/r-basic.yaml'
    script:
        '../../../scripts/bullseye/parse_sites_cts.R'
