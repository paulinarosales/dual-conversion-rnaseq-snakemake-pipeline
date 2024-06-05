
def _input_geneset(wildcards):
    return expand('resources/external/gencode_{release}/{genome}.geneset.tsv', release=GENCODE_RELEASE, genome=GENOME)


rule m6a_overlap_plot:
    input:
        cts_sitesTSV = 'results/site_calling/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_m6A_sites.geneList.tsv',
        ref_geneListTSV = 'resources/ref_sites/High_confidence_genes_plps.tsv',
        ref_ctsTSV = 'resources/ref_sites/A3_plps_sites.tsv',
        genesetTSV = _input_geneset
    output:
       vennPDF = 'results/figures/m6a_sites/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_venn.pdf',
        pairsPDF = 'results/figures/m6a_sites/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_correlation.pdf'
    log:
        'logs/figures/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Chase-time_{chase_time_h}_Bio-rep_{bio_rep}_m6a_ref_overlap.log'
    conda:
        '../../envs/downstream/r-plotter.yaml'
    threads: 12
    resources:
        mem = '12G'
    script:
        '../../scripts/figures/m6a_ref_overlap.R'