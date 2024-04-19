rule fastqc:
    input:
        'resources/fastq_seq/merged/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_R{read}.fq.gz'
    output:
        'results/quality_control/raw_fastq/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_R{read}_fastqc.html',
    params:
        outdir = 'results/quality_control/raw_fastq/{sample_type}_{treatment}_Bio-rep_{bio_rep}',
    conda:
        '../../envs/pre_processing/trim-galore.yaml'
    shell:
        """
            mkdir -p {params.outdir} && \
            fastqc {input} --outdir={params.outdir}
        """