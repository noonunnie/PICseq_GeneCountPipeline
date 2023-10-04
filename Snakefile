# https://github.com/aminnaghdloo/inferCNV/blob/main/Snakefile code

# configfile: "config/config.yaml" # MCF7 cells
configfile: "config/config2.yaml" # SKBR3 cells
samples = list(config['rnaSamples'].values())

rule all:
    input:
        expand('reports/fastqc/{item}_1_fastqc.html', item = samples),
        expand('data/trimmed/{item}_val_1.fq.gz', item = samples),
        expand('data/sorted/{item}.bam', item = samples),
        expand('data/rnaQC/{item}_qc.txt', item = samples),
        expand('data/name_sorted/{item}.bam', item = samples),
        'data/gene_count/features_SKBR3.txt',
        'data/casper/merge_SKBR3.bam'


rule fastqc:
    input:
        r1 = config['rnaReadsDir'] + '/{sample}' + config['read1Suffix'],
        r2 = config['rnaReadsDir'] + '/{sample}' + config['read2Suffix']
    output:
        o1 = 'reports/fastqc/{sample}_1_fastqc.html',
        o2 = 'reports/fastqc/{sample}_2_fastqc.html'
    params:
        outDir = 'reports/fastqc',
    log:
        'logs/fastqc/fastqc_{sample}.log'
    shell:
        "fastqc {input.r1} {input.r2} -o {params.outDir} &> {log}"


rule trimAdapters:
    input:
        r1 = config['rnaReadsDir'] + '/{sample}' + config['read1Suffix'],
        r2 = config['rnaReadsDir'] + '/{sample}' + config['read2Suffix']
    output:
        o1 = 'data/trimmed/{sample}_val_1.fq.gz',
        o2 = 'data/trimmed/{sample}_val_2.fq.gz'
    params:
        minPhred = config['minPhred'],
        minOverlap = config['minOverlap'],
        suffix = '{sample}',
        outDir = 'data/trimmed'
    log:
        'logs/trimAdapters/trimAdapters_{sample}.log'
    shell:
        'trim_galore --paired --quality {params.minPhred} '
        '--stringency {params.minOverlap} --basename {params.suffix} '
        '--output_dir {params.outDir} {input.r1} {input.r2} &> {log}'


rule alignReads:
    input: 
        r1 = 'data/trimmed/{sample}_val_1.fq.gz',
        r2 = 'data/trimmed/{sample}_val_2.fq.gz'
    output:
        sam = temp('data/sorted/{sample}.sam'),
        bam = 'data/sorted/{sample}.bam'
    threads:
        config['nThreads']
    log:
        'logs/hisat2/hisat2_{sample}.log'
    params:
        refIndex = config['refIndex'],
        junction = 'data/junctions/{sample}_hisat2.junction'
    shell:
        'hisat2 --novel-splicesite-outfile {params.junction} '
        '-p {threads} -x {params.refIndex} '
        '-1 {input.r1} -2 {input.r2} -S {output.sam} &> {log}; '
        'samtools sort -@ {threads} -o {output.bam} {output.sam}'


rule rnaQC:
    input:
        'data/sorted/{sample}.bam'
    output:
        bam = temp('data/rnaQC/{sample}_temp.bam'),
        qc = 'data/rnaQC/{sample}_qc.txt'
    params:
        refFlat = config['refFlat'],
        rrna = config['rrna']
    shell:
        'samtools reheader -c \'grep -v ^@PG\' {input} > {output.bam}; '
        'picard CollectRnaSeqMetrics INPUT={output.bam} OUTPUT={output.qc} '
        'REF_FLAT={params.refFlat} RIBOSOMAL_INTERVALS={params.rrna} '
        'STRAND=NONE'


rule sortName:
    input:
        'data/sorted/{sample}.bam'
    output:
        temp('data/name_sorted/{sample}.bam')
    threads:
        config['nThreads']
    log:
        'logs/sortName/sortName_{sample}.log'
    shell:
        'samtools sort -n -@ {threads} -o {output} {input} &> {log}'


rule countGenes:
    input:
        expand('data/name_sorted/{item}.bam', item = samples)
    output:
        'data/gene_count/features_SKBR3.txt'
    params:
        annotation = config['annotation']
    shell:
        'featureCounts -a {params.annotation} -o {output} -t exon -g gene_id '
        '-p {input}'


rule mergeBam:
    input:
        bam = expand('data/sorted/{item}.bam', item = samples),
        bai = expand('data/sorted/{item}.bam', item = samples)
    output:
        'data/casper/merge_SKBR3.bam'
    threads:
        config['nThreads']
    shell:
        'samtools merge -@ {threads} -f -X {output} {input.bam} {input.bai}'


rule extractBAF:
    input:
        'data/casper/merge_SKBR3.bam'
        #'data/sorted/ERR1898567.bam'
    output:
        allele = directory('/scratch2/naghdloo/merge_pileup'),
        snp = 'data/casper/merge_SKBR3.snp'
    params:
        pileupRef = 'data/ref/GRCh38/BAF/genome_pileup',
        list = 'data/ref/GRCh38/BAF/genome.list'
    log:
        step1 = 'logs/casper/extractBAF_step1.log',
        step2 = 'logs/casper/extractBAF_step2.log'
    shell:
        'samtools view {input} | '
        './scripts/BAFExtract -generate_compressed_pileup_per_SAM '
        'stdin {params.list} {output.allele} 50 0; '
        './scripts/BAFExtract -get_SNVs_per_pileup {params.list} '
        '{output.allele} {params.pileupRef} 20 4 0.1 {output.snp} '
