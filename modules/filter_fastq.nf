process merge_fastq {
    label 'filter'
    publishDir(
        path: "${params.output}/filtered_fastq",
        mode: 'copy',
        pattern: ""
    )

    input:
    val(platform)
    val(orientation)
    path(fastq)

    output:
    path('*.merged.fastq.gz'), emit: fastq
    path('*.merged.sample.log'), emit: log

    script:
    """
    ls *.fastq.gz > ${platform}.${orientation}.merged.sample.log
    cat *.fastq.gz > ${platform}.${orientation}.merged.fastq.gz
    """
}

process filter_longread_fastq {
    label 'filter'
    publishDir(
        path: "${params.output}/filtered_fastq",
        mode: 'copy',
        pattern: ""
    )

    input:
    path(fastq)

    output:
    path('*.filtered.fastq'), emit: filtered_fastq

    script:
    """
    zcat ${fastq} | NanoFilt -q 10 -l 150 --headcrop 75 --tailcrop 75 > ${fastq.baseName}.filtered.fastq
    """
}

process filter_shortread_fastq {
    label 'filter'
    publishDir(
        path: "${params.output}/filtered_fastq",
        mode: 'copy',
        pattern: ""
    )

    input:
    tuple val(sample), path(r1_fastq), path(r2_fastq)

    output:
    tuple val(sample), path('*/*.filtered_R1.fastq.gz'), path('*/*.filtered_R2.fastq.gz'), emit: filtered_r1_r2
    tuple val(sample), path('*/*.unpaired_R1.fastq.gz'), path('*/*.unpaired_R2.fastq.gz'), emit: unpaired_r1_r2

    script:
    """
    mkdir -p ${sample}

    java -jar /tools/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 \
    ${r1_fastq} ${r2_fastq} \
    ${sample}/${sample}'.filtered_R1.fastq.gz' \
    ${sample}/${sample}'.unpaired_R1.fastq.gz' \
    ${sample}/${sample}'.filtered_R2.fastq.gz' \
    ${sample}/${sample}'.unpaired_R2.fastq.gz' \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}
