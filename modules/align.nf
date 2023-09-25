process align_nanopore {
    label 'align'
    publishDir(
        path: "${params.output}/alignment/${prefix}/nanopore",
        mode: 'copy',
        pattern: ""
    )

    input:
    val(prefix)
    path(ref_fasta)
    path(nanopore)
    
    output:
    path("nanopore.sortedByCoord.out.bam"), emit: bam
    path("nanopore.sortedByCoord.out.bam.bai"), emit: bai

    script:
    """
        echo "${nanopore}"
        minimap2 -ax map-ont -uf -t ${task.cpus} ${ref_fasta} ${nanopore} | samtools sort -@ ${task.cpus} --verbosity 3 -o nanopore.sortedByCoord.out.bam
        samtools index -@ ${task.cpus} nanopore.sortedByCoord.out.bam
    """
}

process align_illumina {
    label 'align'
    publishDir(
        path: "${params.output}/alignment/${prefix}/illumina",
        mode: 'copy',
        pattern: ""
    )

    input:
    val(prefix)
    path(ref_fasta)
    tuple path(illumina_r1), path(illumina_r2)
    
    output:
    path("illumina.sortedByCoord.out.bam"), emit: bam
    path("illumina.sortedByCoord.out.bam.bai"), emit: bai

    script:
    """
        echo "read 1: ${illumina_r1}"
        echo "read 2: ${illumina_r2}"
        minimap2 -ax sr -t ${task.cpus} ${ref_fasta} ${illumina_r1} ${illumina_r2} | samtools sort -@ ${task.cpus} --verbosity 3 -o illumina.sortedByCoord.out.bam
        samtools index -@ ${task.cpus} illumina.sortedByCoord.out.bam
    """
}

process align_pacbio {
    label 'align'
    publishDir(
        path: "${params.output}/alignment/${prefix}/pacbio",
        mode: 'copy',
        pattern: ""
    )

    input:
    val(prefix)
    path(ref_fasta)
    path(pacbio)
    
    output:
    path("pacbio.sortedByCoord.out.bam"), emit: bam
    path("pacbio.sortedByCoord.out.bam.bai"), emit: bai

    script:
    """
        echo "${pacbio}"
        minimap2 -ax asm20 -t ${task.cpus} ${ref_fasta} ${pacbio} | samtools sort -@ ${task.cpus} --verbosity 3 -o pacbio.sortedByCoord.out.bam
        samtools index -@ ${task.cpus} pacbio.sortedByCoord.out.bam
    """
}

process salmon_index {
    label 'align'
    publishDir(
        path: "${params.output}/alignment/${prefix}/salmon",
        mode: 'copy',
        pattern: ""
    )

    input:
    val(prefix)
    path(ref_fasta)
    
    output:
    path("salmon_index"), emit: index

    script:
    """
    echo $prefix
    
    salmon index -t ${ref_fasta} -i salmon_index -k 31
    """
}

process salmon_illumina {
    label 'align'
    publishDir(
        path: "${params.output}/alignment/${prefix}/salmon",
        mode: 'copy',
        pattern: ""
    )

    input:
    val(prefix)
    tuple path(index), val(sample), path(R1), path(R2)
    
    output:
    tuple path("${sample}"), emit: salmon_quant

    script:
    """
    salmon quant -i ${index} \
        -l a \
        -1 ${R1} -2 ${R2} \
        --writeUnmappedNames \
        -p ${task.cpus} \
        -o ${sample} \
        --writeOrphanLinks --dumpEq > /dev/null 2>&1

    gunzip ${sample}/aux_info/eq_classes.txt.gz
    """
}
