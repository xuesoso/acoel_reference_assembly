process rattle_polish {
    label 'polish'
    publishDir(
        path: "${params.output}/assembly",
        mode: 'copy',
        pattern: ""
    )

    input:
    path(consensi)
    
    output:
    path("*/transcriptome.fasta"), emit: fasta
    path("*/transcriptome.fq"), emit: fastq

    script:
    """
        echo "${consensi}"
        mkdir -p rattle.polished
        rattle polish -i ${consensi} -t ${task.cpus} -o rattle.polished/
        seqtk seq -a rattle.polished/transcriptome.fq > rattle.polished/transcriptome.fasta
    """
}

process pilon {
    label 'polish'
    publishDir(
        path: "${params.output}/reference/${prefix}",
        mode: 'copy',
        pattern: ""
    )

    input:
    val(prefix)
    path(pacbio_bam)
    path(pacbio_bai)
    path(illumina_bam)
    path(illumina_bai)
    path(transcriptome_fasta)
    
    output:
    path("*.fasta"), emit: fasta
    path("*.changes"), emit: changes

    script:
    def mem = "${task.memory}".replaceAll(/[^0-9]/, '')

    """
        echo "${transcriptome_fasta}"
        java -Xmx${mem}G -jar /usr/bin/pilon-1.23.jar \
        --genome ${transcriptome_fasta} --threads ${task.cpus} \
        --pacbio ${pacbio_bam} --frags ${illumina_bam} --output pilon_corrected.reference \
        --verbose --debug --changes --fix all
    """
}

process trim_reference {
    label 'polish'
    publishDir(
        path: "${params.output}/reference/${prefix}",
        mode: 'copy',
        pattern: ""
    )

    input:
    val(prefix)
    path(ref_fasta)
    path(bam)
    path(bai)
    
    output:
    path("*.trimmed.fasta"), emit: fasta

    script:
    """
        python ${workflow.projectDir}/bin/trim_fa.py \
                ${ref_fasta} \
                ${bam} \
                'pilon_corrected.reference.trimmed.fasta' \
                ${task.cpus}
    """
}

process cluster_transcripts {
    // After Salmon clustering, we identify transcripts that belong to the same
    // cluster and pick one with the longest predicted codon length
    label 'filter'
    publishDir(
        path: "${params.output}/reference/${prefix}",
        mode: 'copy',
        pattern: ""
    )

    input:
    val(prefix)
    path(input_ref_fasta)
    path(input_pep)
    path(salmon_cluster)
    
    output:
    path("*.grouped.fasta"), emit: fasta

    script:
    """
        echo ${input_ref_fasta}
        echo ${input_pep}
        echo ${salmon_cluster}

        python ${workflow.projectDir}/bin/cluster_transcripts.py \
                ${input_ref_fasta} \
                ${input_pep} \
                ${salmon_cluster} \
                'pilon_corrected.reference.trimmed.grouped.fasta'
    """
}

