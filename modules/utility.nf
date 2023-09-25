process rev_comp {
    label 'filter'
    publishDir(
        path: "${params.output}/reference/${prefix}",
        mode: 'copy',
        pattern: ""
    )

    input:
    val(prefix)
    path(transcriptome_fasta)
    
    output:
    path("*.reverse_complement.fasta"), emit: fasta

    script:
    """
        echo "${transcriptome_fasta}"
        seqtk seq -r ${transcriptome_fasta} > ${transcriptome_fasta.baseName}.reverse_complement.fasta
    """
}

process busco_score {
    label 'annotate'
    publishDir(
        path: "${params.output}/reference/${prefix}",
        mode: 'copy',
        pattern: ""
    )

    input:
    val(prefix)
    path(transcriptome_fasta)
    
    script:
    """
        echo "${transcriptome_fasta}"
    """
}

