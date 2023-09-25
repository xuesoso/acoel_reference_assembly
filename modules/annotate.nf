process transdecoder {
    label 'annotate'
    publishDir(
        path: "${params.output}/annotate/${prefix}/transdecoder",
        mode: 'copy',
        pattern: "*.{log,pep}"
    )

    input:
    val(prefix)
    path(ref_fasta)

    output:
    path("*.only_longest_ORF.pep"), emit: longest_ORF_pep
    path("*"), emit: others
    
    script:
    """
        echo $ref_fasta
        echo $prefix

        TransDecoder.LongOrfs -t ${ref_fasta} &> transdecoder.longorfs.log
        TransDecoder.Predict --no_refine_starts -t ${ref_fasta} &> transdecoder.predict.log

        /usr/bin/util/get_longest_ORF_per_transcript.pl ${ref_fasta}'.transdecoder.pep' > ${ref_fasta}'.transdecoder.only_longest_ORF.pep'
    """
}

process salmon_grouper {
    container 'combinelab/grouper'
    label 'annotate'
    publishDir(
        path: "${params.output}/annotate/${prefix}/salmon_grouper",
        mode: 'copy',
        pattern: "*.clust"
    )

    input:
    val(prefix)
    path(config_yaml)
    path(salmon_quant)

    output:
    path("mag.flat.clust"), emit: mag_flat_cluster
    path("mag.clust"), emit: mag_cluster
    
    script:
    """
        echo ${config_yaml}

        Grouper --config ${config_yaml}
    """
}
