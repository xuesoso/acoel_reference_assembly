process rattle_assemble {
    label 'assembly'
    publishDir(
        path: "${params.output}/assembly",
        mode: 'copy',
        pattern: ""
    )

    input:
    path(nanopore)
    
    output:
    path("*/clusters.out")

    script:
    """
        echo "${nanopore}"
        mkdir -p rattle.clusters
        rattle cluster -i ${nanopore} -t ${task.cpus} -o rattle.clusters/ --iso --rna --fastq
    """
}

process rattle_correct {
    label 'assembly'
    publishDir(
        path: "${params.output}/assembly",
        mode: 'copy',
        pattern: ""
    )

    input:
    path(nanopore)
    path(clusters)
    
    output:
    path("*/consensi.fq")

    script:
    """
        echo "${nanopore}"
        echo "${clusters}"
        mkdir -p rattle.corrected
        rattle correct -i ${nanopore} -c ${clusters} -t ${task.cpus} -o rattle.corrected/
    """
}

