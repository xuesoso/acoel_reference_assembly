// Set DSL2 syntax
nextflow.enable.dsl = 2

// Import the required modules
include { 
            filter_shortread_fastq as filter_illumina_fastq; \
            filter_longread_fastq as filter_nanopore_fastq; \
            filter_longread_fastq as filter_pacbio_fastq; \
            merge_fastq as merge_illumina_fastq_r1;
            merge_fastq as merge_illumina_fastq_r2;
            } from './modules/filter_fastq.nf'
include { rattle_assemble; rattle_correct } from './modules/rattle.nf'
include { 
            align_nanopore as align_nanopore_round1; \
            align_illumina as align_illumina_round1; \
            align_pacbio as align_pacbio_round1;
            align_nanopore as align_nanopore_round2; \
            align_illumina as align_illumina_round2; \
            align_pacbio as align_pacbio_round2; \
            } from './modules/align.nf'
include { align_illumina as align_illumina_post_round1 } \
          from './modules/align.nf'
include { 
            rattle_polish; \
            pilon as pilon_round1; \
            pilon as pilon_round2; \
            trim_reference; \
            cluster_transcripts \
            } from './modules/polish.nf'
include { salmon_index; salmon_illumina } from './modules/align.nf'
include { 
            transdecoder as transdecoder_round1; \
            transdecoder as transdecoder_round2; \
            salmon_grouper \
            } from './modules/annotate.nf'
include { 
            rev_comp as rev_comp_round1; \
            rev_comp as rev_comp_round2; \
            } from './modules/utility.nf'


// Define the workflow structure
workflow {

    fastq_ch = Channel
                .fromPath(params.input)
                .splitCsv(header: true, sep: ',')
                .map { row -> tuple(params.input_dir + '/' + row.R1,
                                    params.input_dir + '/' + row.R2,
                                    row.sample,
                                    row.platform) }

    illumina_ch = fastq_ch
                    .filter { it[3] == 'illumina' }
                    .map { row -> tuple(row[2], row[0], row[1]) }

    nanopore_ch = fastq_ch
                    .filter { it[3] == 'nanopore' }
                    .map { row -> row[0] }

    pacbio_ch = fastq_ch
                    .filter { it[3] == 'pacbio' }
                    .map { row -> row[0] }

    filter_illumina_ch = filter_illumina_fastq(illumina_ch)

    filter_illumina_ch_r1 = filter_illumina_ch.filtered_r1_r2
                              .map { sample, r1, r2 -> r1 }
                              .collect()

    filter_illumina_ch_r2 = filter_illumina_ch.filtered_r1_r2
                              .map { sample, r1, r2 -> r2 }
                              .collect()

    filter_merge_illumina_r1_ch = merge_illumina_fastq_r1('illumina', 'R1', filter_illumina_ch_r1)
    filter_merge_illumina_r2_ch = merge_illumina_fastq_r2('illumina', 'R2', filter_illumina_ch_r2)
    filter_merge_illumina_ch = filter_merge_illumina_r1_ch.fastq.merge( filter_merge_illumina_r2_ch.fastq )

    filter_nanopore_ch = filter_nanopore_fastq(nanopore_ch)

    filter_pacbio_ch = filter_pacbio_fastq(pacbio_ch)

    rattle_clusters_ch = rattle_assemble(filter_nanopore_ch.filtered_fastq)

    rattle_correct_ch = rattle_correct(filter_nanopore_ch.filtered_fastq, rattle_clusters_ch)

    rattle_polish_ch = rattle_polish(rattle_correct_ch)

    R1_nanopore_align_ch = align_nanopore_round1('round_1',
                        rattle_polish_ch.fasta,
                        filter_nanopore_ch.filtered_fastq)

    R1_illumina_align_ch = align_illumina_round1('round_1',
                        rattle_polish_ch.fasta,
                        filter_merge_illumina_ch)

    R1_pacbio_align_ch = align_pacbio_round1('round_1',
                        rattle_polish_ch.fasta,
                        filter_pacbio_ch.filtered_fastq)

    R1_pilon_ch = pilon_round1('round_1',
                        R1_pacbio_align_ch.bam,
                        R1_pacbio_align_ch.bai,
                        R1_illumina_align_ch.bam,
                        R1_illumina_align_ch.bai,
                        rattle_polish_ch.fasta
                        )

    R1_pilon_rev_ch = rev_comp_round1('round_1', R1_pilon_ch.fasta)

    post_R1_illumina_align_ch = align_illumina_post_round1('post_round_1',
                                                        R1_pilon_ch.fasta,
                                                        filter_merge_illumina_ch)

    R1_pilon_trimmed_ch = trim_reference('round_1',
                                        R1_pilon_ch.fasta,
                                        post_R1_illumina_align_ch.bam,
                                        post_R1_illumina_align_ch.bai,
                                        )
    R1_salmon_index_ch = salmon_index('round_1', R1_pilon_trimmed_ch.fasta)

    R1_salmon_ch = salmon_illumina('round_1',
                                    R1_salmon_index_ch.combine(
                                    filter_illumina_ch.filtered_r1_r2)
                                    )

    R1_transdecoder_ch = transdecoder_round1('round_1',
                                            R1_pilon_trimmed_ch.fasta,
                                            )

    grouper_yaml_ch = Channel.fromPath(params.grouper_config)

    R1_grouper_ch = salmon_grouper('round_1', grouper_yaml_ch, R1_salmon_ch.salmon_quant.collect()) 

    R1_grouped_ref_ch = cluster_transcripts('round_1',
                                        R1_pilon_trimmed_ch.fasta,
                                        R1_transdecoder_ch.longest_ORF_pep,
                                        R1_grouper_ch.mag_cluster) 

    R2_nanopore_align_ch = align_nanopore_round2('round_2',
                        R1_grouped_ref_ch.fasta,
                        filter_nanopore_ch.filtered_fastq)

    R2_illumina_align_ch = align_illumina_round2('round_2',
                        R1_grouped_ref_ch.fasta,
                        filter_merge_illumina_ch)

    R2_pacbio_align_ch = align_pacbio_round2('round_2',
                        R1_grouped_ref_ch.fasta,
                        filter_pacbio_ch.filtered_fastq)

    R2_pilon_ch = pilon_round2('round_2',
                        R2_pacbio_align_ch.bam,
                        R2_pacbio_align_ch.bai,
                        R2_illumina_align_ch.bam,
                        R2_illumina_align_ch.bai,
                        R1_grouped_ref_ch.fasta
                        )

    R2_pilon_rev_ch = rev_comp_round2('round_2', R2_pilon_ch.fasta)

    R2_transdecoder_ch = transdecoder_round2('round_2',
                                            R2_pilon_ch.fasta,
                                            )

    //// separator  ////

    }
