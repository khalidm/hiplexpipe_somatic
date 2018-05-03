'''
Build the pipeline workflow by plumbing the stages together.
'''

from ruffus import Pipeline, suffix, formatter, add_inputs, output_from
from stages import Stages


def make_pipeline(state):
    '''Build the pipeline by constructing stages and connecting them together'''
    # Build an empty pipeline
    pipeline = Pipeline(name='hiplexpipe')
    # Get a list of paths to all the FASTQ files
    fastq_files = state.config.get_option('fastqs')
    # Stages are dependent on the state
    stages = Stages(state)

    # The original FASTQ files
    # This is a dummy stage. It is useful because it makes a node in the
    # pipeline graph, and gives the pipeline an obvious starting point.
    pipeline.originate(
        task_func=stages.original_fastqs,
        name='original_fastqs',
        output=fastq_files)

    # Align paired end reads in FASTQ to the reference producing a BAM file
    pipeline.transform(
        task_func=stages.align_bwa,
        name='align_bwa',
        input=output_from('original_fastqs'),
        # Match the R1 (read 1) FASTQ file and grab the path and sample name.
        # This will be the first input to the stage.
        # Hi-Plex example: OHI031002-P02F04_S318_L001_R1_001.fastq
        # new sample name = OHI031002-P02F04
        filter=formatter(
            '.+/(?P<sample>[a-zA-Z0-9-]+)-(?P<tumor>[TN]+)_(?P<readid>[a-zA-Z0-9-]+)_(?P<lane>[a-zA-Z0-9]+)_R1_(?P<lib>[a-zA-Z0-9-:]+).fastq'),

        # Add one more inputs to the stage:
        #    1. The corresponding R2 FASTQ file
        # Hi-Plex example: OHI031002-P02F04_S318_L001_R2_001.fastq
        add_inputs=add_inputs(
            '{path[0]}/{sample[0]}-{tumor[0]}_{readid[0]}_{lane[0]}_R2_{lib[0]}.fastq'),

        # Add an "extra" argument to the state (beyond the inputs and outputs)
        # which is the sample name. This is needed within the stage for finding out
        # sample specific configuration options
        extras=['{sample[0]}', '{tumor[0]}', '{readid[0]}', '{lane[0]}', '{lib[0]}'],
        # The output file name is the sample name with a .bam extension.
        output='alignments/{sample[0]}/{sample[0]}_{tumor[0]}.bam')

    # Sort the BAM file using Picard
    pipeline.transform(
        task_func=stages.sort_bam_picard,
        name='sort_bam_picard',
        input=output_from('align_bwa'),
        filter=suffix('.bam'),
        output='.sort.bam')

    # High quality and primary alignments
    pipeline.transform(
        task_func=stages.primary_bam,
        name='primary_bam',
        input=output_from('sort_bam_picard'),
        filter=suffix('.sort.bam'),
        output='.primary.bam')

    # index bam file
    pipeline.transform(
        task_func=stages.index_sort_bam_picard,
        name='index_bam',
        input=output_from('primary_bam'),
        filter=suffix('.primary.bam'),
        output='.primary.bam.bai')

    # Clip the primer_seq from BAM File
    (pipeline.transform(
        task_func=stages.clip_bam,
        name='clip_bam',
        input=output_from('primary_bam'),
        filter=suffix('.primary.bam'),
        output='.primary.primerclipped.bam')
        .follows('index_bam'))

    ###### GATK VARIANT CALLING - MuTect2 ######

    # Call somatics variants using MuTect2
    pipeline.transform(
        task_func=stages.call_mutect2_gatk,
        name='call_mutect2_gatk',
        input=output_from('clip_bam'),
        # filter=suffix('.merged.dedup.realn.bam'),
        filter=formatter('.+/(?P<sample>[a-zA-Z0-9-]+)_T.primary.primerclipped.bam'),
        add_inputs=add_inputs(
            '{path[0]}/{sample[0]}_N.primary.primerclipped.bam'),
        # extras=['{sample[0]}'],
        output='variants/mutect2/{sample[0]}.mutect2.vcf')
        # .follows('clip_bam')

    ###### GATK VARIANT CALLING - MuTect2 ######

    # -------- VEP ----------
    # Apply NORM
    (pipeline.transform(
        task_func=stages.apply_vt,
        name='apply_vt',
        input=output_from('call_mutect2_gatk'),
        filter=suffix('.mutect2.vcf'),
        # add_inputs=add_inputs(['variants/ALL.indel_recal', 'variants/ALL.indel_tranches']),
        output='.mutect2.vt.vcf')
        .follows('call_mutect2_gatk'))
    #
    # Apply VEP
    (pipeline.transform(
        task_func=stages.apply_vep,
        name='apply_vep',
        input=output_from('apply_vt'),
        filter=suffix('.mutect2.vt.vcf'),
        # add_inputs=add_inputs(['variants/ALL.indel_recal', 'variants/ALL.indel_tranches']),
        output='.mutect2.vt.vep.vcf')
        .follows('apply_vt'))
    #
    # Apply vcfanno
    (pipeline.transform(
        task_func=stages.apply_vcfanno,
        name='apply_vcfanno',
        input=output_from('apply_vep'),
        filter=suffix('.mutect2.vt.vep.vcf'),
        # add_inputs=add_inputs(['variants/ALL.indel_recal', 'variants/ALL.indel_tranches']),
        output='.mutect2.annotated.vcf')
        .follows('apply_vep'))

    return pipeline
