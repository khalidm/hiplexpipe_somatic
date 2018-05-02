'''
Individual stages of the pipeline implemented as functions from
input files to output files.

The run_stage function knows everything about submitting jobs and, given
the state parameter, has full access to the state of the pipeline, such
as config, options, DRMAA and the logger.
'''

from utils import safe_make_dir
from runner import run_stage
import os
import pysam

# PICARD_JAR = '$PICARD_HOME/lib/picard-1.69.jar'
# PICARD_JAR = '/vlsci/VR0002/kmahmood/Programs/Picard/picard-tools-2.8.3/picard.jar'
PICARD_JAR = '/usr/local/easybuild/software/picard/2.3.0/picard.jar'
SNPEFF_JAR = '/usr/local/easybuild/software/snpEff/4.1d-Java-1.7.0_80/snpEff.jar'

GATK_JAR = '$GATK_HOME/GenomeAnalysisTK.jar'

def java_command(jar_path, mem_in_gb, command_args):
    '''Build a string for running a java command'''
    # Bit of room between Java's max heap memory and what was requested.
    # Allows for other Java memory usage, such as stack.
    java_mem = mem_in_gb - 2
    return 'java -Xmx{mem}g -jar {jar_path} {command_args}'.format(
        jar_path=jar_path, mem=java_mem, command_args=command_args)

def run_java(state, stage, jar_path, mem, args):
    command = java_command(jar_path, mem, args)
    run_stage(state, stage, command)

class Stages(object):
    def __init__(self, state):
        self.state = state
        self.reference = self.get_options('ref_grch37')
        self.dbsnp_hg19 = self.get_options('dbsnp_hg19')
        self.mills_hg19 = self.get_options('mills_hg19')
        self.one_k_g_snps = self.get_options('one_k_g_snps')
        self.one_k_g_indels = self.get_options('one_k_g_indels')
        self.one_k_g_highconf_snps = self.get_options('one_k_g_highconf_snps')
        self.hapmap = self.get_options('hapmap')
        # self.interval_hg19 = self.get_options('exome_bed_hg19')
        # self.CEU_mergeGvcf = self.get_options('CEU_mergeGvcf')
        self.snpeff_conf = self.get_options('snpeff_conf')
        self.bamclipper = self.get_options('bamclipper')
        self.vep_path = self.get_options('vep_path')
        self.vt_path = self.get_options('vt_path')
        self.coord_file = self.get_options('coord_file')
        self.target_bed = self.get_options('target_bed')
        # self.interval_file = self.get_options('interval_file')
        self.primer_file = self.get_options('primer_file')
        self.primer_bedpe_file = self.get_options('primer_bedpe_file')
        self.proportionthresh = self.get_options('proportionthresh')
        self.absthresh = self.get_options('absthresh')
        self.maxvariants = self.get_options('maxvariants')
        # self.fragment_bed = self.get_options('fragment_bed')
        self.annolua = self.get_options('annolua')
        self.anno = self.get_options('anno')
        self.hrfile = self.get_options('hrfile')
        self.other_vep = self.get_options('other_vep')
        self.snpeff_path = self.get_options('snpeff_path')
        self.mutect2_gnomad = self.get_options('mutect2_gnomad')

        # self.GBR_mergeGvcf = self.get_options('GBR_mergeGvcf')
        # self.FIN_mergeGvcf = self.get_options('FIN_mergeGvcf')

    def run_picard(self, stage, args):
        mem = int(self.state.config.get_stage_options(stage, 'mem'))
        return run_java(self.state, stage, PICARD_JAR, mem, args)

    def run_snpeff(self, stage, args):
        mem = int(self.state.config.get_stage_options(stage, 'mem'))
        return run_java(self.state, stage, SNPEFF_JAR, mem, args)

    def run_gatk(self, stage, args):
        mem = int(self.state.config.get_stage_options(stage, 'mem'))
        return run_java(self.state, stage, GATK_JAR, mem, args)

    def get_stage_options(self, stage, *options):
        return self.state.config.get_stage_options(stage, *options)

    def get_options(self, *options):
        return self.state.config.get_options(*options)

    def original_fastqs(self, output):
        '''Original fastq files'''
        # print output
        pass

    def align_bwa(self, inputs, bam_out, sample_id, tumor_id, read_id, lane, lib):
        # def align_bwa(self, inputs, bam_out, sample_id):
        '''Align the paired end fastq files to the reference genome using bwa'''
        fastq_read1_in, fastq_read2_in = inputs
        cores = self.get_stage_options('align_bwa', 'cores')
        safe_make_dir('alignments/{sample}'.format(sample=sample_id))
        read_group = '"@RG\\tID:{readid}\\tSM:{sample}_{tumor_id}_{readid}\\tPU:lib1\\tLN:{lane}\\tPL:Illumina"' \
            .format(readid=read_id, lib=lib, lane=lane, sample=sample_id, tumor_id=tumor_id)
        command = 'bwa mem -M -t {cores} -R {read_group} {reference} {fastq_read1} {fastq_read2} ' \
                  '| samtools view -b -h -o {bam} -' \
                  .format(cores=cores,
                          read_group=read_group,
                          fastq_read1=fastq_read1_in,
                          fastq_read2=fastq_read2_in,
                          reference=self.reference,
                          bam=bam_out)
        run_stage(self.state, 'align_bwa', command)

    # def apply_undr_rover(self, inputs, vcf_output, sample_id, readid):
    #     # def align_bwa(self, inputs, bam_out, sample_id):
    #     '''Apply undr_rover to call variants from paired end fastq files'''
    #     fastq_read1_in, fastq_read2_in = inputs
    #     cores = self.get_stage_options('apply_undr_rover', 'cores')
    #     safe_make_dir('variants/undr_rover')
    #     safe_make_dir('variants/undr_rover/coverdir')
    #     coverfile = "variants/undr_rover/coverdir/" + sample_id + "_" + readid + ".coverage"
    #
    #     command = 'undr_rover --primer_coords {coord_file} ' \
    #               '--primer_sequences {primer_file} ' \
    #               '--reference {reference} ' \
    #               '--out {vcf_output} ' \
    #               '--coverfile {coverfile} ' \
    #               '--proportionthresh {proportionthresh} ' \
    #               '--absthresh {absthresh} ' \
    #               '--max_variants {maxvariants} ' \
    #               '{fastq_read1} {fastq_read2}'.format(
    #                     coord_file=self.coord_file, primer_file=self.primer_file,
    #                     reference=self.reference,
    #                     vcf_output=vcf_output,
    #                     #coverdir=self.coverdir,
    #                     proportionthresh=self.proportionthresh,
    #                     absthresh=self.absthresh,
    #                     maxvariants=self.maxvariants,
    #                     coverfile=coverfile,
    #                     fastq_read1=fastq_read1_in,
    #                     fastq_read2=fastq_read2_in)
    #     run_stage(self.state, 'apply_undr_rover', command)

    def clip_bam(self, bam_in, sorted_bam_out):
        '''Clip the BAM file using Bamclipper'''
        bamclipper_args = '{bamclipper} -b {bam_in} -p {primer_bedpe_file} -n 1'.format(
                          bamclipper=self.bamclipper, bam_in=bam_in, primer_bedpe_file=self.primer_bedpe_file)
        run_stage(self.state, 'clip_bam', bamclipper_args)

    def sort_bam_picard(self, bam_in, sorted_bam_out):
        '''Sort the BAM file using Picard'''
        picard_args = 'SortSam INPUT={bam_in} OUTPUT={sorted_bam_out} ' \
                      'VALIDATION_STRINGENCY=LENIENT SORT_ORDER=coordinate ' \
                      'MAX_RECORDS_IN_RAM=5000000 CREATE_INDEX=True'.format(
                          bam_in=bam_in, sorted_bam_out=sorted_bam_out)
        self.run_picard('sort_bam_picard', picard_args)

    def primary_bam(self, bam_in, sbam_out):
        '''On keep primary alignments in the BAM file using samtools'''
        command = 'samtools view -h -q 1 -f 2 -F 4 -F 8 -F 256 -b ' \
                    '-o {sbam_out} {bam_in}'.format(
                        bam_in=bam_in, sbam_out=sbam_out)
        run_stage(self.state, 'primary_bam', command)

    # index sorted bam file
    def index_sort_bam_picard(self, bam_in, bam_index):
        '''Index sorted bam using samtools'''
        command = 'samtools index {bam_in} {bam_index}'.format(
                          bam_in=bam_in, bam_index=bam_index)
        run_stage(self.state, 'index_sort_bam_picard', command)

    # coverage bam
    def call_mutect2_gatk(self, inputs, vcf_out):
        '''Call somatic variants from using MuTect2'''
        tumor_in, normal_in = inputs
        tumor_samfile = pysam.AlignmentFile(tumor_in, "rb")
        normal_samfile = pysam.AlignmentFile(normal_in, "rb")
        tumor_id = tumor_samfile.header['RG'][0]['SM']
        normal_id = normal_samfile.header['RG'][0]['SM']
        tumor_samfile.close()
        normal_samfile.close()
        # safe_make_dir('variants/mutect2/{sample}'.format(sample=sample_id))
        safe_make_dir('variants/mutect2/')
        command = "gatk Mutect2 -R {reference} " \
            "-I {tumor_in} " \
            "-tumor {tumor_id} " \
            "-I {normal_in} " \
            "-normal {normal_id} " \
            "--germline-resource {mutect2_gnomad} " \
            "--af-of-alleles-not-in-resource 0.00003125 " \
            "-O {out} " \
            "--dont-use-soft-clipped-bases".format(reference=self.reference,
                        tumor_in=tumor_in,
                        normal_in=normal_in,
                        tumor_id=tumor_id,
                        normal_id=normal_id,
                        mutect2_gnomad=self.mutect2_gnomad,
                        out=vcf_out)
        run_stage(self.state, 'call_mutect2_gatk', command)

    # # Merge per lane bam into a single bam per sample
    # def merge_sample_bams(self, bam_files_in, bam_out):
    #     '''Merge per lane bam into a merged bam file'''
    #     bam_files = ' '.join(['INPUT=' + bam for bam in bam_files_in])
    #     picard_args = 'MergeSamFiles {bams_in} OUTPUT={merged_bam_out} ' \
    #                   'VALIDATION_STRINGENCY=LENIENT ' \
    #                   'MAX_RECORDS_IN_RAM=5000000 ASSUME_SORTED=True ' \
    #                   'CREATE_INDEX=True'.format(
    #                       bams_in=bam_files, merged_bam_out=bam_out)
    #     self.run_picard('merge_sample_bams', picard_args)

    # coverage bam interval
    # def target_coverage_bamutil_interval(self, bam_in, coverage_out):
    #     '''Calculate target coverage using bamutil'''
    #     command = 'bam stats --basic --in {bam_in} --regionList {fragment_bed} &> {coverage_out}'.format(
    #                       bam_in=bam_in, fragment_bed = self.fragment_bed, coverage_out=coverage_out)
    #     run_stage(self.state, 'target_coverage_bamutil_interval', command)

    # multicov
    def apply_multicov(self, bam_in, multicov):
        '''Samtools mpileup'''
        # bam_in = bam_in
        bams = ' '.join([bam for bam in bam_in])
        safe_make_dir('coverage')
        command = 'bedtools multicov -bams {bams} -bed {target_bed} > {multicov}'.format(
                          bams=bams, target_bed=self.target_bed, multicov=multicov)
        run_stage(self.state, 'apply_multicov', command)

    # multicov plots
    def apply_multicov_plots(self, bam_in, multicov):
        '''Generate multicov plots'''
        walltime_hours = self.get_stage_options('apply_multicov_plots', 'walltime')
        h, m = walltime_hours.split(':')
        walltime = int(h) * 3600 + int(m) * 60
        # safe_make_dir('variants')
        # command = 'jupyter nbconvert --to html --execute coverage_analysis_main.ipynb ' \
        #             '--ExecutePreprocessor.timeout={walltime}'.format(walltime=walltime)
        command = 'jupyter nbconvert --to html coverage_analysis_main.ipynb ' \
                    '--ExecutePreprocessor.timeout={walltime}'.format(walltime=walltime)
        run_stage(self.state, 'apply_multicov_plots', command)

    # summarize picard
    def apply_summarize_picard(self, input, output):
        '''Summarize picard coverage'''
        input = input
        # bams = ' '.join([bam for bam in bam_in])
        # safe_make_dir('variants')
        command = 'python coverage_summary.py > {output} '.format(
                          output=output)
        run_stage(self.state, 'apply_summarize_picard', command)

    # samtools
    def apply_samtools_mpileup(self, bam_in, mpileup_out_bcf):
        '''Samtools mpileup'''
        # bam_in = bam_in
        bams = ' '.join([bam for bam in bam_in])
        safe_make_dir('variants')
        command = 'samtools mpileup -t DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR -go {mpileup_out_bcf} ' \
                  '-f {reference} {bams}'.format(
                          mpileup_out_bcf=mpileup_out_bcf,reference=self.reference,bams=bams)
        run_stage(self.state, 'apply_samtools_mpileup', command)

    # bcftools
    def apply_bcftools(self, mpileup_in, vcf_out):
        '''Bcftools call variants'''
        mpileup_in = mpileup_in
        # mpileup_in = ' '.join([vcf for vcf in vcf_files_in])
        command = 'bcftools call -vmO v -o {vcf_out} {mpileup_in}'.format(
                          vcf_out=vcf_out,mpileup_in=mpileup_in)
        run_stage(self.state, 'apply_bcftools', command)

    def apply_vt(self, inputs, vcf_out):
        '''Apply NORM'''
        vcf_in = inputs
        cores = self.get_stage_options('apply_vt', 'cores')
        vt_command = "{vt_path} decompose -s {vcf_in} - | {vt_path2} normalize -r {reference} - | " \
                    "{vt_path3} uniq - -o {vcf_out}".format(
                    vt_path=self.vt_path, vcf_in=vcf_in, vt_path2=self.vt_path, reference=self.reference,
                    vt_path3=self.vt_path, vcf_out=vcf_out)
        run_stage(self.state, 'apply_vt', vt_command)

    def apply_vep(self, inputs, vcf_out):
        '''Apply VEP'''
        vcf_in = inputs
        cores = self.get_stage_options('apply_vep', 'cores')
        vep_command = "{vep_path}/variant_effect_predictor.pl --cache --refseq --offline {other_vep} --fasta {reference} " \
                    "-i {vcf_in} --sift b --polyphen b --symbol --numbers --biotype --total_length --hgvs " \
                    "--format vcf -o {vcf_vep} --force_overwrite --vcf " \
                    "--fields Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT," \
                    "Protein_position,BIOTYPE,HGVSc,HGVSp,cDNA_position,CDS_position,HGVSc,HGVSp,cDNA_position,CDS_position,PICK " \
                    "--fork {threads} --flag_pick".format(
                    reference=self.reference, vep_path=self.vep_path, vcf_in=vcf_in, vcf_vep=vcf_out, other_vep=self.other_vep, threads=cores)
        run_stage(self.state, 'apply_vep', vep_command)

    def apply_bcf(self, inputs, vcf_out):
        '''Apply BCF'''
        vcf_in = inputs
        cores = self.get_stage_options('apply_bcf', 'cores')
        command = "bcftools filter -e \"ALT='*'\" {vcf_in} > {vcf_out}".format(cores=cores,
                            vcf_in=vcf_in, vcf_out=vcf_out)
        run_stage(self.state, 'apply_bcf', command)

    def apply_snpeff(self, inputs, vcf_out):
        '''Apply SnpEFF'''
        vcf_in = inputs
        #cores = self.get_stage_options('apply_snpeff', 'cores')  apply_snpeff
        # mem = int(self.state.config.get_stage_options(stage, 'mem'))
        mem = int(self.get_stage_options('apply_snpeff', 'mem')) - 2
        snpeff_command = "java -Xmx{mem}g -jar {snpeff_path} eff -c {snpeff_conf} " \
                    "-canon GRCh37.75 {vcf_in} | bgzip -c > {vcf_out}".format(
                    mem=mem, snpeff_path=self.snpeff_path, snpeff_conf=self.snpeff_conf,
                    vcf_in=vcf_in, vcf_out=vcf_out)
        run_stage(self.state, 'apply_snpeff', snpeff_command)
        #run_snpeff(self.state, 'apply_snpeff', snpeff_command)

    # ORIGINAL WITH RUN_JAVA
    # def apply_snpeff(self, inputs, vcf_out):
    #     '''Apply SnpEFF'''
    #     vcf_in = inputs
    #     #cores = self.get_stage_options('apply_snpeff', 'cores')
    #     snpeff_command = "eff -c {snpeff_conf} -canon GRCh37.75 {vcf_in} > {vcf_out}".format(
    #                 snpeff_conf=self.snpeff_conf, vcf_in=vcf_in, vcf_out=vcf_out)
    #     self.run_snpeff('apply_snpeff', snpeff_command)
    #     #run_snpeff(self.state, 'apply_snpeff', snpeff_command)

    def apply_vcfanno(self, inputs, vcf_out):
        '''Apply anno'''
        vcf_in = inputs
        #cores = self.get_stage_options('apply_snpeff', 'cores')
        anno_command = "./vcfanno_linux64 -lua {annolua} {anno} {vcf_in} > {vcf_out}".format(
                    annolua=self.annolua, anno=self.anno, vcf_in=vcf_in, vcf_out=vcf_out)
        run_stage(self.state, 'apply_vcfanno', anno_command)

    def apply_cat_vcf(self, inputs, vcf_out):
        '''Concatenate and sort undr_rover VCF files for downstream analysis'''
        vcfs = ' '.join([vcf for vcf in inputs])
        # safe_make_dir('variants')
        command = 'vcf-concat {vcfs} | vcf-sort -c | bgzip -c > {vcf_out} '.format(vcfs=vcfs,vcf_out=vcf_out)
        run_stage(self.state, 'apply_cat_vcf', command)

    def apply_tabix(self, input, vcf_out):
        '''bgzip the vcf file in prepartion for bcftools annotation'''
        vcf = input
        command = "tabix -p vcf {vcf}".format(vcf=vcf)
        run_stage(self.state, 'apply_tabix', command)

    def apply_homopolymer_ann(self, inputs, vcf_out):
        '''Apply HomopolymerRun annotation to undr_rover output'''
        vcf_in = inputs
        # safe_make_dir('variants')
        command = "echo \"##INFO=<ID=HRUN,Number=1,Type=String,Description=\"HRun\">\" > header.tmp; "\
                    "bcftools annotate -a {hrfile} -c CHROM,FROM,TO,HRUN " \
                    "-h header.tmp " \
                    "{vcf_in} > {vcf_out}".format(hrfile=self.hrfile,vcf_in=vcf_in,vcf_out=vcf_out)
        run_stage(self.state, 'apply_cat_vcf', command)

    # def apply_cat_vcf(self, inputs, vcf_out):
    #     '''Concatenate and sort undr_rover VCF files for downstream analysis'''
    #     vcfs = ' '.join([vcf for vcf in inputs])
    #     # safe_make_dir('variants')
    #     command = 'vcf-concat {vcfs} | vcf-sort -c | bgzip -c > {vcf_out} '.format(vcfs=vcfs,vcf_out=vcf_out)
    #     run_stage(self.state, 'apply_cat_vcf', command)
