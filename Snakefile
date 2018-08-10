# to qsub snakemake:
#source activate snakemake

#qsub -N tnb_SAMN03979292 -b y -wd \
#/share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/triple_neg_breast_cancer/ \
#-j y -R y -pe smp 32 -V \
#"snakemake -s /share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/triple_neg_breast_cancer/Snakefile --cores 12"


#############################################################
### 0. Set up variables and paths ###
#############################################################

import os
import re
os.system('module load gi/zlib/1.2.8')
os.system('module load phuluu/samtools/1.4')
os.system('module load gi/novosort/precompiled/1.03.08')
os.system('hugfre/HTSeq/0.5.4p3')

# define/create directories:
home_dir = '/share/ScratchGeneral/jamtor/'

project_dir = home_dir + 'projects/hgsoc_repeats/RNA-seq/triple_neg_breast_cancer/'

raw_dir = project_dir + 'raw_files/'

results_dir = project_dir + 'results/'

sra_dir = 'raw_files/sra/'

fastq_dir = 'raw_files/fastq/'
if not os.path.exists(project_dir + '/' + fastq_dir):
	os.makedirs(project_dir + '/' + fastq_dir)

#fastqc_dir = 'results/fastqc/'
#if not os.path.exists(project_dir + '/' + fastqc_dir):
#	os.makedirs(project_dir + '/' + fastqc_dir)

star_GC_dir = 'results/star/GC/'
if not os.path.exists(project_dir + '/' + star_GC_dir):
	os.makedirs(project_dir + '/' + star_GC_dir)

temp_sort_dir = star_GC_dir + 'temp/'
if not os.path.exists(project_dir + '/' + temp_sort_dir):
	os.makedirs(project_dir + '/' + temp_sort_dir)

htseq_dir = 'results/htseq/'
if not os.path.exists(project_dir + '/' + htseq_dir):
	os.makedirs(project_dir + '/' + htseq_dir)

	
SAMPLES = [ re.sub('.sra', '', x) for x in list(os.walk(raw_dir + '/sra'))[0][2] ]

print(' '.join(SAMPLES))


#############################################################
### 1. Extract sras to fastq.gz files and QC ###
#############################################################

# convert all SAMPLE sras to fastqs and move to correct 
# directory:
for s in SAMPLES:
	if os.path.isfile(project_dir + 'raw_files/fastq/' + 
		s + '_pass_1.fastq.gz') and \
		os.path.isfile(project_dir + 'raw_files/fastq/' + 
			s + '_pass_1.fastq.gz'):
		print(project_dir + 'raw_files/fastq/' + 
			s + '_pass_1.fastq.gz and ' +
			project_dir + 'raw_files/fastq/' + s + 
			'_pass_1.fastq.gz already exist')
	else:
		print('Creating ' + project_dir + 'raw_files/fastq/' + 
			s + '.fastq.gz')
		os.system('fastq-dump --outdir ' + project_dir + 
			'raw_files/fastq/' + 
			' --gzip --skip-technical --readids ' + 
			'--read-filter pass --dumpbase --split-3 ' + 
			'--clip ' + sra_dir + '/' + s + '.sra')

rule all:
	input:
		expand(star_GC_dir + '{sample}/Aligned.novosortedByName.out.bam', \
			sample=SAMPLES),
		expand(htseq_dir + '{sample}.gc.htseq.txt', \
			sample=SAMPLES),			
		expand(htseq_dir + '{sample}.custom3.htseq.txt', \
			sample=SAMPLES)

#rule fastqc1:
#	input:
#		fastq_dir + '{sample}_pass_1.fastq.gz'
#	output:
#		fastqc_dir + '{sample}_pass_1_fastqc.html',
#		fastqc_dir + '{sample}_pass_1_fastqc.zip'
#	shell:
#		'fastqc -t 2 -o ' + fastqc_dir + ' {input}'
#
#rule fastqc2:
#	input:
#		fastq_dir + '{sample}_pass_2.fastq.gz'
#	output:
#		fastqc_dir + '{sample}_pass_2_fastqc.html',
#		fastqc_dir + '{sample}_pass_2_fastqc.zip'
#	shell:
#		'fastqc -t 2 -o ' + fastqc_dir + ' {input}'

rule star:
	input:
		fq1 = fastq_dir + '{sample}_pass_1.fastq.gz',
		fq2 = fastq_dir + '{sample}_pass_2.fastq.gz'
	output:
		star_GC_dir + '{sample}/Aligned.out.bam',
		star_GC_dir + '{sample}/Log.final.out',
		star_GC_dir + '{sample}/Chimeric.out.junction',
		star_GC_dir + '{sample}/Chimeric.out.sam'
	threads: 6
	shell:
		'/home/jamtor/local/bin/STAR --runMode alignReads ' +
      	' --readFilesCommand zcat ' +
	    '--genomeDir /home/jamtor/genomes/hg38_ercc/starRef ' +
	    '--outFilterType BySJout ' +
	    '--outSAMattributes NH HI AS NM MD ' +
	    '--outFilterMultimapNmax 999 ' +
	    '--outMultimapperOrder Random ' +
	    '--runRNGseed 666 ' +
	    '--outSAMmultNmax 1 ' +
	    '--outFilterMismatchNmax 999 ' +
	    '--outFilterMismatchNoverReadLmax 0.04 ' +
	    '--alignIntronMin 20 ' +
	    '--alignIntronMax 1500000 ' +
	    '--alignMatesGapMax 1500000 ' +
	    '--alignSJoverhangMin 6 ' +
	    '--alignSJDBoverhangMin 1 ' +
	    '--readFilesIn {input.fq1} {input.fq2} ' +
	    '--outFileNamePrefix ' + project_dir + '/' + 
	    	star_GC_dir + '{wildcards.sample}' +
	    '/ --runThreadN 6 ' +
	    '--outFilterMatchNmin 76 ' +
	  	'--chimSegmentMin 25 ' +
	    '--chimJunctionOverhangMin 25 ' +
	    '--chimScoreMin 0 ' +
	    '--chimScoreDropMax 20 ' +
	    '--chimScoreSeparation 10 ' +
	    '--chimScoreJunctionNonGTAG -1 ' +
	    '--outSAMtype BAM Unsorted'

rule novosort:
	input:
		star_GC_dir + '{sample}/Aligned.out.bam'
	output:
		star_GC_dir + '{sample}/Aligned.novosortedByName.out.bam'
	threads: 6
	shell:
		'/home/jamtor/local/bin/novocraft/novosort -t ' + 
		temp_sort_dir + ' -n -c 6 -m 22G {input} > {output}'


##########################################################################
### 5. Count gc genes from bam using htseq ###
##########################################################################

rule htseq_gc:
	input:
		star_GC_dir + '{sample}/Aligned.novosortedByName.out.bam'
	output:
		htseq_dir + '{sample}.gc.htseq.txt'
	threads: 6
	shell:
		'source activate p2.7env; htseq-count -f bam -i ID -t exon -a 0 ' +
		'-m intersection-strict {input} ' +
		'/home/jamtor/genomes/hg38_ercc/gencode_v24_hg38_annotation.gff ' +
		'>> {output}'


##########################################################################
### 6. Count custom3 genes from bam using htseq ###
##########################################################################

rule htseq_custom3:
	input:
		star_GC_dir + '{sample}/Aligned.novosortedByName.out.bam'
	output:
		htseq_dir + '{sample}.custom3.htseq.txt'
	threads: 6
	shell:
		'source activate p2.7env; htseq-count -f bam -i ID -t exon --stranded=no -a 0 ' +
		'-m intersection-strict {input} ' +
		'/home/jamtor/genomes/repeats/custom3rep.final.gff ' +
		'>> {output}'