#!/bin/bash
# Exit immediately on error
unset MODULERCFILE
set -eu -o pipefail

#-------------------------------------------------------------------------------
# ChipSeq Batch Job Submission Bash script
# Version: 3.1.5-beta
# Created on: 2020-03-08T19:08:54
# Steps:
#   picard_sam_to_fastq: 0 job... skipping
#   trimmomatic: 3 jobs
#   merge_trimmomatic_stats: 1 job
#   bwa_mem_picard_sort_sam: 3 jobs
#   samtools_view_filter: 4 jobs
#   picard_merge_sam_files: 3 jobs
#   picard_mark_duplicates: 4 jobs
#   metrics: 5 jobs
#   homer_make_tag_directory: 3 jobs
#   qc_metrics: 1 job
#   homer_make_ucsc_file: 7 jobs
#   macs2_callpeak: 5 jobs
#   homer_annotate_peaks: 3 jobs
#   homer_find_motifs_genome: 3 jobs
#   annotation_graphs: 1 job
#   ihec_preprocess_files: 6 jobs
#   run_spp: 3 jobs
#   ihec_metrics: 4 jobs
#   multiqc_report: 1 job
#   cram_output: 3 jobs
#   TOTAL: 63 jobs
#-------------------------------------------------------------------------------

OUTPUT_DIR=/home/normie/Documents/mcdevops/genpipes/pipelines/chipseq
JOB_OUTPUT_DIR=$OUTPUT_DIR/job_output
TIMESTAMP=`date +%FT%H.%M.%S`
JOB_LIST=$JOB_OUTPUT_DIR/ChipSeq_job_list_$TIMESTAMP
export CONFIG_FILES="/home/normie/Documents/mcdevops/genpipes/pipelines/chipseq/chipseq.base.ini,/home/normie/Documents/mcdevops/genpipes/pipelines/chipseq/cit.ini,/home/normie/Documents/mcdevops/genpipes/pipelines/chipseq/chipseq.batch.ini"
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR

SEPARATOR_LINE=`seq -s - 80 | sed 's/[0-9]//g'`

#-------------------------------------------------------------------------------
# STEP: trimmomatic
#-------------------------------------------------------------------------------
STEP=trimmomatic
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: trimmomatic.GM12878_chr19_input
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.GM12878_chr19_input
JOB_DONE=job_output/trimmomatic/trimmomatic.GM12878_chr19_input.e17135835584fd4cc26be0a0073c09ec.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'trimmomatic.GM12878_chr19_input.e17135835584fd4cc26be0a0073c09ec.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/GM12878_chr19_input && \
`cat > trim/GM12878_chr19_input/GM12878_chr19_input.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /cvmfs/soft.mugqic/CentOS6/testdata/chipseq/raw_data/chipseq_GM12878_chr19_input.fastq.gz \
  trim/GM12878_chr19_input/GM12878_chr19_input.trim.single.fastq.gz \
  ILLUMINACLIP:trim/GM12878_chr19_input/GM12878_chr19_input.trim.adapters.fa:2:30:15 \
  TRAILING:20 \
  MINLEN:25 \
  2> trim/GM12878_chr19_input/GM12878_chr19_input.trim.log
trimmomatic.GM12878_chr19_input.e17135835584fd4cc26be0a0073c09ec.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# JOB: trimmomatic.GM12878_chr19_CTCF
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.GM12878_chr19_CTCF
JOB_DONE=job_output/trimmomatic/trimmomatic.GM12878_chr19_CTCF.93165316ad08495c2745b4e789aaf3b0.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'trimmomatic.GM12878_chr19_CTCF.93165316ad08495c2745b4e789aaf3b0.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/GM12878_chr19_CTCF && \
`cat > trim/GM12878_chr19_CTCF/GM12878_chr19_CTCF.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /cvmfs/soft.mugqic/CentOS6/testdata/chipseq/raw_data/chipseq_GM12878_chr19_CTCF.fastq.gz \
  trim/GM12878_chr19_CTCF/GM12878_chr19_CTCF.trim.single.fastq.gz \
  ILLUMINACLIP:trim/GM12878_chr19_CTCF/GM12878_chr19_CTCF.trim.adapters.fa:2:30:15 \
  TRAILING:20 \
  MINLEN:25 \
  2> trim/GM12878_chr19_CTCF/GM12878_chr19_CTCF.trim.log
trimmomatic.GM12878_chr19_CTCF.93165316ad08495c2745b4e789aaf3b0.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# JOB: trimmomatic.GM12878_chr19_Rad21
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.GM12878_chr19_Rad21
JOB_DONE=job_output/trimmomatic/trimmomatic.GM12878_chr19_Rad21.0233f173ac848d4f55b8d6814814b842.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'trimmomatic.GM12878_chr19_Rad21.0233f173ac848d4f55b8d6814814b842.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/GM12878_chr19_Rad21 && \
`cat > trim/GM12878_chr19_Rad21/GM12878_chr19_Rad21.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /cvmfs/soft.mugqic/CentOS6/testdata/chipseq/raw_data/chipseq_GM12878_chr19_Rad21.fastq.gz \
  trim/GM12878_chr19_Rad21/GM12878_chr19_Rad21.trim.single.fastq.gz \
  ILLUMINACLIP:trim/GM12878_chr19_Rad21/GM12878_chr19_Rad21.trim.adapters.fa:2:30:15 \
  TRAILING:20 \
  MINLEN:25 \
  2> trim/GM12878_chr19_Rad21/GM12878_chr19_Rad21.trim.log
trimmomatic.GM12878_chr19_Rad21.0233f173ac848d4f55b8d6814814b842.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# STEP: merge_trimmomatic_stats
#-------------------------------------------------------------------------------
STEP=merge_trimmomatic_stats
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: merge_trimmomatic_stats
#-------------------------------------------------------------------------------
JOB_NAME=merge_trimmomatic_stats
JOB_DONE=job_output/merge_trimmomatic_stats/merge_trimmomatic_stats.40c91abbc83f77d68cbc9552d6b56fa6.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'merge_trimmomatic_stats.40c91abbc83f77d68cbc9552d6b56fa6.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/perl/5.22.1 mugqic/pandoc/1.15.2 && \
mkdir -p metrics && \
echo 'Sample	Readset	Raw Single Reads #	Surviving Single Reads #	Surviving Single Reads %' > metrics/trimReadsetTable.tsv && \
grep ^Input trim/GM12878_chr19_input/GM12878_chr19_input.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/GM12878_chr19_input	GM12878_chr19_input	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/GM12878_chr19_CTCF/GM12878_chr19_CTCF.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/GM12878_chr19_CTCF	GM12878_chr19_CTCF	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/GM12878_chr19_Rad21/GM12878_chr19_Rad21.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/GM12878_chr19_Rad21	GM12878_chr19_Rad21	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
cut -f1,3- metrics/trimReadsetTable.tsv | awk -F"	" '{OFS="	"; if (NR==1) {if ($2=="Raw Paired Reads #") {paired=1};print "Sample", "Raw Reads #", "Surviving Reads #", "Surviving %"} else {if (paired) {$2=$2*2; $3=$3*2}; raw[$1]+=$2; surviving[$1]+=$3}}END{for (sample in raw){print sample, raw[sample], surviving[sample], surviving[sample] / raw[sample] * 100}}' \
  > metrics/trimSampleTable.tsv && \
mkdir -p report && \
cp metrics/trimReadsetTable.tsv metrics/trimSampleTable.tsv report/ && \
trim_readset_table_md=`LC_NUMERIC=en_CA awk -F "	" '{OFS="|"; if (NR == 1) {$1 = $1; print $0; print "-----|-----|-----:|-----:|-----:"} else {print $1, $2, sprintf("%\47d", $3), sprintf("%\47d", $4), sprintf("%.1f", $5)}}' metrics/trimReadsetTable.tsv` && \
pandoc \
  /home/normie/Documents/mcdevops/genpipes/bfx/report/Illumina.merge_trimmomatic_stats.md \
  --template /home/normie/Documents/mcdevops/genpipes/bfx/report/Illumina.merge_trimmomatic_stats.md \
  --variable trailing_min_quality=20 \
  --variable min_length=25 \
  --variable read_type=Single \
  --variable trim_readset_table="$trim_readset_table_md" \
  --to markdown \
  > report/Illumina.merge_trimmomatic_stats.md
merge_trimmomatic_stats.40c91abbc83f77d68cbc9552d6b56fa6.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# STEP: bwa_mem_picard_sort_sam
#-------------------------------------------------------------------------------
STEP=bwa_mem_picard_sort_sam
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam.GM12878_chr19_input
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.GM12878_chr19_input
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.GM12878_chr19_input.ee95a874f5f3c3b7db9c9d3b7fc7e470.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'bwa_mem_picard_sort_sam.GM12878_chr19_input.ee95a874f5f3c3b7db9c9d3b7fc7e470.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/GM12878_chr19_input/GM12878_chr19_input && \
bwa mem  \
  -M -t 1 \
  -R '@RG	ID:GM12878_chr19_input	SM:GM12878_chr19_input	LB:1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/bwa_index/Homo_sapiens.hg19.fa \
  trim/GM12878_chr19_input/GM12878_chr19_input.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=/lb/scratch/ -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/lb/scratch/ \
 INPUT=/dev/stdin \
 OUTPUT=alignment/GM12878_chr19_input/GM12878_chr19_input/GM12878_chr19_input.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=135000
bwa_mem_picard_sort_sam.GM12878_chr19_input.ee95a874f5f3c3b7db9c9d3b7fc7e470.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam.GM12878_chr19_CTCF
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.GM12878_chr19_CTCF
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.GM12878_chr19_CTCF.06ef95f9e7b44e7f8afef1a281edb7d3.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'bwa_mem_picard_sort_sam.GM12878_chr19_CTCF.06ef95f9e7b44e7f8afef1a281edb7d3.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/GM12878_chr19_CTCF/GM12878_chr19_CTCF && \
bwa mem  \
  -M -t 1 \
  -R '@RG	ID:GM12878_chr19_CTCF	SM:GM12878_chr19_CTCF	LB:2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/bwa_index/Homo_sapiens.hg19.fa \
  trim/GM12878_chr19_CTCF/GM12878_chr19_CTCF.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=/lb/scratch/ -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/lb/scratch/ \
 INPUT=/dev/stdin \
 OUTPUT=alignment/GM12878_chr19_CTCF/GM12878_chr19_CTCF/GM12878_chr19_CTCF.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=135000
bwa_mem_picard_sort_sam.GM12878_chr19_CTCF.06ef95f9e7b44e7f8afef1a281edb7d3.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam.GM12878_chr19_Rad21
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.GM12878_chr19_Rad21
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.GM12878_chr19_Rad21.ba5e00f2ee96b04ee946dfde0a6054c6.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'bwa_mem_picard_sort_sam.GM12878_chr19_Rad21.ba5e00f2ee96b04ee946dfde0a6054c6.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/GM12878_chr19_Rad21/GM12878_chr19_Rad21 && \
bwa mem  \
  -M -t 1 \
  -R '@RG	ID:GM12878_chr19_Rad21	SM:GM12878_chr19_Rad21	LB:3	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/bwa_index/Homo_sapiens.hg19.fa \
  trim/GM12878_chr19_Rad21/GM12878_chr19_Rad21.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=/lb/scratch/ -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/lb/scratch/ \
 INPUT=/dev/stdin \
 OUTPUT=alignment/GM12878_chr19_Rad21/GM12878_chr19_Rad21/GM12878_chr19_Rad21.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=135000
bwa_mem_picard_sort_sam.GM12878_chr19_Rad21.ba5e00f2ee96b04ee946dfde0a6054c6.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# STEP: samtools_view_filter
#-------------------------------------------------------------------------------
STEP=samtools_view_filter
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter.GM12878_chr19_input
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.GM12878_chr19_input
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.GM12878_chr19_input.4a2db80841bf1afc4640e5d9769679ce.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'samtools_view_filter.GM12878_chr19_input.4a2db80841bf1afc4640e5d9769679ce.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/GM12878_chr19_input/GM12878_chr19_input/GM12878_chr19_input.sorted.bam \
  > alignment/GM12878_chr19_input/GM12878_chr19_input/GM12878_chr19_input.sorted.filtered.bam
samtools_view_filter.GM12878_chr19_input.4a2db80841bf1afc4640e5d9769679ce.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter.GM12878_chr19_CTCF
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.GM12878_chr19_CTCF
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.GM12878_chr19_CTCF.095cc2472648ebbb954471738bbc5ad7.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'samtools_view_filter.GM12878_chr19_CTCF.095cc2472648ebbb954471738bbc5ad7.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/GM12878_chr19_CTCF/GM12878_chr19_CTCF/GM12878_chr19_CTCF.sorted.bam \
  > alignment/GM12878_chr19_CTCF/GM12878_chr19_CTCF/GM12878_chr19_CTCF.sorted.filtered.bam
samtools_view_filter.GM12878_chr19_CTCF.095cc2472648ebbb954471738bbc5ad7.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter.GM12878_chr19_Rad21
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.GM12878_chr19_Rad21
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.GM12878_chr19_Rad21.bf6901d424231c1dffef81d14d3165f9.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'samtools_view_filter.GM12878_chr19_Rad21.bf6901d424231c1dffef81d14d3165f9.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/GM12878_chr19_Rad21/GM12878_chr19_Rad21/GM12878_chr19_Rad21.sorted.bam \
  > alignment/GM12878_chr19_Rad21/GM12878_chr19_Rad21/GM12878_chr19_Rad21.sorted.filtered.bam
samtools_view_filter.GM12878_chr19_Rad21.bf6901d424231c1dffef81d14d3165f9.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_report
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter_report
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter_report.c011bfaadf04f6801111ce6d84c305ae.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'samtools_view_filter_report.c011bfaadf04f6801111ce6d84c305ae.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/pandoc/1.15.2 && \
mkdir -p report && \
pandoc --to=markdown \
  --template /home/normie/Documents/mcdevops/genpipes/bfx/report/ChipSeq.samtools_view_filter.md \
  --variable min_mapq="20" \
  /home/normie/Documents/mcdevops/genpipes/bfx/report/ChipSeq.samtools_view_filter.md \
  > report/ChipSeq.samtools_view_filter.md
samtools_view_filter_report.c011bfaadf04f6801111ce6d84c305ae.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# STEP: picard_merge_sam_files
#-------------------------------------------------------------------------------
STEP=picard_merge_sam_files
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: symlink_readset_sample_bam.GM12878_chr19_input
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.GM12878_chr19_input
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.GM12878_chr19_input.df1f4e3294ce29edfd5eaa992d898508.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'symlink_readset_sample_bam.GM12878_chr19_input.df1f4e3294ce29edfd5eaa992d898508.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
mkdir -p alignment/GM12878_chr19_input && \
ln -s -f GM12878_chr19_input/GM12878_chr19_input.sorted.filtered.bam alignment/GM12878_chr19_input/GM12878_chr19_input.merged.bam
symlink_readset_sample_bam.GM12878_chr19_input.df1f4e3294ce29edfd5eaa992d898508.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# JOB: symlink_readset_sample_bam.GM12878_chr19_CTCF
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.GM12878_chr19_CTCF
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.GM12878_chr19_CTCF.b7e1c933378860fb7a2e2e43a22fa218.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'symlink_readset_sample_bam.GM12878_chr19_CTCF.b7e1c933378860fb7a2e2e43a22fa218.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
mkdir -p alignment/GM12878_chr19_CTCF && \
ln -s -f GM12878_chr19_CTCF/GM12878_chr19_CTCF.sorted.filtered.bam alignment/GM12878_chr19_CTCF/GM12878_chr19_CTCF.merged.bam
symlink_readset_sample_bam.GM12878_chr19_CTCF.b7e1c933378860fb7a2e2e43a22fa218.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# JOB: symlink_readset_sample_bam.GM12878_chr19_Rad21
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.GM12878_chr19_Rad21
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.GM12878_chr19_Rad21.69bc502b8e0a1dc0f5fa6b5f31e6b439.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'symlink_readset_sample_bam.GM12878_chr19_Rad21.69bc502b8e0a1dc0f5fa6b5f31e6b439.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
mkdir -p alignment/GM12878_chr19_Rad21 && \
ln -s -f GM12878_chr19_Rad21/GM12878_chr19_Rad21.sorted.filtered.bam alignment/GM12878_chr19_Rad21/GM12878_chr19_Rad21.merged.bam
symlink_readset_sample_bam.GM12878_chr19_Rad21.69bc502b8e0a1dc0f5fa6b5f31e6b439.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# STEP: picard_mark_duplicates
#-------------------------------------------------------------------------------
STEP=picard_mark_duplicates
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates.GM12878_chr19_input
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.GM12878_chr19_input
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.GM12878_chr19_input.759ba1d0752221000842d310f168a66f.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'picard_mark_duplicates.GM12878_chr19_input.759ba1d0752221000842d310f168a66f.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/lb/scratch/ -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/lb/scratch/ \
 INPUT=alignment/GM12878_chr19_input/GM12878_chr19_input.merged.bam \
 OUTPUT=alignment/GM12878_chr19_input/GM12878_chr19_input.sorted.dup.bam \
 METRICS_FILE=alignment/GM12878_chr19_input/GM12878_chr19_input.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=100000 
picard_mark_duplicates.GM12878_chr19_input.759ba1d0752221000842d310f168a66f.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates.GM12878_chr19_CTCF
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.GM12878_chr19_CTCF
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.GM12878_chr19_CTCF.9d26b919d66225f7a735b92160fd7b3a.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'picard_mark_duplicates.GM12878_chr19_CTCF.9d26b919d66225f7a735b92160fd7b3a.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/lb/scratch/ -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/lb/scratch/ \
 INPUT=alignment/GM12878_chr19_CTCF/GM12878_chr19_CTCF.merged.bam \
 OUTPUT=alignment/GM12878_chr19_CTCF/GM12878_chr19_CTCF.sorted.dup.bam \
 METRICS_FILE=alignment/GM12878_chr19_CTCF/GM12878_chr19_CTCF.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=100000 
picard_mark_duplicates.GM12878_chr19_CTCF.9d26b919d66225f7a735b92160fd7b3a.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates.GM12878_chr19_Rad21
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.GM12878_chr19_Rad21
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.GM12878_chr19_Rad21.067d6bc406ce9d98367f9dea3d3a55a7.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'picard_mark_duplicates.GM12878_chr19_Rad21.067d6bc406ce9d98367f9dea3d3a55a7.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/lb/scratch/ -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/lb/scratch/ \
 INPUT=alignment/GM12878_chr19_Rad21/GM12878_chr19_Rad21.merged.bam \
 OUTPUT=alignment/GM12878_chr19_Rad21/GM12878_chr19_Rad21.sorted.dup.bam \
 METRICS_FILE=alignment/GM12878_chr19_Rad21/GM12878_chr19_Rad21.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=100000 
picard_mark_duplicates.GM12878_chr19_Rad21.067d6bc406ce9d98367f9dea3d3a55a7.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_report
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates_report
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates_report.9d8a9a18a9ff8bc24199a1ff4d4f96fc.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'picard_mark_duplicates_report.9d8a9a18a9ff8bc24199a1ff4d4f96fc.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
mkdir -p report && \
cp \
  /home/normie/Documents/mcdevops/genpipes/bfx/report/ChipSeq.picard_mark_duplicates.md \
  report/ChipSeq.picard_mark_duplicates.md
picard_mark_duplicates_report.9d8a9a18a9ff8bc24199a1ff4d4f96fc.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# STEP: metrics
#-------------------------------------------------------------------------------
STEP=metrics
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: picard_collect_multiple_metrics.GM12878_chr19_input
#-------------------------------------------------------------------------------
JOB_NAME=picard_collect_multiple_metrics.GM12878_chr19_input
JOB_DONE=job_output/metrics/picard_collect_multiple_metrics.GM12878_chr19_input.2b97fe2e2469599585cf84d79d87605b.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'picard_collect_multiple_metrics.GM12878_chr19_input.2b97fe2e2469599585cf84d79d87605b.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 mugqic/R_Bioconductor/3.5.0_3.7 && \
java -Djava.io.tmpdir=/lb/scratch/ -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx4G -jar $PICARD_HOME/picard.jar CollectMultipleMetrics \
 PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT \
 TMP_DIR=/lb/scratch/ \
 REFERENCE_SEQUENCE=/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa \
 INPUT=alignment/GM12878_chr19_input/GM12878_chr19_input.sorted.dup.bam \
 OUTPUT=alignment/GM12878_chr19_input/GM12878_chr19_input.sorted.dup.all.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_collect_multiple_metrics.GM12878_chr19_input.2b97fe2e2469599585cf84d79d87605b.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# JOB: picard_collect_multiple_metrics.GM12878_chr19_CTCF
#-------------------------------------------------------------------------------
JOB_NAME=picard_collect_multiple_metrics.GM12878_chr19_CTCF
JOB_DONE=job_output/metrics/picard_collect_multiple_metrics.GM12878_chr19_CTCF.939483cfb2b050704291205bd5492a9b.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'picard_collect_multiple_metrics.GM12878_chr19_CTCF.939483cfb2b050704291205bd5492a9b.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 mugqic/R_Bioconductor/3.5.0_3.7 && \
java -Djava.io.tmpdir=/lb/scratch/ -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx4G -jar $PICARD_HOME/picard.jar CollectMultipleMetrics \
 PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT \
 TMP_DIR=/lb/scratch/ \
 REFERENCE_SEQUENCE=/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa \
 INPUT=alignment/GM12878_chr19_CTCF/GM12878_chr19_CTCF.sorted.dup.bam \
 OUTPUT=alignment/GM12878_chr19_CTCF/GM12878_chr19_CTCF.sorted.dup.all.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_collect_multiple_metrics.GM12878_chr19_CTCF.939483cfb2b050704291205bd5492a9b.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# JOB: picard_collect_multiple_metrics.GM12878_chr19_Rad21
#-------------------------------------------------------------------------------
JOB_NAME=picard_collect_multiple_metrics.GM12878_chr19_Rad21
JOB_DONE=job_output/metrics/picard_collect_multiple_metrics.GM12878_chr19_Rad21.30f78b3307d8cc8015b545d18b25e7f4.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'picard_collect_multiple_metrics.GM12878_chr19_Rad21.30f78b3307d8cc8015b545d18b25e7f4.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 mugqic/R_Bioconductor/3.5.0_3.7 && \
java -Djava.io.tmpdir=/lb/scratch/ -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx4G -jar $PICARD_HOME/picard.jar CollectMultipleMetrics \
 PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT \
 TMP_DIR=/lb/scratch/ \
 REFERENCE_SEQUENCE=/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa \
 INPUT=alignment/GM12878_chr19_Rad21/GM12878_chr19_Rad21.sorted.dup.bam \
 OUTPUT=alignment/GM12878_chr19_Rad21/GM12878_chr19_Rad21.sorted.dup.all.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_collect_multiple_metrics.GM12878_chr19_Rad21.30f78b3307d8cc8015b545d18b25e7f4.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# JOB: metrics.flagstat
#-------------------------------------------------------------------------------
JOB_NAME=metrics.flagstat
JOB_DONE=job_output/metrics/metrics.flagstat.54d2588a649686fc67293d51f2c2cba8.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'metrics.flagstat.54d2588a649686fc67293d51f2c2cba8.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/samtools/1.3.1 && \
samtools flagstat \
  alignment/GM12878_chr19_input/GM12878_chr19_input.sorted.dup.bam \
  > alignment/GM12878_chr19_input/GM12878_chr19_input.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/GM12878_chr19_CTCF/GM12878_chr19_CTCF.sorted.dup.bam \
  > alignment/GM12878_chr19_CTCF/GM12878_chr19_CTCF.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/GM12878_chr19_Rad21/GM12878_chr19_Rad21.sorted.dup.bam \
  > alignment/GM12878_chr19_Rad21/GM12878_chr19_Rad21.sorted.dup.bam.flagstat
metrics.flagstat.54d2588a649686fc67293d51f2c2cba8.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# JOB: metrics_report
#-------------------------------------------------------------------------------
JOB_NAME=metrics_report
JOB_DONE=job_output/metrics/metrics_report.c80161c14ea84d30f4b76e3c5fc74dd8.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'metrics_report.c80161c14ea84d30f4b76e3c5fc74dd8.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/pandoc/1.15.2 && \
module load mugqic/sambamba/0.7.0 && \
cp /dev/null metrics/SampleMetrics.stats && \
for sample in GM12878_chr19_input GM12878_chr19_CTCF GM12878_chr19_Rad21
do
  flagstat_file=alignment/$sample/$sample.sorted.dup.bam.flagstat
  bam_file=alignment/$sample/$sample.sorted.dup.bam
  supplementarysecondary_alignment=`bc <<< $(grep "secondary" $flagstat_file | sed -e 's/ + [[:digit:]]* secondary.*//')+$(grep "supplementary" $flagstat_file | sed -e 's/ + [[:digit:]]* supplementary.*//')`
  mapped_reads=`bc <<< $(grep "mapped (" $flagstat_file | sed -e 's/ + [[:digit:]]* mapped (.*)//')-$supplementarysecondary_alignment`
  duplicated_reads=`grep "duplicates" $flagstat_file | sed -e 's/ + [[:digit:]]* duplicates$//'`
  duplicated_rate=$(echo "100*$duplicated_reads/$mapped_reads" | bc -l)
  mito_reads=$(sambamba view -c $bam_file chrM)
  mito_rate=$(echo "100*$mito_reads/$mapped_reads" | bc -l)
  echo -e "$sample	$mapped_reads	$duplicated_reads	$duplicated_rate	$mito_reads	$mito_rate" >> metrics/SampleMetrics.stats
done && \
sed -i -e "1 i\Sample	Aligned Filtered Reads #	Duplicate Reads #	Duplicate %	Mitchondrial Reads #	Mitochondrial %" metrics/SampleMetrics.stats && \
mkdir -p report && \
if [[ -f metrics/trimSampleTable.tsv ]]
then
  awk -F "	" 'FNR==NR{trim_line[$1]=$0; surviving[$1]=$3; next}{OFS="	"; if ($1=="Sample") {print trim_line[$1], $2, "Aligned Filtered %", $3, $4, $5, $6} else {print trim_line[$1], $2, $2 / surviving[$1] * 100, $3, $4, $5, $6}}' metrics/trimSampleTable.tsv metrics/SampleMetrics.stats \
  > report/trimMemSampleTable.tsv
else
  cp metrics/SampleMetrics.stats report/trimMemSampleTable.tsv
fi && \
trim_mem_sample_table=`if [[ -f metrics/trimSampleTable.tsv ]] ; then LC_NUMERIC=en_CA awk -F "	" '{OFS="|"; if (NR == 1) {$1 = $1; print $0; print "-----|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:"} else {print $1, sprintf("%\47d", $2), sprintf("%\47d", $3), sprintf("%.1f", $4), sprintf("%\47d", $5), sprintf("%.1f", $6), sprintf("%\47d", $7), sprintf("%.1f", $8), sprintf("%\47d", $9), sprintf("%.1f", $10)}}' report/trimMemSampleTable.tsv ; else LC_NUMERIC=en_CA awk -F "	" '{OFS="|"; if (NR == 1) {$1 = $1; print $0; print "-----|-----:|-----:|-----:|-----:|-----:"} else {print $1, sprintf("%\47d", $2), sprintf("%\47d", $3), sprintf("%.1f", $4), sprintf("%\47d", $5), sprintf("%.1f", $6)}}' report/trimMemSampleTable.tsv ; fi` && \
pandoc --to=markdown \
  --template /home/normie/Documents/mcdevops/genpipes/bfx/report/ChipSeq.metrics.md \
  --variable trim_mem_sample_table="$trim_mem_sample_table" \
  /home/normie/Documents/mcdevops/genpipes/bfx/report/ChipSeq.metrics.md \
  > report/ChipSeq.metrics.md

metrics_report.c80161c14ea84d30f4b76e3c5fc74dd8.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# STEP: homer_make_tag_directory
#-------------------------------------------------------------------------------
STEP=homer_make_tag_directory
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory.GM12878_chr19_input
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.GM12878_chr19_input
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.GM12878_chr19_input.c76ee5271798e3d5a2f1abac1420f693.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'homer_make_tag_directory.GM12878_chr19_input.c76ee5271798e3d5a2f1abac1420f693.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/GM12878_chr19_input \
            alignment/GM12878_chr19_input/GM12878_chr19_input.sorted.dup.bam \
            -genome hg19 \
            -checkGC \
 
homer_make_tag_directory.GM12878_chr19_input.c76ee5271798e3d5a2f1abac1420f693.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory.GM12878_chr19_CTCF
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.GM12878_chr19_CTCF
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.GM12878_chr19_CTCF.96f2ed88162cfce718a31cb84158780c.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'homer_make_tag_directory.GM12878_chr19_CTCF.96f2ed88162cfce718a31cb84158780c.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/GM12878_chr19_CTCF \
            alignment/GM12878_chr19_CTCF/GM12878_chr19_CTCF.sorted.dup.bam \
            -genome hg19 \
            -checkGC \
 
homer_make_tag_directory.GM12878_chr19_CTCF.96f2ed88162cfce718a31cb84158780c.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory.GM12878_chr19_Rad21
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.GM12878_chr19_Rad21
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.GM12878_chr19_Rad21.bde6170523da2d7e9fa19e7368dcecc7.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'homer_make_tag_directory.GM12878_chr19_Rad21.bde6170523da2d7e9fa19e7368dcecc7.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/GM12878_chr19_Rad21 \
            alignment/GM12878_chr19_Rad21/GM12878_chr19_Rad21.sorted.dup.bam \
            -genome hg19 \
            -checkGC \
 
homer_make_tag_directory.GM12878_chr19_Rad21.bde6170523da2d7e9fa19e7368dcecc7.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# STEP: qc_metrics
#-------------------------------------------------------------------------------
STEP=qc_metrics
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: qc_plots_R
#-------------------------------------------------------------------------------
JOB_NAME=qc_plots_R
JOB_DONE=job_output/qc_metrics/qc_plots_R.41437429f9bddc9f5e6282453ed21c1b.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'qc_plots_R.41437429f9bddc9f5e6282453ed21c1b.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/mugqic_tools/2.2.2 mugqic/R_Bioconductor/3.5.0_3.7 && \
mkdir -p graphs && \
Rscript $R_TOOLS/chipSeqGenerateQCMetrics.R \
  ../../../../../../../cvmfs/soft.mugqic/CentOS6/testdata/chipseq/design.chipseq.txt \
  /home/normie/Documents/mcdevops/genpipes/pipelines/chipseq && \
cp /home/normie/Documents/mcdevops/genpipes/bfx/report/ChipSeq.qc_metrics.md report/ChipSeq.qc_metrics.md && \
for sample in GM12878_chr19_input GM12878_chr19_CTCF GM12878_chr19_Rad21
do
  cp --parents graphs/${sample}_QC_Metrics.ps report/
  convert -rotate 90 graphs/${sample}_QC_Metrics.ps report/graphs/${sample}_QC_Metrics.png
  echo -e "----

![QC Metrics for Sample $sample ([download high-res image](graphs/${sample}_QC_Metrics.ps))](graphs/${sample}_QC_Metrics.png)
" \
  >> report/ChipSeq.qc_metrics.md
done
qc_plots_R.41437429f9bddc9f5e6282453ed21c1b.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# STEP: homer_make_ucsc_file
#-------------------------------------------------------------------------------
STEP=homer_make_ucsc_file
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file.GM12878_chr19_input
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.GM12878_chr19_input
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.GM12878_chr19_input.609f59b8c97ed5afe19ea9c5a9ffb8e3.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'homer_make_ucsc_file.GM12878_chr19_input.609f59b8c97ed5afe19ea9c5a9ffb8e3.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 && \
mkdir -p tracks/GM12878_chr19_input && \
makeUCSCfile \
        tags/GM12878_chr19_input > tracks/GM12878_chr19_input/GM12878_chr19_input.ucsc.bedGraph && \
        gzip -c tracks/GM12878_chr19_input/GM12878_chr19_input.ucsc.bedGraph > tracks/GM12878_chr19_input/GM12878_chr19_input.ucsc.bedGraph.gz
homer_make_ucsc_file.GM12878_chr19_input.609f59b8c97ed5afe19ea9c5a9ffb8e3.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_bigWig.GM12878_chr19_input
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file_bigWig.GM12878_chr19_input
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file_bigWig.GM12878_chr19_input.99feeb9509410eef380d1b841320ab54.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'homer_make_ucsc_file_bigWig.GM12878_chr19_input.99feeb9509410eef380d1b841320ab54.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/ucsc/v346 && \
mkdir -p tracks/GM12878_chr19_input/bigWig && \
export TMPDIR=/lb/scratch/ && \
cat tracks/GM12878_chr19_input/GM12878_chr19_input.ucsc.bedGraph | head -n 1 > tracks/GM12878_chr19_input/GM12878_chr19_input.ucsc.bedGraph.head.tmp && \
cat tracks/GM12878_chr19_input/GM12878_chr19_input.ucsc.bedGraph | awk ' NR > 1 ' | sort  --temporary-directory=/lb/scratch/ -k1,1 -k2,2n | \
awk '{if($0 !~ /^[A-W]/) print ""$0; else print $0}' | grep -v "GL\|lambda\|pUC19\|KI\|\KN\|random"  | \
awk '{printf "%s\t%d\t%d\t%4.4g\n", $1,$2,$3,$4}' > tracks/GM12878_chr19_input/GM12878_chr19_input.ucsc.bedGraph.body.tmp && \
cat tracks/GM12878_chr19_input/GM12878_chr19_input.ucsc.bedGraph.head.tmp tracks/GM12878_chr19_input/GM12878_chr19_input.ucsc.bedGraph.body.tmp > tracks/GM12878_chr19_input/GM12878_chr19_input.ucsc.bedGraph.sorted && \
rm tracks/GM12878_chr19_input/GM12878_chr19_input.ucsc.bedGraph.head.tmp tracks/GM12878_chr19_input/GM12878_chr19_input.ucsc.bedGraph.body.tmp && \
bedGraphToBigWig \
  tracks/GM12878_chr19_input/GM12878_chr19_input.ucsc.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  tracks/GM12878_chr19_input/bigWig/GM12878_chr19_input.bw
homer_make_ucsc_file_bigWig.GM12878_chr19_input.99feeb9509410eef380d1b841320ab54.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file.GM12878_chr19_CTCF
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.GM12878_chr19_CTCF
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.GM12878_chr19_CTCF.cfb0143478c57e4b9a6aa47eff4bacfb.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'homer_make_ucsc_file.GM12878_chr19_CTCF.cfb0143478c57e4b9a6aa47eff4bacfb.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 && \
mkdir -p tracks/GM12878_chr19_CTCF && \
makeUCSCfile \
        tags/GM12878_chr19_CTCF > tracks/GM12878_chr19_CTCF/GM12878_chr19_CTCF.ucsc.bedGraph && \
        gzip -c tracks/GM12878_chr19_CTCF/GM12878_chr19_CTCF.ucsc.bedGraph > tracks/GM12878_chr19_CTCF/GM12878_chr19_CTCF.ucsc.bedGraph.gz
homer_make_ucsc_file.GM12878_chr19_CTCF.cfb0143478c57e4b9a6aa47eff4bacfb.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_bigWig.GM12878_chr19_CTCF
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file_bigWig.GM12878_chr19_CTCF
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file_bigWig.GM12878_chr19_CTCF.4d12d26ce31fd036fc9207b7c0be8013.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'homer_make_ucsc_file_bigWig.GM12878_chr19_CTCF.4d12d26ce31fd036fc9207b7c0be8013.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/ucsc/v346 && \
mkdir -p tracks/GM12878_chr19_CTCF/bigWig && \
export TMPDIR=/lb/scratch/ && \
cat tracks/GM12878_chr19_CTCF/GM12878_chr19_CTCF.ucsc.bedGraph | head -n 1 > tracks/GM12878_chr19_CTCF/GM12878_chr19_CTCF.ucsc.bedGraph.head.tmp && \
cat tracks/GM12878_chr19_CTCF/GM12878_chr19_CTCF.ucsc.bedGraph | awk ' NR > 1 ' | sort  --temporary-directory=/lb/scratch/ -k1,1 -k2,2n | \
awk '{if($0 !~ /^[A-W]/) print ""$0; else print $0}' | grep -v "GL\|lambda\|pUC19\|KI\|\KN\|random"  | \
awk '{printf "%s\t%d\t%d\t%4.4g\n", $1,$2,$3,$4}' > tracks/GM12878_chr19_CTCF/GM12878_chr19_CTCF.ucsc.bedGraph.body.tmp && \
cat tracks/GM12878_chr19_CTCF/GM12878_chr19_CTCF.ucsc.bedGraph.head.tmp tracks/GM12878_chr19_CTCF/GM12878_chr19_CTCF.ucsc.bedGraph.body.tmp > tracks/GM12878_chr19_CTCF/GM12878_chr19_CTCF.ucsc.bedGraph.sorted && \
rm tracks/GM12878_chr19_CTCF/GM12878_chr19_CTCF.ucsc.bedGraph.head.tmp tracks/GM12878_chr19_CTCF/GM12878_chr19_CTCF.ucsc.bedGraph.body.tmp && \
bedGraphToBigWig \
  tracks/GM12878_chr19_CTCF/GM12878_chr19_CTCF.ucsc.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  tracks/GM12878_chr19_CTCF/bigWig/GM12878_chr19_CTCF.bw
homer_make_ucsc_file_bigWig.GM12878_chr19_CTCF.4d12d26ce31fd036fc9207b7c0be8013.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file.GM12878_chr19_Rad21
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.GM12878_chr19_Rad21
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.GM12878_chr19_Rad21.5376cfb801be765f1cb7be73beef6010.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'homer_make_ucsc_file.GM12878_chr19_Rad21.5376cfb801be765f1cb7be73beef6010.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 && \
mkdir -p tracks/GM12878_chr19_Rad21 && \
makeUCSCfile \
        tags/GM12878_chr19_Rad21 > tracks/GM12878_chr19_Rad21/GM12878_chr19_Rad21.ucsc.bedGraph && \
        gzip -c tracks/GM12878_chr19_Rad21/GM12878_chr19_Rad21.ucsc.bedGraph > tracks/GM12878_chr19_Rad21/GM12878_chr19_Rad21.ucsc.bedGraph.gz
homer_make_ucsc_file.GM12878_chr19_Rad21.5376cfb801be765f1cb7be73beef6010.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_bigWig.GM12878_chr19_Rad21
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file_bigWig.GM12878_chr19_Rad21
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file_bigWig.GM12878_chr19_Rad21.7f13a88b9795494017eef374ba019292.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'homer_make_ucsc_file_bigWig.GM12878_chr19_Rad21.7f13a88b9795494017eef374ba019292.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/ucsc/v346 && \
mkdir -p tracks/GM12878_chr19_Rad21/bigWig && \
export TMPDIR=/lb/scratch/ && \
cat tracks/GM12878_chr19_Rad21/GM12878_chr19_Rad21.ucsc.bedGraph | head -n 1 > tracks/GM12878_chr19_Rad21/GM12878_chr19_Rad21.ucsc.bedGraph.head.tmp && \
cat tracks/GM12878_chr19_Rad21/GM12878_chr19_Rad21.ucsc.bedGraph | awk ' NR > 1 ' | sort  --temporary-directory=/lb/scratch/ -k1,1 -k2,2n | \
awk '{if($0 !~ /^[A-W]/) print ""$0; else print $0}' | grep -v "GL\|lambda\|pUC19\|KI\|\KN\|random"  | \
awk '{printf "%s\t%d\t%d\t%4.4g\n", $1,$2,$3,$4}' > tracks/GM12878_chr19_Rad21/GM12878_chr19_Rad21.ucsc.bedGraph.body.tmp && \
cat tracks/GM12878_chr19_Rad21/GM12878_chr19_Rad21.ucsc.bedGraph.head.tmp tracks/GM12878_chr19_Rad21/GM12878_chr19_Rad21.ucsc.bedGraph.body.tmp > tracks/GM12878_chr19_Rad21/GM12878_chr19_Rad21.ucsc.bedGraph.sorted && \
rm tracks/GM12878_chr19_Rad21/GM12878_chr19_Rad21.ucsc.bedGraph.head.tmp tracks/GM12878_chr19_Rad21/GM12878_chr19_Rad21.ucsc.bedGraph.body.tmp && \
bedGraphToBigWig \
  tracks/GM12878_chr19_Rad21/GM12878_chr19_Rad21.ucsc.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  tracks/GM12878_chr19_Rad21/bigWig/GM12878_chr19_Rad21.bw
homer_make_ucsc_file_bigWig.GM12878_chr19_Rad21.7f13a88b9795494017eef374ba019292.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_report
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file_report
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file_report.9bb61419202ad8d51bc67fa52720773d.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'homer_make_ucsc_file_report.9bb61419202ad8d51bc67fa52720773d.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
mkdir -p report && \
zip -r report/tracks.zip tracks/*/*.ucsc.bedGraph.gz && \
cp /home/normie/Documents/mcdevops/genpipes/bfx/report/ChipSeq.homer_make_ucsc_file.md report/
homer_make_ucsc_file_report.9bb61419202ad8d51bc67fa52720773d.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# STEP: macs2_callpeak
#-------------------------------------------------------------------------------
STEP=macs2_callpeak
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak.CTCF_Input
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.CTCF_Input
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.CTCF_Input.7d5df0792fb590889ece89a6e75de8e3.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'macs2_callpeak.CTCF_Input.7d5df0792fb590889ece89a6e75de8e3.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/CTCF_Input && \
macs2 callpeak --format BAM --nomodel \
  --tempdir /lb/scratch/ \
  --gsize 2509729011.2 \
  --treatment \
  alignment/GM12878_chr19_CTCF/GM12878_chr19_CTCF.sorted.dup.bam \
  --control \
  alignment/GM12878_chr19_input/GM12878_chr19_input.sorted.dup.bam \
  --name peak_call/CTCF_Input/CTCF_Input \
  >& peak_call/CTCF_Input/CTCF_Input.diag.macs.out
macs2_callpeak.CTCF_Input.7d5df0792fb590889ece89a6e75de8e3.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_bigBed.CTCF_Input
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.CTCF_Input
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.CTCF_Input.1ee61fb47b413e8025e39d7df2592c2f.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'macs2_callpeak_bigBed.CTCF_Input.1ee61fb47b413e8025e39d7df2592c2f.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/CTCF_Input/CTCF_Input_peaks.narrowPeak > peak_call/CTCF_Input/CTCF_Input_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/CTCF_Input/CTCF_Input_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  peak_call/CTCF_Input/CTCF_Input_peaks.narrowPeak.bb
macs2_callpeak_bigBed.CTCF_Input.1ee61fb47b413e8025e39d7df2592c2f.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak.Rad21_Input
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.Rad21_Input
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.Rad21_Input.ab24a644ed068306c57e54f4befa76e4.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'macs2_callpeak.Rad21_Input.ab24a644ed068306c57e54f4befa76e4.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/Rad21_Input && \
macs2 callpeak --format BAM --nomodel \
  --tempdir /lb/scratch/ \
  --gsize 2509729011.2 \
  --treatment \
  alignment/GM12878_chr19_Rad21/GM12878_chr19_Rad21.sorted.dup.bam \
  --control \
  alignment/GM12878_chr19_input/GM12878_chr19_input.sorted.dup.bam \
  --name peak_call/Rad21_Input/Rad21_Input \
  >& peak_call/Rad21_Input/Rad21_Input.diag.macs.out
macs2_callpeak.Rad21_Input.ab24a644ed068306c57e54f4befa76e4.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_bigBed.Rad21_Input
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.Rad21_Input
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.Rad21_Input.1af22083e96b3bf50006fa801f484381.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'macs2_callpeak_bigBed.Rad21_Input.1af22083e96b3bf50006fa801f484381.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/Rad21_Input/Rad21_Input_peaks.narrowPeak > peak_call/Rad21_Input/Rad21_Input_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/Rad21_Input/Rad21_Input_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  peak_call/Rad21_Input/Rad21_Input_peaks.narrowPeak.bb
macs2_callpeak_bigBed.Rad21_Input.1af22083e96b3bf50006fa801f484381.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_report
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_report
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_report.88d19568368a32bbce038010a555c954.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'macs2_callpeak_report.88d19568368a32bbce038010a555c954.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
mkdir -p report && \
cp /home/normie/Documents/mcdevops/genpipes/bfx/report/ChipSeq.macs2_callpeak.md report/ && \
for contrast in CTCF_Input Rad21_Input
do
  cp -a --parents peak_call/$contrast/ report/ && \
  echo -e "* [Peak Calls File for Design $contrast](peak_call/$contrast/${contrast}_peaks.xls)" \
  >> report/ChipSeq.macs2_callpeak.md
done
macs2_callpeak_report.88d19568368a32bbce038010a555c954.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# STEP: homer_annotate_peaks
#-------------------------------------------------------------------------------
STEP=homer_annotate_peaks
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks.CTCF_Input
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.CTCF_Input
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.CTCF_Input.bcba8b663cef7f52108297a4544edfbf.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'homer_annotate_peaks.CTCF_Input.bcba8b663cef7f52108297a4544edfbf.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/mugqic_tools/2.2.2 && \
mkdir -p annotation/CTCF_Input/CTCF_Input && \
annotatePeaks.pl \
    peak_call/CTCF_Input/CTCF_Input_peaks.narrowPeak \
    hg19 \
    -gsize hg19 \
    -cons -CpG \
    -go annotation/CTCF_Input/CTCF_Input \
    -genomeOntology annotation/CTCF_Input/CTCF_Input \
    > annotation/CTCF_Input/CTCF_Input.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/CTCF_Input/CTCF_Input.annotated.csv",
  "annotation/CTCF_Input/CTCF_Input",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.CTCF_Input.bcba8b663cef7f52108297a4544edfbf.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks.Rad21_Input
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.Rad21_Input
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.Rad21_Input.2d8945a7ede8c926145cb860c1d263ba.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'homer_annotate_peaks.Rad21_Input.2d8945a7ede8c926145cb860c1d263ba.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/mugqic_tools/2.2.2 && \
mkdir -p annotation/Rad21_Input/Rad21_Input && \
annotatePeaks.pl \
    peak_call/Rad21_Input/Rad21_Input_peaks.narrowPeak \
    hg19 \
    -gsize hg19 \
    -cons -CpG \
    -go annotation/Rad21_Input/Rad21_Input \
    -genomeOntology annotation/Rad21_Input/Rad21_Input \
    > annotation/Rad21_Input/Rad21_Input.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/Rad21_Input/Rad21_Input.annotated.csv",
  "annotation/Rad21_Input/Rad21_Input",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.Rad21_Input.2d8945a7ede8c926145cb860c1d263ba.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_report
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks_report
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks_report.71c38e7efb27e5b95c9e297fce82ccdf.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'homer_annotate_peaks_report.71c38e7efb27e5b95c9e297fce82ccdf.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
mkdir -p report/annotation/ && \
cp /home/normie/Documents/mcdevops/genpipes/bfx/report/ChipSeq.homer_annotate_peaks.md report && \
for contrast in CTCF_Input Rad21_Input
do
  rsync -avP annotation/$contrast report/annotation/ && \
  echo -e "* [Gene Annotations for Design $contrast](annotation/$contrast/${contrast}.annotated.csv)
* [HOMER Gene Ontology Annotations for Design $contrast](annotation/$contrast/$contrast/geneOntology.html)
* [HOMER Genome Ontology Annotations for Design $contrast](annotation/$contrast/$contrast/GenomeOntology.html)" \
  >> report/ChipSeq.homer_annotate_peaks.md
done
homer_annotate_peaks_report.71c38e7efb27e5b95c9e297fce82ccdf.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# STEP: homer_find_motifs_genome
#-------------------------------------------------------------------------------
STEP=homer_find_motifs_genome
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome.CTCF_Input
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.CTCF_Input
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.CTCF_Input.269fa7a1022424dcb4104ec3af1e9248.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'homer_find_motifs_genome.CTCF_Input.269fa7a1022424dcb4104ec3af1e9248.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/weblogo/3.3 && \
mkdir -p annotation/CTCF_Input/CTCF_Input && \
findMotifsGenome.pl \
  peak_call/CTCF_Input/CTCF_Input_peaks.narrowPeak \
  hg19 \
  annotation/CTCF_Input/CTCF_Input \
  -preparsedDir annotation/CTCF_Input/CTCF_Input/preparsed \
  -p 1
homer_find_motifs_genome.CTCF_Input.269fa7a1022424dcb4104ec3af1e9248.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome.Rad21_Input
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.Rad21_Input
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.Rad21_Input.80d4fdcc2fe1860dd6305f06f3ca5366.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'homer_find_motifs_genome.Rad21_Input.80d4fdcc2fe1860dd6305f06f3ca5366.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/weblogo/3.3 && \
mkdir -p annotation/Rad21_Input/Rad21_Input && \
findMotifsGenome.pl \
  peak_call/Rad21_Input/Rad21_Input_peaks.narrowPeak \
  hg19 \
  annotation/Rad21_Input/Rad21_Input \
  -preparsedDir annotation/Rad21_Input/Rad21_Input/preparsed \
  -p 1
homer_find_motifs_genome.Rad21_Input.80d4fdcc2fe1860dd6305f06f3ca5366.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_report
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome_report
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome_report.47c5b448765201bc91a25ce39be59c9e.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'homer_find_motifs_genome_report.47c5b448765201bc91a25ce39be59c9e.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
    mkdir -p report/annotation/ && \
    cp /home/normie/Documents/mcdevops/genpipes/bfx/report/ChipSeq.homer_find_motifs_genome.md report/ && \
    for contrast in CTCF_Input Rad21_Input
    do
      rsync -avP annotation/$contrast report/annotation/ && \
      echo -e "* [HOMER _De Novo_ Motif Results for Design $contrast](annotation/$contrast/$contrast/homerResults.html)
* [HOMER Known Motif Results for Design $contrast](annotation/$contrast/$contrast/knownResults.html)" \
      >> report/ChipSeq.homer_find_motifs_genome.md
    done
homer_find_motifs_genome_report.47c5b448765201bc91a25ce39be59c9e.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# STEP: annotation_graphs
#-------------------------------------------------------------------------------
STEP=annotation_graphs
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: annotation_graphs
#-------------------------------------------------------------------------------
JOB_NAME=annotation_graphs
JOB_DONE=job_output/annotation_graphs/annotation_graphs.eb11f1cb9fd78e1f08cef1e1f962c084.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'annotation_graphs.eb11f1cb9fd78e1f08cef1e1f962c084.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/mugqic_tools/2.2.2 mugqic/R_Bioconductor/3.5.0_3.7 mugqic/pandoc/1.15.2 && \
mkdir -p graphs && \
Rscript $R_TOOLS/chipSeqgenerateAnnotationGraphs.R \
  ../../../../../../../cvmfs/soft.mugqic/CentOS6/testdata/chipseq/design.chipseq.txt \
  /home/normie/Documents/mcdevops/genpipes/pipelines/chipseq && \
mkdir -p report/annotation/ && \
if [[ -f annotation/peak_stats.csv ]]
then
  cp annotation/peak_stats.csv report/annotation/
peak_stats_table=`LC_NUMERIC=en_CA awk -F "," '{OFS="|"; if (NR == 1) {$1 = $1; print $0; print "-----|-----|-----:|-----:|-----:|-----:|-----:|-----:"} else {print $1, $2,  sprintf("%\47d", $3), $4, sprintf("%\47.1f", $5), sprintf("%\47.1f", $6), sprintf("%\47.1f", $7), sprintf("%\47.1f", $8)}}' annotation/peak_stats.csv`
else
  peak_stats_table=""
fi
pandoc --to=markdown \
  --template /home/normie/Documents/mcdevops/genpipes/bfx/report/ChipSeq.annotation_graphs.md \
  --variable peak_stats_table="$peak_stats_table" \
  --variable proximal_distance="2" \
  --variable distal_distance="10" \
  --variable distance5d_lower="10" \
  --variable distance5d_upper="100" \
  --variable gene_desert_size="100" \
  /home/normie/Documents/mcdevops/genpipes/bfx/report/ChipSeq.annotation_graphs.md \
  > report/ChipSeq.annotation_graphs.md && \
for contrast in CTCF_Input Rad21_Input
do
  cp --parents graphs/${contrast}_Misc_Graphs.ps report/
  convert -rotate 90 graphs/${contrast}_Misc_Graphs.ps report/graphs/${contrast}_Misc_Graphs.png
  echo -e "----

![Annotation Statistics for Design $contrast ([download high-res image](graphs/${contrast}_Misc_Graphs.ps))](graphs/${contrast}_Misc_Graphs.png)
" \
  >> report/ChipSeq.annotation_graphs.md
done
annotation_graphs.eb11f1cb9fd78e1f08cef1e1f962c084.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# STEP: ihec_preprocess_files
#-------------------------------------------------------------------------------
STEP=ihec_preprocess_files
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: ihec_preprocess_symlink.GM12878_chr19_input
#-------------------------------------------------------------------------------
JOB_NAME=ihec_preprocess_symlink.GM12878_chr19_input
JOB_DONE=job_output/ihec_preprocess_files/ihec_preprocess_symlink.GM12878_chr19_input.aabe0bc80e98f20785ffdf1140e5f336.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'ihec_preprocess_symlink.GM12878_chr19_input.aabe0bc80e98f20785ffdf1140e5f336.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
mkdir -p ihec_alignment && \
ln -s -f ../alignment/GM12878_chr19_input/GM12878_chr19_input/GM12878_chr19_input.sorted.bam ihec_alignment/GM12878_chr19_input.merged.bam
ihec_preprocess_symlink.GM12878_chr19_input.aabe0bc80e98f20785ffdf1140e5f336.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# JOB: ihec_preprocess_mark_duplicates.GM12878_chr19_input
#-------------------------------------------------------------------------------
JOB_NAME=ihec_preprocess_mark_duplicates.GM12878_chr19_input
JOB_DONE=job_output/ihec_preprocess_files/ihec_preprocess_mark_duplicates.GM12878_chr19_input.9398c55b6b7cd615e2839e46afffd6a0.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'ihec_preprocess_mark_duplicates.GM12878_chr19_input.9398c55b6b7cd615e2839e46afffd6a0.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
export TMPDIR=/lb/scratch/ && \
java -Djava.io.tmpdir=/lb/scratch/ -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/lb/scratch/ \
 INPUT=ihec_alignment/GM12878_chr19_input.merged.bam \
 OUTPUT=ihec_alignment/GM12878_chr19_input.merged.mdup.bam \
 METRICS_FILE=ihec_alignment/GM12878_chr19_input.merged.mdup.metrics \
 MAX_RECORDS_IN_RAM=100000 
ihec_preprocess_mark_duplicates.GM12878_chr19_input.9398c55b6b7cd615e2839e46afffd6a0.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# JOB: ihec_preprocess_symlink.GM12878_chr19_CTCF
#-------------------------------------------------------------------------------
JOB_NAME=ihec_preprocess_symlink.GM12878_chr19_CTCF
JOB_DONE=job_output/ihec_preprocess_files/ihec_preprocess_symlink.GM12878_chr19_CTCF.0af1d7e6a52d3b6b57ddb26b2d51409c.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'ihec_preprocess_symlink.GM12878_chr19_CTCF.0af1d7e6a52d3b6b57ddb26b2d51409c.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
mkdir -p ihec_alignment && \
ln -s -f ../alignment/GM12878_chr19_CTCF/GM12878_chr19_CTCF/GM12878_chr19_CTCF.sorted.bam ihec_alignment/GM12878_chr19_CTCF.merged.bam
ihec_preprocess_symlink.GM12878_chr19_CTCF.0af1d7e6a52d3b6b57ddb26b2d51409c.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# JOB: ihec_preprocess_mark_duplicates.GM12878_chr19_CTCF
#-------------------------------------------------------------------------------
JOB_NAME=ihec_preprocess_mark_duplicates.GM12878_chr19_CTCF
JOB_DONE=job_output/ihec_preprocess_files/ihec_preprocess_mark_duplicates.GM12878_chr19_CTCF.a61d69b33b73d9c46fe6c2a8c4e7facf.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'ihec_preprocess_mark_duplicates.GM12878_chr19_CTCF.a61d69b33b73d9c46fe6c2a8c4e7facf.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
export TMPDIR=/lb/scratch/ && \
java -Djava.io.tmpdir=/lb/scratch/ -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/lb/scratch/ \
 INPUT=ihec_alignment/GM12878_chr19_CTCF.merged.bam \
 OUTPUT=ihec_alignment/GM12878_chr19_CTCF.merged.mdup.bam \
 METRICS_FILE=ihec_alignment/GM12878_chr19_CTCF.merged.mdup.metrics \
 MAX_RECORDS_IN_RAM=100000 
ihec_preprocess_mark_duplicates.GM12878_chr19_CTCF.a61d69b33b73d9c46fe6c2a8c4e7facf.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# JOB: ihec_preprocess_symlink.GM12878_chr19_Rad21
#-------------------------------------------------------------------------------
JOB_NAME=ihec_preprocess_symlink.GM12878_chr19_Rad21
JOB_DONE=job_output/ihec_preprocess_files/ihec_preprocess_symlink.GM12878_chr19_Rad21.ac280c3896ad23645d9538aa2b342528.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'ihec_preprocess_symlink.GM12878_chr19_Rad21.ac280c3896ad23645d9538aa2b342528.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
mkdir -p ihec_alignment && \
ln -s -f ../alignment/GM12878_chr19_Rad21/GM12878_chr19_Rad21/GM12878_chr19_Rad21.sorted.bam ihec_alignment/GM12878_chr19_Rad21.merged.bam
ihec_preprocess_symlink.GM12878_chr19_Rad21.ac280c3896ad23645d9538aa2b342528.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# JOB: ihec_preprocess_mark_duplicates.GM12878_chr19_Rad21
#-------------------------------------------------------------------------------
JOB_NAME=ihec_preprocess_mark_duplicates.GM12878_chr19_Rad21
JOB_DONE=job_output/ihec_preprocess_files/ihec_preprocess_mark_duplicates.GM12878_chr19_Rad21.b531c519b4a2d90d5aaaf646a10e8640.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'ihec_preprocess_mark_duplicates.GM12878_chr19_Rad21.b531c519b4a2d90d5aaaf646a10e8640.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
export TMPDIR=/lb/scratch/ && \
java -Djava.io.tmpdir=/lb/scratch/ -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/lb/scratch/ \
 INPUT=ihec_alignment/GM12878_chr19_Rad21.merged.bam \
 OUTPUT=ihec_alignment/GM12878_chr19_Rad21.merged.mdup.bam \
 METRICS_FILE=ihec_alignment/GM12878_chr19_Rad21.merged.mdup.metrics \
 MAX_RECORDS_IN_RAM=100000 
ihec_preprocess_mark_duplicates.GM12878_chr19_Rad21.b531c519b4a2d90d5aaaf646a10e8640.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# STEP: run_spp
#-------------------------------------------------------------------------------
STEP=run_spp
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: run_spp.GM12878_chr19_input
#-------------------------------------------------------------------------------
JOB_NAME=run_spp.GM12878_chr19_input
JOB_DONE=job_output/run_spp/run_spp.GM12878_chr19_input.82999c8910dfde33a8d170fa99f0e265.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'run_spp.GM12878_chr19_input.82999c8910dfde33a8d170fa99f0e265.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/samtools/1.3.1 mugqic/mugqic_tools/2.2.2 mugqic/R_Bioconductor/3.5.0_3.7 && \
mkdir -p ihec_metrics && \
Rscript $R_TOOLS/run_spp.R -c=ihec_alignment/GM12878_chr19_input.merged.mdup.bam -savp -out=ihec_metrics/GM12878_chr19_input.crosscor -rf -tmpdir=/lb/scratch/
run_spp.GM12878_chr19_input.82999c8910dfde33a8d170fa99f0e265.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# JOB: run_spp.GM12878_chr19_CTCF
#-------------------------------------------------------------------------------
JOB_NAME=run_spp.GM12878_chr19_CTCF
JOB_DONE=job_output/run_spp/run_spp.GM12878_chr19_CTCF.f6002ce61c3a91b878b9fc299419a359.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'run_spp.GM12878_chr19_CTCF.f6002ce61c3a91b878b9fc299419a359.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/samtools/1.3.1 mugqic/mugqic_tools/2.2.2 mugqic/R_Bioconductor/3.5.0_3.7 && \
mkdir -p ihec_metrics && \
Rscript $R_TOOLS/run_spp.R -c=ihec_alignment/GM12878_chr19_CTCF.merged.mdup.bam -savp -out=ihec_metrics/GM12878_chr19_CTCF.crosscor -rf -tmpdir=/lb/scratch/
run_spp.GM12878_chr19_CTCF.f6002ce61c3a91b878b9fc299419a359.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# JOB: run_spp.GM12878_chr19_Rad21
#-------------------------------------------------------------------------------
JOB_NAME=run_spp.GM12878_chr19_Rad21
JOB_DONE=job_output/run_spp/run_spp.GM12878_chr19_Rad21.66fd80b86c0074942eab730173b5ef7b.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'run_spp.GM12878_chr19_Rad21.66fd80b86c0074942eab730173b5ef7b.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/samtools/1.3.1 mugqic/mugqic_tools/2.2.2 mugqic/R_Bioconductor/3.5.0_3.7 && \
mkdir -p ihec_metrics && \
Rscript $R_TOOLS/run_spp.R -c=ihec_alignment/GM12878_chr19_Rad21.merged.mdup.bam -savp -out=ihec_metrics/GM12878_chr19_Rad21.crosscor -rf -tmpdir=/lb/scratch/
run_spp.GM12878_chr19_Rad21.66fd80b86c0074942eab730173b5ef7b.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# STEP: ihec_metrics
#-------------------------------------------------------------------------------
STEP=ihec_metrics
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: IHEC_chipseq_metrics.GM12878_chr19_Rad21
#-------------------------------------------------------------------------------
JOB_NAME=IHEC_chipseq_metrics.GM12878_chr19_Rad21
JOB_DONE=job_output/ihec_metrics/IHEC_chipseq_metrics.GM12878_chr19_Rad21.479eca519c8650c96d9dfba569a96f77.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'IHEC_chipseq_metrics.GM12878_chr19_Rad21.479eca519c8650c96d9dfba569a96f77.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/mugqic_tools/2.2.2 mugqic/samtools/1.3.1 mugqic/sambamba/0.7.0 mugqic/deepTools/2.5.0.1 && \
mkdir -p ihec_metrics && \
IHEC_chipseq_metrics_max.sh \
    -d ihec_alignment/GM12878_chr19_Rad21.merged.mdup.bam \
    -i ihec_alignment/GM12878_chr19_input.merged.mdup.bam \
    -s GM12878_chr19_Rad21 \
    -j GM12878_chr19_input \
    -t narrow \
    -n 6 \
    -p peak_call/Rad21_Input/Rad21_Input_peaks.narrowPeak \
    -o ihec_metrics \
    -a hg19
IHEC_chipseq_metrics.GM12878_chr19_Rad21.479eca519c8650c96d9dfba569a96f77.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# JOB: IHEC_chipseq_metrics.GM12878_chr19_CTCF
#-------------------------------------------------------------------------------
JOB_NAME=IHEC_chipseq_metrics.GM12878_chr19_CTCF
JOB_DONE=job_output/ihec_metrics/IHEC_chipseq_metrics.GM12878_chr19_CTCF.80c79627471d636ee53b1bdfeeb64a76.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'IHEC_chipseq_metrics.GM12878_chr19_CTCF.80c79627471d636ee53b1bdfeeb64a76.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/mugqic_tools/2.2.2 mugqic/samtools/1.3.1 mugqic/sambamba/0.7.0 mugqic/deepTools/2.5.0.1 && \
mkdir -p ihec_metrics && \
IHEC_chipseq_metrics_max.sh \
    -d ihec_alignment/GM12878_chr19_CTCF.merged.mdup.bam \
    -i ihec_alignment/GM12878_chr19_input.merged.mdup.bam \
    -s GM12878_chr19_CTCF \
    -j GM12878_chr19_input \
    -t narrow \
    -n 6 \
    -p peak_call/CTCF_Input/CTCF_Input_peaks.narrowPeak \
    -o ihec_metrics \
    -a hg19
IHEC_chipseq_metrics.GM12878_chr19_CTCF.80c79627471d636ee53b1bdfeeb64a76.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# JOB: merge_ihec_metrics
#-------------------------------------------------------------------------------
JOB_NAME=merge_ihec_metrics
JOB_DONE=job_output/ihec_metrics/merge_ihec_metrics.000e8c386ab367c978919f619af3026c.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'merge_ihec_metrics.000e8c386ab367c978919f619af3026c.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
cp /dev/null ihec_metrics/IHEC_metrics_AllSamples.tsv && \
for sample in ihec_metrics/IHEC_metrics_chipseq_GM12878_chr19_Rad21.txt ihec_metrics/IHEC_metrics_chipseq_GM12878_chr19_Rad21.txt
do
    header=$(head -n 1 $sample)
    tail -n 1 $sample >> ihec_metrics/IHEC_metrics_AllSamples.tsv
done && \
sed -i -e "1 i\\$header" ihec_metrics/IHEC_metrics_AllSamples.tsv
merge_ihec_metrics.000e8c386ab367c978919f619af3026c.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# JOB: merge_ihec_metrics_report
#-------------------------------------------------------------------------------
JOB_NAME=merge_ihec_metrics_report
JOB_DONE=job_output/ihec_metrics/merge_ihec_metrics_report.eac870efe9af55ff37c6d36e329ef521.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'merge_ihec_metrics_report.eac870efe9af55ff37c6d36e329ef521.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/pandoc/1.15.2 && \
mkdir -p report && \
cp ihec_metrics/IHEC_metrics_AllSamples.tsv report/IHEC_metrics_AllSamples.tsv && \
pandoc --to=markdown \
  --template /home/normie/Documents/mcdevops/genpipes/bfx/report/ChipSeq.ihec_metrics.md \
  --variable ihec_metrics_merged_table="IHEC_metrics_AllSamples.tsv" \
  /home/normie/Documents/mcdevops/genpipes/bfx/report/ChipSeq.ihec_metrics.md \
  > report/ChipSeq.ihec_metrics.md
merge_ihec_metrics_report.eac870efe9af55ff37c6d36e329ef521.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# STEP: multiqc_report
#-------------------------------------------------------------------------------
STEP=multiqc_report
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: multiqc_report
#-------------------------------------------------------------------------------
JOB_NAME=multiqc_report
JOB_DONE=job_output/multiqc_report/multiqc_report.e70cc33c6b0d4316897abea53022654f.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'multiqc_report.e70cc33c6b0d4316897abea53022654f.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/MultiQC/v1.6 && \
export MULTIQC_CONFIG_PATH=/cvmfs/soft.mugqic/CentOS6/software/genpipes/genpipes-3.1.5/bfx/report/multiqc_reports/chipseq_multiqc_config.yaml && \
    multiqc .
multiqc_report.e70cc33c6b0d4316897abea53022654f.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# STEP: cram_output
#-------------------------------------------------------------------------------
STEP=cram_output
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: cram_output.GM12878_chr19_input
#-------------------------------------------------------------------------------
JOB_NAME=cram_output.GM12878_chr19_input
JOB_DONE=job_output/cram_output/cram_output.GM12878_chr19_input.64af3cb825bb5eee4a43d0c13bcbaf74.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'cram_output.GM12878_chr19_input.64af3cb825bb5eee4a43d0c13bcbaf74.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/samtools/1.3.1 && \
samtools view -h -T $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa -C \
  alignment/GM12878_chr19_input/GM12878_chr19_input.sorted.dup.bam \
  > alignment/GM12878_chr19_input/GM12878_chr19_input.sorted.dup.cram
cram_output.GM12878_chr19_input.64af3cb825bb5eee4a43d0c13bcbaf74.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# JOB: cram_output.GM12878_chr19_CTCF
#-------------------------------------------------------------------------------
JOB_NAME=cram_output.GM12878_chr19_CTCF
JOB_DONE=job_output/cram_output/cram_output.GM12878_chr19_CTCF.bc7018b7b6b79669d04f19657ec307a6.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'cram_output.GM12878_chr19_CTCF.bc7018b7b6b79669d04f19657ec307a6.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/samtools/1.3.1 && \
samtools view -h -T $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa -C \
  alignment/GM12878_chr19_CTCF/GM12878_chr19_CTCF.sorted.dup.bam \
  > alignment/GM12878_chr19_CTCF/GM12878_chr19_CTCF.sorted.dup.cram
cram_output.GM12878_chr19_CTCF.bc7018b7b6b79669d04f19657ec307a6.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# JOB: cram_output.GM12878_chr19_Rad21
#-------------------------------------------------------------------------------
JOB_NAME=cram_output.GM12878_chr19_Rad21
JOB_DONE=job_output/cram_output/cram_output.GM12878_chr19_Rad21.2605f2e390f710da81a63f5f674755cc.mugqic.done
COMMAND=$JOB_OUTPUT_DIR/$STEP/${JOB_NAME}_$TIMESTAMP.sh
cat << 'cram_output.GM12878_chr19_Rad21.2605f2e390f710da81a63f5f674755cc.mugqic.done' > $COMMAND
#!/bin/bash
set -eu -o pipefail
module purge && \
module load mugqic/samtools/1.3.1 && \
samtools view -h -T $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa -C \
  alignment/GM12878_chr19_Rad21/GM12878_chr19_Rad21.sorted.dup.bam \
  > alignment/GM12878_chr19_Rad21/GM12878_chr19_Rad21.sorted.dup.cram
cram_output.GM12878_chr19_Rad21.2605f2e390f710da81a63f5f674755cc.mugqic.done
chmod 755 $COMMAND
printf "\n$SEPARATOR_LINE\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \
rm -f $JOB_DONE 
 $COMMAND
MUGQIC_STATE=$?
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE

if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi


#-------------------------------------------------------------------------------
# Call home with pipeline statistics
#-------------------------------------------------------------------------------
LOG_MD5=$(echo $USER-'127.0.0.1-ChipSeq-GM12878_chr19_Rad21.GM12878_chr19_Rad21,GM12878_chr19_CTCF.GM12878_chr19_CTCF,GM12878_chr19_input.GM12878_chr19_input' | md5sum | awk '{ print $1 }')
wget "http://mugqic.hpc.mcgill.ca/cgi-bin/pipeline.cgi?hostname=FOXHOUND&ip=127.0.0.1&pipeline=ChipSeq&steps=picard_sam_to_fastq,trimmomatic,merge_trimmomatic_stats,bwa_mem_picard_sort_sam,samtools_view_filter,picard_merge_sam_files,picard_mark_duplicates,metrics,homer_make_tag_directory,qc_metrics,homer_make_ucsc_file,macs2_callpeak,homer_annotate_peaks,homer_find_motifs_genome,annotation_graphs,ihec_preprocess_files,run_spp,ihec_metrics,multiqc_report,cram_output&samples=3&md5=$LOG_MD5" --quiet --output-document=/dev/null

