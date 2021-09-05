conda env create -n gatk -f gatkcondaenv.yml

#get all the files
#bwa index -p assembled_contigs -a is assembled_contigs.fasta

#to do
#GATK
#R package

projectname="210813_miseq_AR"

inDir="/share/ScratchGeneral/nenbar/projects/Anthony/raw_data/210813_miseq_AR/fastq/fastq"
resultsDir="/share/ScratchGeneral/nenbar/projects/Anthony/results"
mkdir -p $resultsDir

annotationDir="/share/ScratchGeneral/nenbar/projects/Anthony/annotation"
index=$annotationDir/"assembled_contigs"
logDir="/share/ScratchGeneral/nenbar/projects/Anthony/logs"
mkdir -p $logDir

cutadapt_minlen=30

#trimming
trimDir="$resultsDir/"$projectname".trimgalore"
mkdir -p $trimDir
bwaDir="$resultsDir/"$projectname".bwa"
mkdir -p $bwaDir

FWD="CTGGTGATTCCCTGTGACG"
REV="AATGTCTCCACCGCGTTC"
REVrc="GAACGCGGTGGAGACATT"
e1=0.2 #accept 20% error rate in primer sequence
e2=0.2
cutadapt_minlen=30

#load files
files=`ls $inDir/*R1_001.fastq.gz`
files=( $files )

for file in ${files[@]}; do 
  echo $file; 
  #file="/share/ScratchGeneral/nenbar/projects/Anthony/raw_data/210813_miseq_AR/fastq/fastq/7111_S29_L001_R1_001.fastq.gz"

  sampleName=`basename $file | sed s/_L001_R1_001.fastq.gz//`
  echo $sampleName

#  file1=$inDir/$sampleName"_L001_R1_001.fastq.gz"
#  file2=$inDir/$sampleName"_L001_R2_001.fastq.gz"
#  file2rev=$inDir/$sampleName"_L001_R2_001.revcomp.fastq"
#  fileTrim1=$trimDir/$sampleName"_L001_R1_001_val_1.fq.gz"
#  fileTrim2=$trimDir/$sampleName"_L001_R2_001.revcomp_val_2.fq.gz"
  bamFile="$bwaDir/$sampleName.bam"
  sortedBamFile="$bwaDir/$sampleName.sorted"
  samFile="$bwaDir/$sampleName.sam"

  rev_command="vsearch --fastx_revcomp $file2 --fastqout $file2rev"
  gz_command="gzip -f $file2rev"
  trim_line="trim_galore "$file1" "$file2rev".gz -a $FWD -a2 $REVrc --gzip --paired -o $trimDir"
  mapping_line="bwa mem $index $fileTrim1 $fileTrim2 | samtools view -S -b - -o $bwaDir/$sampleName.bam"
  sort_line="samtools sort $bwaDir/$sampleName.bam $sortedBamFile"
  index_line="samtools index $bwaDir/$sampleName.sorted.bam"
  clean_line="rm -f $bwaDir/$sampleName.bam;rm -f $outDir/$sampleName.sam; rm -f "$file2rev".gz"

  qsub -b y -wd $logDir -j y -N rev$sampleName -R y -pe smp 1 -V $rev_command
  qsub -b y -hold_jid rev$sampleName -wd $logDir -j y -N gz$sampleName -R y -pe smp 1 -V $gz_command
  qsub -b y -hold_jid gz$sampleName -wd $logDir -j y -N trim$sampleName -R y -pe smp 1 -V $trim_line
  qsub -b y -hold_jid trim$sampleName -wd $logDir -j y -N map$sampleName -R y -pe smp 1 -V $mapping_line
  qsub -b y -hold_jid map$sampleName -wd $logDir -j y -N sort$sampleName -R y -pe smp 1 -V $sort_line
  qsub -b y -hold_jid sort$sampleName -wd $logDir -j y -N index$sampleName -R y -pe smp 1 -V $index_line
  qsub -b y -hold_jid index$sampleName -wd $logDir -j y -N clean$sampleName -R y -pe smp 1 -V $clean_line

done
#map files


#IGV