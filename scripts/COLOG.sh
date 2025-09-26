#!/bin/bash

set -euC

THREADS=4

PICARD=$CONDA_PREFIX/share/picard-2.10.10/picard.jar
SIMSEQDIR=$CONDA_PREFIX/share/SimSeq
SIMSEQJAR=$CONDA_PREFIX/share/SimSeq/SimSeqNBProject/store/SimSeq.jar
PLOTDEPTHPY=$CONDA_PREFIX/share/Alt_nCov2019_primers/tools/plot_depth.py
TRIMARTICPY=$CONDA_PREFIX/share/Alt_nCov2019_primers/tools/trim_primers/trim_primer_parts.py

FASTAVERTICAL=$CONDA_PREFIX/bin/fasta_to_vertical-type.pl
HVCHANGE=$CONDA_PREFIX/bin/HV-change.pl
ONELINEFASTA=$CONDA_PREFIX/bin/one-line-FASTA.pl
FASTQONELINE=$CONDA_PREFIX/bin/fastq_to_one-line-1.pl
FASTA60FORMAT=$CONDA_PREFIX/bin/fasta-60bp-return-format.pl

ARTICREGIONDB=$CONDA_PREFIX/opt/reference/MN908947_55-29835.fasta
REFPATH=$CONDA_PREFIX/opt/reference/nextstrain/sars-cov-2/wuhan-hu-1/orfs
TRIMPRIMERBED=$CONDA_PREFIX/opt/reference/nCovid-19-primerN7.mod.bed
adapter=$CONDA_PREFIX/opt/reference/Illumina-adapter.fasta

while [ $# -gt 1 ] ; do
  case "$1" in
    -p | --primer  ) TRIMPRIMERBED=$2
                     shift
                     shift ;;
    -t | --threads ) THREADS=$2
                     shift
                     shift ;;
                 * ) echo "ERROR : illegal option"
                     echo "you could overwrite primer bed file full path using -p/--primer option, $TRIMPRIMERBED by default"
                     echo "you could overwrite thread number using -t/--threads option, $THREADS by default"
                     exit ;;
  esac
done

if [ $# -lt 1 ]; then
  echo "ERROR : you need to supply output directory path"
  exit 1
fi

DIR=$1

if [ -s $DIR/result-summary.txt ]; then
  echo "$DIR already finished."
  exit 0
fi

NEXTCLADE_DATA=$REFPATH
REFSEQ=$NEXTCLADE_DATA/reference.fasta

set +e
err=$(bedtools getfasta -fi $REFSEQ -bed $TRIMPRIMERBED 2>&1 > /dev/null)
set -e
if [ ${#err} -gt 0 ];then
  echo "$err"
  exit 1
fi

pwd
echo $DIR
ls -lh $DIR

cd $DIR

ls -lh

echo -e "\nINFO `date`   SARS-CoV-2 analysis start in `hostname`.\n"

export TGZPATH=`pwd`
export SAMPLENAME=`basename $TGZPATH`

cp $REFSEQ ./
REFSEQ=${REFSEQ##*/}

samtools faidx $REFSEQ
bwa index $REFSEQ
formatdb -i $REFSEQ -p F
if [ ! -f ${REFSEQ%.fasta}.dict ]; then
  gatk CreateSequenceDictionary --REFERENCE $REFSEQ --OUTPUT ${REFSEQ%.fasta}.dict
fi
GNMLEN=genome-len.txt
cut -f1,2 ${REFSEQ}.fai > $GNMLEN
TOTALLEN=`seqkit stat -T $REFSEQ | cut -f 5 | tail -n 1`

set +e
fq1=$(ls *_R1*.fastq.gz)
fq2=$(ls *_R2*.fastq.gz)
fastafile=$(ls *.fasta)
set -e

fastqfile1=`echo $fq1 | awk '{print $1}'`
fastqfile2=`echo $fq2 | awk '{print $1}'`


if [ "${fastqfile1}" != "" -a "${fastqfile2}" != "" ]; then

   echo -e "Fastq data input" > fastq-data.txt
   echo -e "\nINFO `date`    Fastq data, read trimming start.\n"

#########################
##### read trimming #####
#########################

   zcat *_R1*.fastq.gz > read1.fastq
   zcat *_R2*.fastq.gz > read2.fastq

   read1fastq=./read1.fastq
   read2fastq=./read2.fastq

   fastq-mcf $adapter -t 0 -l 40 -q 15 $read1fastq $read2fastq -o PE1_header-change-trim-q15-adapter-minLen40.fastq -o PE2_header-change-trim-q15-adapter-minLen40.fastq 

   sickle pe -f PE1_header-change-trim-q15-adapter-minLen40.fastq -r PE2_header-change-trim-q15-adapter-minLen40.fastq -t sanger -l 40 -o R1_trim-sickle.fastq -p R2_trim-sickle.fastq -s single.fastq

   if [ `wc -l R1_trim-sickle.fastq | awk '{print $1}'` -ge 1000000 ]; then
      echo -e "#####  CAUTION! Large read sample (>=250,000 reads/file)! This sample is downsampled by seqkit program.  #####"
      seqkit sample -s 11 -n 250000 R1_trim-sickle.fastq > R1_trim-sickle-sample-by-nubmer.fastq
      seqkit sample -s 11 -n 250000 R2_trim-sickle.fastq > R2_trim-sickle-sample-by-nubmer.fastq
      mv R1_trim-sickle-sample-by-nubmer.fastq R1_trim-sickle.fastq
      mv R2_trim-sickle-sample-by-nubmer.fastq R2_trim-sickle.fastq
   fi

   bwa mem -t 4 $REFSEQ R1_trim-sickle.fastq R2_trim-sickle.fastq | \
      python $TRIMARTICPY $TRIMPRIMERBED trim-ARTIC-primer_read1.fastq trim-ARTIC-primer_read2.fastq

   $FASTQONELINE trim-ARTIC-primer_read1.fastq | awk -F'\t' '{print $1"\t"substr($2,1,length($2)-5)"\t"$3"\t"substr($4,1,length($4)-5)}' > 1.fq
   $FASTQONELINE trim-ARTIC-primer_read2.fastq | awk -F'\t' '{print $1"\t"substr($2,1,length($2)-5)"\t"$3"\t"substr($4,1,length($4)-5)}' > 2.fq
   paste 1.fq 2.fq | awk -F'\t' 'length($2) >= 50 && length($6) >= 50 {print $1"\n"$2"\n"$3"\n"$4 > "trim-ARTIC-primer_read1_over50mer.fastq"} ; length($2) >= 50 && length($6) >= 50 {print $5"\n"$6"\n"$7"\n"$8 > "trim-ARTIC-primer_read2_over50mer.fastq"}' 

   rm -rf PE1_header-change-trim-q15-adapter-minLen40.fastq PE2_header-change-trim-q15-adapter-minLen40.fastq \
       read1.fastq read2.fastq trim-ARTIC-primer_read1.fastq trim-ARTIC-primer_read2.fastq 1.fq 2.fq

   for i in *fastq; do gzip $i ; done

   trimread1fastq=trim-ARTIC-primer_read1_over50mer.fastq.gz
   trimread2fastq=trim-ARTIC-primer_read2_over50mer.fastq.gz

   echo -e "\nINFO `date`    read trimming end.\n"

elif [ "${fastafile}" != "" ]; then

   echo -e "Fasta data input" > fasta-data.txt
   echo -e "\nINFO `date`    Fasta data, simulated read making start.\n"

############################################
##### making simulate read with simseq #####
############################################

N300=`perl -E 'say "N" x 300'`
perl -pse 's/>.+/$N300/g' -- -N300=$N300 *.fasta | $ONELINEFASTA | perl -pe 's/^/>concatenate\n/' | perl -pe 's/-//g' | $FASTA60FORMAT > concatenate.fa

READCOUNT=`perl -pe 's/\n//g' concatenate.fa | wc -m | awk '{OFMT = "%.10f"} {print int($1/3)}'`

java -jar -Xmx2048m $SIMSEQJAR -1 150 -2 150 --error  $SIMSEQDIR/model_bwa_mapping.txt --error2 \
     $SIMSEQDIR/model_bwa_mapping.txt --insert_size 200 --insert_stdev 0 --read_number $READCOUNT --read_prefix \
     read_ --reference concatenate.fa --duplicate_probability 0 --out out.sam

rm concatenate.fa

gatk --java-options "-Xmx2048m" SamToFastq --INPUT out.sam --FASTQ R1_trim-sickle.fastq --SECOND_END_FASTQ R2_trim-sickle.fastq \
     --INCLUDE_NON_PF_READS True --VALIDATION_STRINGENCY SILENT

rm out.sam

for i in *fastq; do gzip $i ; done

trimread1fastq=R1_trim-sickle.fastq.gz
trimread2fastq=R2_trim-sickle.fastq.gz

   echo -e "\nAMiGA-INFO `date`    simulated read making end.\n"

else

echo -e "No fastq or fasta file ..." > error.txt

exit

fi

############################
##### mapping analysis #####
############################

bwa mem -t $THREADS -k 31 -a -M $REFSEQ $trimread1fastq $trimread2fastq -R "@RG\tID:01\tSM:$SAMPLENAME\tPL:Illumina" > out12_mem.sam
gatk --java-options "-XX:ParallelGCThreads=$THREADS" SortSam --INPUT out12_mem.sam --OUTPUT out12_mem.sorted.bam --SORT_ORDER coordinate
gatk --java-options "-XX:ParallelGCThreads=$THREADS" BuildBamIndex --INPUT out12_mem.sorted.bam
# In SARS-CoV-2 genome analysis, the excluding PCR duplicate step is ignored because there is no PCR amplification step in the library construction.
rm out12_mem.sam
gatk3 -XX:ParallelGCThreads=$THREADS -T RealignerTargetCreator -R $REFSEQ -I out12_mem.sorted.bam -o realignment_targets.list
gatk3 -XX:ParallelGCThreads=$THREADS -T IndelRealigner -R $REFSEQ -I out12_mem.sorted.bam -targetIntervals realignment_targets.list -o realigned_reads.bam --filter_bases_not_stored 
mv realigned_reads.bai realigned_reads.bam.bai

gatk --java-options "-Xmx4G -XX:ParallelGCThreads=$THREADS" HaplotypeCaller -R $REFSEQ --emit-ref-confidence GVCF -I realigned_reads.bam -O raw_variants.gvcf
gatk --java-options "-Xmx4G -XX:ParallelGCThreads=$THREADS" GenotypeGVCFs -R $REFSEQ -V raw_variants.gvcf -O raw_variants.vcf
rm raw_variants.gvcf

#################################
##### Indel check (ref seq) #####
#################################

gatk --java-options "-Xmx4G -XX:ParallelGCThreads=$THREADS" SelectVariants -R $REFSEQ -V raw_variants.vcf --select-type-to-include INDEL -O raw_indels.vcf
gatk --java-options "-Xmx4G -XX:ParallelGCThreads=$THREADS" VariantFiltration -R $REFSEQ -V raw_indels.vcf -O GATK-indel-filtered-ref.vcf \
                         -filter "QD < 2.0" --filter-name "QD2"       \
                         -filter "QUAL < 30.0" --filter-name "QUAL30" \
                         -filter "FS > 200.0" --filter-name "FS200"   \
                         -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20"
rm raw_indels.vcf

cp GATK-indel-filtered-ref.vcf GATK-indel-filtered-ref.vcf-before-modi

set +C
awk -F'\t' '$7~/PASS/ && $10~/^1[\/\|]1/ {print $0}' GATK-indel-filtered-ref.vcf-before-modi >> GATK-indel-filtered-ref.vcf
set -C

set +e
diff GATK-indel-filtered-ref.vcf GATK-indel-filtered-ref.vcf-before-modi > GATK-diff-check.txt
set -e

if [ -s GATK-diff-check.txt ]; then
     rm GATK-diff-check.txt
   else
     rm GATK-diff-check.txt \
        GATK-indel-filtered-ref.vcf-before-modi
fi

##################################
##### VarScan SNV extraction #####  
##################################

gatk --java-options "-Xmx4G -XX:ParallelGCThreads=$THREADS" SelectVariants -R $REFSEQ -V raw_variants.vcf --select-type-to-include SNP -O raw_snps.vcf
gatk --java-options "-Xmx4G -XX:ParallelGCThreads=$THREADS" VariantFiltration -R $REFSEQ -V raw_snps.vcf -O GATK-snps-filtered-ref.vcf \
                         -filter "QD < 2.0" --filter-name "QD2"       \
                         -filter "QUAL < 30.0" --filter-name "QUAL30" \
                         -filter "SOR > 4.0" --filter-name "SOR4"     \
                         -filter "FS > 60.0" --filter-name "FS60"     \
                         -filter "MQ < 40.0" --filter-name "MQ40"     \
                         -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
                         -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"
rm raw_variants.vcf raw_variants.gvcf.idx raw_snps.vcf
grep -v "^#" GATK-snps-filtered-ref.vcf | cut -f 1,2 > GATK-snps-filtered-list.txt

samtools mpileup -B -A -f $REFSEQ realigned_reads.bam | varscan mpileup2snp --strand-filter 0 --min-freq-for-hom 0.6 --output-vcf | perl -pe s/Sample1/$SAMPLENAME/ > VarScan-ref.vcf

samtools mpileup -B -A -f $REFSEQ realigned_reads.bam | varscan mpileup2snp --strand-filter 0 --min-freq-for-hom 0.6 > VarScan.default-format.txt
grep "^##" VarScan-ref.vcf > VarScan-ref.vcf.header1
echo '##FORMAT=<ID=CONS,Number=1,Type=Integer,Description="Consensus sequence with VarScan basecalling">' > VarScan-ref.vcf.header2
grep "^#CHROM" VarScan-ref.vcf > VarScan-ref.vcf.header3
grep -v "^#" VarScan-ref.vcf | cut -f 1-10 | perl -pe 's/GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR/GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR:CONS/' > VarScan-ref.vcf.data
cut -f 5 VarScan.default-format.txt | grep -v "Cons:Cov:Reads1:Reads2:Freq:P-value" | perl -pe 's/^/\|/' | awk -F':' '{print $1}' > VarScan.default-format.txt.data
paste VarScan-ref.vcf.data VarScan.default-format.txt.data | perl -pe 's/\t\|/:/' > VarScan.format-merge.txt
cat VarScan-ref.vcf.header1 VarScan-ref.vcf.header2 VarScan-ref.vcf.header3 VarScan.format-merge.txt > VarScan-modified.vcf     # varscanは、75%以上はMixと判定していない
rm VarScan-ref.vcf.header* VarScan-ref.vcf.data VarScan.default-format.txt.data VarScan.format-merge.txt VarScan-ref.vcf VarScan.default-format.txt
igvtools index VarScan-modified.vcf

grep "^#" VarScan-modified.vcf > VarScan-modified-rm-Mix-allele.vcf

while read seq posi ; do \
   awk -F'\t' -v C=$seq -v P=$posi '$1==C && $2==P {print $0}' \
   VarScan-modified.vcf ; done < GATK-snps-filtered-list.txt > SNV-VarScan-GATK.txt

awk -F'\t' '$1~/^#/ || $10~/^1\/1/ {print $0}' SNV-VarScan-GATK.txt >> VarScan-modified-rm-Mix-allele.vcf

igvtools index VarScan-modified-rm-Mix-allele.vcf

rm GATK-snps-filtered-list.txt igv.log SNV-VarScan-GATK.txt GATK-snps-filtered-ref.vc*

############################################################
##### construction of consensus sequence with SNV data #####
#####  for pangolin                                    #####
############################################################

export PATH=/opt/apache/miniconda3/bin:/opt/bio/bin:$PATH

bgzip VarScan-modified-rm-Mix-allele.vcf
tabix -p vcf VarScan-modified-rm-Mix-allele.vcf.gz
export TGZPATH=`pwd`
export SAMPLENAME=`basename $TGZPATH`
cat $REFSEQ | vcf-consensus VarScan-modified-rm-Mix-allele.vcf.gz > vcf-consensus-non-concatenate.fasta
$ONELINEFASTA vcf-consensus-non-concatenate.fasta | grep ">" | perl -pe 's/\%+/\n/g' | grep -v ">" > vcf-consensus-concatenate.sequence
echo -e ">${SAMPLENAME}" > concatenatevcf-header.txt
cat concatenatevcf-header.txt vcf-consensus-concatenate.sequence > vcf-consensus.fasta
rm concatenatevcf-header.txt vcf-consensus-concatenate.sequence vcf-consensus-non-concatenate.fasta
gzip -d VarScan-modified-rm-Mix-allele.vcf.gz
rm VarScan-modified-rm-Mix-allele.vcf.gz.tbi

pangolin -t $THREADS vcf-consensus.fasta

###########################
##### coverage search #####
###########################

bedtools genomecov -d -split -ibam realigned_reads.bam > BED-depth.txt

echo -e 'mean\tmin\tq1\tmedian\tq3\tmax\tiqr\tsstdev\tjarque' > read-depth-stat.txt
awk -F'\t' '$2>=99 && $2<=29796 {print $3}' BED-depth.txt | \
   datamash mean 1 min 1 q1 1 median 1 q3 1 max 1 iqr 1 sstdev 1 jarque 1 >> read-depth-stat.txt
READ_DEPTH_MAX=`awk -F'\t' '$2>=99 && $2<=29796 {print $3}' BED-depth.txt | datamash median 1`

echo -e read depth \(median\)"\t"$READ_DEPTH_MAX > read_depth_info.txt
RAED_DEPTH_CUTOFF=10  # for SARS-CoV-2 mapping analysis
echo -e read depth cutoff"\t"$RAED_DEPTH_CUTOFF >> read_depth_info.txt

awk -F'\t' -v S="${RAED_DEPTH_CUTOFF}" '{ if ($3>=S) print $1"\t"$2"\t"$3"\tPASS_(>=_"S")" ; else print $1"\t"$2"\t"$3"\tFAIL_(<_"S")" }' BED-depth.txt | sort -k1,1 -k2,2n > BED-depth-check.txt

cut -f1,4 BED-depth-check.txt | grep "FAIL" | uniq -c | perl -pe 's/^ +//' | awk '{print $2"\t"$1}' | sort -k 1,1 -k 2,2n > DEPTH_FAIL_POS.txt
sort -k 1,1 -k 2,2n $GNMLEN | perl -pe 's/ /\t/g' > TOTAL_POS_COUNT.txt
echo -e "SeqID\tCovrage (%)\tRef length (bp)\tCovered length (bp)" > Coverage-info.txt
join -1 1 -2 1 -a 1 TOTAL_POS_COUNT.txt DEPTH_FAIL_POS.txt | awk '{a+=$2; b+=$3; print $1"\t"($2-$3)/$2*100"%\t"$2"\t"$2-$3} END {print "Total\t"(a-b)/a*100"%\t"a"\t"a-b}' >> Coverage-info.txt
rm DEPTH_FAIL_POS.txt TOTAL_POS_COUNT.txt
TOTALCOVLEN=`awk -F'\t' 'END{print $4}' Coverage-info.txt`

rm BED-depth.txt

ALLELE_SITE=`awk -F'\t' '$7~/PASS/ {print $(NF)}' VarScan-modified.vcf | awk -F':' -v S="${RAED_DEPTH_CUTOFF}" '{ if ($3>=S) print $7}' | perl -pe 's/\%/ /g' | sort -n | uniq -c | perl -pe 's/^ +//' | awk '$2>=80 {m+=$1} END{print m;}' | awk '{ if ($1>0) print $1 ; else print 0}'`
TOTAL_SITE=`awk -F'\t' '$7~/PASS/ {print $(NF)}' VarScan-modified.vcf | awk -F':' -v S="${RAED_DEPTH_CUTOFF}" '{ if ($3>=S) print $7}' | perl -pe 's/\%/ /g' | sort -n | uniq -c | perl -pe 's/^ +//' | awk '{m+=$1} END{print m;}' | awk '{ if ($1>0) print $1 ; else print 0}' `

echo -e "Total detected sites\t"$TOTAL_SITE > biallele_frequency.txt
echo -e "Allele sites\t"$ALLELE_SITE >> biallele_frequency.txt
echo -e "biallele sites\t"$(($TOTAL_SITE-$ALLELE_SITE)) >> biallele_frequency.txt
echo $ALLELE_SITE $TOTAL_SITE | awk '{ if ($2>0) print "biallele frequency\t"(1-$1/$2)*100"%" ; else print "biallele frequency\t0%" }' >> biallele_frequency.txt
echo -e "biallele frequency in total coverage length\t"`echo $ALLELE_SITE $TOTAL_SITE $TOTALCOVLEN | awk '{ if ($3>0) res=($2-$1)/$3*100; else res="0"; OFMT="%.10f"; print res}'`"%" >> biallele_frequency.txt

java -jar $PICARD CollectInsertSizeMetrics INPUT=realigned_reads.bam OUTPUT=insert_metrics.txt HISTOGRAM_FILE=insert_size_histogram.pdf

rm -f raw_indels.vcf raw_indels.vcf.idx raw_variants.vcf raw_variants.vcf.idx realignment_targets.list metrics.txt

####################
##### assemble #####
####################

if [ "${fastqfile1}" != "" -a "${fastqfile2}" != "" ]; then

   echo -e "\nINFO `date`    Fastq data, assemble start.\n"

skesa --fastq $trimread1fastq,$trimread2fastq  --cores $THREADS --memory 48 > skesa-assembled-contigs-default.fasta

megablast -d $REFSEQ -i skesa-assembled-contigs-default.fasta -F F -e 1.00E-20 -m 8 -o skesa-contigs.megablast.table

awk -F'\t' '{if($9<$10) print $1"\tpositive hit\t"$7"\t"$8"\t"$9"\t"$10; else print $1"\tnegative hit\t"$7"\t"$8"\t"$10"\t"$9}' skesa-contigs.megablast.table | \
awk -F'\t' '$4-$3+1 >= 200 {print $0}' | sort -t$'\t' -n -k5,5 > skesa-contigs.megablast-hit-list.table

awk -F'\t' '{if(min[$1] == "") min[$1]="inf"; \
if(rmin[$1] == "") rmin[$1]="inf"; \
if(min[$1]>$3) min[$1]=$3 ;\
if(max[$1]<$4) max[$1]=$4 ;\
if(strand[$1]<$4) strand[$1]=$2 ;\
if(rmin[$1]>$5) rmin[$1]=$5 ;\
if(rmax[$1]<$6) rmax[$1]=$6} ;\
END{for(key in min) print key"\t"strand[key]"\t"min[key]"\t"max[key]"\t"rmin[key]"\t"rmax[key]}' skesa-contigs.megablast-hit-list.table | sort > extraction-data-skesa-contig.txt

$ONELINEFASTA skesa-assembled-contigs-default.fasta | perl -pe 's/>//g' | perl -pe 's/\%+/\t/g' | sort > skesa-oneline.txt

join -1 1 -2 1 -t "$(printf '\011')" extraction-data-skesa-contig.txt skesa-oneline.txt | sort -t$'\t' -nk5,5 | awk -F'\t' '{print $1"\t"$2"\t"substr($7,$3,$4-$3+1)}' > extraction-data-skesa-contig-trimed-seq.txt

perl -F'\t' -anle '{if (@F[1] =~ /negative/) { $rev = reverse(@F[2]); $rev =~ tr/AaCcGgTtBbDdHhKkMmRrSsVvWwYyNn/TtGgCcAaVvHhDdMmKkYySsBbWwRrNn/;\
 print ">".@F[0]."\n".$rev; } else {  print ">".@F[0]."\n".@F[2]; }}' extraction-data-skesa-contig-trimed-seq.txt > extract-contigs-skesa.fasta

rm skesa-oneline.txt extraction-data-skesa-contig-trimed-seq.txt

megablast -d $REFSEQ -i extract-contigs-skesa.fasta -F F -e 1.00E-20 -m 8 -o extract-contigs-skesa.megablast.table

awk -F'\t' '{if($9<$10) print $1"\tpositive hit\t"$7"\t"$8"\t"$9"\t"$10; else print $1"\tnegative hit\t"$7"\t"$8"\t"$10"\t"$9}' extract-contigs-skesa.megablast.table | \
awk -F'\t' '$4-$3+1 >= 200 {print $0}' | sort -t$'\t' -n -k5,5 > extract-contigs-skesa.megablast-hit-list.table



a5_pipeline.pl --threads=$THREADS $trimread1fastq $trimread2fastq a5-miseq-assemble

cp a5-miseq-assemble.contigs.fasta a5-miseq-assemble.contigs.fasta-before-header-rename
set +C
awk '{print $1}' a5-miseq-assemble.contigs.fasta-before-header-rename > a5-miseq-assemble.contigs.fasta
set -C

megablast -d $REFSEQ -i a5-miseq-assemble.contigs.fasta -F F -e 1.00E-20 -m 8 -o a5-miseq-assemble.contigs.megablast.table

awk -F'\t' '{if($9<$10) print $1"\tpositive hit\t"$7"\t"$8"\t"$9"\t"$10; else print $1"\tnegative hit\t"$7"\t"$8"\t"$10"\t"$9}' a5-miseq-assemble.contigs.megablast.table | \
awk -F'\t' '$4-$3+1 >= 200 {print $0}' | sort -t$'\t' -n -k5,5 > a5-miseq-assemble.contigs.megablast-hit-list.table

awk -F'\t' '{if(min[$1] == "") min[$1]="inf"; \
if(rmin[$1] == "") rmin[$1]="inf"; \
if(min[$1]>$3) min[$1]=$3 ;\
if(max[$1]<$4) max[$1]=$4 ;\
if(strand[$1]<$4) strand[$1]=$2 ;\
if(rmin[$1]>$5) rmin[$1]=$5 ;\
if(rmax[$1]<$6) rmax[$1]=$6} ;\
END{for(key in min) print key"\t"strand[key]"\t"min[key]"\t"max[key]"\t"rmin[key]"\t"rmax[key]}' a5-miseq-assemble.contigs.megablast-hit-list.table | sort > extraction-data-a5miseq-contig.txt

$ONELINEFASTA a5-miseq-assemble.contigs.fasta | perl -pe 's/>//g' | perl -pe 's/\%+/\t/g' | sort > a5miseq-oneline.txt

join -1 1 -2 1 -t "$(printf '\011')" extraction-data-a5miseq-contig.txt a5miseq-oneline.txt | sort -t$'\t' -nk5,5 | awk -F'\t' '{print $1"\t"$2"\t"substr($7,$3,$4-$3+1)}' > extraction-data-a5miseq-contig-trimed-seq.txt

perl -F'\t' -anle '{if (@F[1] =~ /negative/) { $rev = reverse(@F[2]); $rev =~ tr/AaCcGgTtBbDdHhKkMmRrSsVvWwYyNn/TtGgCcAaVvHhDdMmKkYySsBbWwRrNn/;\
 print ">".@F[0]."\n".$rev; } else {  print ">".@F[0]."\n".@F[2]; }}' extraction-data-a5miseq-contig-trimed-seq.txt > extract-contigs-a5miseq.fasta

rm a5miseq-oneline.txt extraction-data-a5miseq-contig-trimed-seq.txt

megablast -d $REFSEQ -i extract-contigs-a5miseq.fasta -F F -e 1.00E-20 -m 8 -o extract-contigs-a5miseq.megablast.table

awk -F'\t' '{if($9<$10) print $1"\tpositive hit\t"$7"\t"$8"\t"$9"\t"$10; else print $1"\tnegative hit\t"$7"\t"$8"\t"$10"\t"$9}' extract-contigs-a5miseq.megablast.table | \
awk -F'\t' '$4-$3+1 >= 200 {print $0}' | sort -t$'\t' -n -k5,5 > extract-contigs-a5miseq.megablast-hit-list.table

rm -rf extract-contigs.txt list.txt a5-miseq-assemble.broken.scaffolds.fasta a5-miseq-assemble.crude.scaffolds.fasta \
       a5-miseq-assemble.ec.fastq.gz a5-miseq-assemble.final.scaffolds.fastq.gz a5-miseq-assemble.final.scaffolds.qvl \
       a5-miseq-assemble.library_1.txt a5-miseq-assemble.library_1.txt.strict a5-miseq-assemble.preproc.libs \
       a5-miseq-assemble.raw1.pe.sort.bam a5-miseq-assemble.raw1.pe.sort.bam.bai a5-miseq-assemble.s* \
       a5-miseq-assemble.tmplibs a5-miseq-assemble.assembly_stats.csv a5-miseq-assemble


echo -e "# query sequence\treference sequence id\tidentity %\talignment length\tnumber of mismatches\tnumber of gap openings\tstart of alignment in query\tend of alignment in query\tstart of alignment in ref\tend of alignment in ref\tE-value\tbit score" > header.txt
for i in *.megablast.table ; do cat header.txt $i > data.txt ; mv data.txt $i ; done

set +C
echo -e "# query sequence\tstrand of query\tstart of alignment in query\tend of alignment in query\tstart of alignment in ref\tend of alignment in ref" > header.txt
set -C
for i in *.megablast-hit-list.table ; do cat header.txt $i > data.txt ; mv data.txt $i ; done

rm header.txt


for i in *fastq; do gzip $i ; done

##############################################
### indel including vcf consensus sequence ###
##############################################

grep "^#" VarScan-modified-rm-Mix-allele.vcf > VarScan-GATK-mix.vcf
set +e
grep -v "^#" VarScan-modified-rm-Mix-allele.vcf > VarScan-GATK-mix-1.txt
set -e

awk -F'\t' '$6>=100 && $7=="PASS" && $NF~/^1[\/\|]/ {print $0}' GATK-indel-filtered-ref.vcf > VarScan-GATK-mix-2.txt

cat VarScan-GATK-mix-1.txt VarScan-GATK-mix-2.txt | sort -k1,1 -nk2,2 >> VarScan-GATK-mix.vcf

bgzip VarScan-GATK-mix.vcf
tabix -p vcf VarScan-GATK-mix.vcf.gz
export TGZPATH=`pwd`
export SAMPLENAME=`basename $TGZPATH`

cat $REFSEQ | bcftools consensus VarScan-GATK-mix.vcf.gz -o vcf-consensus-with-indel-non-concatenate.fasta &> bcftools-consensus-log.txt

set +e
SKIPCOUNT=`grep -c "skipping" bcftools-consensus-log.txt`
set -e

if [ $SKIPCOUNT -eq 0 ]; then

        rm bcftools-consensus-log.txt

else

        rm bcftools-consensus-log.txt
        mkdir retry-mapping-indel
        cd retry-mapping-indel

        mv ../vcf-consensus-with-indel-non-concatenate.fasta ./
        RETRYDB=vcf-consensus-with-indel-non-concatenate.fasta

        bwa index $RETRYDB
        gatk CreateSequenceDictionary --REFERENCE $RETRYDB --OUTPUT ${RETRYDB%.fasta}.dict
        samtools faidx $RETRYDB

        bwa mem -t $THREADS -k 31 -a -M $RETRYDB ../$trimread1fastq ../$trimread2fastq -R "@RG\tID:01\tSM:$SAMPLENAME\tPL:Illumina" > out12_mem.sam
        gatk --java-options "-XX:ParallelGCThreads=$THREADS" SortSam --INPUT out12_mem.sam --OUTPUT out12_mem.sorted.bam --SORT_ORDER coordinate
        gatk --java-options "-XX:ParallelGCThreads=$THREADS" BuildBamIndex --INPUT out12_mem.sorted.bam
        rm out12_mem.sam
        gatk3 -XX:ParallelGCThreads=$THREADS -T RealignerTargetCreator -R $RETRYDB -I out12_mem.sorted.bam -o realignment_targets.list
        gatk3 -XX:ParallelGCThreads=$THREADS -T IndelRealigner -R $RETRYDB -I out12_mem.sorted.bam -targetIntervals realignment_targets.list -o realigned_reads.bam --filter_bases_not_stored
        mv realigned_reads.bai realigned_reads.bam.bai

        gatk --java-options "-Xmx4G -XX:ParallelGCThreads=$THREADS" HaplotypeCaller -R $RETRYDB --emit-ref-confidence GVCF -I realigned_reads.bam -O raw_variants.gvcf --read-filter ReadStrandFilter --keep-reverse-strand-only false
        gatk --java-options "-Xmx4G -XX:ParallelGCThreads=$THREADS" GenotypeGVCFs -R $RETRYDB -V raw_variants.gvcf -O raw_variants.vcf
        rm raw_variants.gvcf

        gatk --java-options "-Xmx4G -XX:ParallelGCThreads=$THREADS" SelectVariants -R $RETRYDB -V raw_variants.vcf --select-type-to-include INDEL -O raw_indels.vcf
        gatk --java-options "-Xmx4G -XX:ParallelGCThreads=$THREADS" VariantFiltration -R $RETRYDB -V raw_indels.vcf -O GATK-indel-filtered-ref.vcf \
                                 -filter "QD < 2.0" --filter-name "QD2"       \
                                 -filter "QUAL < 30.0" --filter-name "QUAL30" \
                                 -filter "FS > 200.0" --filter-name "FS200"   \
                                 -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20"

        cp GATK-indel-filtered-ref.vcf GATK-indel-filtered-ref.vcf-before-modi

        set +C
        grep "^#" GATK-indel-filtered-ref.vcf-before-modi > GATK-indel-filtered-ref.vcf
        set -C
        awk -F'\t' '$7~/PASS/ && $10~/^1[\/\|]1/ {print $0}' GATK-indel-filtered-ref.vcf-before-modi >> GATK-indel-filtered-ref.vcf
        set +e
        diff GATK-indel-filtered-ref.vcf GATK-indel-filtered-ref.vcf-before-modi > GATK-diff-check.txt
        set -e

        if [ -s GATK-diff-check.txt ]; then
             rm GATK-diff-check.txt
           else
             rm GATK-diff-check.txt \
                GATK-indel-filtered-ref.vcf-before-modi
        fi

        gatk --java-options "-Xmx4G -XX:ParallelGCThreads=$THREADS" SelectVariants -R $RETRYDB -V raw_variants.vcf --select-type-to-include SNP -O raw_snps.vcf
        gatk --java-options "-Xmx4G -XX:ParallelGCThreads=$THREADS" VariantFiltration -R $RETRYDB -V raw_snps.vcf -O GATK-snps-filtered-ref.vcf \
                                 -filter "QD < 2.0" --filter-name "QD2"       \
                                 -filter "QUAL < 30.0" --filter-name "QUAL30" \
                                 -filter "SOR > 4.0" --filter-name "SOR4"     \
                                 -filter "FS > 60.0" --filter-name "FS60"     \
                                 -filter "MQ < 40.0" --filter-name "MQ40"     \
                                 -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
                                 -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"

        rm raw_variants.vcf raw_variants.gvcf.idx raw_snps.vcf
        grep -v "^#" GATK-snps-filtered-ref.vcf | cut -f 1,2 > GATK-snps-filtered-list.txt
        samtools mpileup -B -A -f $RETRYDB realigned_reads.bam | varscan mpileup2snp --strand-filter 0 --min-freq-for-hom 0.6 --output-vcf | perl -pe s/Sample1/$SAMPLENAME/ > VarScan-ref.vcf
        samtools mpileup -B -A -f $RETRYDB realigned_reads.bam | varscan mpileup2snp --strand-filter 0 --min-freq-for-hom 0.6 > VarScan.default-format.txt
        grep "^##" VarScan-ref.vcf > VarScan-ref.vcf.header1
        echo '##FORMAT=<ID=CONS,Number=1,Type=Integer,Description="Consensus sequence with VarScan basecalling">' > VarScan-ref.vcf.header2
        grep "^#CHROM" VarScan-ref.vcf > VarScan-ref.vcf.header3
        grep -v "^#" VarScan-ref.vcf | cut -f 1-10 | perl -pe 's/GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR/GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR:CONS/' > VarScan-ref.vcf.data
        cut -f 5 VarScan.default-format.txt | grep -v "Cons:Cov:Reads1:Reads2:Freq:P-value" | perl -pe 's/^/\|/' | awk -F':' '{print $1}' > VarScan.default-format.txt.data
        paste VarScan-ref.vcf.data VarScan.default-format.txt.data | perl -pe 's/\t\|/:/' > VarScan.format-merge.txt
        cat VarScan-ref.vcf.header1 VarScan-ref.vcf.header2 VarScan-ref.vcf.header3 VarScan.format-merge.txt > VarScan-modified.vcf     # varscanは、75%以上はMixと判定していない
        rm VarScan-ref.vcf.header* VarScan-ref.vcf.data VarScan.default-format.txt.data VarScan.format-merge.txt VarScan-ref.vcf VarScan.default-format.txt
        igvtools index VarScan-modified.vcf

        grep "^#" VarScan-modified.vcf > VarScan-modified-rm-Mix-allele.vcf

        while read seq posi ; do \
           awk -F'\t' -v C=$seq -v P=$posi '$1==C && $2==P {print $0}' \
           VarScan-modified.vcf ; done < GATK-snps-filtered-list.txt > SNV-VarScan-GATK.txt

        awk -F'\t' '$1~/^#/ || $10~/^1\/1/ {print $0}' SNV-VarScan-GATK.txt >> VarScan-modified-rm-Mix-allele.vcf

        igvtools index VarScan-modified-rm-Mix-allele.vcf

        rm GATK-snps-filtered-list.txt igv.log SNV-VarScan-GATK.txt GATK-snps-filtered-ref.vc*

        grep "^#" GATK-indel-filtered-ref.vcf > GATK-indel-VarScan-snp.vcf
        set +e
        grep -v "^#" GATK-indel-filtered-ref.vcf > GATK-indel-VarScan-snp-1.vcf
        grep -v "^#" VarScan-modified-rm-Mix-allele.vcf >> GATK-indel-VarScan-snp-1.vcf
        set -e
        LANG=C sort -k2,2 -n GATK-indel-VarScan-snp-1.vcf >> GATK-indel-VarScan-snp.vcf
        rm GATK-indel-VarScan-snp-1.vcf

        bgzip GATK-indel-VarScan-snp.vcf
        tabix -p vcf GATK-indel-VarScan-snp.vcf.gz

        cat $RETRYDB | bcftools consensus GATK-indel-VarScan-snp.vcf.gz -o vcf-consensus-with-indel-non-concatenate-retry.fasta &> bcftools-consensus-log.txt

        mv vcf-consensus-with-indel-non-concatenate-retry.fasta ../vcf-consensus-with-indel-non-concatenate.fasta

        cd ../

        rm -rf retry-mapping-indel

fi

$ONELINEFASTA vcf-consensus-with-indel-non-concatenate.fasta | grep ">" | perl -pe 's/\%+/\n/g' | grep -v ">" > vcf-consensus-with-indel-concatenate.sequence
echo -e ">${SAMPLENAME}" > concatenatevcf-header.txt
cat concatenatevcf-header.txt vcf-consensus-with-indel-concatenate.sequence > vcf-consensus-with-indel.fasta
rm concatenatevcf-header.txt vcf-consensus-with-indel-concatenate.sequence \
   vcf-consensus-with-indel-non-concatenate.fasta VarScan-GATK-mix.vcf.gz.tbi \
   VarScan-GATK-mix.vcf.gz VarScan-GATK-mix-1.txt VarScan-GATK-mix-2.txt

############################
### large deletion check ###
############################

DB=vcf-consensus-with-indel.fasta

formatdb -p F -o F -i $DB

### skesa contig curation ###

megablast -d $DB -i extract-contigs-skesa.fasta -F F -e 1.00E-20 -m 8 -o extract-contigs-skesa-ref-vcf.megablast.table

awk -F'\t' '{if($9<$10) print $1"\tpositive hit\t"$7"\t"$8"\t"$9"\t"$10; else print $1"\tnegative hit\t"$7"\t"$8"\t"$10"\t"$9}' \
   extract-contigs-skesa-ref-vcf.megablast.table | awk -F'\t' '$4-$3+1 >= 200 {print $0}' | \
   sort -t$'\t' -n -k5,5 > extract-contigs-skesa-ref-vcf.megablast-hit-list.table


### A5miseq contig curation ###

megablast -d $DB -i extract-contigs-a5miseq.fasta -F F -e 1.00E-20 -m 8 -o extract-contigs-a5miseq-ref-vcf.megablast.table

awk -F'\t' '{if($9<$10) print $1"\tpositive hit\t"$7"\t"$8"\t"$9"\t"$10; else print $1"\tnegative hit\t"$7"\t"$8"\t"$10"\t"$9}' \
   extract-contigs-a5miseq-ref-vcf.megablast.table | awk -F'\t' '$4-$3+1 >= 200 {print $0}' | \
   sort -t$'\t' -n -k5,5 > extract-contigs-a5miseq-ref-vcf.megablast-hit-list.table


### large gap region extraction ###

BLASTDATA=extract-contigs-skesa-ref-vcf.megablast-hit-list.table
awk -F'\t' '{print $6+1}' $BLASTDATA > skesa-gap-left.txt
awk -F'\t' 'NR>1 {print $5-1}' $BLASTDATA  > skesa-gap-right.txt
paste skesa-gap-left.txt skesa-gap-right.txt | awk -F'\t' '$1<$2 {print $1"\t"$2}' > skesa-gap-list.txt
rm skesa-gap-left.txt skesa-gap-right.txt

BLASTDATA=extract-contigs-a5miseq-ref-vcf.megablast-hit-list.table
awk -F'\t' '{print $6+1}' $BLASTDATA > a5miseq-gap-left.txt
awk -F'\t' 'NR>1 {print $5-1}' $BLASTDATA  > a5miseq-gap-right.txt
paste a5miseq-gap-left.txt a5miseq-gap-right.txt | awk -F'\t' '$1<$2 {print $1"\t"$2}' > a5miseq-gap-list.txt
rm a5miseq-gap-left.txt a5miseq-gap-right.txt

export TGZPATH=`pwd`
export SAMPLENAME=`basename $TGZPATH`

cat a5miseq-gap-list.txt skesa-gap-list.txt | sort | \
   uniq -c | grep "    2 " | awk -F' +' -v SEQNAME=${SAMPLENAME} '{print SEQNAME"\t"$3}' > large-deletion.txt

rm vcf-consensus-with-indel.fasta.n* *-ref-vcf.megablast* *-gap-list.txt formatdb.log

########################
### mapping analysis ###
########################

DB=vcf-consensus-with-indel.fasta
trimread1fastq=trim-ARTIC-primer_read1_over50mer.fastq.gz
trimread2fastq=trim-ARTIC-primer_read2_over50mer.fastq.gz

export TGZPATH=`pwd`
export SAMPLENAME=`basename $TGZPATH`

bwa index $DB
gatk CreateSequenceDictionary --REFERENCE $DB --OUTPUT ${DB%.fasta}.dict
samtools faidx $DB

bwa mem -t 4 -k 31 -a -M $DB $trimread1fastq $trimread2fastq -R "@RG\tID:01\tSM:$SAMPLENAME\tPL:Illumina" > out12_mem_ref_vcf-cons-indel.sam
gatk --java-options "-XX:ParallelGCThreads=$THREADS" SortSam --INPUT out12_mem_ref_vcf-cons-indel.sam --OUTPUT out12_mem_ref_vcf-cons-indel.sorted.bam --SORT_ORDER coordinate
gatk --java-options "-XX:ParallelGCThreads=$THREADS" BuildBamIndex --INPUT out12_mem_ref_vcf-cons-indel.sorted.bam
rm out12_mem_ref_vcf-cons-indel.sam
gatk3 -XX:ParallelGCThreads=$THREADS -T RealignerTargetCreator -R $DB -I out12_mem_ref_vcf-cons-indel.sorted.bam -o realignment_targets_ref_vcf-cons-indel.list
gatk3 -XX:ParallelGCThreads=$THREADS -T IndelRealigner -R $DB -I out12_mem_ref_vcf-cons-indel.sorted.bam -targetIntervals realignment_targets_ref_vcf-cons-indel.list -o realigned_reads_ref_vcf-cons-indel.bam --filter_bases_not_stored 
mv realigned_reads_ref_vcf-cons-indel.bai realigned_reads_ref_vcf-cons-indel.bam.bai


#################################################################
### masking for low coverage region and trim 5' and 3' region ###
#################################################################

bedtools genomecov -d -split -ibam realigned_reads_ref_vcf-cons-indel.bam > BED-depth_ref_vcf-cons-indel.txt

RAED_DEPTH_CUTOFF=10  # for SARS-CoV-2 mapping analysis

awk -F'\t' -v S="${RAED_DEPTH_CUTOFF}" '{ if ($3>=S) print $1"\t"$2"\t"$3"\tPASS_(>=_"S")" ; else print $1"\t"$2"\t"$3"\tFAIL_(<_"S")" }' BED-depth_ref_vcf-cons-indel.txt | sort -k1,1 -k2,2n > BED-depth-check_ref_vcf-cons-indel.txt

$FASTAVERTICAL  vcf-consensus-with-indel.fasta | \
grep -v ">" > vcf-consensus-with-indel-vertical.fasta


paste BED-depth-check_ref_vcf-cons-indel.txt vcf-consensus-with-indel-vertical.fasta | \
   awk -F'\t' '{if($4~/PASS/ && $2>=55 && $2<=29835) {print $0"\t"$5} \
   else if($4~/FAIL/ && $2>=55 && $2<=29835) {print $0"\tN"}}' > BED-depth-check_ref_vcf-cons-indel-trim-with-seq.txt

while read seq start end ; do \
   awk -F'\t' -v C=$seq -v S=$start -v E=$end '$1==C && $2>=S && $2<=E {print $0}' \
   BED-depth-check_ref_vcf-cons-indel-trim-with-seq.txt ; done < large-deletion.txt > large-deletion-BED.txt

echo -e ">$SAMPLENAME" > vcf-consensus-with-indel-trim.fasta
diff BED-depth-check_ref_vcf-cons-indel-trim-with-seq.txt large-deletion-BED.txt | \
   awk -F'\t' '$1~/^< / {print $6}' | perl -pe 's/\n//g' | perl -pe 's/$/\n/g' | perl -pe 's/N+$//'\
   | perl -pe 's/N+.{0,4}A+$//' >> vcf-consensus-with-indel-trim.fasta

     ### curation of 5' and 3' region for ARTIC data ###

     cat $ARTICREGIONDB vcf-consensus-with-indel-trim.fasta > MN908947_55-29835-vs-vcf-consensus-with-indel-trim.fasta
     mafft MN908947_55-29835-vs-vcf-consensus-with-indel-trim.fasta > MN908947_55-29835-vs-vcf-consensus-with-indel-trim.aln

     ALNSTART=`$ONELINEFASTA MN908947_55-29835-vs-vcf-consensus-with-indel-trim.aln | perl -pe 's/\%+/\t/g' | \
               awk -F'\t' '$1~/MN908947_55-29835/ {if(match($2,/^-+/)) \
               {num = substr($2, RSTART, RLENGTH); print length(num)} else {print "0"}}'`

     ALNEND=`$ONELINEFASTA MN908947_55-29835-vs-vcf-consensus-with-indel-trim.aln | perl -pe 's/\%+/\t/g' | \
               awk -F'\t' '$1~/MN908947_55-29835/ {if(match($2,/-+$/)) \
               {num = substr($2, RSTART, RLENGTH); print length(num)} else {print "0"}}'`

     $ONELINEFASTA vcf-consensus-with-indel-trim.fasta | grep ">" | perl -pe 's/\%+/\t/g' | \
        awk -F'\t' -v S=$ALNSTART -v E=$ALNEND '{print $1"\n"substr($2, S+1, length($2)-S-E)}' > vcf-consensus-with-indel-trim-ARTIC-region.fasta

     rm MN908947_55-29835-vs-vcf-consensus-with-indel-trim.*

     mv vcf-consensus-with-indel-trim-ARTIC-region.fasta vcf-consensus-with-indel-trim.fasta

rm -f vcf-consensus-with-indel.fast* vcf-consensus-with-indel.dict *_vcf-cons-indel* vcf-consensus-with-indel-vertical.fasta \
   BED-depth-check_ref_vcf-cons-indel-trim-with-seq.txt large-deletion-BED.txt 5-and-3-end-check-ARTIC-region.txt


############################################
### final complete genome sequence check ###
############################################

LENCUTOFF=28000

if [ `grep -c ">" extract-contigs-skesa.fasta` -eq 1 ] && \
   [ `grep -v ">" extract-contigs-skesa.fasta | perl -pe 's/\n//g' | wc -c` -ge $LENCUTOFF ]; then

     REFCONTIGLEN=`grep -v ">" vcf-consensus-with-indel-trim.fasta | perl -pe 's/\n//g' | wc -c`
     QUECONTIGLEN=`grep -v ">" extract-contigs-skesa.fasta | perl -pe 's/\n//g' | wc -c`
     
     bl2seq -j vcf-consensus-with-indel-trim.fasta -i extract-contigs-skesa.fasta \
        -F F -p blastn -e 1.00E-10 -D 1 | grep -v "^#" > 5-and-3-end-check.txt

     QUE5ENDPOSI=`cut -f 7,8 5-and-3-end-check.txt | perl -pe 's/\t/\n/g' | sort -n | head -n 1`
     QUE3ENDPOSI=`cut -f 7,8 5-and-3-end-check.txt | perl -pe 's/\t/\n/g' | sort -n | tail -n 1`

     REF5ENDPOSI=`cut -f 9,10 5-and-3-end-check.txt | perl -pe 's/\t/\n/g' | sort -n | head -n 1`
     REF3ENDPOSI=`cut -f 9,10 5-and-3-end-check.txt | perl -pe 's/\t/\n/g' | sort -n | tail -n 1`

     echo ">$SAMPLENAME" > skesa-complete.fasta
     grep -v ">" vcf-consensus-with-indel-trim.fasta | awk -v s=$REF5ENDPOSI '{print substr($0, 1, s-1)}' >> skesa-complete.fasta
     QUEEXTACTLEN=`echo $(($QUE3ENDPOSI - $QUE5ENDPOSI + 1 ))`
     grep -v ">" extract-contigs-skesa.fasta | awk -v s=$QUE5ENDPOSI -v l=$QUEEXTACTLEN '{print substr($0, s, l)}' >> skesa-complete.fasta
     REFEXTACTLEN=`echo $(($REFCONTIGLEN - $REF3ENDPOSI ))`
     grep -v ">" vcf-consensus-with-indel-trim.fasta | awk -v s=$REF3ENDPOSI -v l=$REFEXTACTLEN '{print substr($0, s+1, l)}' >> skesa-complete.fasta
     $ONELINEFASTA skesa-complete.fasta | grep ">" | perl -pe 's/\%+/\n/g' > $SAMPLENAME.skesa-complete.fasta
     rm skesa-complete.fasta
     
     rm 5-and-3-end-check.txt

elif [ `grep -c ">" extract-contigs-skesa.fasta` -le 6 ] && \
     [ `grep -c "_Circ" extract-contigs-skesa.fasta` -eq 0 ] && \
     [ `grep -v ">" extract-contigs-skesa.fasta | perl -pe 's/\n//g' | wc -c` -ge $(($LENCUTOFF-2000)) ]; then

    #########################################################################
    ### skesa contig scaffold and gap filling with vcf consensus sequence ###
    #########################################################################

    echo "final complete genome sequence check: skesa contig scaffold and gap filling with vcf consensus sequence"

    ragtag.py scaffold vcf-consensus-with-indel-trim.fasta extract-contigs-skesa.fasta -o ragtag_output_scaffold

    cd ragtag_output_scaffold
    tgsgapcloser --scaff ragtag.scaffold.fasta --reads ../vcf-consensus-with-indel-trim.fasta --ne --min_match 50 --output skesa-ragtag-scaffold-gapfill-vcf  # min_matchを下げることで、バラバラでもいける。
    perl -pe 's/_RagTag//g' skesa-ragtag-scaffold-gapfill-vcf.scaff_seqs | $ONELINEFASTA | grep ">" | perl -pe 's/\%+/\n/g' > extract-contigs-skesa-gapfill-vcf-consensus.fasta

    REFCONTIGLEN=`grep -v ">" ../vcf-consensus-with-indel-trim.fasta | perl -pe 's/\n//g' | wc -c`
    QUECONTIGLEN=`grep -v ">" extract-contigs-skesa-gapfill-vcf-consensus.fasta | perl -pe 's/\n//g' | wc -c`

    bl2seq -j ../vcf-consensus-with-indel-trim.fasta -i extract-contigs-skesa-gapfill-vcf-consensus.fasta \
       -F F -p blastn -e 1.00E-10 -D 1 | grep -v "^#" > 5-and-3-end-check.txt

    QUE5ENDPOSI=`cut -f 7,8 5-and-3-end-check.txt | perl -pe 's/\t/\n/g' | sort -n | head -n 1`
    QUE3ENDPOSI=`cut -f 7,8 5-and-3-end-check.txt | perl -pe 's/\t/\n/g' | sort -n | tail -n 1`

    REF5ENDPOSI=`cut -f 9,10 5-and-3-end-check.txt | perl -pe 's/\t/\n/g' | sort -n | head -n 1`
    REF3ENDPOSI=`cut -f 9,10 5-and-3-end-check.txt | perl -pe 's/\t/\n/g' | sort -n | tail -n 1`

    echo ">$SAMPLENAME" > skesa-vcf-gapfill.fasta
    grep -v ">" ../vcf-consensus-with-indel-trim.fasta | awk -v s=$REF5ENDPOSI '{print substr($0, 1, s-1)}' >> skesa-vcf-gapfill.fasta
    QUEEXTACTLEN=`echo $(($QUE3ENDPOSI - $QUE5ENDPOSI + 1 ))`
    grep -v ">" extract-contigs-skesa-gapfill-vcf-consensus.fasta | awk -v s=$QUE5ENDPOSI -v l=$QUEEXTACTLEN '{print substr($0, s, l)}' >> skesa-vcf-gapfill.fasta
    REFEXTACTLEN=`echo $(($REFCONTIGLEN - $REF3ENDPOSI ))`
    grep -v ">" ../vcf-consensus-with-indel-trim.fasta | awk -v s=$REF3ENDPOSI -v l=$REFEXTACTLEN '{print substr($0, s+1, l)}' >> skesa-vcf-gapfill.fasta
    $ONELINEFASTA skesa-vcf-gapfill.fasta | grep ">" | perl -pe 's/\%+/\n/g' > skesa-vcf-gapfill-2.fasta && mv skesa-vcf-gapfill-2.fasta skesa-vcf-gapfill.fasta

    nextclade run --input-dataset=${NEXTCLADE_DATA} \
                              --output-tsv=Nextclade-result-insert-check.tsv -j $THREADS skesa-vcf-gapfill.fasta

    rm -f *.aligned.fasta *.errors.csv *.gene*fasta *.insertions.csv

    $HVCHANGE Nextclade-result-insert-check.tsv | \
    egrep 'clade|substitutions|deletions|insertions|aaSubstitutions|aaDeletions|total|errors|Nextclade_pango' \
       > Nextclade-result-data-list.txt
    TOTALINSCHECK=`awk -F'\t' '$1=="totalInsertions" {print $2}' Nextclade-result-data-list.txt`
    rm Nextclade-result-data-list.txt

    cd ../
         
      COUNTZERO=0
      INSTHRESHOLD=20
      if [ `grep -v ">" ./ragtag_output_scaffold/skesa-vcf-gapfill.fasta | grep -c "N"` -eq $COUNTZERO ] && \
         [ `grep -v ">" ./ragtag_output_scaffold/skesa-vcf-gapfill.fasta | perl -pe 's/\n//g' | wc -c` -ge $LENCUTOFF ] && \
         [ $TOTALINSCHECK -le $INSTHRESHOLD ] ; then

        cat ./ragtag_output_scaffold/skesa-vcf-gapfill.fasta > $SAMPLENAME.vcf-consensus-complete.fasta
 
       elif [ `grep -v ">" ./ragtag_output_scaffold/skesa-vcf-gapfill.fasta | perl -pe 's/\n//g' | wc -c` -ge $LENCUTOFF ] && \
            [ $TOTALINSCHECK -le $INSTHRESHOLD ] ; then

        cat ./ragtag_output_scaffold/skesa-vcf-gapfill.fasta > $SAMPLENAME.draft.fasta

       elif [ `grep -v ">" vcf-consensus-with-indel-trim.fasta | grep -c "N"` -eq 0 ] && \
            [ `grep -v ">" vcf-consensus-with-indel-trim.fasta | perl -pe 's/\n//g' | wc -c` -ge $LENCUTOFF ]; then

         cat vcf-consensus-with-indel-trim.fasta > $SAMPLENAME.vcf-consensus-complete.fasta

       else

         cat vcf-consensus-with-indel-trim.fasta > $SAMPLENAME.draft.fasta

      fi

else

      echo "final complete genome sequence check: only vcf consensus sequence, many skesa contigs"

      if [ `grep -v ">" vcf-consensus-with-indel-trim.fasta | grep -c "N"` -eq 0 ] && \
         [ `grep -v ">" vcf-consensus-with-indel-trim.fasta | perl -pe 's/\n//g' | wc -c` -ge $LENCUTOFF ]; then
           cat vcf-consensus-with-indel-trim.fasta > $SAMPLENAME.vcf-consensus-complete.fasta
      else
           cat vcf-consensus-with-indel-trim.fasta > $SAMPLENAME.draft.fasta
      fi

fi

   echo -e "\nINFO `date`    assemble end.\n"

elif [ -f fasta-data.txt ]; then

echo -e "Input-data is fasta format. Skipped assmble step." > fasta-data-info.txt

else

echo -e "No fastq or fasta file ..." > error2.txt

fi


##############################################################
##### Nextclade and pangolin analysis using final contig #####
##############################################################

if [ -f fastq-data.txt -o -f fasta-data.txt ]; then
   rm -f Nextclade-result.tsv
   echo "Start Nextclade. "`date`
   # Create Dummy entry.
   echo -e \
"seqName\tclade\tqc.overallScore\tqc.overallStatus\ttotalGaps\ttotalInsertions\ttotalMissing\ttotalMutations\ttotalNonACGTNs\
\ttotalPcrPrimerChanges\tsubstitutions\tdeletions\tinsertions\tmissing\tnonACGTNs\tpcrPrimerChanges\taaSubstitutions\ttotalAminoacidSubstitutions\
\taaDeletions\ttotalAminoacidDeletions\talignmentEnd\talignmentScore\talignmentStart\tqc.missingData.missingDataThreshold\tqc.missingData.score\
\tqc.missingData.status\tqc.missingData.totalMissing\tqc.mixedSites.mixedSitesThreshold\tqc.mixedSites.score\tqc.mixedSites.status\
\tqc.mixedSites.totalMixedSites\tqc.privateMutations.cutoff\tqc.privateMutations.excess\tqc.privateMutations.score\tqc.privateMutations.status\
\tqc.privateMutations.total\tqc.snpClusters.clusteredSNPs\tqc.snpClusters.score\tqc.snpClusters.status\tqc.snpClusters.totalSNPs\terrors\n\
Analysis error, Try again.\tAnalysis error, Try again.\tAnalysis error, Try again.\tAnalysis error, Try again.\tAnalysis error, Try again.\
\tAnalysis error, Try again.\tAnalysis error, Try again.\tAnalysis error, Try again.\tAnalysis error, Try again.\tAnalysis error, Try again.\
\tAnalysis error, Try again.\tAnalysis error, Try again.\tAnalysis error, Try again.\tAnalysis error, Try again.\tAnalysis error, Try again.\
\tAnalysis error, Try again.\tAnalysis error, Try again.\tAnalysis error, Try again.\tAnalysis error, Try again.\tAnalysis error, Try again.\
\tAnalysis error, Try again.\tAnalysis error, Try again.\tAnalysis error, Try again.\tAnalysis error, Try again.\tAnalysis error, Try again.\
\tAnalysis error, Try again.\tAnalysis error, Try again.\tAnalysis error, Try again.\tAnalysis error, Try again.\tAnalysis error, Try again.\
\tAnalysis error, Try again.\tAnalysis error, Try again.\tAnalysis error, Try again.\tAnalysis error, Try again.\tAnalysis error, Try again.\
\tAnalysis error, Try again.\tAnalysis error, Try again.\tAnalysis error, Try again.\tAnalysis error, Try again.\tAnalysis error, Try again.\tAnalysis error, Try again." > Nextclade-result.tsv

   nextclade run --input-dataset=${NEXTCLADE_DATA} \
                              --output-tsv=Nextclade-result.tsv -j $THREADS $SAMPLENAME*.fasta

   NEXTCLADEVERSION=`nextclade --version | perl -pe 's/nextclade //g'`
   rm -f $SAMPLENAME*.aligned.fasta $SAMPLENAME*.errors.csv $SAMPLENAME*.gene*fasta $SAMPLENAME*.insertions.csv
   mv lineage_report.csv lineage_report_vcf-mapping-version.csv
   pangolin -t $THREADS $SAMPLENAME*.fasta
else

   echo -e "No fastq or fasta file ... Stop Nextclade analysis" >> error.txt
   NEXTCLADEVERSION="ERROR"
   exit

fi


##################
### plot detph ###
##################

export TGZPATH=`pwd`
export SAMPLENAME=`basename $TGZPATH`

export PATH=/opt/bio/bin:$PATH

ln -s realigned_reads.bam ${SAMPLENAME}.bam

$PLOTDEPTHPY \
-i ${SAMPLENAME}.bam \
-p $TRIMPRIMERBED \
-r $REFSEQ \
-s --only_primer_mismatch -o plot_depth.pdf

unlink ${SAMPLENAME}.bam

#############################################
##### json file and result summary data #####
#############################################

TGZPATH=`pwd`
SAMPLENAME=`basename $TGZPATH`

READ_DEPTH_MEDIAN=`awk -F'\t' 'END {print $4}' read-depth-stat.txt`
COVERAGE=`awk -F'\t' 'END{print $2}' Coverage-info.txt`
TOAL_SNV=`awk -F'\t' 'NR==1 {print $2}' biallele_frequency.txt`
SNV=`awk -F'\t' 'NR==2 {print $2}' biallele_frequency.txt`
MIX_SNV=`awk -F'\t' 'NR==3 {print $2}' biallele_frequency.txt`
MIX_SNV_RATIO_1=`awk -F'\t' 'NR==4 {print $2}' biallele_frequency.txt`
MIX_SNV_RATIO_2=`awk -F'\t' 'NR==5 {print $2}' biallele_frequency.txt`



if [ "${fastqfile1}" != "" -a "${fastqfile2}" != "" ]; then
  INPUT_DATA="raw read data"
  COUNT=`zcat *R1_001.fastq.gz | wc -l` ; TOTAL_READ=`echo -e "$(( $COUNT /2 ))"`
  COUNT=`zcat R1_trim-sickle.fastq.gz | wc -l` ; TRIMMED_READ=`echo -e "$(( $COUNT /2 ))"`
  COUNT=`zcat trim-ARTIC-primer_read1_over50mer.fastq.gz | wc -l` ; MAPPED_READ=`echo -e "$(( $COUNT /2 ))"`
  TRIMMED_READ_RATIO=`echo "scale=2; 100 * $TRIMMED_READ / $TOTAL_READ " | bc` 
  MAPPED_READ_RATIO=`echo "scale=2; 100 * $MAPPED_READ / $TRIMMED_READ " | bc`
  SKESA_LEN=`grep -v ">" extract-contigs-skesa.fasta | perl -pe 's/\n//' | wc -m`
  SKESA_NUM_CONTIGS=`grep -c ">" extract-contigs-skesa.fasta`
  A5_LEN=`grep -v ">" extract-contigs-a5miseq.fasta | perl -pe 's/\n//' | wc -m`
  A5_NUM_CONTIGS=`grep -c ">" extract-contigs-a5miseq.fasta`
  INPUT_LEN="0"
  INPUT_NUM_CONTIGS="0"
   DATA1=$SAMPLENAME.skesa-complete.fasta
   DATA2=$SAMPLENAME.a5miseq-complete.fasta
   DATA3=$SAMPLENAME.vcf-consensus-complete.fasta
   DATA4=$SAMPLENAME.draft.fasta
    if [ -f ${DATA1} ]; then
        DATA_LEN=`grep -v ">" $DATA1 | perl -pe 's/\n//' | wc -m`
        FINAL_CONTIG_LEN=$DATA_LEN
        FINAL_CONTIG_METHOD="de novo assemble (SKESA v. 2.3.0)"
        FINAL_CONTIG_STATUS="complete"
    elif [ -f ${DATA2} ]; then
        DATA_LEN=`grep -v ">" $DATA2 | perl -pe 's/\n//' | wc -m`
        FINAL_CONTIG_LEN=$DATA_LEN
        FINAL_CONTIG_METHOD="de novo assemble (A5-miseq v. 20140604)"
        FINAL_CONTIG_STATUS="complete"
    elif [ -f ${DATA3} ]; then
        DATA_LEN=`grep -v ">" $DATA3 | perl -pe 's/\n//' | wc -m`
        FINAL_CONTIG_LEN=$DATA_LEN
        FINAL_CONTIG_METHOD="read mapping (bwa v. 0.7.13-r1126, picard v. 2.10.10, GATK v. 3.8.0, VarScan v. 2.4.3, VCFtools v. 0.1.12b, SKESA v. 2.3.0, A5-miseq v. 20140604)"
        FINAL_CONTIG_STATUS="complete"
    elif [ -f ${DATA4} ]; then
        DATA_LEN=`grep -v ">" $DATA4 | perl -pe 's/\n//' | wc -m`
        FINAL_CONTIG_LEN=$DATA_LEN
        FINAL_CONTIG_METHOD="read mapping (bwa v. 0.7.13-r1126, picard v. 2.10.10, GATK v. 3.8.0, VarScan v. 2.4.3, VCFtools v. 0.1.12b, SKESA v. 2.3.0, A5-miseq v. 20140604)"
        FINAL_CONTIG_STATUS="draft"
    fi
else
  INPUT_DATA="contig data"
  TOTAL_READ="0"
  TRIMMED_READ="0"
  MAPPED_READ="0"
  TRIMMED_READ_RATIO="0"
  MAPPED_READ_RATIO="0"
  SKESA_LEN="0"
  SKESA_NUM_CONTIGS="0"
  A5_LEN="0"
  A5_NUM_CONTIGS="0"
  FINAL_CONTIG_LEN=`grep -v ">" $SAMPLENAME.fasta | perl -pe 's/\n//' | wc -m`
  FINAL_NUM_CONTIG=`grep -c ">" $SAMPLENAME.fasta`
  FINAL_CONTIG_METHOD="(input data)"
  FINAL_CONTIG_STATUS="(input data)"
  INPUT_LEN=$FINAL_CONTIG_LEN
  INPUT_NUM_CONTIGS=$FINAL_NUM_CONTIG
fi

MEDIAN_INSERT_SIZE=`grep -A1 "MEDIAN_INSERT_SIZE" insert_metrics.txt | awk -F'\t' 'END {print $1}'`
MEAN_INSERT_SIZE=`grep -A1 "MEDIAN_INSERT_SIZE" insert_metrics.txt | awk -F'\t' 'END {print $5}'`
PANGOLIN_PROBABILITY=`awk -F',' 'END {print $3}' lineage_report.csv` # abolition of pangolin data, chose nextclade lineage data
READ_DEPTH_MEAN=`awk -F'\t' 'END {print $1}' read-depth-stat.txt`
READ_DEPTH_MIN=`awk -F'\t' 'END {print $2}' read-depth-stat.txt`
READ_DEPTH_Q1=`awk -F'\t' 'END {print $3}' read-depth-stat.txt`
READ_DEPTH_MEDIAN=`awk -F'\t' 'END {print $4}' read-depth-stat.txt`
READ_DEPTH_Q3=`awk -F'\t' 'END {print $5}' read-depth-stat.txt`
READ_DEPTH_MAX=`awk -F'\t' 'END {print $6}' read-depth-stat.txt`
READ_DEPTH_IQR=`awk -F'\t' 'END {print $7}' read-depth-stat.txt`
READ_DEPTH_SSTDEV=`awk -F'\t' 'END {print $8}' read-depth-stat.txt`

$HVCHANGE Nextclade-result.tsv | \
   egrep 'clade|substitutions|deletions|insertions|aaSubstitutions|aaDeletions|total|errors|Nextclade_pango' \
   > Nextclade-result-data-list.txt
NEXTCLADECLADE=`awk -F'\t' '$1=="clade" {print $2}' Nextclade-result-data-list.txt`
NEXTCLADESUBST=`awk -F'\t' '$1=="substitutions" {print $2}' Nextclade-result-data-list.txt`
NEXTCLADEDELET=`awk -F'\t' '$1=="deletions" {print $2}' Nextclade-result-data-list.txt`
NEXTCLADEINSER=`awk -F'\t' '$1=="insertions" {print $2}' Nextclade-result-data-list.txt`
NEXTCLADEAASUB=`awk -F'\t' '$1=="aaSubstitutions" {print $2}' Nextclade-result-data-list.txt`
NEXTCLADEAADEL=`awk -F'\t' '$1=="aaDeletions" {print $2}' Nextclade-result-data-list.txt`
NEXTCLADETOTALGAP=`awk -F'\t' '$1=="totalDeletions" {print $2}' Nextclade-result-data-list.txt`
NEXTCLADETOTALINS=`awk -F'\t' '$1=="totalInsertions" {print $2}' Nextclade-result-data-list.txt`
NEXTCLADETOTALMIS=`awk -F'\t' '$1=="totalMissing" {print $2}' Nextclade-result-data-list.txt`
NEXTCLADETOTALMUT=`awk -F'\t' '$1=="totalSubstitutions" {print $2}' Nextclade-result-data-list.txt`
NEXTCLADETOTALNON=`awk -F'\t' '$1=="totalNonACGTNs" {print $2}' Nextclade-result-data-list.txt`
NEXTCLADEAASUBCNT=`awk -F'\t' '$1=="totalAminoacidSubstitutions" {print $2}' Nextclade-result-data-list.txt`
NEXTCLADEAADELCNT=`awk -F'\t' '$1=="totalAminoacidDeletions" {print $2}' Nextclade-result-data-list.txt`
PANGOLIN_LINEAGE=`awk -F'\t' '$1=="Nextclade_pango" {print $2}' Nextclade-result-data-list.txt`
PANGOLIN_LINEAGE=`if [ ${#PANGOLIN_LINEAGE} != 0 ]; then echo -e ${PANGOLIN_LINEAGE}; else echo -e "Unassingned"; fi` # 2023-01-27 add
NEXTCLADEERROR=`awk -F'\t' '$1=="errors" {print $2}' Nextclade-result-data-list.txt | perl -pe 's/"//g'`
NEXTCLADECLADE=`if [ ${#NEXTCLADECLADE} != 0 ]; then echo -e ${NEXTCLADECLADE}; \
                elif [ ${#NEXTCLADEERROR} != 0 ]; then echo -e ${NEXTCLADEERROR}; \
                else echo -e "No data"; fi`
rm Nextclade-result-data-list.txt

CHECK_NO_READ=`tail -n 1 read-depth-stat.txt | cut -f 1`
  if [ ${CHECK_NO_READ} == "mean" ]; then
        READ_DEPTH_MEDIAN="0"
        COVERAGE="0%"
        TOAL_SNV="0"
        SNV="0"
        MIX_SNV="0"
        MIX_SNV_RATIO_1="0"
        MIX_SNV_RATIO_2="0%"
        FINAL_CONTIG_STATUS="No data"
        FINAL_CONTIG_METHOD="No data"
        FINAL_CONTIG_LEN="0"
        SKESA_LEN="0"
        SKESA_NUM_CONTIGS="0"
        A5_LEN="0"
        A5_NUM_CONTIGS="0"
        INPUT_LEN="0"
        INPUT_NUM_CONTIGS="0"
        MEDIAN_INSERT_SIZE="0"
        MEAN_INSERT_SIZE="0"
        READ_DEPTH_MEAN="0"
        READ_DEPTH_MIN="0"
        READ_DEPTH_Q1="0"
        READ_DEPTH_MEDIAN="0"
        READ_DEPTH_Q3="0"
        READ_DEPTH_MAX="0"
        READ_DEPTH_IQR="0"
        READ_DEPTH_SSTDEV="0"
        NEXTCLADECLADE="No data"
        NEXTCLADEVERSION="ERROR"
        NEXTCLADESUBST="No data"
        NEXTCLADEDELET="No data"
        NEXTCLADEINSER="No data"
        NEXTCLADEAASUB="No data"
        NEXTCLADEAADEL="No data"
        NEXTCLADETOTALGAP="No data"
        NEXTCLADETOTALINS="No data"
        NEXTCLADETOTALMIS="No data"
        NEXTCLADETOTALMUT="No data"
        NEXTCLADETOTALNON="No data"
        NEXTCLADEAASUBCNT="No data"
        NEXTCLADEAADELCNT="No data"
        rm *vcf-consensus-complete.fasta
  fi

set +u
echo -e "{
  \"sequences\" : [
    {
      \"@project_name\" : \"$SAMPLENAME\",
      \"input_data\" : \"$INPUT_DATA\",
      \"read_depth_median\" : \"$READ_DEPTH_MEDIAN\",
      \"coverage\" : \"$COVERAGE\",
      \"total_detected_snv_sites\" : \"$TOAL_SNV\",
      \"allele_sites\" : \"$SNV\",
      \"mix_allele_sites\" : \"$MIX_SNV\",
      \"mix_allele_frequency\" : \"$MIX_SNV_RATIO_1\",
      \"mix_allele_frequency_in_total_coverage_length\" : \"$MIX_SNV_RATIO_2\",
      \"contig_status\" : \"$FINAL_CONTIG_STATUS\",
      \"assembly_method\" : \"$FINAL_CONTIG_METHOD\",
      \"total_contig_length\" : \"$FINAL_CONTIG_LEN\",
      \"skesa_total_contig_length\" : \"$SKESA_LEN\",
      \"num_skesa_contig\" : \"$SKESA_NUM_CONTIGS\",
      \"a5miseq_total_contig_length\" : \"$A5_LEN\",
      \"num_a5miseq_contig\" : \"$A5_NUM_CONTIGS\",
      \"input_contig_length\" : \"$INPUT_LEN\",
      \"num_input_contig\" : \"$INPUT_NUM_CONTIGS\",
      \"total_read\" : \"$TOTAL_READ\",
      \"trimmed_read\" : \"$TRIMMED_READ\",
      \"mapped_read\" : \"$MAPPED_READ\",
      \"trimmed_read_per\" : \"$TRIMMED_READ_RATIO\",
      \"mapped_read_per\" : \"$MAPPED_READ_RATIO\",
      \"median_insert_size\" : \"$MEDIAN_INSERT_SIZE\",
      \"mean_insert_size\" : \"$MEAN_INSERT_SIZE\",
      \"lineage_pangolin\" : \"$PANGOLIN_LINEAGE\",
      \"conflict_pangolin\" : \"$PANGOLIN_PROBABILITY\",
      \"read_depth_mean\" : \"$READ_DEPTH_MEAN\",
      \"read_depth_min\" : \"$READ_DEPTH_MIN\",
      \"read_depth_q1\" : \"$READ_DEPTH_Q1\",
      \"read_depth_q3\" : \"$READ_DEPTH_Q3\",
      \"read_depth_max\" : \"$READ_DEPTH_MAX\",
      \"read_depth_iqr\" : \"$READ_DEPTH_IQR\",
      \"read_depth_sstdev\" : \"$READ_DEPTH_SSTDEV\",
      \"clade\" : \"$NEXTCLADECLADE\",
      \"nextcladeversion\" : \"$NEXTCLADEVERSION\",
      \"total_gaps\" : \"$NEXTCLADETOTALGAP\",
      \"total_insertions\" : \"$NEXTCLADETOTALINS\",
      \"total_missing\" : \"$NEXTCLADETOTALMIS\",
      \"total_mutations\" : \"$NEXTCLADETOTALMUT\",
      \"total_non_acgtns\" : \"$NEXTCLADETOTALNON\",
      \"substitutions\" : \"$NEXTCLADESUBST\",
      \"deletions\" : \"$NEXTCLADEDELET\",
      \"insertions\" : \"$NEXTCLADEINSER\",
      \"total_amino_acid_substitutions\" : \"$NEXTCLADEAASUBCNT\",
      \"total_amino_acid_deletions\" : \"$NEXTCLADEAADELCNT\",
      \"aa_substitutions\" : \"$NEXTCLADEAASUB\",
      \"aa_deletions\" : \"$NEXTCLADEAADEL\"     
    }
  ]
}" > $SAMPLENAME.json

echo -e "project_name : $SAMPLENAME
Input data : $INPUT_DATA
Read depth (median) : $READ_DEPTH_MEDIAN
Coverage (%) (ref length:$TOTALLEN) : $COVERAGE
Total detected SNV sites : $TOAL_SNV
Allele sites (≥ 80%) : $SNV
Mix allele sites (< 80%) : $MIX_SNV
Mix allele frequency : $MIX_SNV_RATIO_1
Mix allele frequency in total coverage length : $MIX_SNV_RATIO_2
Contig status : $FINAL_CONTIG_STATUS
Assembly method : $FINAL_CONTIG_METHOD
Total contig length (bp) : $FINAL_CONTIG_LEN
SKESA total contig length (bp) : $SKESA_LEN
Num. SKESA contig : $SKESA_NUM_CONTIGS
A5miseq total contig length (bp) : $A5_LEN
Num. A5miseq contig : $A5_NUM_CONTIGS
Input contig length (bp) : $INPUT_LEN
Num. input contig : $INPUT_NUM_CONTIGS
Total read : $TOTAL_READ
Trimmed read : $TRIMMED_READ
Mapped read : $MAPPED_READ
Trimmed read (%) : $TRIMMED_READ_RATIO
Mapped read (%) : $MAPPED_READ_RATIO
Median insert size : $MEDIAN_INSERT_SIZE
Mean insert size : $MEAN_INSERT_SIZE
Lineage (pangolin) : $PANGOLIN_LINEAGE
Conflict (pangolin) : $PANGOLIN_PROBABILITY
Read depth (mean) : $READ_DEPTH_MEAN
Read depth (min) : $READ_DEPTH_MIN
Read depth (q1) : $READ_DEPTH_Q1
Read depth (median) : $READ_DEPTH_MEDIAN
Read depth (q3) : $READ_DEPTH_Q3
Read depth (max) : $READ_DEPTH_MAX
Read depth (iqr) : $READ_DEPTH_IQR
Read depth (sstdev) : $READ_DEPTH_SSTDEV
Clade : $NEXTCLADECLADE
NextCladeVersion : $NEXTCLADEVERSION
Total Gaps : $NEXTCLADETOTALGAP
Total Insertions : $NEXTCLADETOTALINS
Total Missing : $NEXTCLADETOTALMIS
Total Mutations : $NEXTCLADETOTALMUT
Total NonACGTNs : $NEXTCLADETOTALNON
Substitutions : $NEXTCLADESUBST
Deletions : $NEXTCLADEDELET
Insertions : $NEXTCLADEINSER
Total Amino acid Substitutions : $NEXTCLADEAASUBCNT
Total Amino acid Deletions : $NEXTCLADEAADELCNT
AA substitutions : $NEXTCLADEAASUB
AA deletions : $NEXTCLADEAADEL" > result-summary.txt
set -u
 
#########################
##### Download data #####
#########################

mkdir $SAMPLENAME-result
cp result-summary.txt $SAMPLENAME-result

mkdir $SAMPLENAME-result/mapping-data-and-mutation-site
cp VarScan-modified.vcf $SAMPLENAME-result/mapping-data-and-mutation-site
cp VarScan-modified-rm-Mix-allele.vcf $SAMPLENAME-result/mapping-data-and-mutation-site
egrep -v 'GATKCommandLine|reference=file' GATK-indel-filtered-ref.vcf > $SAMPLENAME-result/mapping-data-and-mutation-site/GATK-indel-filtered-ref.vcf
cp realigned_reads.bam $SAMPLENAME-result/mapping-data-and-mutation-site
cp realigned_reads.bam.bai $SAMPLENAME-result/mapping-data-and-mutation-site
cp BED-depth-check.txt $SAMPLENAME-result/mapping-data-and-mutation-site
cp vcf-consensus.fasta $SAMPLENAME-result/mapping-data-and-mutation-site
set +e
cp vcf-consensus-with-indel-trim.fasta $SAMPLENAME-result/mapping-data-and-mutation-site
cp large-deletion.txt $SAMPLENAME-result/mapping-data-and-mutation-site
set -e
cp $REFSEQ $SAMPLENAME-result/mapping-data-and-mutation-site/ref-seq_MN908947.fasta

set +e
mkdir $SAMPLENAME-result/contigs
cp extract-contigs-a5miseq.fasta $SAMPLENAME-result/contigs
cp extract-contigs-skesa.fasta $SAMPLENAME-result/contigs
cp a5-miseq-assemble.contigs.fasta $SAMPLENAME-result/contigs
cp skesa-assembled-contigs-default.fasta $SAMPLENAME-result/contigs
cp *table $SAMPLENAME-result/contigs

mv $SAMPLENAME-result/contigs $SAMPLENAME-result/contigs-
mkdir -p $SAMPLENAME-result/contigs/source_before_modified
mv $SAMPLENAME-result/contigs-/* $SAMPLENAME-result/contigs/source_before_modified
rmdir $SAMPLENAME-result/contigs-
mv $SAMPLENAME-result/contigs/source_before_modified/extr* $SAMPLENAME-result/contigs
set -e

cp $SAMPLENAME*.fasta $SAMPLENAME-result/

mkdir $SAMPLENAME-result/other
cp lineage_report.csv $SAMPLENAME-result/other
cp insert_size_histogram.pdf $SAMPLENAME-result/other
cp plot_depth.pdf $SAMPLENAME-result/other
cp Nextclade-result.tsv $SAMPLENAME-result/other

zip -rm $SAMPLENAME-result.zip  $SAMPLENAME-result/

#####################################
##### remove intermediate files #####
#####################################

echo -e "Total read : $TOTAL_READ
Trimmed read : $TRIMMED_READ
Mapped read : $MAPPED_READ" > read-count.txt

set +e
rm trim-ARTIC-primer_* R[12]_trim-sickle.fastq.gz out12_mem.sorted* single.fastq.gz
set -e


echo -e "\nINFO `date`   SARS-CoV-2 analysis finish.\n"
