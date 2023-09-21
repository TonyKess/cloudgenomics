# Setting up

First we use azcopy to move all of our genomic data to a directory 
```
azcopy login
azcopy copy <reads folder of raw genomic data for project> ./
```
Next we make a file of Ids to pass through the software pipeline, these map to each individual we sequenced for R1 and R2 reads on an Illumina sequencer

```
cd reads
ls *fastq.gz | sed 's/\_..fastq.gz//' | sort | uniq > ../SRIds
```
Now we use [screen](https://www.gnu.org/software/screen/) to ensure processes run even when we are disconnected from the virtual machine 
```
screen -S genotyping
```

We add some variables to our environment
```
gatk37=/datadisk0/software/GenomeAnalysisTK.jar
genome=<genome>
export gatk37
export genome
```
# Data preprocessing
Activate our trimming environment and carry out adapter trimming
```
mamba activate fastp

cat Ids | \
parallel -j 40 \
  'fastp -i reads/{}_R1.fastq.gz \
  -I reads/{}_R2.fastq.gz \
  -o trim/{}_R1.trimmed.fastq.gz \
  -O trim/{}_R2.trimmed.fastq.gz \
  --thread 4 '
mamba deactivate
```
Now we carry out alignment, first indexing our genome
```
mamba activate align
mkdir align

bwa-mem2 index $genome
samtools faidx $genome
```
Then aligning our trimmed reads to the indexed genome

```
while read ind;
  do echo $ind\.bam ;
  RGID=$(echo $ind |  sed 's/i5.*/i5/') ;
  SMID=$(echo $ind | sed 's/NS.*i5.//') ;
  LBID=$(echo $ind | sed 's/.UDP.*//');
  bwa-mem2 mem \
  -t 200 \
  -R "@RG\tID:$RGID\tSM:$SMID\tLB:$LBID" \
  $genome \
  trim/$ind\_R1.trimmed.fastq.gz  trim/$ind\_R2.trimmed.fastq.gz | \
  samtools sort -o align/$ind\.sorted.bam -T $ind -@ 100 -m 20G ;
  done < Ids
```
Now we conduct deduplication using GATK

```
mamba activate gatk4
mkdir /datadisk0/gatktemp

cat Ids | \
  parallel --tmpdir /datadisk0/gatktemp --jobs 40 \
  'gatk --java-options "-Xmx3G" \
   MarkDuplicates \
   I=align/{}.sorted.bam \
   O=align/{}.deDup.bam M=align/{}_deDupMetrics.txt \
   REMOVE_DUPLICATES=true \
   TMP_DIR=/datadisk0/gatktemp '
```
We now prepare for indel realignment, the final step in preprocessing of our genomic data.
First we make a sequence dictionary
```
gatk CreateSequenceDictionary -R $genome
mamba deactivate
```

Index our dedpuplicated reads:
```
mamba activate align
cat Ids | \
  parallel --jobs 180  'samtools index align/{}.deDup.bam '

mamba deactivate
```

And carry out the two steps of indel realignment

```
mamba activate gatk37

cat Ids | parallel --tmpdir /datadisk0/gatktemp \
--jobs 180 ' java -Xmx5g  \
-Djava.io.tmpdir=/datadisk0/gatktemp \
-jar $gatk37 \
-T RealignerTargetCreator \
-R $genome \
-I align/{}.deDup.bam \
-o align/{}.intervals'


cat Ids | parallel --tmpdir /datadisk0/gatktemp \
  --jobs 180 ' java -Xmx5g  \
  -Djava.io.tmpdir=/datadisk0/gatktemp \
  -jar $gatk37 \
  -T IndelRealigner \
  -R $genome \
  -I align/{}.deDup.bam \
  -targetIntervals align/{}.intervals \
  -o align/{}.realigned.bam '

mamba deactivate
```

# Genotype calling and genotype likelihoods

This section assumes we have made a chromosome list containing chromosomes of interest. These can be obtained from the genome.fasta.fai file that was output by samtools faidx, or manually selected in a text editor. We also ened a list of our bam files that have been deduplicated and aligned

```
ls align/*real*bam > bamlist.tsv
```
We then reactivate our alignment environment that contains ANGSD, which we can use for genotyping. bcftools, freebayes, and GATK are also options - this pipeline uses ANGSD due to contingency. We had lcWGS data first, and so ANGSD was previously a big part of our toolset. 

```
mamba activate align
ls

cat Chroms.tsv | parallel --jobs 25 ' angsd -nThreads 8 \
-dobcf 1 \
-gl 1 \
-dopost 1 \
-dogeno 5 \
-domajorminor 1 \
-doMaf 1 \
-docounts 1 \
-dumpCounts 2 \
-doQsDist 1 \
-minMapQ 30 \
-minQ 30 \
-minInd 650 -setMinDepth 1000 \
-SNP_pval 2e-6 \
-uniqueOnly 1 \
-only_proper_pairs 1 \
-remove_bads 1 \
-minMaf 0.05 \
-r {}: \
-out angsd_out/{} \
-bam bamlist.tsv '
```

We can now take our output bcf files and convert to vcf
```
ls angsd_out/*bcf | sed 's/.bcf//' | parallel 'bcftools view {}.bcf -O z -o {}.vcf.gz '
```
And index them for any downstream analyses we are interested in

```
ls angsd_out/*bcf | sed 's/.bcf//' | parallel 'bcftools view {}.bcf -O z -o {}.vcf.gz '
ls angsd_out/*vcf.gz | parallel ' tabix {} '
```
From here we can work directly with the variable-coverage samples, or we can do imputation.

