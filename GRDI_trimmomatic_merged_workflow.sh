#!/bin/bash
#$ -S /bin/bash
#$ -m eas -M galen.guo@canada.ca
#$ -o /isilon/ottawa-rdc/users/shared/chenw_lab/galen/temp_files
#$ -e /isilon/ottawa-rdc/users/shared/chenw_lab/galen/temp_files

source ~/miniconda3/bin/activate
conda activate anvio-7

#directories

GALEN=/isilon/ottawa-rdc/users/shared/chenw_lab/galen
WORK=$GALEN/GRDI
RAW=$WORK/raw_combined
origin=$WORK/raw_data
output=$WORK/output
input=$WORK/input

TRIM=$WORK/output/trimmomatic_merged

MEGAHIT2016=$output/megahit_2016
MEGAHIT2017=$output/megahit_2017
MEGAHIT2018=$output/megahit_2018

MAP2016=$output/mapping_2016
MAP2017=$output/mapping_2017
MAP2018=$output/mapping_2018

ANVIO2016=$output/anvio2016
ANVIO2017=$output/anvio2017
ANVIO2018=$output/anvio2018

####################################################################
### QC of fastqc pre-trim
####################################################################


cd $RAW
mkdir fastqc
fastqc $RAW/*fq.gz -t 10 -o $RAW/fastqc/ -f fastq
multiqc $RAW/fastqc/ -o $RAW/fastqc/

###################################################################
### removal of adapter using trimmomatic
###################################################################

echo Trimming
mkdir $output/trimmomatic_merged


for i in `awk '{print $1}' $RAW/sample_name`;
do
echo $i
    file1=$input/$i_R1.fq.gz
    file2=$input/$i_R2.fq.gz
	java -jar $GALEN/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 20 $file1 $file2 $TRIM/i.pair1.fq.gz $TRIM/$i.unpair1.fq.gz $TRIM/$i.pair2.fq.gz $TRIM/$i.unpair2.fq.gz ILLUMINACLIP:$GALEN/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 HEADCROP:6 CROP:120
done


echo Trimming_done

####################################################################
### QC of fastqc post-trim
####################################################################
echo QC_of_Trimming

mkdir $TRIM/fastqc
fastqc -t 20 $TRIM/*fq.gz -o $TRIM/fastqc/
multiqc $TRIM/fastqc/*.pair*

echo QC_of_Trimming_done

## move files to their respective "year" folder

mkdir $TRIM/GRDI_2016
mkdir $TRIM/GRDI_2017
mkdir $TRIM/GRDI_2018

mv $TRIM/*2016* $TRIM/GRDI_2016/
mv $TRIM/*2017* $TRIM/GRDI_2017/
mv $TRIM/*2018* $TRIM/GRDI_2018/

####################################################################
### co-assembly megahit
### 
####################################################################


# 2016
echo co-assembling 2016
 
F=$(echo $(ls $TRIM/GRDI_2016/*.pair1.fq.gz)  | sed "s/ /,/g")
R=$(echo $(ls $TRIM/GRDI_2016/*.pair2.fq.gz)  | sed "s/ /,/g")
S=$(echo $(ls $TRIM/GRDI_2016/*.unpair*)  | sed "s/ /,/g")

MEGAHIT2016=$output/megahit_2016

megahit -1 $F -2 $R -r $S -o $MEGAHIT2016 -t 50 --min-contig-len 1000 -m 0.85 --continue

echo co-assembling_done

cp $MEGAHIT2016/final.contigs.fa $MEGAHIT2016/contigs2016.fa


# 2017
echo co-assembling 2017

F=$(echo $(ls $TRIM/GRDI_2017/*.pair1.fq.gz)  | sed "s/ /,/g")
R=$(echo $(ls $TRIM/GRDI_2017/*.pair2.fq.gz)  | sed "s/ /,/g")
S=$(echo $(ls $TRIM/GRDI_2017/*.unpair*)  | sed "s/ /,/g")

MEGAHIT2017=$output/megahit_2017

megahit -1 $F -2 $R -r $S -o $MEGAHIT2017 -t 50 --min-contig-len 1000 -m 0.85 --continue

echo co-assembling_done

cp $MEGAHIT2017/final.contigs.fa $MEGAHIT2017/contigs2017.fa


echo co-assembling 2018

F=$(echo $(ls $TRIM/GRDI_2018/*.pair1.fq.gz)  | sed "s/ /,/g")
R=$(echo $(ls $TRIM/GRDI_2018/*.pair2.fq.gz)  | sed "s/ /,/g")
S=$(echo $(ls $TRIM/GRDI_2018/*.unpair*)  | sed "s/ /,/g")

MEGAHIT2018=$output/megahit_2018

megahit -1 $F -2 $R -r $S -o $MEGAHIT2018 -t 50 --min-contig-len 1000 -m 0.85 --continue

echo co-assembling_done

cp $MEGAHIT2018/final.contigs.fa $MEGAHIT2018/contigs2018.fa



####################################################################
### loading into anvio SEPERATED BY YEAR!
### 
####################################################################

#2016

ANVIO2016=$output/anvio2016

mkdir $ANVIO2016

anvi-script-reformat-fasta $MEGAHIT2016/contigs2016.fa -o $ANVIO2016/contigs2016.fa -l 1000 –report-file --simplify-names

###create anvio database
GRDI	

anvi-gen-contigs-database -f $ANVIO2016/contigs2016.fa  -o $ANVIO2016/contigs2016.db -n GRDI2016 -T 25

anvi-run-hmms -c $ANVIO2016/contigs2016.db -T 25

#run once:  anvi-setup-kegg-kofams --kegg-data-dir $GALEN/KOFAM
#run once: anvi-setup-ncbi-cogs --cog-data-dir $GALEN/COG --num-threads 25

anvi-run-kegg-kofams -c $ANVIO2016/contigs2016.db -T 25 --hmmer-program hmmsearch --kegg-data-dir $GALEN/KOFAM --keep-all-hits --just-do-it

anvi-run-ncbi-cogs -c $ANVIO2016/contigs2016.db --cog-data-dir $GALEN/COG --num-threads 20 --search-with diamond --temporary-dir-path $GALEN/temp


anvi-get-sequences-for-gene-calls -c $ANVIO2016/contigs2016.db \
                                    --get-aa-sequences \
                                    -o $ANVIO2016/aa_seq_2016.fa

anvi-display-contigs-stats $ANVIO2016/contigs2016.db  --report-as-text  --output-file $ANVIO2016/contigs2016_post-hmm-cogs_kofams.txt


#2017

ANVIO2017=$output/anvio2017

mkdir $ANVIO2017

anvi-script-reformat-fasta $MEGAHIT2017/contigs2017.fa -o $ANVIO2017/contigs2017.fa -l 0 --simplify-names

###create anvio database
GRDI	

anvi-gen-contigs-database -f $ANVIO2017/contigs2017.fa  -o $ANVIO2017/contigs2017.db -n GRDI2017 -T 25

anvi-run-hmms -c $ANVIO2017/contigs2017.db -T 25

anvi-run-kegg-kofams -c $ANVIO2017/contigs2017.db -T 25 --hmmer-program hmmsearch --keep-all-hits --kegg-data-dir $GALEN/KOFAM --just-do-it

anvi-run-ncbi-cogs -c $ANVIO2017/contigs2017.db --cog-data-dir $GALEN/COG --num-threads 25 --search-with diamond --temporary-dir-path $ANVIO/temp

anvi-get-sequences-for-gene-calls -c $ANVIO2017/contigs2017.db \
                                    --get-aa-sequences \
                                    -o $ANVIO2017/aa_seq_2017.fa

anvi-display-contigs-stats $ANVIO2017/contigs2017.db  --report-as-text  --output-file $ANVIO2017/contigs2017_post-hmm-cogs_kofams.txt


#2018

ANVIO2018=$output/anvio2018

mkdir $ANVIO2018

anvi-script-reformat-fasta $MEGAHIT2018/contigs2018.fa -o $ANVIO2018/contigs2018.fa -l 0 --simplify-names

###create anvio database
GRDI	

anvi-gen-contigs-database -f $ANVIO2018/contigs2018.fa  -o $ANVIO2018/contigs2018.db -n GRDI2018 -T 25

anvi-run-hmms -c $ANVIO2018/contigs2018.db -T 25

anvi-run-kegg-kofams -c $ANVIO2018/contigs2018.db -T 25 --hmmer-program hmmsearch --keep-all-hits --kegg-data-dir $GALEN/KOFAM --just-do-it 

anvi-run-ncbi-cogs -c $ANVIO2018/contigs2018.db --cog-data-dir $GALEN/COG --num-threads 20 --search-with diamond --temporary-dir-path $ANVIO/temp

anvi-get-sequences-for-gene-calls -c $ANVIO2018/contigs2018.db \
                                    --get-aa-sequences \
                                    -o $ANVIO2018/aa_seq_2018.fa

anvi-display-contigs-stats $ANVIO2018/contigs2018.db  --report-as-text  --output-file $ANVIO2018/contigs2018_post-hmm-cogs_kofams.txt



####################################################################
### mapping with bowtie2
### 
####################################################################

echo mapping_2016

## create mapping directory

MAP2016=$output/mapping_2016
mkdir $MAP2016
cd $MAP2016

## building mapping files

bowtie2-build $ANVIO2016/contigs2016.fa $MAP2016/map2016_contigs --threads 25

### mapping read to contigs

cd $TRIM


for sample in `awk '{print $1}' $TRIM/sample_name`;
do
R1="${sample}.pair1.fq.gz"
R2="${sample}.pair2.fq.gz"
sam="${sample}.sam" 
bowtie2 -x $MAP2016/map2016_contigs -1 GRDI_2016/$R1 -2 GRDI_2016/$R2 -S $MAP2016/$sam --threads 25
done


conda deactivate

## echo convert sam to bam

# transfering sample name file to mapping folder.

cp $RAW/sample_name $MAP2016

for sample in `awk '{print $1}' $MAP2016/sample_name`;
do
sam=".sam"
bam=".bam"
samtools view -S -b $MAP2016/${sample}.sam > $MAP2016/${sample}.bam
done


conda activate anvio-7


## bam to anvio_bam profiling

for sample in `awk '{print $1}' $MAP2016/sample_name`;
do
anvi-init-bam $MAP2016/${sample}.bam -o $MAP2016/${sample}_anvi.bam
done


echo mapping_done




echo mapping_2017

## create mapping directory

MAP2017=$output/mapping_2017
mkdir $MAP2017
cd $MAP2017



## building mapping files

bowtie2-build $ANVIO2017/contigs2017.fa $MAP2017/map2017_contigs --threads 25

### mapping read to contigs

cd $TRIM


for sample in `awk '{print $1}' $TRIM/sample_name`;
do
R1="${sample}.pair1.fq.gz"
R2="${sample}.pair2.fq.gz"
sam="${sample}.sam" 
bowtie2 -x $MAP2017/map2017_contigs -1 GRDI_2017/$R1 -2 GRDI_2017/$R2 -S $MAP2017/$sam --threads 25
done


conda deactivate

## echo convert sam to bam

# transfering sample name file to mapping folder.

cp $RAW/sample_name $MAP2017

for sample in `awk '{print $1}' $MAP2017/sample_name`;
do
sam=".sam"
bam=".bam"
samtools view -S -b $MAP2017/${sample}.sam > $MAP2017/${sample}.bam
done


conda activate anvio-7


## bam to anvio_bam profiling

for sample in `awk '{print $1}' $MAP2017/sample_name`;
do
anvi-init-bam $MAP2017/${sample}.bam -o $MAP2017/${sample}_anvi.bam
done


echo mapping_done




echo mapping_2018

## create mapping directory

MAP2018=$output/mapping_2018
mkdir $MAP2018
cd $MAP2018



## building mapping files

bowtie2-build $ANVIO2018/contigs2018.fa $MAP2018/map2018_contigs --threads 25

cd $TRIM


for sample in `awk '{print $1}' $TRIM/sample_name`;
do
R1="${sample}.pair1.fq.gz"
R2="${sample}.pair2.fq.gz"
sam="${sample}.sam" 
bowtie2 -x $MAP2018/map2018_contigs -1 GRDI_2018/$R1 -2 GRDI_2018/$R2 -S $MAP2018/$sam --threads 25
done


conda deactivate

## echo convert sam to bam

# transfering sample name file to mapping folder.

cp $RAW/sample_name $MAP2018

for sample in `awk '{print $1}' $MAP2018/sample_name`;
do
sam=".sam"
bam=".bam"
samtools view -S -b $MAP2018/${sample}.sam > $MAP2018/${sample}.bam
done


conda activate anvio-7


## bam to anvio_bam profiling

for sample in `awk '{print $1}' $MAP2018/sample_name`;
do
anvi-init-bam $MAP2018/${sample}.bam -o $MAP2018/${sample}_anvi.bam
done


echo mapping_done



####################################################################
### Profiling cont'd very long. bam --> anvio
### PROFILE DONT LIKE SAMPLE NAME WITH "-" 
####################################################################

### change bam file name from "-" to "_"


cd $MAP2016
find . -depth -name '*-*' -exec rename '-' '_' {} +


## ### do everything in one file, using more cpu. (5 days)

mkdir $ANVIO2016/profile
cp $TRIM/sample_name $ANVIO2016/profile


file_ext="_anvi.bam"
for sample in `awk '{print $1}' $ANVIO2016/profile/sample_name`;
do
echo "$sample$file_ext"
anvi-profile -i $MAP2016/"$sample$file_ext" -c $ANVIO2016/contigs2016.db --output-dir $ANVIO2016/profile/$sample --sample-name $sample -T 100 --min-contig-length 1000
done

## merge profile

mkdir $ANVIO2016/profile_merged_2016
anvi-merge  $ANVIO2016/profile/*/PROFILE.db -o $ANVIO2016/profile_merged_2016 -c $ANVIO2016/contigs2016.db -S GRDI2016 -W


cd $MAP2017
find . -depth -name '*-*' -exec rename '-' '_' {} +


## do everything in one file, using more cpu.

mkdir $ANVIO2017/profile
cp $TRIM/sample_name $ANVIO2017/profile


file_ext="_anvi.bam"
for sample in `awk '{print $1}' $ANVIO2017/profile/sample_name`;
do
echo "$sample$file_ext"
anvi-profile -i $MAP2017/"$sample$file_ext" -c $ANVIO2017/contigs2017.db --output-dir $ANVIO2017/profile/$sample --sample-name $sample -T 100 --min-contig-length 1000
done

## merge profile
mkdir $ANVIO2017/profile_merged_2017
anvi-merge  $ANVIO2017/profile/*/PROFILE.db -o $ANVIO2017/profile_merged_2017 -c $ANVIO2017/contigs2017.db -S GRDI2017 -W


cd $MAP2018
find . -depth -name '*-*' -exec rename '-' '_' {} +


## ### do everything in one file, using more cpu. (5 days)

mkdir $ANVIO2018/profile
cp $TRIM/sample_name $ANVIO2018/profile


file_ext="_anvi.bam"
for sample in `awk '{print $1}' $ANVIO2018/profile/sample_name`;
do
echo "$sample$file_ext"
anvi-profile -i $MAP2018/"$sample$file_ext" -c $ANVIO2018/contigs2018.db --output-dir $ANVIO2018/profile/$sample --sample-name $sample -T 100 --min-contig-length 1000
done

## merge profile
mkdir $ANVIO2018/profile_merged_2018
anvi-merge  $ANVIO2018/profile/*/PROFILE.db -o $ANVIO2018/profile_merged_2018 -c $ANVIO2018/contigs2018.db -S GRDI2018 -W


####################################################################
### BINNING! Finally
### Concoct
####################################################################

anvi-cluster-contigs -p $ANVIO2016/profile_merged_2016/PROFILE.db -c $ANVIO2016/contigs2016.db -C concoct_2016 --driver concoct -T 50 --just-do-it
anvi-cluster-contigs -p $ANVIO2017/profile_merged_2017/PROFILE.db -c $ANVIO2017/contigs2017.db -C concoct_2017 --driver concoct -T 50 --just-do-it
anvi-cluster-contigs -p $ANVIO2018/profile_merged_2018/PROFILE.db -c $ANVIO2018/contigs2018.db -C concoct_2018 --driver concoct -T 50 --just-do-it


####################################################################
### metabat2
####################################################################


anvi-cluster-contigs -p $ANVIO2016/profile_merged_2016/PROFILE.db -c $ANVIO2016/contigs2016.db  -C metabat2_2016 --driver metabat2 -T 50 --just-do-it
anvi-cluster-contigs -p $ANVIO2016/profile_merged_2017/PROFILE.db -c $ANVIO2016/contigs2017.db  -C metabat2_2017 --driver metabat2 -T 50 --just-do-it
anvi-cluster-contigs -p $ANVIO2016/profile_merged_2018/PROFILE.db -c $ANVIO2016/contigs2018.db  -C metabat2_2018 --driver metabat2 -T 50 --just-do-it

####################################################################
### maxbin2
####################################################################

anvi-cluster-contigs -p $ANVIO2016/profile_merged_2016/PROFILE.db -c $ANVIO2016/contigs2016.db  -C maxbin2_2016 --driver maxbin2 -T 50 --just-do-it
anvi-cluster-contigs -p $ANVIO2016/profile_merged_2017/PROFILE.db -c $ANVIO2016/contigs2017.db  -C maxbin2_2017 --driver maxbin2 -T 50 --just-do-it
anvi-cluster-contigs -p $ANVIO2016/profile_merged_2018/PROFILE.db -c $ANVIO2016/contigs2018.db  -C maxbin2_2018 --driver maxbin2 -T 50 --just-do-it


####################################################################
### dastool
####################################################################


anvi-cluster-contigs -p $ANVIO2016/profile_merged_2016/PROFILE.db -c $ANVIO2016/contigs2016.db -S concoct_2016,metabat2_2016 --driver dastool --search-engine usearch -C dastool_2016 -T 50 --just-do-it
anvi-cluster-contigs -p $ANVIO2017/profile_merged_2017/PROFILE.db -c $ANVIO2017/contigs2017.db -S concoct_2017,metabat2_2017 --driver dastool --search-engine usearch -C dastool_2017 -T 50 --just-do-it
anvi-cluster-contigs -p $ANVIO2018/profile_merged_2018/PROFILE.db -c $ANVIO2018/contigs2018.db -S concoct_2018,metabat2_2018 --driver dastool -C dastool_2018 -T 50 --just-do-it


####################################################################
### summarize!
####################################################################


anvi-summarize -p $ANVIO2016/profile_merged_2016/PROFILE.db -c $ANVIO2016/contigs2016.db -o $ANVIO2016/sample_summary_dastool -C dastool_2016 --init-gene-coverages
anvi-summarize -p $ANVIO2017/profile_merged_2017/PROFILE.db -c $ANVIO2017/contigs2017.db -o $ANVIO2017/sample_summary_dastool -C dastool_2017 --init-gene-coverages
anvi-summarize -p $ANVIO2018/profile_merged_2018/PROFILE.db -c $ANVIO2018/contigs2018.db -o $ANVIO2018/sample_summary_dastool -C dastool_2018 --init-gene-coverages

####################################################################
### refine!!
####################################################################

anvi-refine ...

####################################################################
### estimate metabolism (KEGG)
####################################################################

### rename bin and remove bin that has < 50% completion
## import collection in txt
## manually remove unwanted bin, and rename bins 
cd $ANVIO2016/
anvi-export-collection -C dastool_2016 -p $ANVIO2016/profile_merged_2016/PROFILE.db 

anvi-import-collection $ANVIO2016/collection-dastool_2016.txt -C FINAL -p $ANVIO2016/profile_merged_2016/PROFILE.db -c $ANVIO2016/contigs2016.db --bins-info $ANVIO2016/collection-dastool_2016-info.txt

cd $ANVIO2017/
anvi-export-collection -C dastool_2017 -p $ANVIO2017/profile_merged_2017/PROFILE.db 

anvi-import-collection $ANVIO2017/collection-dastool_2017.txt -C FINAL -p $ANVIO2017/profile_merged_2017/PROFILE.db -c $ANVIO2017/contigs2017.db --bins-info $ANVIO2017/collection-dastool_2017-info.txt

cd $ANVIO2018/
anvi-export-collection -C dastool_2018 -p $ANVIO2018/profile_merged_2018/PROFILE.db 

anvi-import-collection $ANVIO2018/collection-dastool_2018.txt -C FINAL -p $ANVIO2018/profile_merged_2018/PROFILE.db -c $ANVIO2018/contigs2018.db --bins-info $ANVIO2018/collection-dastool_2018-info.txt

# estimate metabolism (KEGG) of collection, must run anvi-setup-kegg-kofams first


anvi-estimate-metabolism -c $ANVIO2016/contigs2016.db -p $ANVIO2016/profile_merged_2016/PROFILE.db  -C FINAL --kegg-data-dir $GALEN/KOFAM/
anvi-estimate-metabolism -c $ANVIO2017/contigs2017.db -p $ANVIO2017/profile_merged_2017/PROFILE.db  -C FINAL --kegg-data-dir $GALEN/KOFAM/
anvi-estimate-metabolism -c $ANVIO2018/contigs2018.db -p $ANVIO2018/profile_merged_2018/PROFILE.db  -C FINAL --kegg-data-dir $GALEN/KOFAM/

####################################################################
### estimate taxonomy of bins (using GTDB)
####################################################################
conda install -c bioconda diamond
anvi-setup-scg-taxonomy #needs diamond (conda install diamond)

anvi-run-scg-taxonomy -c $ANVIO2016/contigs2016.db -T 100
anvi-run-scg-taxonomy -c $ANVIO2017/contigs2017.db -T 100
anvi-run-scg-taxonomy -c $ANVIO2018/contigs2018.db -T 100



cd $ANVIO2016/
anvi-estimate-scg-taxonomy -c $ANVIO2016/contigs2016.db -C FINAL -p $ANVIO2016/profile_merged_2016/PROFILE.db -T 100
cd $ANVIO2017/
anvi-estimate-scg-taxonomy -c $ANVIO2017/contigs2017.db -C FINAL -p $ANVIO2017/profile_merged_2017/PROFILE.db -T 100
cd $ANVIO2018/
anvi-estimate-scg-taxonomy -c $ANVIO2018/contigs2018.db -C FINAL -p $ANVIO2018/profile_merged_2018/PROFILE.db -T 100

####################################################################
### summarize! again..
####################################################################

anvi-summarize -p $ANVIO2016/profile_merged_2016/PROFILE.db -c $ANVIO2016/contigs2016.db -o $ANVIO2016/sample_summary_FINAL -C FINAL --init-gene-coverages
anvi-summarize -p $ANVIO2017/profile_merged_2017/PROFILE.db -c $ANVIO2017/contigs2017.db -o $ANVIO2017/sample_summary_FINAL -C FINAL --init-gene-coverages
anvi-summarize -p $ANVIO2018/profile_merged_2018/PROFILE.db -c $ANVIO2018/contigs2018.db -o $ANVIO2018/sample_summary_FINAL -C FINAL --init-gene-coverages


####################################################################
### checkm
####################################################################

conda activate 

checkm lineage_wf $ANVIO2016/sample_summary_FINAL/bin_by_bin/bins_fa $ANVIO2016/sample_summary_FINAL//checkm -x fa -t 25
checkm lineage_wf $ANVIO2017/sample_summary_FINAL/bin_by_bin/bins_fa $ANVIO2017/sample_summary_FINAL/checkm -x fa -t 25
checkm lineage_wf $ANVIO2018/sample_summary_FINAL/bin_by_bin/bins_fa $ANVIO2018/sample_summary_FINAL/checkm -x fa -t 25

#ssu finder?

checkm ssu_finder $ANVIO2016/contigs2016.fa $ANVIO2016/sample_summary_FINAL/bin_by_bin/bins_fa $ANVIO2016/sample_summary_FINAL/ssu_finder -x fa
checkm ssu_finder $ANVIO2017/contigs2017.fa $ANVIO2017/sample_summary_FINAL/bin_by_bin/bins_fa $ANVIO2017/sample_summary_FINAL/ssu_finder -x fa
checkm ssu_finder $ANVIO2018/contigs2018.fa $ANVIO2018/sample_summary_FINAL/bin_by_bin/bins_fa $ANVIO2018/sample_summary_FINAL/ssu_finder -x fa


####################################################################
### gtdb-tk
####################################################################
conda activate gtdb-tk

GTDBTK_DATA_PATH="/isilon/ottawa-rdc/users/shared/chenw_lab/galen/GTDB_release95/"

gtdbtk de_novo_wf --genome_dir $ANVIO2016/sample_summary_dastool_refine/bin_by_bin/bins_fa --bacteria -x fa --outgroup_taxon p__Chlamydiae --out_dir $ANVIO2016/sample_summary_dastool_refine/gtdb_output --cpus 20
gtdbtk ani_rep --genome_dir $ANVIO2016/sample_summary_dastool_refine/bin_by_bin/bins_fa --out_dir $ANVIO2016/sample_summary_dastool_refine/gtdb_output -x fa --cpus 20
gtdbtk classify_wf --genome_dir $ANVIO2016/sample_summary_dastool_refine/bin_by_bin/bins_fa --out_dir $ANVIO2016/sample_summary_dastool_refine/gtdb_output --outgroup_taxon p__Chlamydiae -x fa --cpus 20

gtdbtk de_novo_wf --genome_dir $ANVIO2017/sample_summary_dastool_refine/bin_by_bin/bins_fa --bacteria -x fa --outgroup_taxon p__Chlamydiae --out_dir $ANVIO2017/sample_summary_dastool_refine/gtdb_output --cpus 20
gtdbtk ani_rep --genome_dir $ANVIO2017/sample_summary_dastool_refine/bin_by_bin/bins_fa --out_dir $ANVIO2017/sample_summary_dastool_refine/gtdb_output -x fa --cpus 20
gtdbtk classify_wf --genome_dir $ANVIO2017/sample_summary_dastool_refine/bin_by_bin/bins_fa --out_dir $ANVIO2017/sample_summary_dastool_refine/gtdb_output --outgroup_taxon p__Chlamydiae -x fa --cpus 20

gtdbtk de_novo_wf --genome_dir $ANVIO2018/sample_summary_dastool_refine/bin_by_bin/bins_fa --bacteria -x fa --outgroup_taxon p__Chlamydiae --out_dir $ANVIO2018/sample_summary_dastool_refine/gtdb_output --cpus 20
gtdbtk ani_rep --genome_dir $ANVIO2018/sample_summary_dastool_refine/bin_by_bin/bins_fa --out_dir $ANVIO2018/sample_summary_dastool_refine/gtdb_output -x fa --cpus 20
gtdbtk classify_wf --genome_dir $ANVIO2018/sample_summary_dastool_refine/bin_by_bin/bins_fa --out_dir $ANVIO2018/sample_summary_dastool_refine/gtdb_output --outgroup_taxon p__Chlamydiae -x fa --cpus 20

####################################################################
### DRAM
####################################################################

conda activate DRAM

cd $ANVIO2016/sample_summary_FINAL/bin_by_bin/bins_fa
DRAM.py annotate -i '*.fa' -o $ANVIO2016/sample_summary_FINAL/DRAM --threads 100

### Distilling (summarizing all data into a very nice graph)

DRAM.py distill -i $ANVIO2016/sample_summary_FINAL/DRAM/annotations.tsv -o $ANVIO2016/sample_summary_FINAL/DRAM/genome_summaries --trna_path $ANVIO2016/sample_summary_FINAL/DRAM/trnas.tsv --rrna_path $ANVIO2016/sample_summary_FINAL/DRAM/rrnas.tsv


conda activate DRAM


cd $ANVIO2017/sample_summary_FINAL/bin_by_bin/bins_fa
DRAM.py annotate -i '*.fa' -o $ANVIO2017/sample_summary_FINAL/DRAM --threads 100

### Distilling (summarizing all data into a very nice graph)

DRAM.py distill -i $ANVIO2017/sample_summary_FINAL/DRAM/annotations.tsv -o $ANVIO2017/sample_summary_FINAL/DRAM/genome_summaries --trna_path $ANVIO2017/sample_summary_FINAL/DRAM/trnas.tsv --rrna_path $ANVIO2017/sample_summary_FINAL/DRAM/rrnas.tsv



conda activate DRAM

cd $ANVIO2018/sample_summary_FINAL/bin_by_bin/bins_fa
DRAM.py annotate -i '*.fa' -o $ANVIO2018/sample_summary_FINAL/DRAM --threads 100

### Distilling (summarizing all data into a very nice graph)

DRAM.py distill -i $ANVIO2018/sample_summary_FINAL/DRAM/annotations.tsv -o $ANVIO2018/sample_summary_FINAL/DRAM/genome_summaries --trna_path $ANVIO2018/sample_summary_FINAL/DRAM/trnas.tsv --rrna_path $ANVIO2018/sample_summary_FINAL/DRAM/rrnas.tsv









#################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################
###################################### code not used, but kept just in case   #############################################################################################################################
#################################################################################################################################################################################################################################
### GENE FUNCTION CALLS
### interproscan on amino-acid seq from anvio
##################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################

## split amino acid file into multiple smaller fasta files (5000 seq per file)

mkdir $ANVIO2016/interproscan
mkdir $ANVIO2016/interproscan/iprs

perl /isilon/ottawa-rdc/users/shared/chenw_lab/galen/split_fasta.pl -i $ANVIO2016/aa_seq_2016.fa -o $ANVIO2016/interproscan/iprs -n 5000

## Create sample list of all file created do manually

ls $ANVIO2016/interproscan/iprs/ > $ANVIO2016/interproscan/iprs/iprs_list


for sample in `awk '{print $1}' $ANVIO2016/interproscan/iprs/iprs_list`;
do
/isilon/ottawa-rdc/users/shared/chenw_lab/galen/interproscan-5.48-83.0/interproscan.sh -i $ANVIO2016/interproscan/iprs/${sample} -f tsv \
	-d $ANVIO2016/interproscan/ \
	--tempdir $WORK/temp/ \
	--disable-precalc \
	-appl Pfam,PIRSF,SUPERFAMILY,TIGRFAM,PANTHER \
	--iprlookup \
	--goterms \
	--pathways
done

### create new folder to house concatenate of all smaller fxnal annotation output (the script dont like a folder with too much clutter, im guessing)

cat $ANVIO2016/interpro/*.tsv > $ANVIO2016/interpro/iprs_output/all_iprs.tsv

### script to clean up and allow import to anvio

## create iprs2anvio.sh file found here: https://github.com/xvazquezc/stuff/blob/master/iprs2anvio.sh

/isilon/ottawa-rdc/users/shared/chenw_lab/galen/interproscan-5.48-83.0/iprs2anvio.sh -i $ANVIO2016/interproscan/all_iprs.tsv -o $ANVIO/interproscan/function2016 -g -p -r


### importing functional annotation to anvio

anvi-import-functions -c $ANVIO2016/contigs2016.db -i $ANVIO2016/interproscan/function2016_iprs2anvio.tsv


####################################################################
### TAXONOMY CALLS
### centrifuge  on amino-acid seq from anvio
####################################################################
mkdir $ANVIO/centrifuge


centrifuge -f -x $CENTRIFUGE_BASE/p+h+v/p+h+v $ANVIO/amino-acid-sequences.fa -S $ANVIO/centrifuge/centrifuge_hits.tsv

## Make sure there is two files in the work directory ($ANVIO/centrifuge/)

anvi-import-taxonomy-for-genes -c $ANVIO/contigs.db -i $ANVIO/centrifuge/centrifuge_report.tsv $ANVIO/centrifuge/centrifuge_hits.tsv -p centrifuge

### adding Single copy gene into database

# run this first, only once.

anvi-setup-scg-databases


anvi-run-scg-taxonomy -c $ANVIO/contigs.db -T 10

anvi-estimate-genome-taxonomy \ 
				-c $ANVIO/contigs.db \
                              	--metagenome-mode \
				--output-file $ANVIO/scg_est_output.txt

# to run after profiling!

anvi-estimate-genome-taxonomy -c $ANVIO/contigs.db \
                              -p $ANVIO/profile_merged/PROFILE.db \
                              --metagenome-mode \
                              --compute-scg-coverages


