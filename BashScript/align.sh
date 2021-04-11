#!/bin/bash
file1="$1"
file2="$2"
hdmilength="$3"
whitelists="$4"
prefix="$5"
starpath="$6"
seqtkpath="$7"
geneIndex="$8"
file1_basename="$(basename -s .fastq.gz $file1)"
file1_dir="$(dirname "${file1}")"
file1_final_postfix="_modified.fastq"
file1_final="$file1_dir/$file1_basename$file1_final_postfix"
file1_final_gz="$file1_final.gz"

cd $seqtkpath
echo 'sucess'
./seqtk trimfq -q 0 -l $hdmilength $file1 > ./file1_trim.fastq
pigz -p 8 ./file1_trim.fastq
#paste <(zcat ./file1_trim.fastq.gz) <(zcat $file2) | perl -lane 'if ( $. % 4 == 1 ) { print "$F[0] $F[1]"; } elsif ( $. % 4 == 3 ) { print "+"; } else { print substr($F[0],0,$hdmilength).substr($F[1],0,9).substr($F[0],50); }' > $file1_final
if [  "$hdmilength" -eq 30 ]
then 
	paste <(zcat ./file1_trim.fastq.gz) <(zcat $file2) | perl -lane 'if ( $. % 4 == 1 ) { print "$F[0] $F[1]"; } elsif ( $. % 4 == 3 ) { print "+"; } else { print substr($F[0],0,30).substr($F[1],0,9).substr($F[0],50); }' > $file1_final
else 
	paste <(zcat ./file1_trim.fastq.gz) <(zcat $file2) | perl -lane 'if ( $. % 4 == 1 ) { print "$F[0] $F[1]"; } elsif ( $. % 4 == 3 ) { print "+"; } else { print substr($F[0],0,20).substr($F[1],0,9).substr($F[0],50); }' > $file1_final


#paste <(zcat $file1) <(zcat $file2) | perl -lane 'if ( $. % 4 == 1 ) { print "$F[0] $F[1]"; } elsif ( $. % 4 == 3 ) { print "+"; } else { print substr($F[0],0,30).substr($F[1],0,9).substr($F[0],50); }' > $file1_final
#paste <(zcat $file1_trim) <(zcat $file2) | perl -lane 'if ( $. % 4 == 1 ) { print "$F[0] $F[1]"; } elsif ( $. % 4 == 3 ) { print "+"; } else { print substr($F[0],0,30).substr($F[1],0,9).substr($F[0],50); }' > $file1_final_fq

pigz -p 8 $file1_final
rm ./file1_trim.fastq.gz
rndstart=`expr 1 + $hdmilength`
cd $starpath
./STAR    --genomeDir  $geneIndex \
          --readFilesIn  $file2 $file1_final_gz  \
          --outSAMtype BAM SortedByCoordinate  \
          --readFilesCommand zcat \
          --runDirPerm All_RWX \
          --outFileNamePrefix $prefix  \
          --soloType CB_UMI_Simple \
          --soloCBstart 1 --soloCBlen $hdmilength \
          --soloUMIstart $rndstart --soloUMIlen 9 \
          --soloCBwhitelist $whitelists \
          --runThreadN 6 \
          --soloBarcodeReadLength 0 \
          --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
          --outFilterScoreMinOverLread 0 \
          --outFilterMatchNminOverLread 0 \
          --clip3pAdapterSeq AAAAAAAAAA \
          --clip3pAdapterMMp 0.1 \
          --soloFeatures Gene GeneFull SJ Velocyto \
          --limitOutSJcollapsed 1000000 \
          --soloCellFilter None 


