#!/bin/bash
miseq="$1"
hiseq="$2"
hdmilength="$3"
miseq_pos="./spatialcoordinates.txt"
whitelists="./whitelist.txt"


echo " Extract HDMI sequences, position, and STARsolo Whitelists"
zcat $miseq | sed -n '1~4s/:/ /gp' | cut -d ' ' -f 4-7 >  ./pos-MiSeq-temp.txt

zcat $miseq | perl -lane 'print $_ if ( $. % 4 == 2 )'  | cut -c 1-20  > ./HDMIs-MiSeq-temp.txt
cat ./HDMIs-MiSeq-temp.txt | rev | tr ACGTN TGCAN > ./HDMIs-MiSeq-temp-rev.txt
paste  ./HDMIs-MiSeq-temp-rev.txt ./pos-MiSeq-temp.txt | column -s $'\t' -t > ./MiSeq-temp-revHDMIs-pos.txt
#remove duplicated HDMIs
awk '!seen[$1]++' ./MiSeq-temp-revHDMIs-pos.txt  > $miseq_pos
##bottom tiles for MiSeq
if [  "$hdmilength" -eq 30 ]
then
	zcat $hiseq | perl -lane 'print $_ if ( $. % 4 == 2 )'  | cut -c 1-30 > ./HDMI_SeqScope_2nd.txt
	zcat $miseq | perl -lane 'print $_ if ( $. % 4 == 2 )'  | cut -c 3-32  > ./HDMIs-MiSeq-temp.txt
	cat ./HDMIs-MiSeq-temp.txt | rev | tr ACGTN TGCAN > ./HDMIs-MiSeq-temp-rev.txt
	paste  ./HDMIs-MiSeq-temp-rev.txt ./pos-MiSeq-temp.txt | column -s $'\t' -t > ./MiSeq-temp-revHDMIs-pos.txt
	#remove duplicated HDMIs
	awk '!seen[$1]++' ./MiSeq-temp-revHDMIs-pos.txt  > $miseq_pos
	cat $miseq_pos | awk '{ print $1 }' > $whitelists
else 
	zcat $hiseq | perl -lane 'print $_ if ( $. % 4 == 2 )'  | cut -c 1-20 > ./HDMI_SeqScope_2nd.txt
	zcat $miseq | perl -lane 'print $_ if ( $. % 4 == 2 )'  | cut -c 1-20  > ./HDMIs-MiSeq-temp.txt
	cat ./HDMIs-MiSeq-temp.txt | rev | tr ACGTN TGCAN > ./HDMIs-MiSeq-temp-rev.txt
	paste  ./HDMIs-MiSeq-temp-rev.txt ./pos-MiSeq-temp.txt | column -s $'\t' -t > ./MiSeq-temp-revHDMIs-pos.txt
	#remove duplicated HDMIs
	awk '!seen[$1]++' ./MiSeq-temp-revHDMIs-pos.txt  > $miseq_pos
	cat $miseq_pos | awk '{ if ($3 > 2100) { print $1 } }' > $whitelists

fi

rm ./pos-MiSeq-temp.txt
rm ./HDMIs-MiSeq-temp.txt
rm ./MiSeq-temp-revHDMIs-pos.txt
