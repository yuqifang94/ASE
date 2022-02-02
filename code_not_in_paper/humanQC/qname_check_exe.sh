#!/bin/bash
# Input
f="$1"
type="$2"
echo "Reading in bam file $f"
fn=`echo $f|sed 's/\/.*\///g'|sed 's/\.bam//g'`
samtools view -@24 "$f"| cut -f 1-2 > $fn.txt
dat=$fn.txt
echo "starting scan $f"
echo "setting is $type"
if [ "$type" == "SE" ];then
	echo "start qname check with SE"
	cut -f 1 $dat| sort | uniq -c | awk -v limit=1 '{if($1 > limit){print $2}else{print $2>"'$fn'.correct"}}'
	echo "start FLAG check with SE"
	awk '{if ($2 != 0 && $2 !=16 ){print $1}else{print $1>"'$fn'.flag_check"}}' $dat
		
elif [ "$type" == "PE" ];then
	echo "start qname check with PE"
	cut -f 1 $dat| sort | uniq -c | awk -v limit=2 '{if($1 != limit){print $2}else{print $2>"'$fn'.correct"}}'
	echo "start FLAG check with PE"
	awk '{a[$1]+=$2}END{for (i in a) print i,a[i]>"'$fn'.flag"}' $dat
	awk '{if ($2 !=246){print $1}else{print $1>"'$fn'.flag_check"}}' $fn.flag	
	
else
	echo "WRONG type input"

fi
rm $fn.txt
#rm $fn.correct
#rm $fn.flag
echo "Finish scan $f"
