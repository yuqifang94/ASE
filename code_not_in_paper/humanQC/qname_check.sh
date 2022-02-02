#!/bin/bash
echo "HUES64"
for f in /scratch/groups/afeinbe2/shared/CpelAsm/data/HUES64/bam/*genome[12]*bam ;do
	echo $f
	qname_sub.sh $f PE
done


echo -e "H9"
for f in /scratch/groups/afeinbe2/shared/CpelAsm/data/H9/bam/*genome[12]*bam ;do
	echo "$f"
	qname_sub.sh $f SE
done

echo -e "skin03"
for f in /scratch/groups/afeinbe2/shared/CpelAsm/data/skin03/bam/*genome[12]*bam ;do
	echo "$f"
	qname_sub.sh $f PE
done

echo -e "HuFGM02"
for f in /scratch/groups/afeinbe2/shared/CpelAsm/data/HuFGM02/bam/*genome[12]*bam ;do
	echo "$f"
	qname_sub.sh $f PE
done

echo -e "STL001"
for f in /scratch/groups/afeinbe2/shared/CpelAsm/data/STL001/bam/*genome[12]*bam ;do
	echo "$f"
	qname_sub.sh $f SE
done

echo -e "STL002"
for f in /scratch/groups/afeinbe2/shared/CpelAsm/data/STL002/bam/*genome[12]*bam ;do
	echo "$f"
	qname_sub.sh $f SE
done

echo -e "STL003"
for f in /scratch/groups/afeinbe2/shared/CpelAsm/data/STL003/bam/*genome[12]*bam ;do
	echo "$f"
	qname_sub.sh $f SE
done

echo -e "STL011"
for f in /scratch/groups/afeinbe2/shared/CpelAsm/data/STL011/bam/*genome[12]*bam ;do
	echo "$f"
	qname_sub.sh $f SE
done
echo -e "H1"
for f in /scratch/groups/afeinbe2/shared/CpelAsm/data/H1/bam/*genome[12]*bam ;do
	echo "$f"
	qname_sub_NSF.sh $f PE
done
echo -e "112"
for f in /scratch/groups/afeinbe2/shared/CpelAsm/data/112/bam/*genome[12]*bam ;do
	echo "$f"
	qname_sub.sh $f PE
done
echo -e "150"
for f in /scratch/groups/afeinbe2/shared/CpelAsm/data/150/bam/*genome[12]*bam ;do
	echo "$f"
	qname_sub.sh $f PE
done
echo -e "149"
for f in /scratch/groups/afeinbe2/shared/CpelAsm/data/149/bam/*genome[12]*bam ;do
	echo "$f"
	qname_sub.sh $f PE
done
echo -e "HuFGM02"
for f in /scratch/groups/afeinbe2/shared/CpelAsm/data/HuFGM02/bam/*genome[12]*bam ;do
	echo "$f"
	qname_sub.sh $f PE
done



