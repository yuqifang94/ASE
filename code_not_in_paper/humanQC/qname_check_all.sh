#!/bin/bash
echo "HUES64"
for f in /scratch/groups/afeinbe2/shared/JuliASM/data/HUES64/bam/*genome[12]*bam ;do
	echo $f
	qname_check_exe.sh $f PE
done


echo -e "H9"
for f in /scratch/groups/afeinbe2/shared/JuliASM/data/H9/bam/*genome[12]*bam ;do
	echo "$f"
	qname_check_exe.sh $f SE
done

echo -e "skin03"
for f in /scratch/groups/afeinbe2/shared/JuliASM/data/skin03/bam/*genome[12]*bam ;do
	echo "$f"
	qname_check_exe.sh $f PE
done

echo -e "HuFGM02"
for f in /scratch/groups/afeinbe2/shared/JuliASM/data/HuFGM02/bam/*genome[12]*bam ;do
	echo "$f"
	qname_check_exe.sh $f PE
done

echo -e "STL001"
for f in /scratch/groups/afeinbe2/shared/JuliASM/data/STL001/bam/*genome[12]*bam ;do
	echo "$f"
	qname_check_exe.sh $f SE
done

echo -e "STL002"
for f in /scratch/groups/afeinbe2/shared/JuliASM/data/STL002/bam/*genome[12]*bam ;do
	echo "$f"
	qname_check_exe.sh $f SE
done

echo -e "STL003"
for f in /scratch/groups/afeinbe2/shared/JuliASM/data/STL003/bam/*genome[12]*bam ;do
	echo "$f"
	qname_check_exe.sh $f SE
done

echo -e "STL011"
for f in /scratch/groups/afeinbe2/shared/JuliASM/data/STL011/bam/*genome[12]*bam ;do
	echo "$f"
	qname_check_exe.sh $f SE
done
