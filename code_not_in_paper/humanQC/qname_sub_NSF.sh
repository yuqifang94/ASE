#!/bin/bash
f=$1
type="$2"
fn=`echo $f|sed 's/\/.*\///g'|sed 's/\.bam//g'`
sbatch <<EOJ
#!/bin/bash
#SBATCH --job-name=${fn}_qn
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G
#SBATCH --time=00:20:00
#SBATCH -p shared,parallel
#SBATCH --chdir /scratch/users/yfang27@jhu.edu/yfang/NSF/qname_check
#SBATCH -o /scratch/users/yfang27@jhu.edu/yfang/NSF/qname_check/log/${fn}.out
#SBATCH -e /scratch/users/yfang27@jhu.edu/yfang/NSF/qname_check/log/${fn}.err
qname_check_exe.sh ${f} ${type}
EOJ
sleep 1
