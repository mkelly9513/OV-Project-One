for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
  name=`echo py_chr$i`
  jobfile=${name}.job.sh

  cat <<EOF > ${jobfile}
#!/bin/bash
#SBATCH --job-name pytestg.${name}
#SBATCH --mem 50g
#SBATCH --time 14-0
#SBATCH --partition allnodes
#SBATCH --output ${name}_pytestg.job.out
#SBATCH --error ${name}_.pytestg.job.err


python OVLP_CNV_Whole_Genome.py --chrom chr$i --cnumber $i

EOF

  chmod +x ${jobfile}
  sbatch --export=ALL ${jobfile}

done

