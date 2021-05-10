for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
  name=`echo py_chr$i`
  jobfile=${name}.job.sh

  cat <<EOF > ${jobfile}
#!/bin/bash
#SBATCH --job-name pytest.${name}
#SBATCH --mem 50g
#SBATCH --partition allnodes
#SBATCH --output ${name}_pytest.job.out
#SBATCH --error ${name}_.pytest.job.err


python SEOVLP_CNV.py --chrom chr$i --cnumber $i

EOF

  chmod +x ${jobfile}
  sbatch --export=ALL ${jobfile}

done

