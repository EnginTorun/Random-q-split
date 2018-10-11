#!/bin/bash
#
cat>ph.sh<<EOF
#!/bin/bash -l
# Submission script for Iris
#SBATCH -J q${1}_$2
#SBATCH -n 2
#SBATCH --ntasks-per-node=2
#SBATCH --ntasks-per-socket=2
#SBATCH --time=0-08:00:00
#SBATCH -p batch
#SBATCH --qos=qos-batch
module load swenv/default-env/devel
module load phys/Yambo/4.2.4-intel-2018a-QuantumESPRESSO-6.1 # QE 6.1 is loaded with yambo by this module
srun -n \$SLURM_NTASKS ph.x -inp $3
EOF
cp ph.sh q$1 ; cd q$1 ; sbatch ph.sh ; cd ..
