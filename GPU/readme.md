Check out here: 

https://www.vasp.at/forum/viewtopic.php?p=22828&hilit=OMP_NUM_THREADS#p22828



Curently trying with the bashfile

#!/bin/bash
#SBATCH --job-name=mos2_scf
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -p beta-gpu
#SBATCH --mem=200g
#SBATCH --qos gpu_access
#SBATCH --gres=gpu:2
#SBATCH -t 1-

module add vasp/6.3.2

export OMP_NUM_THREADS=2
export MKL_THREADING_LAYER=INTEL

mpirun -np 2 vasp_std



https://enccs.github.io/vasp-best-practices/compile/

they try 
---
https://help.rc.ufl.edu/doc/Slurm_and_GPU_Use

#!/bin/bash
#SBATCH --job-name=vasptest
#SBATCH --output=vasp.out
#SBATCH --error=vasp.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=email@ufl.edu
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=8
#SBATCH --ntasks-per-socket=4
#SBATCH --mem-per-cpu=7000mb
#SBATCH --distribution=cyclic:cyclic
#SBATCH --partition=gpu
#SBATCH --gres=gpu:geforce:4
#SBATCH --time=00:30:00

echo "Date      = $(date)"
echo "host      = $(hostname -s)"
echo "Directory = $(pwd)"

module purge
module load cuda/10.0.130  intel/2018  openmpi/4.0.0 vasp/5.4.4

T1=$(date +%s)
srun --mpi=pmix_v3 vasp_gpu
T2=$(date +%s)

ELAPSED=$((T2 - T1))
echo "Elapsed Time = $ELAPSED"
