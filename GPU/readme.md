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
