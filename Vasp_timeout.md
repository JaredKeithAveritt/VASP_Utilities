#RUN a script after timeout

#!/bin/bash
#SBATCH --job-name=GA-opt
#SBATCH -N 2
#SBATCH -p debug_queue
#SBATCH --ntasks-per-node=44
#SBATCH --time=04:00:00
module purge
module load hdf5/1.12.1-openmpi
module list
export PATH=$PATH:/nas/longleaf/home/manojw/VASP/vasp.6.3.2/bin
# Function to handle job resubmission upon receiving SIGTERM
handle_timeout() {
    echo "Job timed out. Run vasp_timeout script..."
    python /nas/longleaf/home/jarkeith/vasp_timeout.py "$(pwd)"
}
# Trap SIGTERM signal and call handle_timeout function
trap 'handle_timeout' SIGTERM mpirun vasp_std
# Main Command
mpirun vasp_std
# If script completes successfully, cancel trap
trap - SIGTERM



# rename CONTCAR to POSCAR and Submit SLURM job
import shutil
import subprocess
import sys


def copy_and_submit(current_directory):
    # Define the source and destination file paths
    src_file = f"{current_directory}/CONTCAR"
    dest_file = f"{current_directory}/POSCAR"

    try:
        # Copy CONTCAR to POSCAR
        shutil.copy(src_file, dest_file)
        print(f"Copied {src_file} to {dest_file}")

        # Submit the SLURM job
        submit_command = ["sbatch", "batch.sh"]
        submission_result = subprocess.run(submit_command, cwd=current_directory, capture_output=True, text=True)

        # Print the result of the submission
        print(f"SLURM submission result: {submission_result.stdout}")
    except Exception as e:
        print(f"An error occurred: {e}")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        current_directory = sys.argv[1]
        copy_and_submit(current_directory)
    else:
        print("Directory path was not passed as an argument.")

