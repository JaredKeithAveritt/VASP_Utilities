#this code takes inputs of 
#num_atoms = the number of atoms in your system 
#cores_per_node =the number of cores per node for the partition you want to use

import math

num_atoms = 177
cores_per_node = 40


def calc_ncore(num_atoms,cores_per_node):
    num_nodes = math.ceil(num_atoms / cores_per_node)
    total_num_cores = num_nodes * cores_per_node
    sqrt_total_num_cores = int(math.sqrt(total_num_cores))
    # Check each integer starting from 2 up to sqrt_total_num_cores
    for test_ncore in range(2, sqrt_total_num_cores+1):
        if cores_per_node % test_ncore == 0 and test_ncore % 2 == 0:
            ncore=test_ncore
    return ncore, num_nodes,total_num_cores, cores_per_node

ncore, num_nodes , total_num_cores, cores_per_node=calc_ncore(num_atoms,cores_per_node)

print('add this line to the INCAR file: ')
print('NCORE = ', ncore)
print('and for the .sh file set: ')
print('#SBATCH -N ', num_nodes)
print('#SBATCH --ntasks-per-node= ', cores_per_node)
