import math

num_atoms = #set the number of atoms here
cores_per_node = # set the number of cores per node for the partition you want to use for the calculation


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
print('and for the .sh file set: )
print('#SBATCH -N ', ncore)
print('#SBATCH --ntasks-per-node= ', ncore)
