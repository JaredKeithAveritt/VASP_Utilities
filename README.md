https://www.materialscloud.org/learn/sections/VNL7RL/a-gentle-introduction-to-dft-calculations-april-2020

# Band structure calculations using VASP and vaspkit

## Step 1. Optimize Structures using VASP

create a folder named `opt`

create an `INCAR` file
```fortran
### Relaxation 
 ENCUT = 500
 SIGMA = 0.010000
 PREC = Accurate
 IBRION = 2
 ISMEAR = 0
 NSW = 1000
 IVDW = 12                             
 LREAL = Auto
 EDIFF = 1.0e-7
 EDIFFG = -0.01
```
create a `KPOINTS` file (example)
```fortran
K-Spacing Value to Generate K-Mesh: 0.030
0
Monkhorst-Pack
  15  15   1
0.0  0.0  0.0
```

generate POTCAR and POSCAR

note: increase kpoints untill energy convergence
note: N N 1 for graphene
note: N M ? for mos2 ? 

## Step 2. Convert POSCAR to standardized primitive cell

In the `opt` folder containing the `CONTCAR` 

```bash
cp POSCAR POSCAR.opt
rm POSCAR 
cp CONTCAR POSCAR
```
open vaspkit and select option `602`
```bash
rm POSCAR
```
you should have no `POSCAR` file, but instead two files: `POSCAR.opt` and `PRIMCELL.vasp`

an additional file that is generated `SYMMETRY` 


## Step 3. Calculate CHGCAR file by performing a single-point self-consistent calculation

Create a new folder called `scf` and copy the `PRIMCELL.vasp` file into that folder as `POSCAR` using:

```bash
cp ../opt/PRIMCELL.vasp POSCAR
```

copy `POTCAR` and `batch.sh` files from `opt` folder to `scf` folder

```bash
cp ../opt/POTCAR .
cp ../opt/batch.sh .
```

Create `INCAR`
```fortran
Global Parameters
NCORE  =  4
ISTART =  1            #(Read existing wavefunction, if there)
ISPIN  =  1            #(Non-Spin polarised DFT)
LREAL  = .FALSE.       #(Projection operators: automatic)
ENCUT  =  500        #(Cut-off energy for plane wave basis set, in eV)
PREC   =  Accurate   #(Precision level: Normal or Accurate, set Accurate when perform structure lattice relaxation calculation)
LWAVE  = .TRUE.        #(Write WAVECAR or not)
LCHARG = .TRUE.        #(Write CHGCAR or not)
ADDGRID= .TRUE.        #(Increase grid, helps GGA convergence)
IVDW = 12 
# LVTOT  = .TRUE.      #(Write total electrostatic potential into LOCPOT or not)
# LVHAR  = .TRUE.      #(Write ionic + Hartree electrostatic potential into LOCPOT or not)
# NELECT =             #(No. of electrons: charged cells, be careful)
# LPLANE = .TRUE.      #(Real space distribution, supercells)
# NWRITE = 2           #(Medium-level output)
# KPAR   = 2           #(Divides k-grid into separate groups)
# NGXF    = 300        #(FFT grid mesh density for nice charge/potential plots)
# NGYF    = 300        #(FFT grid mesh density for nice charge/potential plots)
# NGZF    = 300        #(FFT grid mesh density for nice charge/potential plots)
 
Static Calculation
ISMEAR =  0            #(gaussian smearing method)
SIGMA  =  0.05         #(please check the width of the smearing)
LORBIT =  11           #(PAW radii for projected DOS)
NEDOS  =  2001         #(DOSCAR points)
NELM   =  60           #(Max electronic SCF steps)
EDIFF  =  1E-08        #(SCF energy convergence, in eV)
```

Create `KPATH.in` file  using VASPKIT option `302` and copy/save it as `KPOINTS` file, 

```bash
cp KPATH.in KPOINTS
```
it will look like this (this can be double checked and visualized in [Link: seeK-path](https://www.materialscloud.org/work/tools/seekpath):

```plaintext
K-Path Generated by VASPKIT.
   20
Line-Mode
Reciprocal
   0.0000000000   0.0000000000   0.0000000000     GAMMA          
   0.5000000000   0.0000000000   0.0000000000     X              
 
   0.5000000000   0.0000000000   0.0000000000     X              
   0.5000000000   0.5000000000   0.0000000000     S              
 
   0.5000000000   0.5000000000   0.0000000000     S              
   0.0000000000   0.5000000000   0.0000000000     Y              
 
   0.0000000000   0.5000000000   0.0000000000     Y              
   0.0000000000   0.0000000000   0.0000000000     GAMMA    
```

## Step 4. Calculate DOS

NOTE: make sure that `PRECISION` and `ENCUT` are the same in the `INCAR` file
Copy `CHGCAR` file from previous step to new folder called `dos`
Copy `WAVECAR` file from `scf` folder to `dos` folder

```bash
cp ../scf/POTCAR .
cp ../scf/POSCAR .
cp ../scf/batch.sh .
cp ../scf/CHGCAR .
cp ../scf/WAVECAR .
```


Create `INCAR`
```fortran
 ISTART = 2
 PREC = Accurate
 LREAL = Auto
 ENCUT = 500
 EDIFF = 1.0e-7
 IVDW = 12    
 ISMEAR = 0
 SIGMA = 0.010000
 ICHARG = 11
```
## Step 4. Generate DOS plots
note: must not be done using terminus, use the default terminal and have -X11 forwarding turned on. 

in `dos` folder, run vaspkit `211` to get the band structure. a `band.png` figure will result. Additional files including `` , `` , ``  are also generated.

`BAND.dat` and `BAND_RECORMATTED.dat` can be used in ORIGIN.
in `BAND_RECORMATTED.dat` , the first column is the length of the K-path in units of 1/Angstrom, the following columns are the energy of each bands.

in `KLABELS`, contains the positions of high symmetry points on band structure figures
