Certainly! Here is the provided content in markdown format:

```markdown
## Step 1. Optimize Structures using INCAR

```fortran
Global Parameters (SR) Standard Relaxation, geometry optimization
ISTART =  1            (Read existing wavefunction, if there)
ISPIN  =  1            (Non-Spin polarised DFT)
# ICHARG =  11         (Non-self-consistent: GGA/LDA band structures)
LREAL  = .FALSE.       (Projection operators: automatic)
ENCUT  =  400        (Cut-off energy for plane wave basis set, in eV)
# PREC   =  Accurate   (Precision level: Normal or Accurate, set Accurate when perform structure lattice relaxation calculation)
LWAVE  = .FALSE.        (Write WAVECAR or not)
LCHARG = .FALSE.        (Write CHGCAR or not)
ADDGRID= .TRUE.        (Increase grid, helps GGA convergence)
# LVTOT  = .TRUE.      (Write total electrostatic potential into LOCPOT or not)
# LVHAR  = .TRUE.      (Write ionic + Hartree electrostatic potential into LOCPOT or not)
# NELECT =             (No. of electrons: charged cells, be careful)
# LPLANE = .TRUE.      (Real space distribution, supercells)
# NWRITE = 2           (Medium-level output)
# KPAR   = 2           (Divides k-grid into separate groups)
# NGXF    = 300        (FFT grid mesh density for nice charge/potential plots)
# NGYF    = 300        (FFT grid mesh density for nice charge/potential plots)
# NGZF    = 300        (FFT grid mesh density for nice charge/potential plots)
```

## Step 2. Convert POSCAR to standardized primitive cell

In the folder containing the POSCAR open vaspkit and select option 602
Create a new folder and copy the PRIMCELL.vasp file into that folder as POSCAR using:
```bash
cp ../PRIMCELL.vasp POSCAR
```

## Step 3. Calculate CHGCAR file by performing a single-point self-consistent calculation

Starting with POSCAR files above and POTCAR and batch.sh files
Create INCAR

```fortran
Global Parameters, (ST) Static-Calculation , single-point self-consistent calculation
ISTART =  1            (Read existing wavefunction, if there)
ISPIN  =  1            (Non-Spin polarised DFT)
# ICHARG =  11         (Non-self-consistent: GGA/LDA band structures)
LREAL  = .FALSE.       (Projection operators: automatic)
ENCUT  =  400        (Cut-off energy for plane wave basis set, in eV)
PREC   =  Accurate   (Precision level: Normal or Accurate, set Accurate when perform structure lattice relaxation calculation)
LWAVE  = .TRUE.        (Write WAVECAR or not)
LCHARG = .TRUE.        (Write CHGCAR or not)
ADDGRID= .TRUE.        (Increase grid, helps GGA convergence)
# LVTOT  = .TRUE.      (Write total electrostatic potential into LOCPOT or not)
# LVHAR  = .TRUE.      (Write ionic + Hartree electrostatic potential into LOCPOT or not)
# NELECT =             (No. of electrons: charged cells, be careful)
# LPLANE = .TRUE.      (Real space distribution, supercells)
# NWRITE = 2           (Medium-level output)
# KPAR   = 2           (Divides k-grid into separate groups)
# NGXF    = 300        (FFT grid mesh density for nice charge/potential plots)
# NGYF    = 300        (FFT grid mesh density for nice charge/potential plots)
# NGZF    = 300        (FFT grid mesh density for nice charge/potential plots)

Static Calculation
ISMEAR =  0            (Gaussian smearing, metals:1)
SIGMA  =  0.05         (Smearing value in eV, metals:0.2)
LORBIT =  11           (PAW radii for projected DOS)
NEDOS  =  2001         (DOSCAR points)
NELM   =  60           (Max electronic SCF steps)
EDIFF  =  1E-08        (SCF energy convergence, in eV)
```

C. Create KPOINTS File, for graphene/MoS2 that looks like

```plaintext
Comment Line
0
Monkhorst-Pack
5	5	5
0.0	0.0	0.0
```

## Step 5. Calculate DOS

A. Move CHGCAR file from previous step to new folder

B. Create INCAR

```fortran
ICHARG = 11 #read in chgcar file
LWAVE = .FALSE. # do not write WAVECAR
LCHARG = .FALSE.
LVTOT = .FALSE.
LVHAR = .FALSE.
LELF = .FALSE.
LORBIT = 11
NEDOS = 1000
NCORE = 10
##### SCF #####
ENCUT = 400
ISMEAR = 0
SIGMA = 0.05
EDIFF = 1E-6
NELMIN = 5
NELM = 300
GGA = PE
LREAL = .FALSE.
PREC = Accurate
```

C. Create KPATH.in file  using VASPKIT and copy/save it as KPOINTS file, it will look like this:

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
```

If you have any more requests or if there's anything else I can help you with, feel free to let me know!