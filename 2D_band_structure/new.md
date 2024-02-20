Below are examples of what your `INCAR` files might look like for each step in a typical VASP band structure calculation process. Keep in mind that specific values (like `ENCUT`, `EDIFF`, etc.) should be tailored to your particular system and the level of precision required for your study.

### 1. Structural Optimization (Relaxation):

This `INCAR` file is for the initial step where you relax the atomic positions and the lattice parameters to minimize the energy of the system.

```plaintext
SYSTEM = Your_system_name  ! Name of your system
ISTART = 0                 ! Start from scratch
ENCUT = 520                ! Plane-wave cutoff energy in eV, adjust as needed
EDIFF = 1E-5               ! Energy convergence criterion
IBRION = 2                 ! Use the conjugate gradient algorithm
ISIF = 3                   ! Relax internal coordinates and cell volume
NSW = 100                  ! Number of steps for ionic relaxation
NWRITE = 2                 ! Reduce output verbosity
```

### 2. Self-Consistent Field (SCF) Calculation:

After relaxation, you perform an SCF calculation to get a converged charge density. The `INCAR` for this might look like:

```plaintext
SYSTEM = Your_system_name
ISTART = 0                 ! Start from scratch or use 1 to read the previous wavefunction
ENCUT = 520                ! Keep consistent with the relaxation step
EDIFF = 1E-6               ! Tighter energy convergence criterion
ISMEAR = 0                 ! Gaussian smearing, suitable for metals; for insulators, might use ISMEAR = -5
SIGMA = 0.05               ! Width of smearing in eV
LCHARG = .TRUE.            ! Create CHGCAR
LWAVE = .TRUE.             ! Create WAVECAR, needed for non-SCF step
```

### 3. Non-SCF Calculation for Band Structure:

For the band structure calculation (where you don't update the charge density), your `INCAR` would change like this:

```plaintext
SYSTEM = Your_system_name
ISTART = 1                 ! Start from existing wavefunctions
ICHARG = 11                ! Use existing charge density without updating
ENCUT = 520                ! Consistent with previous steps
EDIFF = 1E-6               ! Same convergence criteria as SCF step
ISMEAR = 0                 ! Same as SCF; adjust if your system is an insulator
SIGMA = 0.05               ! Smearing width, same as SCF
LCHARG = .FALSE.           ! Do not write CHGCAR to save disk space
LWAVE = .FALSE.            ! Do not write WAVECAR to save disk space (optional, might be needed for further analysis)
NBANDS = X                 ! Set to a suitable value, often higher than in SCF step
```

For the band structure calculation, ensure you're using the `KPOINTS` file that specifies the path through the Brillouin zone connecting high-symmetry points, rather than a grid used in the SCF calculation.

Remember:

- The `SYSTEM` tag is just a label and can be any string that helps you identify the calculation.
- `NBANDS` (number of bands) should be set to a suitable value, especially if you're interested in states above the Fermi level; sometimes it's increased from the default value to ensure all relevant bands are included.
- `ENCUT`, `EDIFF`, `ISMEAR`, and `SIGMA` should be consistent with your previous calculations for consistency, but you can adjust these based on what is typical for your system and the level of precision required.
- The `INCAR` settings can vary based on the specifics of your system (e.g., metal vs. insulator) and your computational goals. Always review and adjust parameters based on your system's needs and the best practices in the literature.
