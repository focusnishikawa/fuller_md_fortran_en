# fuller_md_fortran_en — Fullerene Crystal NPT Molecular Dynamics (Fortran 95)

NPT molecular dynamics simulation codes for C60/C70/C72/C74/C76/C84 fullerene crystals (Fortran 95).
Three force field models (LJ rigid-body, Molecular Mechanics, AIREBO) with
Serial / OpenMP parallelization modes selectable at compile time.
Ported from the [C++ version fuller_md](https://github.com/focusnishikawa/fuller_md).

## Directory Structure

```
fuller_md_fortran_en/
├── README.md              ← This file
├── FullereneLib/          ← Fullerene molecular coordinate data (.cc1)
│   ├── C60-76/            ←   C60(Ih), C70(D5h), C72(D6d), C74(D3h), C76(D2,Td)
│   └── C84/               ←   C84 No.01–No.24 (24 isomers)
├── src/                   ← Source code & scripts
│   ├── Build_fuller.sh    ←   Build script
│   ├── fuller_LJ_npt_md_core_serial.f90          [1] LJ rigid-body core (Serial only)
│   ├── fuller_LJ_npt_md_core_serial_omp_acc.f90   [2] LJ rigid-body core (Serial/OMP)
│   ├── fuller_LJ_npt_md_serial_omp_acc.f90        [3] LJ rigid-body full (Serial/OMP)
│   ├── fuller_LJ_npt_mmmd_serial_omp_acc.f90      [4] Molecular mechanics full (Serial/OMP)
│   └── fuller_airebo_npt_md_serial_omp_acc.f90    [5] AIREBO full (Serial/OMP)
└── bin/                   ← Executables (created during build)
```

## Source Files

### Core Versions [1][2] — Learning & Benchmarking

Parameters fixed in source code (T=300K, P=0GPa, dt=1fs, 1000 steps).
Only argument is `nc` (cell size). No restart or OVITO output.

| # | File | Description | Parallelization |
|---|------|-------------|-----------------|
| 1 | `fuller_LJ_npt_md_core_serial.f90` | LJ rigid-body Serial only. With parallelization guide comments | Serial |
| 2 | `fuller_LJ_npt_md_core_serial_omp_acc.f90` | LJ rigid-body 2-mode unified. `!$OMP` directive switching | Serial/OMP |

### Full Versions [3][4][5] — Production Runs

All runtime options supported. Restart save/resume, OVITO XYZ output, stop control (abort.md/stop.md).

| # | File | Force Field | Default dt |
|---|------|-------------|-----------|
| 3 | `fuller_LJ_npt_md_serial_omp_acc.f90` | LJ rigid-body intermolecular | 1.0 fs |
| 4 | `fuller_LJ_npt_mmmd_serial_omp_acc.f90` | Molecular mechanics (Bond+Angle+Dihedral+Improper+LJ) | 0.1 fs |
| 5 | `fuller_airebo_npt_md_serial_omp_acc.f90` | AIREBO (REBO-II + LJ) | 0.5 fs |

## Build

### Prerequisites

- **Serial/OpenMP**: gfortran (GCC Fortran compiler)
- macOS: `brew install gcc`

### Build Script

```bash
cd fuller_md_fortran_en

# Build Serial + OpenMP (default)
src/Build_fuller.sh

# Serial only
src/Build_fuller.sh serial

# OpenMP only
src/Build_fuller.sh omp

# All modes
src/Build_fuller.sh all

# Clean build artifacts
src/Build_fuller.sh clean
```

Override compiler with `FC` environment variable:

```bash
FC=/opt/gcc/bin/gfortran src/Build_fuller.sh
```

### Generated Executables

| Executable | Source | Mode |
|-----------|--------|------|
| `fuller_LJ_core_serial_pure` | [1] | Serial |
| `fuller_LJ_core_serial` | [2] | Serial |
| `fuller_LJ_core_omp` | [2] | OpenMP |
| `fuller_LJ_npt_md_serial` | [3] | Serial |
| `fuller_LJ_npt_md_omp` | [3] | OpenMP |
| `fuller_LJ_npt_mmmd_serial` | [4] | Serial |
| `fuller_LJ_npt_mmmd_omp` | [4] | OpenMP |
| `fuller_airebo_npt_md_serial` | [5] | Serial |
| `fuller_airebo_npt_md_omp` | [5] | OpenMP |

## Usage

### Core Version

```bash
cd fuller_md_fortran_en

# Default (3x3x3, N=108 molecules, 1000 steps)
bin/fuller_LJ_core_serial_pure

# Specify cell size (5x5x5, N=500 molecules)
bin/fuller_LJ_core_omp 5
bin/fuller_LJ_core_omp --cell=5
```

### Full Version — Basic

```bash
# LJ rigid-body (default: C60, FCC 3x3x3, 298K, 0GPa, 10000 steps)
bin/fuller_LJ_npt_md_serial

# Specify temperature, pressure, steps
bin/fuller_LJ_npt_md_omp --temp=500 --pres=1.0 --step=50000

# Cold start (4K) + warmup + production
bin/fuller_LJ_npt_md_serial --coldstart=2000 --warmup=3000 --step=20000
```

### Full Version — OVITO Output

```bash
bin/fuller_LJ_npt_md_omp --step=10000 --ovito=100
bin/fuller_LJ_npt_mmmd_omp --step=20000 --ovito=200
bin/fuller_airebo_npt_md_omp --step=10000 --ovito=100
```

### Full Version — Restart

```bash
bin/fuller_LJ_npt_md_serial --step=50000 --restart=5000
bin/fuller_LJ_npt_md_serial --resfile=restart_LJ_serial_00025000.rst
```

### Full Version — Stop Control

```bash
mkdir abort.md    # Immediate stop (saves restart if enabled)
mkdir stop.md     # Stop at next restart checkpoint
```

### Full Version — Help

```bash
bin/fuller_LJ_npt_md_serial --help
bin/fuller_LJ_npt_mmmd_serial --help
bin/fuller_airebo_npt_md_serial --help
```

## Runtime Options (Full Versions)

| Option | Description | Default |
|--------|-------------|---------|
| `--help` | Show help | — |
| `--fullerene=<name>` | Fullerene species | C60 |
| `--crystal=<fcc\|hcp\|bcc>` | Crystal structure | fcc |
| `--cell=<nc>` | Unit cell repeats | 3 |
| `--temp=<K>` | Target temperature [K] | 298.0 |
| `--pres=<GPa>` | Target pressure [GPa] | 0.0 |
| `--step=<N>` | Production steps | 10000 |
| `--dt=<fs>` | Time step [fs] | Force field dependent |
| `--init_scale=<s>` | Lattice scale factor | 1.0 |
| `--seed=<n>` | Random seed | 42 |
| `--coldstart=<N>` | Cold-start steps at 4K | 0 |
| `--warmup=<N>` | Warmup ramp steps 4K→T | 0 |
| `--from=<step>` | Averaging start step | Auto (3/4 point) |
| `--to=<step>` | Averaging end step | nsteps |
| `--mon=<N>` | Monitor print interval | Auto |
| `--ovito=<N>` | OVITO XYZ output interval (0=off) | 0 |
| `--restart=<N>` | Restart save interval (0=off) | 0 |
| `--resfile=<path>` | Resume from restart file | — |
| `--libdir=<path>` | Fullerene library path | FullereneLib |

### Molecular Mechanics [4] Additional Options

| Option | Description | Default |
|--------|-------------|---------|
| `--ff_kb=<kcal/mol>` | Bond stretch force constant | 469.0 |
| `--ff_kth=<kcal/mol>` | Angle bend force constant | 63.0 |
| `--ff_v2=<kcal/mol>` | Dihedral force constant | 14.5 |
| `--ff_kimp=<kcal/mol>` | Improper dihedral force constant | 15.0 |

## Physical Model

- **Ensemble**: NPT (constant temperature & pressure)
- **Thermostat**: Nose-Hoover chain
- **Barostat**: Parrinello-Rahman
- **Time integration**: Velocity-Verlet
- **Periodic boundaries**: 3D (triclinic cell)
- **Neighbor list**: Symmetric full list (no Newton's 3rd law)

### LJ Rigid-body [1][2][3]

- Intermolecular: Lennard-Jones (C-C, sigma=3.4A, epsilon=2.63meV)
- Rigid-body rotation: Quaternion representation
- DOF: COM translation(3) + rotation(3) per molecule

### Molecular Mechanics [4]

- Intramolecular: Bond stretching + Angle bending + Dihedral + Improper
- Intermolecular: LJ + Coulomb (if charges present)
- All-atom DOF

### AIREBO [5]

- REBO-II (Brenner 2002): Covalent bonds (bond order potential)
- LJ (Stuart 2000): Intermolecular van der Waals
- All-atom DOF

## Units

| Quantity | Unit |
|----------|------|
| Distance | A (Angstrom) |
| Mass | amu (atomic mass unit) |
| Energy | eV (electron volt) |
| Time | fs (femtosecond) |
| Temperature | K (Kelvin) |
| Pressure | GPa (Gigapascal) |

## Porting Notes (C++ to Fortran 95)

- 0-based arrays → 1-based arrays (all index computations adjusted)
- `std::mt19937` → Fortran `RANDOM_NUMBER` + Box-Muller method
- `std::chrono` → `CPU_TIME`
- `std::clamp` → `clamp_val` function
- `#ifdef _OPENMP` → `!$OMP` sentinel (ignored by non-OpenMP compilers)
- Flat 1D array layout preserved
- MODULE structure containing all subroutines

## Platforms

- macOS (Homebrew GCC / gfortran) — Serial/OpenMP
- Linux (gfortran) — Serial/OpenMP

## License

This project is released under the [BSD 3-Clause License](LICENSE).

## Other Languages & Versions

| Repository | Language | Description |
|-----------|----------|-------------|
| [fuller_md](https://github.com/focusnishikawa/fuller_md) | C++ (Japanese) | Original |
| [fuller_md_en](https://github.com/focusnishikawa/fuller_md_en) | C++ (English) | C++ English |
| [fuller_md_Julia](https://github.com/focusnishikawa/fuller_md_Julia) | Julia (English) | Julia English |
| [fuller_md_Julia_ja](https://github.com/focusnishikawa/fuller_md_Julia_ja) | Julia (Japanese) | Julia Japanese |
| [fuller_md_fortran](https://github.com/focusnishikawa/fuller_md_fortran) | Fortran 95 (Japanese) | Fortran Japanese |
| [fuller_md_fortran_en](https://github.com/focusnishikawa/fuller_md_fortran_en) | Fortran 95 (English) | This repository |
