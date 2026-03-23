# fuller_md_fortran_en — Fullerene Crystal NPT Molecular Dynamics (Fortran 95)

NPT molecular dynamics simulation code for C60 fullerene crystals (Fortran 95).
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
│   └── fuller_LJ_npt_md_core_serial.f90   [1] LJ rigid-body core (Serial)
└── bin/                   ← Executables (created during build)
```

## Source Files

### Core Version [1] — Learning & Benchmarking

Parameters fixed in source code (T=300K, P=0GPa, dt=1fs, 1000 steps).
Only argument is `nc` (cell size).

| # | File | Description |
|---|------|-------------|
| 1 | `fuller_LJ_npt_md_core_serial.f90` | LJ rigid-body Serial. Complete port from C++ |

## Build

### Prerequisites

- **gfortran** (GCC Fortran compiler)
- macOS: `brew install gcc`

### Build Commands

```bash
cd fuller_md_fortran_en
src/Build_fuller.sh
```

Override compiler with `FC` environment variable:

```bash
FC=/opt/gcc/bin/gfortran src/Build_fuller.sh
```

### Generated Executables

| Executable | Source | Mode |
|-----------|--------|------|
| `fuller_LJ_core_serial` | [1] | Serial |

## Usage

```bash
cd fuller_md_fortran_en

# Default (3x3x3, N=108 molecules, 1000 steps)
bin/fuller_LJ_core_serial

# Specify cell size (5x5x5, N=500 molecules)
bin/fuller_LJ_core_serial 5
bin/fuller_LJ_core_serial --cell=5
```

## Physical Model

- **Ensemble**: NPT (constant temperature & pressure)
- **Thermostat**: Nose-Hoover chain
- **Barostat**: Parrinello-Rahman
- **Time integration**: Velocity-Verlet
- **Periodic boundaries**: 3D (triclinic cell)
- **Rigid-body rotation**: Quaternion representation
- **Intermolecular potential**: Lennard-Jones (C-C, sigma=3.4A, epsilon=2.63meV)
- **Neighbor list**: Symmetric full list (no Newton's 3rd law)

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

- 0-based arrays to 1-based arrays (all index computations adjusted)
- `std::mt19937` replaced with Fortran `RANDOM_NUMBER` + Box-Muller method
- `std::chrono` replaced with `CPU_TIME`
- `std::clamp` replaced with `clamp_val` function
- Flat 1D array layout preserved (for future OpenACC compatibility)
- MODULE structure containing all subroutines

## Platforms

- macOS (Homebrew GCC / gfortran)
- Linux (gfortran)

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
