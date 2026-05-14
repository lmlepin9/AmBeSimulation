# AmBe Source Simulation

Geant4/ROOT simulation for AmBe-style neutron sources with optional
multithreading, MPI execution, source-casing selection, and ROOT output for
spectrum and emerging-particle studies.

## Dependencies
- Geant4
- ROOT
- OpenMPI or other MPI interface
- [bhamnuclear/Geant4-MPI](https://github.com/bhamnuclear/Geant4-MPI/)

MPI can be disabled in CMake and main. In any case, it can be chosen at runtime whether to run with or without MPI via the `AMBE_MPI` environment flag: set this either `TRUE` or `FALSE`.

The executable also accepts `PHYSLIST=<name>` to use a Geant4 reference physics
list instead of the default local `PhysicsList`.

## Build

From the repository root:

```bash
mkdir -p build
cd build
cmake -DG4mpi_DIR=$PWD/../Geant4-MPI/install/lib/G4mpi-11.3.1 ..
make
```

Replace `11.3.1` with the installed Geant4/G4mpi version.

The CMake build copies the `macros/` and `AmBeData/` directories into `build/`,
so either source-tree or build-tree macro paths can be used depending on where
the executable is run.

## Docker

A Docker workflow is available for reproducible builds with ROOT `6.30.06`,
Geant4 `11.3.1`, OpenMPI, and the pinned `Geant4-MPI` revision.

```bash
docker build -t ambe-sim .
docker run --rm -it \
  --user "$(id -u):$(id -g)" \
  -v "$PWD:/workspace/AmBeSimulation" \
  ambe-sim
```

See `DOCKER.md` for custom non-MPI and MPI examples.

## MPI USAGE
Use with [https://github.com/bhamnuclear/Geant4-MPI/](https://github.com/bhamnuclear/Geant4-MPI/), here available in the `Geant4-MPI` folder.

Compile Geant4 MPI first using `Compile.sh` or otherwise, then compile this simulation:

```bash
mkdir build
cd build
cmake -DG4mpi_DIR=$PWD/../Geant4-MPI/install/lib/G4mpi-11.3.1 ..
make
```

Run with:

```bash
AMBE_MPI=true mpirun --mca accelerator rocm -n 2 ambe -t 8 macros/shortlong.mac
```

or:

```bash
AMBE_MPI=true mpirun -n 2 ambe -t 8 macros/shortlong.mac
```


## Current Status
1. Scoring has been implemented correctly
2. Full Multithreading and MPI support
3. `241AmBe`, `239PuBe`, uranium spontaneous fission (`USF`), and monoenergetic neutron modes are available from the CLI. Other sources have only been partially implemented.
4. Source casing options include the cylindrical approximation, X3 casing, and an N02 stainless-steel capsule approximation.
5. World/material modes include air/default placement, no water, water world, and D2O world.
6. Optional output can store initial neutrons, secondary gamma scoring, azimuthal surface scoring, fission-fragment scoring, neutron tracking diagnostics, and detailed emerging-particle records.


## Usage
Various options have been implemented using CLI arguments. It is also possible to use macros, but there is no full support for these and they have not been fully tested. Some macro commands can replace CLI commands, but some macro commands cannot replace CLI commands.

Note that CLI commands are more stable and do not cause crashing when run once the simulation is initialised.

### Cli commands
Options list, as implemented in `main.cc`:

- `-h`, `--help` -> print the built-in usage message and exit.
- `-i`, `--interactive` -> enable the interactive Geant4 UI. This is enabled by default for non-MPI runs with no arguments.
- `-t`, `--threads <N>` -> specify the number of Geant4 worker threads. If using MPI, this is the number of threads allocated to each rank. Default: `1`.
- `-nw`, `--nowater <MODE>` -> select the world/water configuration: `0` air/default placement, `1` no water, `2` water world, `3` D2O world.
- `-r`, `--isotope <MODE>` -> select the neutron source: `-1` default `241Am`, `0` `241Am`, `1` `239Pu`, `2` `USF`. Monoenergetic neutrons can be requested with `SingleX.YZ`, where `X.YZ` is the energy in MeV, for example `Single2.45`.
- `-cs`, `--casing <MODE>` -> select the casing geometry: `0` cylindrical approximation, `1` X3 casing, `2` N02 stainless-steel capsule approximation.
- `-aom`, `--amO2MassMicrogram <UG>` -> set a fixed AmO2/PuO2 oxide mass in micrograms while preserving the default Be/oxide mass ratio. Values `<= 0` use the default active-material density and composition.
- `-in`, `--initialNeutrons` -> store initial neutrons from the primary generator.
- `-ffs`, `--fissionFragmentsScore` -> enable scoring of fission fragments.
- `-nt`, `--neutronTracking` -> enable neutron tracking diagnostics/scoring.
- `-sg`, `--scoregamma` -> add secondary gamma scoring to the spectrum tree.
- `-ass`, `--azimuthalSurfaceScoring` -> score the azimuthal distribution of outgoing neutrons.
- `-sp`, `--saveEmerging` -> save detailed emerging-particle data.

Any unrecognized non-MPI argument is treated as a Geant4 macro file and executed
with `/control/execute`.

### Output

The standard ROOT output is written as spectrum files named like
`AmBe-Spectrum-N<rank>.root`, with per-thread files merged by the master thread
where applicable.

When `-sp`/`--saveEmerging` is enabled, the simulation also writes
`*-EmergingParticles-N<rank>.root` files containing an `EmergingParticles` TTree.
The tree includes rank, thread, event id, track id, parent id, PDG code, vertex,
momentum, and process information for particles emerging from the source
container.

### Examples

Run a short non-MPI job with four threads:

```bash
AMBE_MPI=false ./ambe -t 4 macros/short.mac
```

Run the N02 casing model with fixed oxide mass and emerging-particle output:

```bash
AMBE_MPI=false ./ambe -t 4 -cs 2 -aom 33 -sp macros/short.mac
```

Run a monoenergetic neutron source:

```bash
AMBE_MPI=false ./ambe -r Single2.45 macros/short.mac
```
