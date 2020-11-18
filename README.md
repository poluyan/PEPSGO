# PEPSGO

Use your custom global optimization algorithm to predict the three-dimensional structure of a peptide from its amino acid sequence. Using multivariate quantile function the continuous search space formed by Rosetta fragments and Dunbrack rotamer library. The search space is unit hypercube. 

Results using Adaptive Differential Evolution. Native (green) vs predicted (cyan).

![Alt text](./pics.png)

## Requirements

To compile from source, you will need:
 * Some Linux distro (Ubuntu/Arch Linux/etc)
 * C++ compiler with C++11 support (current Makefile is GCC-based)
 * [Header-only library mveqf](https://github.com/poluyan/mveqf)
 * OpenMP for multithreading
 * [The Rosetta software suite](https://www.rosettacommons.org/software)

## Installation

 * Get the [Rosetta weekly source package](https://www.rosettacommons.org/software) and build it with SCons.
 * Download the source code:
```
git clone https://github.com/poluyan/PEPSGO
cd PEPSGO/
```
 * Modify the CMakeLists.txt. Set the path to the Rosetta build with the kernel and gcc versions which were used in Rosetta build:
```
set( ROSETTAMAINPATH "/ssdwork/psp_builds/rosetta_src_2020.37.61417_bundle/main" )
set( ROSETTALINUXVER "5.8" )
set( ROSETTACPPCOMP "clang" )
set( ROSETTACPPVER "10.0" )
set( MVEQFDIR "${CMAKE_SOURCE_DIR}" )
```
 * Then compile PEPSGO from source with `cmake .` and `make` commands. This will generate library file `libpepsgo.so` and demos presented in `bin` directory.
 * The file `run_testX.sh` in each test process the following `input` directory and get the proper input for the program.

## Input files

 * `sequence.fasta` Peptide sequence in FASTA format.
 * `fragments.Nmers` Fragment file formatted in Rosetta-style (like the output from fragment picker). Fragments can be any length.
 * `native.data` File with a path to the native structure.
 * `prediction.horiz` Secondary structure prediction in PSIPRED HFORMAT.

## Usage

The `demos/testX` presents the example of using PEPSGO. 

## License

The PEPSGO is distributed under Apache License 2.0 and it is open-source software. Feel free to make a copy and modify the source code, but keep the copyright notice and license intact. The PEPSGO is Rosetta-based, please read [Rosetta Licensing Information](https://www.rosettacommons.org/software) before using.
