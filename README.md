# PEPSGO

PEPtide Structure prediction by Global Optimization of a potential energy function. Use your custom global optimization algorithm to predict the three-dimensional structure of a peptide from its amino acid sequence. It uses Rosetta Dunbrack rotamer library

## Requirements

To compile from source, you will need:
 * Some Unix distro (Ubuntu/Arch Linux/etc)
 * C++ compiler with C++11 support (current Makefile is GCC-based)
 * [The Rosetta software suite](https://www.rosettacommons.org/software)

## License

The PEPSGO is distributed under Apache License 2.0 and it is open-source software. The PEPSGO is Rosetta-based, please read [Rosetta Licensing Information](https://www.rosettacommons.org/software) before using.

## Some results

Using Adaptive Differential Evolution. Native (green) vs predicted (cyan)

![Alt text](./pic.png)