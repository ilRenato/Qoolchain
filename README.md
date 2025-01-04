# Qoolchain

Qoolchain is a QUBO preprocessing toolchain designed to reduce problem size and improve the performance of optimization solvers. It is compatible with any QUBO-compliant solver and, in particular, it is optimized for the **Grover Adaptive Search** (**GAS**) hybrid quantum-classical algorithm.
The toolchain is implemented in C++ and compiled into a Python extension module using **Cython**. The repository also contains a collection of tests and benchmarks to validate and compare its performance with state-of-the-art implementations.

The contents of this repository, along with the Qoolchain objectives and results are described in the following paper. Please cite this paper if you use this repository and find it helpful.

**Qoolchain: A QUBO Preprocessing Toolchain for Enhancing Quantum Optimization**  
Authors: Giacomo Orlandi, Deborah Volpe, Mariagrazia Graziano, Giovanna Turvani 
Published in: *Advanced Quantum Technologies, 2024*  
DOI: [https://doi.org/10.1002/qute.202400384](https://doi.org/10.1002/qute.202400384) 

---

## Table of Contents
1. [Overview](#overview)
2. [Directory Structure](#directory-structure)
3. [Attribution](#attribution)
4. [Dependencies](#dependencies)
5. [Building the Python Module](#building-the-python-module)
6. [Usage](#usage)
7. [Running Tests and Benchmarks](#running-tests-and-benchmarks)
8. [License](#license)

---

## Overview

This project provides:
1. A **C++ implementation** compiled into a Python extension module using **Cython**.
2. A suite of **tests and benchmarks** to validate correctness and measure performance.

The goal of this project is to combine the speed of C++ with the usability of Python, allowing users to call highly optimized C++ functions seamlessly from Python.

---

## Directory Structure

The repository is organized as follows:

.
├── Benchmark_tests
│   ├── bounds_utils.py
│   ├── CompareResult.py
│   ├── create_Q_matrix.py
│   ├── Dwave_Toolchain
│   │   └── …
│   ├── GraphColoring
│   │   └── …
│   ├── MaxClique
│   │   └── …
│   ├── MaxCut
│   │   ├── Gset
│   │   │   └── …
│   │   └── …
│   ├── MinimumVertexCover
│   │   └── …
│   ├── NumberPartitioning
│   │   └── …
│   └── ToolchainTestScripts.py
├── LICENSE
├── README.md
├── requirements.txt
├── structure.txt
└── Toolchain_implementation
    ├── FlowFunctions.cpp
    ├── FlowFunctions.h
    ├── GraphFunctions.cpp
    ├── GraphFunctions.h
    ├── ImplicationNetwork.cpp
    ├── ImplicationNetwork.h
    ├── PreprocessingToolchain.cpp
    ├── PreprocessingToolchain.pyx
	├── PreQubo.cpp
    ├── PreQubo.h
    ├── ResidualNetwork.cpp
    ├── ResidualNetwork.h
	└── setup.py
	
A more detailed directory structure is reported in `structure.txt`.
	
---

## Attribution

This repository includes code adapted from the [D-Wave preprocessing toolchain](https://github.com/dwavesystems/dwave-preprocessing) under the Apache-2.0 license.

The files and folders in the `Benchmark_tests/Dwave_Toolchain` directory are based on work from the Original Repository. Changes have been made to the code to make comparisons with this project.
The only modification required was in the file `Benchmark_tests/Dwave_Toolchain/dwave/preprocessing/composites/fix_variables.py`, to extract the fixed variables (or persistencies) for comparison between those found by Qoolchain and by the Dwave toolchain.
Specifically, we modified a single line of code as follows:
```bash
# Original line
return SampleSet.from_future(sampleset, _hook)
# Modified line
return SampleSet.from_future(sampleset, _hook), fixed_variables
```

Original repository:
- Name: dwave-preprocessing
- URL: https://github.com/dwavesystems/dwave-preprocessing
- License: Apache-2.0

The tests conducted on the graph coloring problem (in the `Benchmark_tests/GraphColoring` folder) use the FullIns and Myciel benchmarks from [https://mat.tepper.cmu.edu/COLOR02/](https://mat.tepper.cmu.edu/COLOR02/).
The Max-cut Gset tests (in the `Benchmark_tests/MaxCut/Gset` folder) use benchmarks from [https://web.stanford.edu/~yyye/yyye/Gset/](https://web.stanford.edu/~yyye/yyye/Gset/).

---

## Dependencies

Before building and running the project, ensure you have **Python**==3.9.21 installed with the packages in `requirements.txt` and **GCC** for compiling the C++ code.

---

## Building the Python module

To build the C++ implementation into a Python extension module run the setup script in the `Toolchain_implementation` folder as follows:
```bash
python setup.py build_ext --inplace
```

---

## Usage

Once the Python extension module has been built, Qoolchain can be imported in a Python file by adding the module path to the Python’s module search path, and finally just using:
```python
import QuboPreprocessing
```

---

## Running Tests and Benchmarks

To reproduce the benchmark tests, the `Benchmark_tests/Dwave_Toolchain` folder contains a slightly modified version of the preprocessing toolchain provided by D-Wave, as described in the [Attribution](#attribution) section.
The modifications were necessary to extract variables identified as persistencies without solving the optimization problem. The original toolchain, when installed and imported, solves the optimization problem directly after preprocessing. However, for comparative purposes (e.g., measuring reduction percentages in the initial QUBO problem size), it was needed to obtain the persistencies data structure.
In order to use it and compare results with Qoolchain, compile the Python extension module with the following command:
```bash
python setup.py build_ext --inplace
```
Once compiled, this modified version can be used by all DwavePreprocessing.py files inside each benchmark folder to run tests on the D-Wave toolchain.
NOTE: This procedure and the Python package versions (specified below) were tested for the article and are compatible with Windows only. On Linux, a compatibility issue arises when Cythonizing D-Wave Pyrex files due to the ConstNumeric definition in cyutilities.pxd from the dimod package.
The following package versions were used to ensure reproducibility of the results:
- Cython 3.0.10
- dimod	0.12.14
- numpy	1.26.4
- setuptools 60.2.0

The version of the D-Wave toolchain included in this repository is configured to work with these package versions to reproduce the same results of the published paper.
However, if you wish to run the latest version of the toolchain, only a single line of code needs to be modified. You can make this change directly in the file mentioned in the [Attribution](#attribution) section.

Therefore, to run the latest version of the toolchain also in a Linux machine:
- clone the D-wave repository ([https://github.com/dwavesystems/dwave-preprocessing](https://github.com/dwavesystems/dwave-preprocessing)),
- update the fix_variables.py file by applying the modification described above.
- build the toolchain by launching the setup.py as above.

---

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.