# treams

The package `treams` provides a framework to simplify computations of the electromagnetic scattering of waves at finite and at periodic arrangements of particles based on the T-matrix method.

## Installation

### Installation using pip

To install the package with `pip`, use:

```bash
pip install treams
```
If you're using the system-wide installed version of Python, you might consider the --user option.


For further information, please visit: [treams GitHub Repository](https://github.com/tfp-photonics/treams).

## Description
The example here is a demonstration of integrating  the T-matrix format project requirements into  `treams`. It is possible to load a `tmat.h5` file into the native T-matrix class, which includes reading the necessary data such as modes and embedding, and automatic unit conversion. Next, it is simple to compute optical quantities of the scatterer or an arrangement of scatterers. 
We calculate the coupling between two reference TiO2 cylinders with a radius of 250 nm and a height of 300 nm in air separated horizontally by 600 nm. The resulting average extinction cross-sections for the individual cylinder and the two cylinders are plotted.