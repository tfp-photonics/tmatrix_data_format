# T-Matrix Data Format

The **T-Matrix Data Format** is an HDF5-based format for storing T-matrices used in electromagnetic scattering, together with the metadata required for their unambiguous interpretation and exchange between different softwares.

## Reference publication

The format and its conventions are described in

[N. Asadova *et al.*, T-matrix representation of optical scattering response:
Suggestion for a data format, J. Quant. Spectrosc. Radiat. Transfer **333**,
109310 (2025).](https://doi.org/10.1016/j.jqsrt.2024.109310)

An updated version of the manuscript is available on [arXiv:2408.10727](https://arxiv.org/abs/2408.10727), which includes corrections of several typographical errors in the published version.

## Format

A `.tmat.h5` file stores the T-matrix together with the metadata required for its interpretation and reuse. The basic organization of the HDF5 file is:

```text
hdf5 file
├── name, description, keywords, storage_format_version
├── tmatrix
├── modes
│   ├── l
│   ├── m
│   └── polarization
├── frequency
│   └── unit
├── embedding
│   ├── name, description, keywords
│   ├── relative_permittivity
│   └── relative_permeability
├── scatterer
│   ├── material
│   │   ├── name, description, keywords
│   │   ├── relative_permittivity
│   │   └── relative_permeability
│   └── geometry
│       └── name, description, keywords, shape, unit
└── computation
    ├── software, method, name, description, keywords
    ├── files
    ├── mesh
    └── method_parameters
```

The complete specification of the groups, datasets, attributes, and optional entries is given in the reference publication. Example files following the format are available in [`example_data/`](example_data/).

## VSWF convention

The T-matrix relates the incident- and scattered-field expansion coefficients according to

```math
\mathbf{p} = \mathbf{T}\mathbf{a}.
```

The associated Legendre functions are defined as

```math
P_l^m(x)
=
\frac{(-1)^m}{2^l l!}
(1-x^2)^{m/2}
\frac{\mathrm{d}^{\,l+m}}{\mathrm{d}x^{\,l+m}}
(x^2-1)^l .
```

The spherical harmonics are

```math
Y_{lm}(\theta,\phi)
=
\sqrt{
\frac{2l+1}{4\pi}
\frac{(l-m)!}{(l+m)!}
}
P_l^m(\cos\theta)
e^{im\phi}.
```

The vector spherical harmonic is defined as

```math
\mathbf{X}_{lm}(\theta,\phi)
=
\frac{1}{\sqrt{l(l+1)}}
\mathbf{L}Y_{lm}(\theta,\phi),
\qquad
\mathbf{L}
=
\frac{\mathbf{r}\times\nabla}{i}.
```

The vector spherical wave functions are

```math
\mathbf{M}_{lm}^{(n)}
=
z_l^{(n)}(kr)\mathbf{X}_{lm},
```

and

```math
\mathbf{N}_{lm}^{(n)}
=
\frac{\nabla}{k}\times\mathbf{M}_{lm}^{(n)}.
```

For regular and outgoing waves,

```math
z_l^{(1)}(x)=j_l(x),
\qquad
z_l^{(3)}(x)=h_l^{(1)}(x).
```

The adopted time dependence is

```math
e^{-i\omega t}.
```

The complete definitions and conventions are given in Appendices B and C of the reference publication.



## Validation

[`utilities/validate.py`](utilities/validate.py) checks T-matrix files for compliance with the data format.

In addition to structural and metadata checks, the validator can test physical properties declared in the file, including reciprocity, losslessness/passivity, and symmetry relations where applicable.

For example:

```bash
python utilities/validate.py example_data/example.tmat.h5
```

Validation of the file format and applicable physical constraints does not replace numerical convergence tests or verification of the underlying T-matrix calculation.


## Implementations and examples

This repository provides implementations and examples accompanying the published data format. T-matrices can be generated using several methods and code implementations, including ADDA, COMSOL, JCMsuite, Meep, nanoBEM, ONELAB, and SMARTIES.

The repository also contains examples for importing and using standardized T-matrix files with multiscattering codes such as TERMS and treams.

See the corresponding directories for solver-specific implementations and examples.  Where available, additional instructions are provided in the README files within those directories.

## Reference data

[`reference_data/`](reference_data/) contains a reference T-matrix for checking the normalization and conventions of independent implementations. The reference system is an asymmetric arrangement of four spheres designed to break spatial symmetries and thus provide non-zero T-matrix elements for comparison.

[`example_data/`](example_data/) contains representative T-matrix files for different scatterers and can be used for testing readers, converters, visualization tools, and other software using the format.


## T-matrix database

T-matrix files following this data format can be submitted to the
[Daphona T-matrix Portal](https://tmatrix.scc.kit.edu/). The portal provides a searchable database of standardized T-matrix data for their distribution and reuse.

When using data from or submitting data to the database, you can cite:

[N. Asadova, K. Boussaoud, J. Meyer, F. Tristram, and C. Rockstuhl, *T-matrix database to promote information-driven research in nanophotonics*, Opt. Mater. Express **16**, 1534–1550 (2026).](https://doi.org/10.1364/OME.592265)

## Contributing

Contributions are welcome, including support for additional T-matrix solvers etc.

If you would like to add support for the T-Matrix Data Format to another numerical or analytical method, the utilities in [`tmatrix_tools`](tmatrix_tools/) provide reusable functionality for creating compliant `.tmat.h5` files and can be used as a basis for a new implementation. 

New implementations should include a check demonstrating consistency with the normalization and conventions of the data format. The reference tetrahedron provided in reference_data/ can be used as a benchmark for this purpose, although other suitable validation approaches are also possible.

Please also provide an example script that computes a T-matrix with the respective method and stores it in the standard format. After review, the example can be included in this repository alongside the existing solver-specific implementations.
