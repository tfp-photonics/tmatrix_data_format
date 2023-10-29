# T-Matrix Computation with Comsol

Compute the T-Matrix of arbitrarily shaped objects and using bi-isotropic media. The
embedding medium must be reciprocal. If the geometry allows it, for example
rotationally symmetric objects, computations in a reduced number of dimensions are
possible.

## Warning

Comsol uses the time evolution $`\exp(i \omega t)`$. This means that all values,
especially the material property parameters, must be complex conjugated with respect to
the opposite convention $`\exp(-i \omega t)`$ that is often used in physics.

## File format

The models are distributed using java code for two reasons. First, it simplifies
portability between Comsol versions due to the relatively stable API. Second, text-based
files are better suited for git.

## Comsol Version Requirement

Speaking of compatibility, there is a major difference between Comsol 5.4 and Comsol 5.5
regarding this code: the associated Legendre polynomials are not implemented in the
former one. I was not able to implement them using obvious ways, namely recursion or a
series representation. Instead, I implemented them manually up to degree 5, such that
(3D) T-matrix computations are possible up to degree 4. If you want to use them,
uncomment the corresponding lines in the java files. For Comsol 5.5 and above no such
restrictions apply.

# General Usage

Compile the files using `comsolcompile foo.java` on Windows and
`comsol compile foo.java` on Linux. Then, run the file `foo.class` with Comsol, for
example with `comsolbatch -inputfile foo.class` or `comsol batch -inputfile foo.class`,
respectively. Alternatively, if you want to adjust the model in the GUI, comment the
last lines in the java files to avoid running the studies immediately after opening the
class file. You can then run the studies manually in the GUI.

Define your object, your materials (see warning), and the parameters for the
decomposition. The latter includes the modes that are used in the calculation and
`r_decomp` the radius of the decomposition. The sphere (or cylinder for cylindrical
T-matrices) defined by this radius must completely enclose the object of interest.
Additionally, define the radius of the domain `r_domain >= r_decomp` and the thickness
of the perfectly matched layers (PML) `d_pml` around the domain. These three parameters
define the first object in the geometry, which uses the first material for its
properties. The first object and material are later used for the embedding material in
the decomposition of the scattered fields.

In the case of cylindrical T-matrices in a 3D domain or 2D axisymmetric domain, one has
to additionally specify the periodicity `az`.

As computation object, all of the examples use either a sphere or an infinitely extended
cylinder defined by `r_obj`. Change the geometry or the parameters to your liking.

In the first study the actual solutions of the scattering problems are calculated. A
frequency sweep is assumed per default. Adjust it to your needs. You can define
additionally other parametric sweeps. Run the first study manually, if needed.

The second study is only used as an auxiliary for the post-processing. It does not solve
anything, but evaluates the necessary integrals and computes the T-matrix entries. Make
sure to mirror the frequency sweep of the first study and also all parameter sweeps, you
added to the first study. Run this second study. Then you can get the results via
*Global Evaluation* of the second study in the *Results* section, returning the
variables `ap` and `am`.

## Summary

0. Compile model
1. Define geometry
2. Define materials (remember complex conjugation)
3. Define decomposition parameters
4. Define frequency and other parameter sweeps
5. Run both studies
6. Calculate derived values
7. Export table

# T-Matrix

The usual T-Matrix uses vector spherical wave (VSW) functions as a basis set. To include
chiral materials, we use the VSWs in a formulation of well-defined helicity. For
non-reciprocal media, the eigensolutions get much more complicated, so they cannot be
included in a simple way (for the embedding medium). The same is true for anisotropy, so
the embedding medium must be homogeneous, isotropic and reciprocal. An additional
constraint is, that the object must have a finite size.

There are two java files for the T-matrix calculation, `tmatrix.java` and
`tmatrix_axisym.java` for general three-dimensional objects and for axisymmetric
objects, respectively. If applicable, the latter is more accurate, since the geometry
that needs to be solved is only two-dimensional and the diagonality of the T-matrix with
respect to `m` is enforced.

## Math

We start by defining the vector spherical harmonics

```math
\boldsymbol X_{lm}(\theta, \varphi)
= \frac{1}{\sqrt{l(l+1)}} \boldsymbol L Y_{lm}(\theta, \varphi) \\
= \underbrace{\mathrm i \sqrt{\frac{(2l+1)}{4\pi l (l+1)} \frac{(l-m)!}{(l+m)!}}}_{\gamma_{lm}}
\mathrm e^{\mathrm i m \varphi}
\bigg( \boldsymbol{\hat{\theta}}
\underbrace{\frac{\mathrm i m}{\sin\theta} P_l^m(\cos\theta)}_{\mathrm i \pi_{lm}(\theta)}
- \boldsymbol{\hat{\varphi}}
\underbrace{\frac{\partial}{\partial\theta} P_l^m(\cos\theta)}_{\tau_{lm}(\theta)}
\bigg)
```

and get the vector spherical wave functions with

```math
\boldsymbol M_{lm}^{(n)}(kr, \theta, \varphi)
= \boldsymbol X_{lm}(\theta, \varphi) z_l^{(n)}(kr) \\
\boldsymbol N_{lm}^{(n)}(kr, \theta, \varphi)
= \frac{\nabla}{k} \times \boldsymbol M_{lm}^{(n)}{lm}(kr, \theta, \varphi) \\
\boldsymbol A_{lmp}^{(n)}(k_pr, \theta, \varphi)
= \frac{1}{\sqrt{2}}
\left(
\boldsymbol N_{lm}^{(n)}(k_pr, \theta, \varphi)
+ p \boldsymbol M_{lm}^{(n)}(k_pr, \theta, \varphi)
\right)\,.
```

The latter are modes of well-defined helicity, where we introduced the wave numbers
$`k_\pm = k_0 (\sqrt{\epsilon_r \mu_r} \pm \kappa)\,`$. The incident wave uses spherical
Bessel functions of the first kind ($`n = 1`$). The scattered wave is expressed with
spherical Hankel functions of the first kind ($`n = 1`$).

We can decompose the scattered field as

```math
\boldsymbol E_{\text{sca}}(\boldsymbol r)
= \sum_{l,m,p} a_{lmp} \boldsymbol A_{lmp}^{(3)}(k_pr, \theta, \varphi)\,.
```

We can project onto the different modes by using

```math
\int \mathrm d\Omega
\boldsymbol X_{lm}^\ast(\theta, \varphi)
\boldsymbol M_{lm}^{(n)}(kr, \theta, \varphi)
= z_l^{(n)}(kr) \\
\int \mathrm d\Omega
\boldsymbol r \times \boldsymbol X_{lm}^\ast(\theta, \varphi)
\boldsymbol N_{lm}^{(n)}(kr, \theta, \varphi)
= \frac{z_l^{(n)}(kr)}{kr} + \frac{\partial }{\partial(kr)} z_l^{(n)}(kr) \\
\int \mathrm d\Omega
\boldsymbol X_{lm}^\ast(\theta, \varphi)
\boldsymbol N_{lm}^{(n)}(kr, \theta, \varphi)
= 0 \\
\int \mathrm d\Omega
\boldsymbol r \times \boldsymbol X_{lm}^\ast(\theta, \varphi)
\boldsymbol M_{lm}^{(n)}(kr, \theta, \varphi)
= 0
```

such that the coefficients $`a_{lm+}`$ and $`a_{lm-}`$ can be computed by solving

```math
\frac{1}{\sqrt{2}}
\begin{pmatrix}
h_l^{(1)}(k_+ r) & -h_l^{(1)}(k_- r) \\
\psi_l^{(1)}(k_+ r) & \psi_l^{(1)}(k_- r) \\
\end{pmatrix}
\begin{pmatrix}
a_{lm,+} \\
a_{lm,-}
\end{pmatrix}
=
\begin{pmatrix}
\int \mathrm d\Omega \boldsymbol X_{lm}^\ast(\theta, \varphi) \boldsymbol E_{\text{sca}}(\boldsymbol r) \\
\int \mathrm d\Omega \boldsymbol r \times \boldsymbol X_{lm}^\ast(\theta, \varphi) \boldsymbol E_{\text{sca}}(\boldsymbol r)
\end{pmatrix}
```
where $`\psi_l^{(n)}(x) = \frac{h_l^{(n)}(x)}{x} + \frac{\partial}{\partial x} h_l^{(n)}(x)\,`$.

# Cylindrical T-Matrix

Two-dimensional (or cylindrical) T-matrices use vector cylindrical wave (VCW) functions
as a basis set. All remarks regarding bi-anisotropy for the T-matrix also apply here.
The difference is, that the VCWs allow the objects to be infinitely extended in the
z-direction, either being uniform or periodic along this axis.

There are three java files for the calculation, `tmatrixc.java`, `tmatrixc_axisym.java`,
and `tmatrixc_uni.java`. The first two of them are for periodic objects. Equivalently to
the T-matrix case, the object can either have a general shape or be axisymmetric. The
last file is for infinitely extended uniform objects. This has, like the axisymmetric
case, only a two-dimensional computation domain. Again, axisymmetry enforces diagonality
with respect to `m`.

In the periodic cases only `kz` values that differ by a multiple of the reciprocal
lattice vector are considered in the model, the number of reciprocal lattice vectors
included in each direction is `n_kz`. For uniform cases, different `kz` values do not
couple.

## Math

The cylindrical vector waves are defined as

```math
\boldsymbol M_{mk_z}^{(n)}(\rho, \varphi, z)
= \mathrm e^{\mathrm i m \varphi + \mathrm i k_z z}
\left[
\mathrm i m \frac{Z_m^{(n)}(k_\rho \rho)}{k_\rho \rho} \boldsymbol{\hat\rho}
- Z_m^{(n)\prime}(k_\rho \rho) \boldsymbol{\hat\varphi}
\right] \\
\boldsymbol N_{mk_z}^{(n)}(\rho, \varphi, z)
= \frac{\nabla}{k} \times \boldsymbol M_{mk_z}^{(n)}(\rho, \varphi, z ) \\
\boldsymbol A_{m k_z p}^{(n)}(\rho, \varphi, z)
= \frac{1}{\sqrt{2}} \left(\boldsymbol N_{mk_z}^{(n)}(\rho, \varphi, z )
+ p \boldsymbol M_{mk_z}^{(n)}(\rho, \varphi, z ) \right)\,.
```

where $`k_\rho = \sqrt{k^2 - k_z^2}`$ and $`Z_m^{(n)}`$ are the Bessel or Hankel
functions. The projection onto different modes in the case of chiral media, and
therefore different $`k`$ and $`k_\rho`$ values, is not as straightforward as in the
cylindrical case. Similar to the spherical case we expand the scattered wave using
Hankel functions. We use the integrals

```math
I_1
= \frac{1}{2\pi a_z}
\int_0^{2\pi} \mathrm d \varphi
\int_0^{a_z} \mathrm d_z
\mathrm e^{\mathrm i m \varphi + \mathrm i k_z z}
\boldsymbol{\hat z} \boldsymbol E_{\text{sca}}(\boldsymbol r)
= \frac{1}{\sqrt 2}
\left(
\frac{k_{\rho+} H_m^{(1)}(k_{\rho+} \rho)}{k_{\rho+} \rho} a_{mk_z +}
+ \frac{k_{\rho-} H_m^{(1)}(k_{\rho-} \rho)}{k_{\rho-} \rho} a_{mk_z -}
\right) \\
I_2
= \frac{1}{2\pi a_z}
\int_0^{2\pi} \mathrm d \varphi
\int_0^{a_z} \mathrm d_z
\mathrm e^{\mathrm i m \varphi + \mathrm i k_z z}
\left(
-\mathrm i m \frac{H_m^{(1)}(k_{\rho+} \rho)}{k_{\rho+} \rho} \boldsymbol{\hat\rho}
- H_m^{(1)\prime}(k_{\rho+} \rho) \boldsymbol{\hat\varphi}\right)
\boldsymbol E_{\text{sca}}(\boldsymbol r) \\
= \frac{1}{\sqrt 2}
\left[
H_{m+1}^{(1)}(k_{\rho+} \rho) H_{m-1}^{(1)}(k_{\rho+} \rho) a_{mk_z +}
- \frac{H_{m+1}^{(1)}(k_{\rho+} \rho) H_{m-1}^{(1)}(k_{\rho-} \rho)(k_- - k_z)
+ H_{m+1}^{(1)}(k_{\rho-} \rho) H_{m-1}^{(1)}(k_{\rho+} \rho)(k_- + k_z)}{2k_-}a_{mk_z -}
\right]
```

to get a system of linear equations to determine $`a_{mk_z\pm}`$. The choice of these
integrals is somewhat arbitrary, e.g., for the second integral one could use equally
well the values for negative polarization. In the case of $`k_z=0`$ these integrals are
similar to the spherical case integral in the sense, that they separate
$`\boldsymbol M`$ and $`\boldsymbol N`$ modes.

# Closing remarks

The provided files can also be used as a starting point for any bi-isotropic calculation
and are not necessarily restricted to the computation of T-matrix coefficients.

# Possible enhancements

* Add bi-anisotropy for the object. Anisotropic permittivity and permeability may work
  out of the box.
* Calculate the cylindrical T-matrix for objects sandwiched between PECs
