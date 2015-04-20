# f(J) models code  [![Build Status](https://drone.io/github.com/lposti/fJmodels/status.png)](https://drone.io/github.com/lposti/fJmodels/latest) 

Code for generating action-based distribution function models for spheroids,
as described in [Posti et al. (2015)](http://adsabs.harvard.edu/abs/2015MNRAS.447.3060P).

## Features

The code generates an axisymmetric self-consistent model by specifying a Distribution Function (DF) which is a double-power law function of the action integrals. The model's density (and higher moments) are computed integrating the DF, then the self-consistent potential is found by iterative procedure (convergence criteria are currently in the form of distance from the potential of the previous iteration and on the tensor Virial theorem). The action integrals are computed in the Staeckel Fudge ([Binney 2012](http://adsabs.harvard.edu/abs/2012MNRAS.426.1324B)) approximation.

- the models have a double-power law density profile: the central and outer slopes in the density profile are determined mainly by the DF's power law slope in the small- and large-action regimes. The code has currently *hardwired* such DF's parameters such that the final profile follows that of classical models, such as Hernquist, Isochrone and NFW.
- the code stores the coefficients of the multipole expansion of e.g., density, potential and velocity dispersion distributions in output files at every iterations. Such output files can be managed easily with the Python tool [pyfJmod](https://github.com/lposti/pyfJmod).
- the models can currently be
  - self-consistent single-component
  - single-component with external potential (to be fully implemented...)
  - self-consistent two-components
- it uses a  22-points *tanh-sinh quadrature* (see [Takahasi & Mori 1973](http://www.ems-ph.org/journals/show_abstract.php?issn=0034-5318&vol=9&iss=3&rank=12)) three-dimensional integrator to optimize a (typically) highly peaked DF at the centre. Different choices for the integration rule (e.g., Gauss-Legendre) and for the number of points used are implemented and can be switched on.
- the code benefits from a multi-threaded OpenMP implementation, vector-register instructions (SSE/AVX) and compiler's optimization flags.

## Installation

Simply download the package as a `.tar.gz` file (see Releases) or clone it (*recommended*)
```
git clone https://github.com/lposti/fJmodels
```

Then `cd` to the unpacked/cloned directory and
```
./configure && make
```

If the installation was successful, an executable named `fJmodels` was created in the current directory.

## Usage

The code can be launched as follows:
```
./fJmodels <input_filename>
```
where the argument is the filename of the parameter file and it is optional: if not specified the code assumes that the input parameter file is `param.txt`.
The outputs will be written in the directory named `models`, which must exist.

### Parameter file

The code can currently compute both 1-component and 2-component models, for which all the parameters are specified in the parameter file.
An example is given as the file `param.txt`.

All the DFs are of the following general form:
![alt text][DF]

[DF]: doc/imgs/DF.png "Distribution Function"
and A, B are slopes defined by the model type desired (see [Posti et al. 2015](http://adsabs.harvard.edu/abs/2015MNRAS.447.3060P)).
The above DF is an even function in the azimuthal velocity (since it depends only on the absolute value of the vertical angular momentum ![alt text] (doc/imgs/Jphi.png "vertical angular momentum")), so the models do not roatate. We parametrize rotation adding an odd component in ![alt text] (doc/imgs/Jphi.png "vertical angular momentum"):
![alt text] (doc/imgs/rot.png)

The hyberbolic tangent is an odd function and the parameter `\chi` controls the steepness of the model's rotation curve. If `\chi` is greater than zero, the odd term is added to the DF with k=0.5 (currently *hardwired* for a maximally rotating configuration), else it isn't. 

Follows a schematic descriptions of the implemented parameters that can be set in the file to specify the model.

| Parameter Name | Description: mandatory for the first component, optional for the second | Default Value |
|:-------------- |:-----------------------------------------------------------------------:| -------------:|
| `model / 2:model` | Model type: currently `Hernquist`, `Isochrone`, `NFW` | `Hernquist` | 
| `h(J):dphi / 2:h(J):dphi` | DF parameter `\delta_\phi` for h(**J**) | 0.55 |
| `h(J):dz / 2:h(J):dz` | DF parameter `\delta_z` for h(**J**) | 0.55 |
| `g(J):dphi / 2:g(J):dphi` | DF parameter `\delta_\phi` for g(**J**) | 1.0 |
| `g(J):dz / 2:g(J):dz` | DF parameter `\delta_z` for g(**J**) | 1.0 |
| `chi / 2:chi` | Controls the steepness of the rotation curve. If greater than 0 an odd part is added to the DF | 0.0 |
| `mass / 2:mass` | Component's mass (parameter `M_0` in [Posti et al. 2015](http://adsabs.harvard.edu/abs/2015MNRAS.447.3060P)) | 1.0 |
| `r0 / 2:r0` | Component's scale radius (defines parameter `J_0=sqrt(GM_0*r_0)` as in [Posti et al. 2015](http://adsabs.harvard.edu/abs/2015MNRAS.447.3060P)) | 1.0 |
| `q / 2:q` | Flattening of the initial guess potential | 1.0 |
| `itermax` | Defines the number of iterations computed | 5 | 

The default configuration generates a nearly-isotropic, spherical, non-rotating, Hernquist-like model. 

### Output

Use [pyfJmod](https://github.com/lposti/pyfJmod) package for analysis and plotting f(J) models data.

## References

Main paper: [Posti et al. (2015)](http://adsabs.harvard.edu/abs/2015MNRAS.447.3060P)

Staeckel Fudge: [Binney (2012)](http://adsabs.harvard.edu/abs/2012MNRAS.426.1324B)
