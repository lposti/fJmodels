# f(J) models code  [![Build Status](https://drone.io/github.com/lposti/fJmodels/status.png)](https://drone.io/github.com/lposti/fJmodels/latest) 

Code for generating action-based distribution function models for spheroids,
as described in [Posti et al. (2014)](http://adsabs.harvard.edu/abs/2014arXiv1411.7897P).

## Installation

Simply download the package as a `.tar.gz` file (see Releases) or clone it
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
and A, B are slopes defined by the model type desired (see [Posti et al. (2015)](http://adsabs.harvard.edu/abs/2015MNRAS.447.3060P)).
The above DF is an even function in the azimuthal velocity (since it depends only on the absolute value of the vertical angular momentum $J_\phi$), so the models do not roatate. We parametrize rotation adding an odd component in $J_\phi$:
$$
f_{\rm tot}({\bf J}) = (1-k)f({\bf J}) + k\tanh\left(\frac{\chi J_\phi}{J_0}\right)f({\bf J})
$$


One component is mandatory (and also all its parameters are), the second is optional.
- `itermax [optional]` defines the number of iterations computed. It is defaulted to 5.
- `model / 2:model [mandatory / optional]` 1- and 2-component model type: currently it can be `Hernquist`, `Isochrone`, `NFW`
- `h(J):dphi / 2:h(J):dphi [mandatory / optional]` 1- and 2-component model DF parameter
- `h(J):dz / 2:h(J):dz [mandatory / optional]` 1- and 2-component model DF parameter
- `g(J):dphi / 2:g(J):dphi [mandatory / optional]` 1- and 2-component model DF parameter
- `g(J):dz / 2:g(J):dz [mandatory / optional]` 1- and 2-component model DF parameter
- `chi / 2:chi [mandatory / optional]` 1- and 2-component model DF parameter controlling the steepness of the rotation curve (0 is non-rotating).
- `mass / 2:mass [mandatory / optional]` 1- and 2-component model mass (parameter `M0`)
- `r0 / 2:r0 [mandatory / optional]` 1- and 2-component model scale radius: together with `M0` defines `J0=sqrt(GM0*r0)`
- `q / 2:q [mandatory / optional]` 1- and 2-component model initial flattening of the guess potential

## References

Main paper: [Posti et al. (2015)](http://adsabs.harvard.edu/abs/2015MNRAS.447.3060P)

Staeckel Fudge: [Binney (2012)](http://adsabs.harvard.edu/abs/2012MNRAS.426.1324B)
