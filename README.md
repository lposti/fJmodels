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
and A, B are slopes defined by the model type desired (see [Posti et al. 2015](http://adsabs.harvard.edu/abs/2015MNRAS.447.3060P)).
The above DF is an even function in the azimuthal velocity (since it depends only on the absolute value of the vertical angular momentum ![alt text] (doc/imgs/Jphi.png "vertical angular momentum")), so the models do not roatate. We parametrize rotation adding an odd component in ![alt text] (doc/imgs/Jphi.png "vertical angular momentum"):
![alt text] (doc/imgs/rot.png)

The hyberbolic tangent is an odd function and the parameter `\chi` controls the steepness of the model's rotation curve.

| Name | Description: mandatory for the first component, optional for the second | Default Value |
|:---- |:-----------------------------------------------------------------------:| -------------:|
| `model / 2:model` | Model type: currently `Hernquist`, `Isochrone`, `NFW` | `Hernquist` | 
| `h(J):dphi / 2:h(J):dphi` | DF parameter `\delta_\phi` for h(**J**) | 0.5 |
| `h(J):dz / 2:h(J):dz` | DF parameter `\delta_z` for h(**J**) | 0.5 |
| `g(J):dphi / 2:g(J):dphi` | DF parameter `\delta_\phi` for g(**J**) | 1.0 |
| `g(J):dz / 2:g(J):dz` | DF parameter `\delta_z` for g(**J**) | 1.0 |
| `chi / 2:chi` | Controls the steepness of the rotation curve | 1.0 |
| `mass / 2:mass` | Component's mass (parameter `M0` in [Posti et al. 2015](http://adsabs.harvard.edu/abs/2015MNRAS.447.3060P)) | 1.0 |
| `r0 / 2:r0` | Component's scale radius (defines parameter `J0=sqrt(GM0*r0)` as in [Posti et al. 2015](http://adsabs.harvard.edu/abs/2015MNRAS.447.3060P)) | 1.0 |
| `q / 2:q` | Flattening of the initial guess potential | 1.0 |
|:----- |:-------------------------------------:| --------:|
| `itermax` | Defines the number of iterations computed | 5 | 


### Output

Use [pyfJmod](https://github.com/lposti/pyfJmod) package for analysis and plotting f(J) models data.

## References

Main paper: [Posti et al. (2015)](http://adsabs.harvard.edu/abs/2015MNRAS.447.3060P)

Staeckel Fudge: [Binney (2012)](http://adsabs.harvard.edu/abs/2012MNRAS.426.1324B)
