# pyOptSparse - PYthon OPTimization (Sparse) Framework

[![Build Status](https://travis-ci.com/mdolab/pyoptsparse.svg?branch=master)](https://travis-ci.com/mdolab/pyoptsparse)
[![Coverage Status](https://coveralls.io/repos/github/mdolab/pyoptsparse/badge.svg?branch=master)](https://coveralls.io/github/mdolab/pyoptsparse?branch=master)
[![Documentation Status](https://readthedocs.com/projects/mdolab-pyoptsparse/badge/?version=latest)](https://mdolab-pyoptsparse.readthedocs-hosted.com/en/latest/?badge=latest)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

pyOptSparse is an object-oriented framework for formulating and solving nonlinear constrained optimization problems in an efficient, reusable, and portable manner.
It is a fork of pyOpt that uses sparse matrices throughout the code to more efficiently handle large-scale optimization problems.
Many optimization techniques can be used in pyOptSparse, including both gradient-based and gradient-free methods.
A visualization tool called OptView also comes packaged with pyOptSparse, which shows the optimization history through an interactive GUI.
An example output from OptView is shown below.

![Example](doc/OptView.png)

## Optimizer Support
pyOptSparse provides Python interfaces for a number of optimizers.
ALPSO, CONMIN, IPOPT, NLPQLP, NSGA2, PSQP, SLSQP, ParOpt and SNOPT are currently tested and supported.
NOMAD interface is also provided, but it is not tested nor supported.

We do not provide the source code for SNOPT and NLPQLP, due to their restrictive license requirements.
Please contact the authors of the respective optimizers if you wish to obtain them.
Furthermore, ParOpt and IPOPT are available as a open source package but must be installed separately.
Please see the documentation page of each optimizer for purchase and installation instructions.

## Documentation
Please see the [documentation](https://mdolab-pyoptsparse.readthedocs-hosted.com/) for installation details and API documentation.

To locally build the documentation, enter the `doc` folder and enter `make html` in terminal.
You can then view the built documentation in the `_build` folder.

## Testing
Testing is done with the `testflo` package developed by the openMDAO team, which can be installed via `pip install testflo`.
To run the tests, simply type `testflo .` in the root directory.

## Citation
A pyOptSparse journal paper does not exist, instead please cite pyOpt and the authors of the respective optimization
algorithms in any publication for which you find it useful.
For more background, theory, and figures, see the [pyOpt journal article](http://mdolab.engin.umich.edu/sites/default/files/pyOpt.pdf).

Perez, R. E., Jansen, P. W., and Martins, J. R. R. A., “pyOpt: A Python-Based Object-Oriented Framework for Nonlinear
Constrained Optimization,” Structural and Multidisciplinary Optimization, Vol. 45, No. 1, January 2012, pp. 101–118.
doi:10.1007/s00158-011-0666-3.

```
@article{Perez2012a,
	Author = {Ruben E. Perez and Peter W. Jansen and Joaquim R. R. A. Martins},
	Doi = {10.1007/s00158-011-0666-3},
	Journal = {Structural and Multidisciplinary Optimization},
	Month = {January},
	Number = {1},
	Pages = {101--118},
	Title = {{pyOpt}: A {Python}-Based Object-Oriented Framework for Nonlinear Constrained Optimization},
	Volume = {45},
	Year = {2012},
	Annote = {10.1007/s00158-011-0666-3}}
```

## License
pyOptSparse is licensed under the GNU Lesser General Public License.
See `LICENSE` for the full license.

## Copyright
Copyright (c) 2011 University of Toronto\
Copyright (c) 2014 University of Michigan\
Additional copyright (c) 2014 Gaetan K. W. Kenway, Ruben Perez, Charles A. Mader, and\
Joaquim R. R. A. Martins\
All rights reserved.
