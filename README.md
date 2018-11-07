## What is `doa-tools`

A set of MATLAB functions for direction-of-arrival (DOA) estimation related
applications, including basic array designs, various DOA estimators, and tools
to compute performance bounds. It serves as a small toolbox for
[my research](http://research.wmz.ninja/research.html)
related to array signal processing.

This repository also contains MATLAB codes that produce figures similar to
those in my papers (may not be exactly the same due to the randomness of Monte
Carlo simulations). You can browse them [here](examples/experiments).

**If somehow you do not have access to MATLAB, you can use the Python version [here](https://github.com/morriswmz/doatools.py).**

## Features

* Several array design and difference coarray related functions.
* Commonly used DOA estimators including MVDR beamformer, MUSIC,
  root-MUSIC, ESPRIT, etc.
* Sparsity-based DOA estimator.
* Functions to compute the [Cramér-Rao bounds](https://en.wikipedia.org/wiki/Cram%C3%A9r%E2%80%93Rao_bound).
* Functions to compute the asymptotic covariance matrix of the estimation errors
  of MUSIC estimators (including difference coarray based).
* Functions to visualize the estimation results.
* Several useful utility functions, including a simple progress bar to
  display the simulation progress.

## Getting started

Run `startup.m` to add `doa-tools` to your MATLAB search path.

Refer to the examples [here](examples/).

## License

The source code is released under the [MIT](LICENSE.md) license.

## Citation

If you find my code helpful. You are welcomed to cite my papers
[here](http://research.wmz.ninja/research.html).
