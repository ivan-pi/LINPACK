# LINPACK

```
LINPACK is a collection of Fortran subroutines that analyze and
solve linear equations and linear least-squares probles.  The
package solves linear systems whose matrices are general, banded,
symmetric indefinite, symmetric positive definite, triangular,
and tridiagonal square.  In addition, the package computes
the QR and singular value decompositions of rectangular matrices
and applies them to least-squares problems.  LINPACK uses
column-oriented algorithms to increase efficiency by preserving
locality of reference.

LINPACK was designed for supercomputers in use in the 1970s and
early 1980s.  LINPACK has been largely superceded by LAPACK
which has been designed to run efficiently on shared-memory, vector
supercomputers.
```

**Original site at Netlib:** http://www.netlib.org/linpack/

**Users' Guide:**

> Dongarra, J. J., Moler, C. B., Bunch, J. R., & Stewart, G. W. (1979). *LINPACK users' guide*. Society for Industrial and Applied Mathematics. https://doi.org/10.1137/1.9781611971811

## fpm usage

To use LINPACK within your `fpm` project, add the following to your package manifest file (`fpm.toml`):

```toml
[dependencies]
LINPACK = { git = "https://github.com/ivan-pi/LINPACK.git" }
```
