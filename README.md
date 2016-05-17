# B2AC - Ellipse fitting in Python and C

|  Branch       | Build status     |
| :------------ | ---------------: |
| `master`      | [![Build Status](https://travis-ci.org/hbldh/b2ac.svg?branch=master)](https://travis-ci.org/hbldh/b2ac) |
| `develop`     | [![Build Status](https://travis-ci.org/hbldh/b2ac.svg?branch=develop)](https://travis-ci.org/hbldh/b2ac) |

Python and C implementations of an ellipse fitting algorithm in double and fixed point precision. 

The implementations here were never meant for production, but were written as an examination of how much
speedup could be attained from implementing ellipse fitting methods in fixed point and at what precision cost.

The ellipse fitting method in \[2\] (an improved version of the one in \[1\]) was implemented.

## Installation

Install simply by calling:

    pip install git+https://www.github.com/hbldh/b2ac
    
Numpy needs to be installed prior to `b2ac`, since it compiles the C extension using Numpy headers.

## Testing

Test with nosetests:

    nosetests tests

## Usage

An ellipse can be fitted using Python methods:

```python
import numpy as np
import b2ac.preprocess
import b2ac.fit
import b2ac.conversion

points = np.array([[3, 0], [0, 5], [3, 10], [6, 5]])

# Optional removal of mean value from points to obtain better
# fitting performance, especially in integer precision. 
points, x_mean, y_mean = b2ac.preprocess.remove_mean_values(points)

# Fit using NumPy methods in double precision.
conic_numpy = b2ac.fit.fit_improved_B2AC_numpy(points)
# Fit using own written methods in double precision.
conic_double = b2ac.fit.fit_improved_B2AC_double(points)
# Fit using own written methods in 64-bit integer precision.
conic_int = b2ac.fit.fit_improved_B2AC_int(points)

# Convert from conic coefficient form to general ellipse form.
general_form_numpy = b2ac.conversion.conic_to_general_1(conic_numpy)
general_form_numpy[0][0] += x_mean
general_form_numpy[0][1] += y_mean

general_form_double = b2ac.conversion.conic_to_general_1(conic_double)
general_form_double[0][0] += x_mean
general_form_double[0][1] += y_mean

general_form_int = b2ac.conversion.conic_to_general_int(conic_int)
general_form_int[0][0] += x_mean
general_form_int[0][1] += y_mean

```

> The mathematical notation described below follows the one used in \[2\].

### Solution strategy

First, the scatter matrix **S** is calculated according to regular specifications,
and the **M** matrix as well. The inverse of **S<sub>3</sub>** is calculated through 
definition of [3x3 matrix inverse](http://mathworld.wolfram.com/MatrixInverse.html).

A QR algorithm with largest value shift, using Givens rotations for both 
tridiagonalization and actual QR steps, (See \[3\]) is then used
to find eigenvalues of the matrix **M**. With these eigenvalues, we
apply an inverse iteration method for calculating the 
corresponding eigenvectors.

The algorithms returns the conic coefficients defining one fitted ellipse.
These can then be transformed to general form: center point, 
axes lengths and rotation angle.

#### Special considerations for integer version

> The integer version uses 64 bit integers!

The calculations of **S** and **M** has special handling. The fact that
this method uses a lot of squared values can easily lead to overflow, especially if
32 bit integers are used. 

The major thing to remember here is that when calculating inverses,
the division of the determinant is postponed, thus using the determinant as
a scale factor during calculations. This scale factor can then be removed 
first when the sought eigenvector has been found. 

#### Limitations for fixed point implementation

For integer versions, it is of great importance to first remove the 
mean values from the points, thus making all values much smaller. This will
drastically improve the fitting results.

Another important aspect is the number of points used for estimating an
ellipse. Using `>100` points, the matrix values in the estimation also grows
and the possibility of overflow is increased. Subsampling of points is recommended.


## References

1.  [Fitzgibbon, A., Pilu, M., & Fisher, R.B.,
    "Direct least square fitting of ellipses," Pattern Analysis and Machine Intelligence,
    IEEE Transactions on , vol.21, no.5, pp.476-480, May 1999]
    (http://research.microsoft.com/pubs/67845/ellipse-pami.pdf)
    
2.  [Halir, R., & Flusser, J., "Numerically stable direct least squares fitting of ellipses." 
     Proc. 6th International Conference in Central Europe on Computer Graphics and Visualization. 
     WSCG. Vol. 98. 1998](http://autotrace.sourceforge.net/WSCG98.pdf)

3.  Golub, G.H. & Van Loan, C.F., 2012. Matrix Computations (4th Ed.). 
    Johns Hopkins University Press, Baltimore, MD, USA.
