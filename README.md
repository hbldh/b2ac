# B2AC - Ellipse fitting in Python and C

|  Branch       | Build status     |
| :------------ | ---------------: |
| `master`      | [![Build Status](https://travis-ci.org/hbldh/b2ac.svg?branch=master)](https://travis-ci.org/hbldh/b2ac) |
| `develop`     | [![Build Status](https://travis-ci.org/hbldh/b2ac.svg?branch=develop)](https://travis-ci.org/hbldh/b2ac) |

Python and C implementations of an ellipse fitting algorithm in double and fixed point precision. 

The implementations here were never meant for production, but were written as an examination of how much
speedup could be attained from implementing ellipse fitting methods in fixed point and at what precision cost.

The ellipse fitting method in \[2\] (an improved version of the one in (1)) was implemented.

### Solution strategy for double precision

TBD. 

### Solution strategy for fixed point precision

TBD.

#### Limitations for fixed point implementation.

TBD.

## Installation

Install simply by calling:

    pip install https://github.com/hbldh/ellipse-fitting
    
Numpy needs to be installed prior to `b2ac`, since it compiles the C extension using Numpy headers.

## Testing

Tests will be eventually be implemented.

## Usage

An ellipse can be fitted using Python methods:

```python
import numpy as np
import b2ac.fit

points = np.array([[3, 0], [0, 5], [3, 10], [6, 5]])

# Fit using NumPy methods in double precision.
ellipse_numpy = b2ac.fit.fit_improved_B2AC_numpy(points)
# Fit using own written methods in double precision.
ellipse_double = b2ac.fit.fit_improved_B2AC_double(points)
# Fit using own written methods in 64-bit integer precision.
ellipse_int = b2ac.fit.fit_improved_B2AC_int(points)

```

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