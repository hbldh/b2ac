# B2AC - Ellipse fitting in Python and C

Python and C implementations of an ellipse fitting algorithm in double and fixed point precision. 

The implementations here were never meant for production, but were written as an examination of how much
speedup could be attained from implementing ellipse fitting methods in fixed point and at what precision cost.

The ellipse fitting method in (2) (an improved version of the one in (1)) was implemented.

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

```python
import b2ac.ellipse

ellipse = b2ac.ellipse.def fit_improved_B2AC_double(points)

```

## References

1.  [Fitzgibbon, A.; Pilu, M.; Fisher, R.B.,
    "Direct least square fitting of ellipses," Pattern Analysis and Machine Intelligence,
    IEEE Transactions on , vol.21, no.5, pp.476-480, May 1999]
    (http://research.microsoft.com/pubs/67845/ellipse-pami.pdf)
    
2.  [Hal?r, R., & Flusser, J.; "Numerically stable direct least squares fitting of ellipses." 
     Proc. 6th International Conference in Central Europe on Computer Graphics and Visualization. 
     WSCG. Vol. 98. 1998](http://autotrace.sourceforge.net/WSCG98.pdf)

3.  Golub, G.H. & Van Loan, C.F.; 2012. Matrix Computations (4th Ed.). 
    Johns Hopkins University Press, Baltimore, MD, USA.