# py_pcha
**Fast Python implementation of Archetypal Analysis using Principle Convex Hull Analysis (PCHA).**

*From the source article* [[1][0]]:

"Archetypal analysis (AA) proposed by Cutler and Breiman (1994) [[2][1]] estimates the principal convex hull (PCH) 
of a data set. As such AA favors features that constitute representative ‘corners’ of the data, i.e., distinct 
aspects or **archetypes**."

All code contained in this package was originally written in Matlab. The Matlab package is available [here][2].
Matlab package also handles sparse- and kernel matrices. 

Matlab implementation by: Morten Mørup.
Python implementation by: Ulf Aslak Jensen.

## Install:

```
pip install py_pcha
# or
easy_install py_pcha
```

## Example use:

```python
import numpy as np
from py_pcha.PCHA import PCHA

dimensions = 15
examples = 100
X = np.random.random((dimensions, examples))

XC, S, C, SSE, varexpl = PCHA(X, noc=3, delta=0.1)

print "   # Arc 1     # Arc 2     # Arc 3\n", XC
   # Arc 1     # Arc 2     # Arc 3
[[ 0.32588061  0.3940908   0.71705364]
 [ 0.69790165  0.50729565  0.34076419]
 [ 0.79184963  0.43616783  0.22377323]
 [ 0.36865992  0.51199461  0.68595464]
 [ 0.55887694  0.46533484  0.54946409]
 [ 0.29774011  0.90728239  0.26895903]
 [ 0.33116078  0.87118458  0.26744578]
 [ 0.65678325  0.3104401   0.56770064]
 [ 0.37132093  0.32720999  0.76015795]
 [ 0.31707091  0.44002078  0.81080826]
 [ 0.87002607  0.24002814  0.40317367]
 [ 0.33147574  0.48692694  0.72084014]
 [ 0.2591176   0.81004636  0.34852488]
 [ 0.79427686  0.49692525  0.28712657]
 [ 0.39198509  0.50703908  0.67609915]]
```
	
**Notice:** PCHA takes a 2D-array of shape (dimensions, examples). The same shape applies to any output from the function.
Therefore, the archetypes contained in returned matrix `XC` will be the column vectors. 



[0]: https://scholar.google.com/citations?view_op=view_citation&hl=en&user=RQonsgMAAAAJ&citation_for_view=RQonsgMAAAAJ:J_g5lzvAfSwC
[1]: https://scholar.google.co.il/citations?view_op=view_citation&hl=en&user=9x63d4gAAAAJ&citation_for_view=9x63d4gAAAAJ:u-x6o8ySG0sC
[2]: http://www.mortenmorup.dk/MMhomepageUpdated_files/Page327.htm
