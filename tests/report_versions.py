import numpy as np
import qutip as qu
import sympy
import scipy

print("numpy version:", np.__version__)
print("scipy version:", scipy.__version__) 

try:
    print("sympy version:", sympy.__version__)
except:
    print("sympy not installed")

try:
    print("qutip version:", qu.__version__)
except:
    print("qutip not installed")
