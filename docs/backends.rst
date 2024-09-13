.. _Backends:

Backends
========

In order to benefit as much as possible from the vast range of existing python functionality, this simulation
library supports four separate backends. Each backend uses a distinct data type for quantum objects 
(kets, bras, operators, vectorized operators and superoperators) and provides an elementary set of functions that construct and
operate on these objects. Backends are specified with the :code:`method` keyword argument. 

.. warning::
   Although we tried to minimize deviations between backends some small differences may remain. Further,
   some advanced functionality is not supported for the sympy backend. 

Backend Overview
----------------

.. list-table:: Backends
   :widths: 25 25 50
   :header-rows: 1

   * - Backend Name
     - Data Type of Quantum Objects
     - Benefits / Remarks
   * - Numpy [HMW+20]_
     - Custom data type :code:`QNdarray` based on :code:`numpy.ndarray` 
     - Utilizes numpy instead of specific libraries; optimizes computational performance and portability.
   * - Qutip [JNN20]_
     - :code:`qutip.Qobj` 
     - Whenever possible, qutip functionality is being used. Ideal for users that are familiar with qutip, allows to smoothly extend the SimOS capabilities with qutip functionality. 
   * - Scipy Sparse  [VGO+20]_
     - Custom data type :code:`Qcsc_matrix` based on :code:`scipy.sparse.csc_matrix`
     - The main data type for all quantum objects is a scipy sparse matrix, thereby reducing storage space requirements and accelerating calculation for sparse systems.
   * - Sympy [MSP+17]_
     - Cutom data type :code:`QSMatrix` based on  :code:`sympy.matrix`
     - Uses sympy functionality and therefore allows to obtain symbolic expressions. Note that some functionality of the library does not support this backend. In general, the sympy backend  is not suited for extensive time dynamics simulations but a didactic tool allowing to obtain an intuitive understanding of a systems spin physics.
.. list-table::
   

Quantum Object Data Type
------------------------

The backend defines the data type of all quantum objects, i.e. kets, bras, operators , vectorized operators 
and superoperators. For numpy, scipy and sympy, the quantum objects are custom data types. 
They are based on native data types of the underlying python packages that were subclassed to 
enable the most important functionality of qutips :code:`qutip.Qobj`. 
Important features of these quantum object data type classes are:

#. The python product operator :code:`*` is overloaded such that a matrix multiplication is performed instead of elementwise multiplication. 
#. An attribute :code:`Q.dims` stores the structure of the quantum objects Hilbert space.
#. In analogy to qutip, a set of additional functions directly operates on instances of the quantum object class. Note that these functions can also be accessed as :code:`f(Q, *args, **kwargs)` instead of :code:`Q.f(*args, **kwargs)`.

.. list-table:: Functions Operating on Quantum Object Class
   :widths: 25 25 50
   :header-rows: 1

   * - Function
     - Command
     - Description
   * - Dagger
     - :code:`Q.dag()`
     - returns the adjoint of the quantum object, output is again a quantum object
   * - Conjugate
     - :code:`Q.conj()`  
     - returns the conjugate of the quantum object, output is again a quantum object
   * - Transpose
     - :code:`Q.trans()`
     - returns the transpose of the quantum object, output is again a quantum object
   * - Trace
     - :code:`Q.tr()`
     - returns the trace of the quantum object, output is a scalar
   * - Diagonal
     - :code:`Q.diag()`  
     - returns the diagonal of the quantum object, output is a numpy.ndarray 
   * - Unit
     - :code:`Q.unit()`
     - returns the normalized quantum object, output is again a quantum object
   * - Conjugate
     - :code:`Q.transform(U)`  
     - performs a basis transformation defined by the transformation matrix U, output is again a quantum object
   * - Matrix Exponential
     - :code:`Q.expm()`
     - returns the matrix exponential of the quantum object, output is again a quantum object
   * - Partial Trace
     - :code:`Q.ptrace(i)`  
     - performs a partial trace, keeping only the selected ith subsystem 
   * - Eigenstates
     - :code:`Q.eigenstates()`  
     - returns eigenvalues and eigenvectors of the quantum object
.. list-table::


To transform an object to a quantum object of any backend the :code:`tidyup(Q)` function can be used. To extract the data of a quantum object as a standard :code:`numpy.ndarray()` the :code:`data(Q)` function can be applied.

.. py:currentmodule:: simos.spinoperations
.. automodule:: simos.spinoperations
   :members: tidyup
   :undoc-members:
.. py:currentmodule:: simos.spinoperations
.. automodule:: simos.spinoperations
   :members: data
   :undoc-members:

.. _QMMethods:

Manipulating Quantum Objects
----------------------------

A series of functions for construction and basic manipulation of quantum objects is provided 
for all backends. For functions that construct a quantum object from scratch, the backend is specified with the  :code:`method` keyword argument.
For all other functions, the backend is automatically detected via the data type of the input argument and type conversion to a native quantum data type of the backend is performed before the actual operation is executed.



Hilbert Space Methods
^^^^^^^^^^^^^^^^^^^^^

.. py:currentmodule:: simos.qmatrixmethods
.. automodule:: simos.qmatrixmethods
   :members: isbra, isket, isoper, ket2dm, dm2ket, ket, bra, expect
   :undoc-members:


Liouville Space Methods
^^^^^^^^^^^^^^^^^^^^^^^

.. py:currentmodule:: simos.qmatrixmethods
.. automodule:: simos.qmatrixmethods
   :members: operator_to_vector, vector_to_operator, issuper, spost, spre, lindbladian, liouvillian, applySuperoperator
   :undoc-members:


Fokker-Planck Methods
^^^^^^^^^^^^^^^^^^^^^

.. py:currentmodule:: simos.qmatrixmethods
.. automodule:: simos.qmatrixmethods
   :members: focker_to_dm, dm_to_focker 
   :undoc-members:


.. _BuildYourOwn:

Additional Methods for Creating Operators
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. py:currentmodule:: simos.qmatrixmethods
.. automodule:: simos.qmatrixmethods
   :members: tensor, directsum, block_diagonal, jmat, diags, identity
   :undoc-members:




.. rubric:: References

.. [VGO+20] Pauli Virtanen, Ralf Gommers, Travis E. Oliphant, Matt Haberland, Tyler Reddy, David Cournapeau, Evgeni Burovski, Pearu Peterson, Warren Weckesser, Jonathan Bright, Stéfan J. van der Walt, Matthew Brett, Joshua Wilson, K. Jarrod Millman, Nikolay Mayorov, Andrew R. J. Nelson, Eric Jones, Robert Kern, Eric Larson, CJ Carey, İlhan Polat, Yu Feng, Eric W. Moore, Jake VanderPlas, Denis Laxalde, Josef Perktold, Robert Cimrman, Ian Henriksen, E.A. Quintero, Charles R Harris, Anne M. Archibald, Antônio H. Ribeiro, Fabian Pedregosa, Paul van Mulbregt, and SciPy 1.0 Contributors. SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python. Nature Methods, 17(3), 261-272 (2020).

.. [HMW+20] Harris, C.R., Millman, K.J., van der Walt, S.J. et al. Array programming with NumPy. Nature 585, 357–362 (2020).

.. [JNN20] Johansson, J.R., Nation, P.D. and Nori, F. QuTiP : An open-source Python framework for the dynamics of open quantum systems. Comp. Phys. Comm. 183, 1760-1772 (2012).

.. [MSP+17] Meurer A, Smith CP, Paprocki M, Čertík O, Kirpichev SB, Rocklin M, Kumar A, Ivanov S, Moore JK, Singh S, Rathnayake T, Vig S, Granger BE, Muller RP, Bonazzi F, Gupta H, Vats S, Johansson F, Pedregosa F, Curry MJ, Terrel AR, Roučka Š, Saboo A, Fernando I, Kulal S, Cimrman R, Scopatz A. SymPy: symbolic computing in Python. Peer J Computer Science, 3(103) (2017).


