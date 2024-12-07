.. _System:

System Initialisation
=====================

The core functionality of this library is to provide an easy way to construct arbitrarily complicated systems of spins and electronic levels. 
Below you find a hands-on introduction on how to initialise an open quantum system in python with SimOS.

.. note:
   If you do not want to formally initialise a system and still want to work with common 
   spin operators, you might use :ref:`this <BuildYourOwn>` set of functions for construction of spin operators, 
   their tensor products and direct sum as well as further utilities. 


System Construction
-------------------

Arbitrarily complicated systems of spins and electronic levels can be constructed using only two mathematical operations, (i) a tensor product and (ii) a direct sum, to combine the Hilbert spaces of the individual system components. If a tensor product is used, the overall size of the composite Hilbert space is the product of the dimensions of the individual components. If a direct sum is used, the overall Hilbert space size results as the sum of the individual dimensions. 

.. figure:: /img/core_basicoperations.png
   :width: 400px
   :align: center

   Composite Hilbert spaces of spins or electronic levels can be combined via tensor products (left) and direct sums (right). 

The central part of this library is the class :class:`System` which is used to define the system structure and which holds a complete set of operators after construction.
The class constructor is called with a :code:`systemarray` that is essentially a recipe on how to build the quantum system. It specifies all members of the system and how 
their individual Hilbert spaces are combined.


Each member of the spin system is defined as a python dictionary with keys

#. **name**: The name of member:
#. **val**:  Multiplicity of the level, i.e. val = 0.5 for spin 1/2, val = 1 for spin 1 and so on or val = 0.0 for a single electronic level. 

These dictionaries are then combined in the :code:`systemarray` as a series of nested lists and tuples. Lists indicate that Hilbert spaces are combined using tensor products 
while tuples indicate combination with a direct sum. The system class instance (typically in our examples denoted as :code:`s`) holds identity operators 
for all members of the quantum system (e.g., :code:`s.Aid` for a member with name :code:`A`) as well as  x, y, z, lowering, raising and projection operators for all spins (e.g.,
:code:`s.Ax, s.Ay, s.Az, s.Aplus, s.Aminus, s.Ap[0.5]` and :code:`s.Ap[-0.5]` for a spin 1/2 with name :code:`A`). 

.. warning::
   Names of system members must only contain alphanumerical letters.
   Special characters are reserved and used internally to distinguish native from non-native system members.

.. note::
   Additional keys can be used to specify further parameters of system members, such as the spin type, isotope, relaxation times, spatial positions etc. Some functions of SimOS
   require that these parameters are specified.
   Keys can be read, added, deleted or modified after system construction with the :class:`get_properties`, :class:`set_properties` and :class:`del_properties` functions of the system. 


As a first example, we generate a system of two spins S=1/2 (e.g. two electron spins),
whose composite Hilbert space size is 2x2 = 4. 

.. figure:: /img/core_2spinsystem.png
   :width: 450px
   :align: center

   Constructing a system of two coupled spins. 

.. code-block:: python

   import simos as sos
   S = {'name': 'S', 'val': 0.5}
   I = {'name': 'I', 'val': 0.5}
   system_array = [S,I]
   s = sos.System(system_array)

The spin operators in the joint Hilbert space are stored as class attributes,
e.g. the z-operators of the spins are obtained as :code:`s.Sz` and :code:`s.Iz`.
Another very simple example is a pair of electronic levels, 
e.g. an electronic ground and an electronic excited state, with no spin.
The composite Hilbert space size is 1+1=2.

.. figure:: /img/core_2levelsystem.png
   :width: 450px
   :align: center

   Constructing a system of two electronic levels. 

.. code-block:: python

   GS= {'name': 'GS', 'val': 0}
   ES = {'name': 'ES', 'val': 0}
   system_array = (GS,ES)
   s = sos.System(system_array)

The identity operators of the levels in the combined Hilbert space are obtained as
:code:`s.GSid` and :code:`s.ESid`. More complicated examples, e.g. the NV center in diamond and 
spin-correlated radical pairs are showcased in our examples section.



Basis Transformations
---------------------
Upon system construction, the spin operators are initialized
in the Zeeman basis, spanned by the magnetic quantum numbers of the individual spin
members. To include alternative basis sets and transform operators between various bases,
the :class:`System` class provides users with two specific methods.

Coupled Product States of Spins
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Coupled product states of pairs or groups of spins are useful representations for systems 
in which spin-spin interactions dominate the system Hamiltonian.
In SimOS, they can be constructed with the :code:`add_ghostspin` method of the :class:`System`. 
The simplest example, the coupled representation of two coupled spin 1/2 particles as a singlet ( :math:`S_{tot} = \frac{1}{2}-\frac{1}{2} = 0` ) and
triplet ( :math:`S_{tot} = \frac{1}{2} + \frac{1}{2} = 1` ) is illustrated 
below. Here, spins :code:`S` and :code:`I` are coupled to obtain the singlet and triplet  
ghostspins  :code:`C\_1` and :code:`C\_3`.

.. figure:: /img/core_couplespins.png
   :width: 450px
   :align: center

.. code-block:: python

   S  = {"name": "S", "val": 1/2}
   I  = {"name": "I", "val": 1/2}
   s = sos.System([S,I])
   s.add_ghostspin("C", ["S", "I"])

The operators of the ghostspins are generated in full analogy to the native spins 
of the system. For example, the projector onto the the :math:`m_S=0` level of the triplet
is obtained as :code:`C\_3p[0]`. The matrix representations of these operators are
still formulated in the Zeeman basis. However, the :code:`add_ghostspin` method also 
constructs the transformation matrices to transform operators to or from the coupled basis and 
stores them as attributes of the :code:`System`. In our simple example, the transformation matrices 
can be assessed as :code:`s.toC` and :code:`s.fromC` and allow to transform any system operator to or from the Zeeman 
to the coupled Singlet-Triplet basis.

.. code-block:: python

   # Transform z Operators of Singlet and Triplet from
   # Zeeman to  Singlet-Triplet Basis
   C1z_st = s.C_1z.transform(s.toC) 
   C3z_st = s.C_3z.transform(s.toC)

The spin operators of the ghost-spins are named according to the following conventions:

#. The name that the user has specified, e.g. in our case C.
#. An underscore.
#. An integer or a series of integers specifying the spin multiplicity.
   If more than two spins are coupled, the same multiplicity occurs multiple times for the coupled spin. 
   To distinguish them, the multiplicities of their coupling history are included in their name.
#. The actual name of the operator, e.g. z for a z operator, id for the identity.

If :code:`add_ghostspin` is called with the option :code:`return_names = True` the names of the ghost spins (i.e. parts 1-3 of the above list) are returned.
This can be helpful if larger groups of spins are coupled. 


Generic Basis Transformations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


A completely user-defined basis may be defined using a second method, :code:`add_basis`, 
by providing (i) a transformation matrix from the Zeeman basis to the new basis and (ii) a name 
for the basis and a list of names for all basis states. Here, we illustrate this method for the example
of a pair of electronic levels. The method stores identity operators of all new
basis states as well as transformation matrices for back and forth conversion between the Zeeman and the alternate basis. 

.. figure:: /img/core_couplelevels.png
   :width: 450px
   :align: center

.. code-block:: python

   GS= {'name': 'GS', 'val': 0}
   ES = {'name': 'ES', 'val': 0}
   system_array = (GS,ES)
   s = sos.System(system_array)
   T= _np.array([[0.5, 0.5], [0.5, -0.5]])
   s.add_basis(T,  "A", ["g", "e"])

   # Transform one of these states to the new basis.
   gid = s.A_gid.transform(s.toA)

The nomenclature of the alternate basis levels has the following convention:

#. The name that the user has specified for the alternate basis. 
#. An underscore.
#. The name that the user has specified for the new basis state.



Subspaces
---------

The extraction of sub-parts of a multipartite quantum system is an essential 
task during simulation routines and is supported in SimOS. For composite systems whose
combined Hilbert space was spanned with a tensor product,
the subsystem is extracted with a partial trace.
If a composite system was obtained with a direct sum, subsystems are directly obtained as subsets of the full system operator
by the means of projection.
For more complicated multipartite systems, our algorithm extracts the desired subsystem with up to three steps in a top-down approach.
This routine does not remove members unless they are fully separable from the desired subsystem. 

Let us first consider the simple system of two spins S=1/2 that we have already encountered before. The combined Hilbert space is separable 
into the individual Hilbert spaces since it was constructed with a tensor product.

.. figure:: /img/core_extractptrace.png
   :width: 450px
   :align: center

.. code-block:: python

   op = s.Sz 
   op_subS, op_subI, reverse = sos.subsystem(s, op, "S")


Next, we consider our system of two electronic levels. Here the extraction is performed with a projection
operation since the two members were combined with a direct sum.

.. figure:: /img/core_extractproject.png
   :width: 450px
   :align: center

.. code-block:: python

   op = s.GSid 
   op_subGS, op_subES, reverse = sos.subsystem(s, op, "A")

In addition to the operator in the desired subspace (e.g. :code:`op_subS` and :code:`op_subGS`), 
the :code:`subsystem` routine also returns the subspace for the remaining system (e.g. :code:`op_subI` and :code:`op_subES`).
If the extraction was performed in multiple steps, a list of operators is returned for the remaining system. Each entry of the list
is the subspace which was traced out in a single step of the exctraction procedure.
Further, the routine returns a third argument which is essentially an instruction that can be used
to again construct the full system from density matrices of the individual components. Users may
utilize the  :code:`reverse_subsystem` for this purpose. 



Syntax Reference
----------------

.. py:currentmodule:: simos.core
.. automodule:: simos.core
   :members: System, spinops, subsystem, reverse_subsystem
   :undoc-members:

