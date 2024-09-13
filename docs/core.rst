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

Arbitrarily complicated systems of spins and electronic levels can be constructed using only two mathematical operations, (i) a tensor product and (ii) a direct sum, to combine the Hilbert spaces of the indidual system components. If a tensor product is used, the overall size of the composite Hilbert space is the product of the dimensions of the individual components. If a direct sum is used, the overall Hilbert space size results as the sum of the individual dimensions. 

.. figure:: /img/core_basicoperations.png
   :width: 400px
   :align: center

   Composite Hilbert spaces of spins or electronic levels can be combined via tensor products (left) and direct sums (right). 

The central part of this library is the class :class:`System` which is used to define the system structure and which holds a complete set of operators after construction.
The class constructor is called with a :code:`systemarray` that is essentially a recipe on how to build the quantum system. It specifies all members of the system and how 
their individual Hilbert spaces are combined.


Each member of the spin system is defined as a python dictionary with keys

#. 'name': The name of member:
#. 'val': Multiplicity of the level, i.e. val = 0.5 for spin 1/2, val = 1 for spin 1 and so on or val = 0.0 for a single electronic level. 

These dictionaries are then combined in the :code:`systemarray` as a series of nested lists and tuples. Lists indicate that Hilbert spaces are combined using tensor products 
while tuples indicate combination with a direct sum.

.. warning::
   Names of system members must only contain alphanumerical letters.
   Special characters are reserved and used internally to distinguish native from non-native system members.

.. note::
   Additional keys can be used to specify further parameters of system members, such as the spin type, isotope, relaxation times, spatial positions etc. Some functions of SimOS
   require that these parameters are specified.
   Keys can be read, added, deleted or modified after system construction with the :class:`get_properties`, :class:`set_properties` and :class:`del_properties` functions of the system. 


Two Spin 1/2
^^^^^^^^^^^^

As a first example, we generate a system of two spins S=1/2 (e.g. two electron spins),
whose composite Hilbert space size is 2x2 = 4. 

.. code-block:: python

   import simos as sos
   A = {'name': 'A', 'val': 0.5}
   B = {'name': 'B', 'val': 0.5}
   system_array = [A,B]
   s = sos.System(system_array)

The spin operators in the joint Hilbert space are stored as class attributes,
e.g. the z-operators of the spins are obtained as :code:`s.Az` and :code:`s.Bz`.


Pair of Electronic Ground and Excited State
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Another very simple example is a pair of electronic levels, 
e.g. an electronic ground and an electronic excited state, with no spin.
The composite Hilbert space size is 1+1=2.

.. code-block:: python

   GS= {'name': 'GS', 'val': 0}
   ES = {'name': 'ES', 'val': 0}
   system_array = [GS,ES]
   s = sos.System(system_array)

The NV Center in Diamond
^^^^^^^^^^^^^^^^^^^^^^^^

The NV center in diamond is a prototypical example for ODMR experiments and a promising platform for quantum science and information applications.
A simplified model for a negatively charged NV center encompasses three electronic levels, an electron spin S=1 and a nuclear spin of the nitrogen S=1/2 (for synthetic NV centers). 
We will construct the Hilbert space of this system in three steps:

We start by only considering the electronic levels of the NV centers negative charge state, i.e. a pair of electronic ground (GS) and excited states (ES) as well as a metastable shelving state (SS). 
The combined Hilbert state is simply the direct sum of the individual, one-dimensional Hilbert spaces.

.. code-block:: python

   GS  = {"name": "GS", "val": 0}
   ES  = {"name": "ES", "val": 0}
   SS  = {"name": "SS", "val": 0}
   NV = sos.System((GS, ES, SS))

.. figure:: /img/core_NVstep1.png
   :width: 260px
   :align: center

In a next step  we add the electron spin of the NV center, which results in a fine structure for the GS and ES electronic states.

.. code-block:: python

   S  = {"name": "S", "val": 1}
   NVfine = sos.System(([(GS, ES), S], SS))

.. figure:: /img/core_NVstep2.png
   :width: 475px
   :align: center

Further we include the hyperfine structure induced by coupling to the NV centers nuclear spin.

.. code-block:: python

   N  = {"name": "N", "val": 1/2}
   NVhf = sos.System([([(GS, ES), S], SS), N])

.. figure:: /img/core_NVstep3.png
   :width: 800px
   :align: center


The class instance :code:`NVhf` now holds identity operators for all members (e.g. :code:`NVhf.GSid`) as well as x, y, z,
lowering and raising spin operators for the electron and the nitrogens nuclear spin (e.g. :code:`NVhf.Sz`). Further, projection operators are available for all spin sublevels (e.g.
:code:`NVhf.Sp[0]` for projection onto the mS = 0 state of the NV centers electron spin). 


Basis Transformations
---------------------

Coupled Product States of Spins
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Upon system construction, the spin operators are initialized in the Zeeman basis, spanned by the magnetic quantum numbers of the individual spin members.
However, :class:`System` provides a method :code:`add_ghostspin` for construction of coupled product states of pairs or groups of spins. 
The simplest example is a pair of two strongly coupled spin S=1/2 (for example, a spin-correlated radical pair), that are conveniently described in a singlet-triplet basis. 

.. figure:: /img/core_st.png
   :width: 350px
   :align: center

In SimOS, we can easily access the singlet and triplet states by calling the :code:`add_ghostspin` method after system initialisation.

.. code-block:: python

   A  = {"name": "A", "val": 1/2}
   B  = {"name": "B", "val": 1/2}
   s = sos.System([A,B])
   s.add_ghostspin("C", ["A", "B"])

The system now also contains spin operators of the so-called ghost-spin C, for example the
z spin operators of the Singlet (S=1) and Triplet (S=3) states can now be accessed as :code:`s.C_1z` and :code:`s.C_3z`.
In addition, a transformation matrix for the spin coupling has been created ( :code:`s.toC`) allowing for basis transformations of operators from the zeeman basis to the singlet-triplet basis.

.. code-block:: python
   # Transform z Operators of Singlet and Triplet in Singlet-Triplet Basis
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

.. note::
   The add-basis functionality is covered in our example notebooks on photophysics of NV centers (quantum-mechanical model). Have a look there if you are unsure on
   how to include this functionality in your code.


Generic Basis Transformations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The  :code:`add_basis` method of the :class:`System` class further allows to specify 
alternate bases for arbitrary linear combinations of individual system levels.
In this case, the identity operators of the new basis states are generated and again, 
the transformation matrices for basis interconversion are generated. 

The nomenclature of the alternate basis levels has the following convention:

#. The name that the user has specified for the alternate basis. 
#. An underscore.
#. The name that the user has specified for the new basis state.

.. note::
   The add-basis functionality is covered in our example notebooks on photophysics of NV centers (quantum-mechanical model). Have a look there if you are unsure on
   how to include this functionality in your code.


Subspaces
---------

Sometimes it can be very useful to extract subsystems of a full system.
If the composite Hilbert space of a quantum system was constructed with tensor products,
the individual subsystems can be separated with a partial trace.
If the subsystems were combined with a direct sum, no well-defined mathematical operation exists
to separate them, however one can simply perform a projection operation onto the subsystem and drop all other
dimensions of the Hilbert space.

The SimOS library offers a function :code:`subsystem` that extracts subsystems with a partial trace 
whenever possible and performs the projection-exctration in all other cases.
We will demonstrate the different cases and usage of :code:`subsystem`  with two simple examples.

Let us first consider the simple system of two spins S=1/2 that we have already encountered before. The combined Hilbert space is separable into the individual Hilbert spaces since it was constructed with a tensor product.
In the code below we take the z-operator of spin A in the composite Hilbert space and extract the subspace of spin A. 
Subspace requires a 

.. code-block:: python

   op = s.Az 
   op_subA, _ = sos.subsystem(s, op, "A")


Again, we turn to our model NV center to demonstrate an example which is not separable into subsystems. Here, our combined Hilbert space was not only constructed with tensor products since the shelving state
was added with a direct sum. We can still use subsystem to extract the subsystem of the
excited electronic state. Here, we take the z-operator of the electronic spin of the NV center and project it onto the subspace
of the excited electronic state.

.. code-block:: python

   op = NVhf.Sz
   op_subES, _ = sos.subsystem(NVhf, op, "ES")

In addition to the operator in the Hilbert subspace, the :code:`subsystem` routine also returns a dictionary which specifies how to project an operator again from the subspace onto the full system space 
if used as an input for the  :code:`reverse_subsystem` method.  


.. note::
   The subsystem functionality is utilized in some of our example notebooks, for example the introdcution to SCRP spin dynamics. Have a look there if you are unsure on how to
   incorporate this method into your code.


Syntax Reference
----------------

.. py:currentmodule:: simos.core
.. automodule:: simos.core
   :members: System, spinops, subsystem, reverse_subsystem
   :undoc-members:

