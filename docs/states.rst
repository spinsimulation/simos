.. _states:

States
======

States of (ensembles) of closed and open quantum systems can be expressed in various formalisms 
and spaces (e.g. Hilbert, Liouville or higher dimensional).
In SimOS we support four distinct state representations as summarized in the table below.
In analogy to Qutips Quantum Object, each quantum object in SimOS has a :code:`dims` attribute,  which indicates
the dimensionality and separable substructure of the associated vector space. 
Here, we illustrate the quantum dimensionality of the state objects for a composite two-member system
with individual Hilbert space dimensions :math:`n` and :math:`m` . For the Fokker-Planck space state representation,
:math:`x` indicates the dimensionality of the classical subspace. 

.. list-table:: States
   :widths: 10 10 10 30
   :header-rows: 1

   * - State Type 
     - Structured Dimensionality (2 Spin System Example)
     - Shape
     - Remarks
   * - Ket :math:`\ket{\psi}`
     - :code:`[[n, m], [1,1]]` 
     - :math:`(1, n \cdot m)`
     - A pure state in Hilbert space (state vector, unit vector in Hilbert space). 
   * - Density Matrix :math:`\rho`
     - :code:`[[n,m], [n,m]]` 
     - :math:`(n \cdot m , n \cdot m)`
     - A pure or mixed state in Hilbert space.
   * - Vectorized Density Matrix :math:`\vec{\rho}`
     - :code:`[ [[n,m], [n,m]], [1]]` 
     - :math:`( (n \cdot m)^2 , 1)`
     - Density matrix liouville space representation, allowing for direct multiplication with a superoperator.
   * - Fokker Planck Vectorized Density Matrix :math:`\vec{\rho}_{\mathrm{FP}}`
     - :code:`[ [x, (n \cdot m)^2], [1]]`
     - :math:`( (n \cdot m)^2 \cdot x , 1)`
     - Density matrix fokker planck space representation.  
.. list-table::
   
|

The figure below further illustrates the characteristics of the Hilbert and Liouville
space quantum state representations and their interconversion.
To transform kets and density matrices, users can utilize the :code:`ket2dm()` and :code:`dm2ket()` methods. 
Density matrices can be interconverted with their vecotrized Liouville space representations with the
:code:`dm2vec()` and :code:`vec2dm()` methods.
Although not shown here explicity, conversion to the Fokker-Planck space representation is further
supported in the :code:`dm2fp()` and :code:`fp2dm` methods.  

.. figure:: /img/states_schematic.png
   :width: 650px
   :align: center


|
Initialising States
-------------------

SimOS provides a series of methods to facilitate initialisation of states as documented below.

Pure States
^^^^^^^^^^^

To facilitate the initialisation of pure states of multipartite quantum systems, SimOS provides 
the :code:`state()` method. The scheme below illustrates initialisation of a pure state
for the example of a simple 2-spin system. The desired state is defined as a string which contains all
member names and, in the case of spins, magnetic sublevels in square parentheses.  

.. figure:: /img/states_example.png
   :width: 650px
   :align: center

|

The obtained state vectors can succesively be transformed to density matrices with the :code:`ket2dm()` method; mixed states can be prepared as weighted superpositions of "pure" 
density matrices.

.. warning::
  We strongly recommend to explicitly specify the spin-state of all spin-members (i.e. members with :code:`val` :math:`> 1/2`) when creating an initial state with the :code:`state()` method. 
  If the spin-state of a relevant spin member (i.e. a spin that is coupled to one or multiple of the selected electronic levels) is not explicitly specified, the created state is a coherent superposition in the Zeeman basis.
  The :code:`state()` routine does currently not  permit omitting spins (for performance reasons),  however, this functionality might be deprecated in future releases.  

Thermal and Polarized States
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Magnetic resonance experiments often assume an initial state with thermal or otherwise well-defined non-equilibrium polarizations of individual
system members. SimOS provides a method :code:`thermal_dm()` that returns the thermal density matrix for a given temperature and system Hamiltonian. The method :code:`pol_spin()` returns the 
density matrix for a spin :code:`S= 1/2` with a well-defined polarization, i.e.

.. math::
  \begin{pmatrix}
  \frac{1+p}{2} & 0 \\
  0 & \frac{1-p}{2}
  \end{pmatrix}.



Syntax Reference
----------------


.. py:currentmodule:: simos.states
.. automodule:: simos.states
   :members: state, thermal_state, pol_spin, state_product