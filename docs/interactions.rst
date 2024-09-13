
.. _Coherent:

Coherent Interactions
=====================

The coherent dynamics of spin ensembles are mostly governed by interactions of spins with magnetic and electric fields
and pairwise spin-spin interactions. The following functions may be used to generate the respective Hamiltonians.

Of course you do not have to utilize these functions to introduce coherent interactions - especially for simple,
isotropic interactions it is very easy to just define them "by hand" using the operators of your spin system.
For example, a Zeeman interaction of a spin can simply be introduced as

.. code-block:: python

   S  = {"name": "S", "val": 1/2}
   s = sos.System([S])
   Hz = sos.ye * s.Sz

This way any interaction that is not covered by our helper-functions can also be introduced (for example, a zero-field-splitting).

.. note::
   SimOS provides you with a set of :ref:`physical constants <Constants>` as well as a 
   complete set of :ref:`gyromagnetic ratios <GyroRatios>`.


Spin-Field Interactions
-----------------------

The interaction of nuclear or electronic spins with surrounding magnetic or electric fields is described by a Hamiltonian of the form  

.. math::
   \hat{H} =  \vec{B}^T \cdot \mathbf{\gamma} \cdot \vec{\hat{S}}

where :math:`\vec{B}^T` is the  field vector and :math:`\mathbf{\gamma}` is a tensor that specifies the potentially anisotropic coupling of the spin magnetic dipole (or: electric quarupole, in the case of electric fields) moment to the magnetic (electric) field.
For isotropic spin-field interactions, this simplifies to 

.. math::
   \hat{H} =  \gamma \vec{B}^T \cdot \vec{\hat{S}}

with a scalar coupling constant :math:`\gamma`. SimOS provides two methods for spin-field interactions:
The method :code:`field_interaction` allows to specify generic spin-field interactions for individual spins.
The :code:`auto_zeeman_interaction` automatically calculates the isotropic Zeeman interactions for an entire spin system.


.. py:currentmodule:: simos.coherent
.. automodule:: simos.coherent
   :members: field_interaction
   :undoc-members:

.. py:currentmodule:: simos.coherent
.. automodule:: simos.coherent
   :members: auto_zeeman_interaction
   :undoc-members:


Spin-Spin Interactions
-----------------------

Any interaction between a pair of spins is of the form  

.. math::
   \hat{H} =  \vec{\hat{S}_1}^T \cdot \mathbf{A} \cdot \vec{\hat{S}_2}

where :math:`\mathbf{A}` is a tensor that specifies the strength of the potentially anisotropic coupling. 

Generic Spin-Coupling
^^^^^^^^^^^^^^^^^^^^^

SimOS provides a method :code:`spinspin_coupling` for generic spin-spin interactions. 

.. py:currentmodule:: simos.coherent
.. automodule:: simos.coherent
   :members: spinspin_coupling
   :undoc-members:

Dipolar Coupling
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The dipolar coupling between pairs of nuclear or electronic spins is a through-space magnetic interaction of their associated magnetic dipoles. 
The dipolar coupling Hamiltonian

.. math::
   \hat{H}_{dip} = -\frac{\mu_0 \hbar \gamma_1 \gamma_2}{4 \pi}\frac{1}{r^3} \left(  3(\vec{S_1}\cdot \vec{r_u})(\vec{S_2}\cdot \vec{r_u}) - \vec{S_1}\cdot \vec{S_2}\right)   \\ 
   =  \vec{S_1}^T \cdot \mathbf{D} \cdot \vec{S_2}

can be calculated from the gyromagnetic ratios :math:`\gamma_1, \gamma_2` of the spins and the distance
and orientation of their connecting vector :math:`r` with respect to the external magnetic field. The resulting dipolar coupling tensor
:math:`\mathbf{D}` is a 3x3 matrix. Alternatively, the dipolar coupling may be expressed as a "dipolar alphabet" with components:

.. math::
   \mathcal{A} =  -\frac{\mu_0 \hbar \gamma_1 \gamma_2}{4 \pi}\frac{1}{r^3}  (3 \cos{\theta}^2 - 1 ) S_{1,z}S_{2,z}
.. math::   
   \mathcal{B} =  \frac{\mu_0 \hbar \gamma_1 \gamma_2}{4 \pi}\frac{1}{r^3} (3 \cos{\theta}^2 - 1 )  (S_{1,+}S_{2,-}+S_{1,-}S_{2,+})
.. math::   
   \mathcal{C} =  -\frac{\mu_0 \hbar \gamma_1 \gamma_2}{4 \pi}\frac{1}{r^3} \frac{3}{2}\sin{\theta}\cos{\theta} e^{-1i\phi} (S_{1,x}S_{2,z}+S_{1,z}S_{2,x} + i (S_{1,y}S_{2,z}+S_{1,z}S_{2,y}))
.. math::   
   \mathcal{D} =  -\frac{\mu_0 \hbar \gamma_1 \gamma_2}{4 \pi}\frac{1}{r^3}\frac{3}{2}\sin{\theta}\cos{\theta} e^{1i\phi} (S_{1,x}S_{2,z}+S_{1,z}S_{2,x} + i (S_{1,y}S_{2,z}+S_{1,z}S_{2,y})) 
.. math::   
   \mathcal{E} =  -\frac{\mu_0 \hbar \gamma_1 \gamma_2}{4 \pi}\frac{1}{r^3} \frac{3}{4}\sin{\theta}^2e^{-2i\phi}  (S_{1,x}S_{2,x} - S_{1,y}S_{2,y} + i (S_{1,x}S_{2,y}+S_{1,y}S_{2,x})) 
.. math::   
   \mathcal{F} =  -\frac{\mu_0 \hbar \gamma_1 \gamma_2}{4 \pi}\frac{1}{r^3} \frac{3}{4}\sin{\theta}^2e^{2i\phi}  (S_{1,x}S_{2,x} - S_{1,y}S_{2,y} - i (S_{1,x}S_{2,y}+S_{1,y}S_{2,x})) 

SimOS provides two methods for purely dipolar couplings. The :code:`dipolar_spatial` method only returns the spatial part, i.e. the dipolar coupling tensor, of the interaction while 
the :code:`dipolar_coupling` method returns the complete dipolar coupling Hamiltonian. 

.. py:currentmodule:: simos.coherent
.. automodule:: simos.coherent
   :members: dipolar_spatial
   :undoc-members:

.. py:currentmodule:: simos.coherent
.. automodule:: simos.coherent
   :members: dipolar_coupling
   :undoc-members:

.. py:currentmodule:: simos.coherent
.. automodule:: simos.coherent
   :members: auto_pairwise_coupling
   :undoc-members:


.. _Incoherent:

Incoherent Interactions
=======================

Three common sources of incoherent dynamics in systems featuring ODMR are:

#. Incoherent optical excitation and decay of electronic transitions.
#. Dissipative interaction with a quantum mechanical
   environment (bath).
#. Stochastic modulations of (classical) system parameters (e.g. rotational diffusion in a liquid environment, flow in a field gradient, static field drifts).   

Below we describe how to introduce these incoherent dynamics in SimOS.


Optical Transitions
-------------------

Optical excitation and decay events are characterized by collapse operators
of the type :math:`\ket{m}\bra{n}` for pairs of electronic levels :math:`m, n` and classical transition rates 
that are readily available for many systems.

In SimOS you can simply include incoherent level transitions in three steps:

#. Specify all transitions and their rates in a dictionary.
#. Call the :code:`tidyup_ratedict` method for your dictionary. This will ensure that the rates have been 
   entered in the correct format.
#. Generate the collapse operators with the :code:`transition_operators` method.


The rate dictionary is set up in the following way:

#. The keys of the dictionary are strings that specify between which levels a transition occurs. The "source" and "sink" levels, i.e. the starting and the end points of the transition,
   are single or multiple member names of the system, separated by commatas. If an excitation only occurs from or to a specific sublevel of a member, the sublevel is specified in square brackets behind the members' name.
   The direction of the transition is indicated with arrows between the source and sink names. 
#. The values of the keys specify the transition rates. 

Let us showcase this flow and how to set up a rate dictionary with a simple example - the Nitrogen-Vacancy center in diamond. As shown in the
previous section of this documentation, we can construct a minimal model of a negatively charged NV center as:

.. code-block:: python

   GS  = {"name": "GS", "val": 0}
   ES  = {"name": "ES", "val": 0}
   SS  = {"name": "SS", "val": 0}
   S  = {"name": "S", "val": 1}
   NV = sos.System(([(GS, ES), S], SS))

The scheme below shows excitation and decay paths of this minimal model, involving excitation from the ground to the excited state
and succesive radiative and non-radiative decay.  

.. figure:: /img/interactions_NVsimplescheme.png
   :width: 475px
   :align: center

We now want to construct the collapse operators for the laser excitation. We consider off-resonant, i.e. incoherent, excitation from the electronic
ground to the electronic excited state of the NV center. Note that the excitation is spin conserving, i.e.
the spin state is preserved during the optical excitation. We therefore have to introduce separate
excitation processes for all spin sublevels. 

.. code-block:: python

   # Step 1: Prepare dictionary 
   laser_rates = {}
   laser_rates["GS,S[0]->ES,S[0]"] = 10e6
   laser_rates["GS,S[1]->ES,S[1]"] = 10e6
   laser_rates["GS,S[-1]->ES,S[-1]"] = 10e6
   # Step 2: Tidy-up dictionary
   laser_rates = sos.tidyup_ratedict(NV, laser_rates)
   # Step 3: Generate the collapse operators.
   c_ops = sos.transition_operators(NV, laser_rates)

The returned :code:`c_ops` is a list of collapse operators for the three excitation channels that can be used to generate
a Liouvillian superoperator.  The collapse operators for the decay channels can be constructed in an analogous fashion.
You can find the full code for the whole simulation of this example in our examples section.

.. py:currentmodule:: simos.incoherent
.. automodule:: simos.incoherent
   :members: tidyup_ratedict, transition_operators
   :undoc-members:


Relaxation
----------

Both, the interaction with a quantum mechanical bath as well as stochastic modulation of (classical)
system parameters, induces a relaxation of coherences and non-equilibrium populations of the system.

The construction of suitable collapse operators is usually quite copmlicated. 
In both cases, construction of collapse operators requires knowledge about the time
correlation functions of the bath and the spectral density of the system. 
The underlying theory was originally developed for the case of a true quantum mechanical bath 
and later adapted for the semi-classical case. 

SimOS provides a method :code:`relaxation_operators` that generates collapse operators
for longitudinal and transverse spin relaxation based on empirical relaxation rates. 

.. py:currentmodule:: simos.incoherent
.. automodule:: simos.incoherent
   :members: relaxation_operators
   :undoc-members:

If relaxation is induced by stoachstic modulation of (classical) system parameters, you can further
use our Fokker-Planck module to explicitly model these dynamics instead of approximating them
with a set of collapse operators.