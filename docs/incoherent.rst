.. _Incoherent:

Incoherent Interactions
=======================

Three common sources of incoherent dynamics in systems featuring ODMR are:

#. Incoherent optical excitation and decay of electronic transitions.
#. Dissipative interaction with a quantum mechanical
   environment (bath).
#. Stochastic modulations of (classical) system parameters (e.g. rotational diffusion in a liquid environment, flow in a field gradient, static field drifts).   

Stochastic modulations can be included with our :ref:`Fokker-Planck <FokkerPlanck>` module. A description of SimOS functionality to incorporate optical excitation/decay and empirical
spin relaxation is described in the following.

.. figure:: /img/interactions_incoherentscheme.png
   :width: 500px
   :align: center

|

Optical Transitions
-------------------

Optical excitation and decay events are characterized by collapse operators
of the type :math:`\ket{m}\bra{n}` for pairs of electronic levels :math:`m, n` and classical transition rates 
that are readily available for many systems.

SimOS allows to incorporate incoherent level transitions in three steps:

* Specify all transitions and their rates in a dictionary.
* Utilize the :code:`tidyup_ratedict()` method to tidy up the dictionary; i.e. to ensure that the rates were entered in a correct format. 
* Utilize the :code:`transition_operators()` method to generate the collapse operators.


The rate dictionary is set up in the following way:

* The keys of the dictionary are strings that specify between which levels a transition occurs The "source" and "sink" levels, i.e. the starting and the end points of the transition, are single or multiple member names of the system, separated by commatas. If an excitation only occurs from or to a specific sublevel of a member, the sublevel is specified in square brackets behind the members' name. The direction of the transition is indicated with arrows between the source and sink names. 
* The values of the keys specify the transition rates. 

Let us showcase this flow and how to set up a rate dictionary with a simple example - the Nitrogen-Vacancy center in diamond, which is introduced extensively in a :ref:`separate section <NV>`. We can construct a minimal model 
of a negatively charged NV center as:

.. code-block:: python

   GS  = {"name": "GS", "val": 0}
   ES  = {"name": "ES", "val": 0}
   SS  = {"name": "SS", "val": 0}
   S  = {"name": "S", "val": 1}
   NV = sos.System(([(GS, ES), S], SS))


We now want to construct the collapse operators for the laser excitation, i.e. off-resonant excitation from the electronic
ground to the electronic excited state of the NV center. Note that the excitation is spin conserving, i.e.
the spin state is preserved during the optical excitation. We hence must introduce separate
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


Spin Relaxation
---------------

Both, the interaction with a quantum mechanical bath as well as stochastic modulation of (classical)
system parameter, induces a relaxation of coherences and non-equilibrium populations of the quantum system.

The construction of suitable collapse operators is usually non-trivial and requiires  
detailed knowledge about the time correlation functions of the bath and the spectral density of the system. 
The underlying theory was originally developed for the case of a true quantum mechanical bath 
and later adapted for the semi-classical case. 

SimOS provides a method :code:`relaxation_operators` that generates collapse operators
for longitudinal and transverse spin relaxation based on empirical relaxation rates in the Lindblad formalism. 
If relaxation is induced by stoachstic modulation of (classical) system parameters, you may find our 
:ref:`Fokker-Planck <FokkerPlanck>` module helpful which allows to explicitly model these dynamics instead of approximating them
with a set of collapse operators.



Syntax Reference
----------------

.. py:currentmodule:: simos.incoherent
.. automodule:: simos.incoherent
   :members: tidyup_ratedict, transition_operators, relaxation_operators
   :undoc-members:
