Examples
========

In our Virtual Lab, we have put together a series of example jupyter notebooks that provide a hands-on introduction
to SimOS syntax and spin dynamics simulations in general. 
If you are new to spin dynamics simulations, you might find our For-Beginners section useful.
Further examples feature two protoypical systems with spin-dependent fluorescence, the nitrogen-vacancy
center in diamond and spin-correlated radical-pairs. 
Below you find an overview of all examples and links to their jupyter notebooks. 


For Beginners
-------------

Familiarize yourself with spin dynamics simulations using our example notebook `Simple_NMR_experiments.ipynb <https://simos.kherb.io/virtual/lab/index.html?path=examples%2FSimple_NMR_experiments.ipynb>`_. 
It covers the most essential syntax of the SimOS library and introduces basic concepts, e.g. rabi osicllations, inversion recovery or a DEER experiment. 
There is also a symbolic version, `Simple_NMR_experiments_symbolic.ipynb <https://simos.kherb.io/virtual/lab/index.html?path=examples%2FSimple_NMR_experiments_symbolic.ipynb>`_ that showcases the sympy backend
of SimOS and enables a more intuitive understanding of the described experiments.


The NV Center in Diamond
------------------------

Color centers in diamond and silicon carbide have been studied for years and demonstrated significant potential in nanoscale sensing of magnetic and electric fields and as long-lived
quantum memories. The nitrogen vacancy center is the most well studied point defect in diamond and the research focus of the developers of the SimOS library. A series of jupyter notebooks
cover some essential photophysics and ODMR experiments of the NV center.

Photophysics
^^^^^^^^^^^^

The electronic ground state of the negatively charged NV center is a spin triplet (S= 1) whose spin state can be optically prepared and read-out using green light illumination. 
This optical addressability enables detection of single NV centers with high sensitivity by means of optically detected magnetic resonance.

In the example notebook `NV.ipynb <https://simos.kherb.io/virtual/lab/index.html?path=examples%2FNV.ipynb>`_ we explore
the photophysics of the NV center and theoretically assess optical spin-polarization and spin-dependent fluorescence. 
Starting from a pseudo-classical rate model we will gradually evolve to a fully quantum mechanical model that accounts for
effects of magnetic and electric fields as well as temperature.

Nanoscale NMR
^^^^^^^^^^^^^

The NV center can detect the free induction decay of single nuclear spins. 
In the example notebook `NV_weak_measurements.ipynb <https://simos.kherb.io/virtual/lab/index.html?path=examples%2FNV_weak_measurements.ipynb>`_, we theoretically assess tracking the precession
of individual nucelar spins with a single NV center.

Spin-Correlated Radical Pairs
-----------------------------

Photogenerated spin-correlated radical pairs are entangled electron spin pairs formed in well-defined spin states. 
Spin chemistry of SCRPs plays important roles in biology, e.g., in solar energy harvesting by photosynthetic reaction centers, 
and has been implicated in the leading hypothesis for magnetoreception via the radical pair mechanism  in blue-light sensitive cryptochromes. 
Synthetic molecules mimicking these biologically relevant species have also garnered recent interest for chemical approaches to quantum sensing, e.g.of magnetic fields,
and quantum information science. 
The following jupyter notebooks provide an introduction to the spin dynamics of radical pair systems.


Tutorial on SCRPs: Initial States, Hamiltonians and Basic Dynamics with Sympy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the example notebook 
`SCRP_basics.ipynb <https://simos.kherb.io/virtual/lab/index.html?path=examples%2FSCRP_basics.ipynb>`_. 
, we explore the basic spin dynamics of spin-correlated radical pairs in a symbolic framework utilizing the
sympy backend of SimOS.  


Tutorial on Magnetic Field Effects 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Photogenerated spin-correlated radical pairs (SCRPs) typically occur as short lived intermediates 
in cascades of biochemical transformations.
Although these reactions often involve free energy changes of several keV, 
the chemical faith of SCRPs can be highly susceptible to magnetic fields as low as only a few $\mu$- Tesla.  
A prominent example is magnetoreception in migratory birds which is most likely mediated by
blue-light photoreceptors in the birds retina.

In the example notebook `SCRP_magnetic_field_effects.ipynb <https://simos.kherb.io/virtual/lab/index.html?path=examples%2FSCRP_magnetic_field_effects.ipynb>`_.
we explain how the chemical faith of SCRPS can be sensitive to even the smallest of magnetic fields and perform coherent simulations 
on a simple example SCRP system.


Research Example : Optical Detection of 2J Resonance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the example notebook `SCRP_optical_detection_MFE.ipynb <https://simos.kherb.io/virtual/lab/index.html?path=examples%2FSCRP_optical_detection_MFE.ipynb>`_.
we explore an advanced system and reproduce experimental findings of recent work by Buck et al.  
Here, the 2J resonance in a strongly coupled radical pair was detected via fluorescence lifetime. 
We qualitatively reproduce these findings by simulating a quantum master equation in Lindblad form 
taking into account both coherent spin dynamics and incoherent optical excitation and decay events.


