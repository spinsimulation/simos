A Primer on Quantum Mechanics
=============================

Here, we try to cover the most essential concepts of quantum mechanics that are required to understand
the basics of spin dynamics simulations. They are neither excessive nor complete and do not replace
a course on quantum mechanics or the study of standard textbook literature or review articles [CCH24]_ [MD20]_.


Quantum Mechanics in SimOS
--------------------------

The following scheme provides an overview on the different types of quantum objects and the
differential equations that describe their time dynamics. 

.. figure:: /img/primer_overview.png
   :width: 700px
   :align: center

Closed quantum systems, represented by state vectors or density operators, evolve
in a coherent manner under the system Hamiltonian, as decribed by the Schrödinger and  
von Neumann equations. All quantum objects are vectors or matrices in an N-dimensional Hilbert space.
SimOS provides excessive functionality to construct Hilbert spaces of composite quantum systems and 
helper functions for formulations of Hamiltonians that are common in spin dynamics.
In addition, we provide a series of standard matrix operations that are common in quantum mechanics.

Open quantum systems on the other hand, exchange energy with their environment, and this incoherent 
interaction interferes with their coherent time evolution. Accurate simulations of their time dynamics
requires a quantum master equation (QME). Here, our focus lies on the QME in Lindblad form, 
for which we provide various engines for computationally efficient time propagation. The combined
coherent and incoherent dynamics, described by Hamiltonians and collapse operators,
are formulated as superoperators in Liouville space and the resulting Liouvillian acts on vectorized
density operators. 

If you are familiar with quantum mechanics. just continue from here and:

#. Familiarize yourself with our  :ref:`backends <Backends>` and implementations of basic :ref:`matrix operations <QMMethods>`.
#. Discover our :ref:`hands-on introduction <System>` on how to construct composite Hilbert spaces in SimOS.
#. Check out helper methods to onstruct :ref:`coherent <Coherent>` and :ref:`incoherent <incoherent>` interactions.
#. Simulate :ref:`time evolution <TimePropagation>` of Liovuille-von-Neumann, Liouville Master equations or explore our variable
   module for simulations of spatial dynamics (i.e. Stochastic Schrödinger/Liouville equations).


Closed Quantum Systems
----------------------

The Hilbert Space
^^^^^^^^^^^^^^^^^

*Postulate:* For any isolated physical system there exists a vector space with an inner product. 
The system is completely described by a state vector, that is a unit vector in this 
Hilbert space (or: state space).

In the following we will refer to these state vectors that describe pure states in Hilbert space as kets
:math:`|\psi\rangle` and their adjoints (i.e. conjugate transposes) as bras :math:`\langle\psi|`.  The state of a single spin 1/2, for example
a single proton or electron, is adequately described by a state vector in a two-dimensional Hilbert space, i.e.
:math:`\ket{\psi} = c_1\ket{\phi_1} + c_2\ket{\phi_2}`. A common basis for this Hilbert space is the Zeeman
basis - simply speaking we interprete  :math:`\ket{\phi_1}` as spin-up and :math:`\ket{\phi_2}` as spin-down.

The state vector representation of the closed quantum system is, however, already unsufficient if we start 
to consider ensembles of a system or time-averages of multiple repetitions of the same experiment. If you are
not working with single spins but rather with millions of copies thereof  - as it is the case for all conventional
spin resonance - a state vector can no longer describe your system. It does not allow to describe the
probability distribution of various state vectors that comes naturally with the ensemble. 
To describe these mixed states, described by a set of pure state vectors with probabilities
:math:`\{ |\psi_i\rangle , p_i \}`, we instead utilize density operators (or: density matrices) that we construct 
from the state vector distribution as

.. math::
   \rho =  \sum_i p_i  |\psi_i\rangle \langle\psi_i|

The density matrix must have trace :math:`\text{Tr}(\rho) = 1` at any time and has to be 
positive semi-definite. For the density matrix of a single pure state, 
i.e. :math:`\rho = \ket{\psi} \bra{\psi}`, we further find :math:`\text{Tr}(\rho^2) = 1`. 
while for a mixed state, :math:`\text{Tr}(\rho^2) < 1`.

Once we know how to descibe our system, a natural next question is: How do we perform a measurement on it?
In spin dynamics, this question could translate into determining the spin polarization or
single-quantum coherence of a system, quantities that are directly related to our observables (e.g. a free-induction decay
or a fluorescence signal). 

*Postulate*: All possible measurements in a quantum system are described by a Hermitian operator or observable.

We can calculate the expected outcome of a measurement defined by an operator :math:`\hat{O}` on a state 
as

.. math::
   \langle \hat{O} \rangle = Tr[\hat{O}\rho]

In the case of a spin polarization, this operator would be the z-operator of the spin; in the case of single-quantum
coherence, we would choose a raising or lowering operator for :math:`\hat{O}`. 
A special operator is the Hamiltonian :math:`H`, whose expectation values are the system energies. 


The Schrödinger and the Liouville-von-Neumann Equations 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We have now learned how to describe states of closed quantum systems and how to perform measurements on them, but
of course these states are not static. They evolve under the system Hamiltonian in a coherent manner. 

*Postulate:* The coherent time evolution of a  pure state of a closed quantum system, i.e. a state vector, is described by the Schrödinger equation

.. math:: \frac{\partial}{\partial t} \ket{\psi(t)} = -\frac{i}{\hbar} H \ket{\psi(t)} = -i \mathcal{H} \ket{\psi(t)}.

Here we have introduced the very common notation :math:`\mathcal{H} = \frac{H}{\hbar}`, thus all Hamiltonians will be expressed in units of angular frequency.
For a density matrix :math:`\rho`, we can reformulate the Schrödinger equation to obain the von-Neumann equation

.. math::
    \frac{\partial}{\partial t} \rho(t) = -i [\mathcal{H}, \rho(t)].

In the simplest case, the Hamiltonian is time-independent. In this case, the solution of the Schrödinger equation is given by its propagator

.. math::
    U(t) = \exp(-i \mathcal{H} t).

and the time evolution of state vectors and density matrices is obtained as 

.. math::
    \ket{\psi(t)} = U(t) \ket{\psi(0)}, \quad \rho(t) = U(t) \rho(0) U^\dagger(t),

respectively.

Thus, for time independent Hamiltonians, computing the dynamics of a  closed quantum system is equivalent to computing 
the propagator of the Hamiltonian, i.e. the matrix exponential of the
quantum system. 
However, for most systems, the Hamiltonian is time-dependent. In this case, we can compute the time evolution 
using a time-ordered exponential. Esentially, we assume that the Hamiltonian is piecewise constant for small
enough time intervals. We compute the time evolution for each time interval and obtain the overall evolution operator
by multiplying the time evolution operators of the individual time intervals. Of course, this makes the simulation
computationally expensive - we now have to compute a series instead of a single propagator. 


Open Quantum Systems
--------------------

When dealing with open quantum systems, the situation is significantly more complex. 
The system is not isolated, but interacts with an environment and an accurate simulation
of their dynamics requires a quantum master equation that allows for incoherent time dynamics in addition
to the coherent system evolution. 

Typically, the environment is treated by introducting a so-called ancilla system.
The combined system of the open system and the ancilla is then closed and can be described by
a coherent time evolution.  At the end of the calculation, the ancilla is traced out,
using the partial trace operator and the reduced density matrix of the open system is obtained.

Lindblad-Form of the Quantum Master Equation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There exist multiple formulations of the quantum master equation, that differ slightly in their properties and 
applicability. Here we utilize the Lindblad form of the quantum master equation

.. math::
    \frac{\partial}{\partial t} \rho(t) = -i [\mathcal{H}, \rho(t)] + \sum_k \left( L_k \rho(t) L_k^\dagger - \frac{1}{2} \{L_k^\dagger L_k, \rho(t)\} \right).

where :math:`L_k` are collapse (or: jump) operators that describe the non-coherent part of the time-evolution resulting from the interaction of the
system with the environment (i.e. relaxation processes, incoherent excitations due to photons, or other
processes that lead to a change in the state of the system.)
Derviation of the Lindblad form involves several assumptions, most importantly that the dynamics are Markovian
and the system is weakly damped  such that the density matrix of the system remains separable from the baths degrees 
of freedom at all times.
However and very importantly, the Lindblad form of the quantum master equation guarantees a completely positive and trace preserving map (CPTP).
This means that the density matrix will always remain positive semi-definite and have trace 1 during evolution - i.e.
it will always remain a valid density matrix and the solution of the QME will always be physical.  

To obtain a generator of the Lindblad form of the quantum master equation, we have to introduce the so-called Liouville space. This vector space is of dimensionality :math:`d^2`, where :math:`d` is the dimensionality of the Hilbert space of the system. The Liouville space is spanned by the tensor product of the basis vectors of the Hilbert space. Operators on this space are so-called superoperators. The density matrix of the system is then represented by a vector in the Liouville space.

.. math::
    \vec{\rho} = \mathrm{vec}(\rho).

Technically, this vectorization is achieved by stacking the columns of the density matrix on top of each other. The Linblad equation in the Liouville space is then given by

.. math::
    \frac{\partial}{\partial t} \vec{\rho}(t) = \mathcal{L} \vec{\rho}(t) = (\mathcal{H} + \mathcal{G}) \vec{\rho}(t).

Here, :math:`\mathcal{L}` is the Liouvillian that takes the role of the Hamiltonian in the Schrödinger equation. The Liouvillian is the sum of the Hamiltonian part :math:`\mathcal{H}` and the Lindblad generator :math:`\mathcal{G}`. The Hamiltonian part is given by

.. math::
    \mathcal{H} = -i \mathcal{H} \otimes \mathbb{1} + i \mathbb{1} \otimes \mathcal{H}. 

The Lindblad generator is given by

.. math::
    \mathcal{G} = \sum_k L_k \otimes L_k^\dagger - \frac{1}{2} \left( \mathbb{1}\otimes L_k^\dagger L_k - L_k^\dagger L_k \otimes \mathbb{1} \right).

and the time evolution of the density matrix is obtained as

.. math::
    \vec{\rho}(t) = \exp(\mathcal{L} t) \vec{\rho}(0).

For time-dependent Liouvillians, we can apply the same trick as for time-dependent Hamiltonians: We assume that the Liouvillian is piecewise constant and then compute the time evolution for each time interval. The time evolution operator is then obtained by multiplying the time evolution operators of the individual time intervals.



Alternate Formulations of the Quantum Master Equation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The mathematical gurantee of the complete positivity and trace preservation of the Lindblad form 
of the quantum master equation is a very important advantage compared to other formulations of the quantum
master equation. 
Still, alternate formulations of the QME are commonly used in spin dynamics simulations.
For example, conventional NMR and EPR relaxation theory heavily relies on Bloch-Redfield-Wangness (BRW) theory,that
is very related to the Lindblad formulation and mainly differs in the absence of a final secular approximation.
Unfortunately, BRW does neither ensure a phyiscal solution nor relaxation to a finite temperature equilibrium state unless further corrections are applied. 



.. rubric:: References

.. [CCH24] Campaioli, F., Cole, J.H. and Hapuarachchi, H. Quantum Master Equations: Tips and Tricks for Quantum Optics, Quantum Computing, and Beyond. PRX Quantum 5, 020202 (2024).
.. [MD20] Manzano, D.A short introduction to the Lindblad master equation. AIP Advances 10, 025106 (2020).

