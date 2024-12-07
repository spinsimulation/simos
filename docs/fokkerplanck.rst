.. _FokkerPlanck:


Simulation of Stochastic Spatial Dynamics and Parameter Averages
================================================================

Although the Schrödinger, Liouville, and Lindblad equations are powerful tools to simulate the time-evolution of quantum systems, they are not well suited to incorporate stochastic dynamics.
Stochastic contributions can arise for ensembles of quantum systems, (multiple system copies or measurement repetitions) if the Hamiltonian or Liouvillian is not uniform
across the ensemble. In this case, individual ensemble members evolve along different trajectories and the ensemble aver-
age differs from the result for a single member.

.. figure:: /img/random_walk.gif
   :width: 500px
   :align: center


If these trajectories are not dynamically intertwined (i.e. a system copy does always stay on a single trajectroy as shown above) the individual systems
can be evolved separately with subsequent averaging of their results. An example for such a "static" stochastic situation is a powder average. Here, we interpret
the ensemble as multiple system copies with different spatial orientations relative to a magnetic field and thus varying contributions of anisotropic magnetic interactions (e.g., anisotropic chemical
shifts or dipolar couplings). 

However, the trajectories can be intertwined if system jumps occur throughout the time evolution:

.. figure:: /img/random_walk2.gif
   :width: 500px
   :align: center

The ensemble members can then no longer be evolved in an independent manner and simulating the ensemble dynamics becomes non-trivial. A common example for this "dynamic" stochastic situation
are spatial dynamics in magnetic resonance experiments, for example magic angle spinning or rotational diffusion.


The Fokker-Planck formalism
---------------------------

Recently, Kuprov and coworkers advocated for the Fokker–Planck formalism as a universal and elegant approach to include
stochastic contributions in magnetic resonance simulations [IK16]_.  SimOS uses this framework in a generalized manner, simplifying the incorporation
of arbitrary stochastic dynamics in simulation routines. 

The stochastic dynamics is formulated in a higher dimensional Fokker–Planck space, obtained as the tensor product of the original system space
(i.e., the Hilbert or Liouville space of the multipartite system) and a classical state space. The latter discretizes the stochastically varying, classical conditions 
(e.g., molecule orientation, rotor phase, field value of fluctuating fields) and its dimensionality determines how accurate the stochastic dynamics
are being captured.

.. figure:: /img/fokkerplanck_space.png
   :width: 300px
   :align: center

|

The equation of motion in Fokker--Planck space results as 


.. math::
   \frac{\partial}{\partial t} \vec{\rho}_{\mathrm{FP}} = -i \mathcal{F} \vec{\rho}_{\mathrm{FP}} = -i \mathcal{Q} \vec{\rho}_{\mathrm{FP}} + \omega \left( \left[ \frac{\partial}{\partial \zeta}\right]^n \otimes \mathbb{1} \right) \vec{\rho}_{\mathrm{FP}} 


where :math:`\mathcal{F}` is the Fokker--Planck superoperator and and :math:`\vec{\rho}_{\mathrm{FP}}` is a vectorized, Fokker--Planck space representation of the ensemble state.
The Fokker--Planck superoperator :math:`\mathcal{F}` has two contributions: 

* a block diagonal part :math:`\mathcal{Q}` where every block holds the Hamiltonian or Liouvillian of a distinct classical condition


.. figure:: /img/fokkerplanck_hamiltonian.png
   :width: 200px
   :align: center

|

*  an off-diagonal part :math:`\omega \left( \left[ \frac{\partial}{\partial \zeta}\right] \otimes \mathbb{1} \right)` that exchanges populations between the classical sub-spaces at a frequency :math:`\omega`. The Fourier differentiation matrix :math:`\left[ \frac{\partial}{\partial \zeta}\right]^n` with respect to the classical coordinate :math:`\zeta` may be of order :math:`n=0,1,2`.
   If :math:`n=0`, no exchange of population occurs between the classical subspace; if :math:`n=1`, the dynamics are linear; if :math:`n=2`, diffusive dynamics are introduced.   

.. figure:: /img/fokkerplanck_mixer.png
   :width: 300px
   :align: center


|

The Fokker--Planck space representation of the ensemble state  :math:`\vec{\rho}_{\mathrm{FP}}` is obtained by concatenating the state vectors or density operators
of individual ensemble members (i.e. one state per classical parameter value) to a single vector.

.. figure:: /img/fokkerplanck_states.png
   :width: 450px
   :align: center

|

Implementation in SimOS
-----------------------
.. Attention:: 
   The Fokker-Planck routines currently only support the `qutip` backend.


In SimOS, we provide a specialized object, :class:`StochasticLiouvilleParameters`, to construct all permutations of the classical state space variables.
It is initialized similar to a conventional python dictionary. Each dictionary key represents an independent stochastic parameter (e.g. an angle or a position coordinate)
for which users can specify 

* :code:`values` A list of values that discretizes the stochastic parameter. The amount of provided values sets the dimension of the classical subspace and thus has a strong influence on the computational cost of the Fokker-Planck simulation. 
* :code:`dynamics` A dictionary that specifies the characteristics of the stochastic dynamics as :code:`{ n : w }` where :code:`n` is the order of the Fourier diffrentiation matrix and :code:`w` is the frequency of the dynamics.
* :code:`weights` A list of weights for the individual stochastic parameter values. 

The code-example below illustrates how to implement two independent stochastic variables. The parameters :code:`a` and :code:`b` both take three values, 1,2,3 and 10, 20,30, respectively. The parameter :code:`b` is
related to diffusive dynamics as the key of the dynamics-dictionary is 2. 

.. code-block:: python

   params = sos.StochasticLiouvilleParameters()
   params['a'].values = [1,2,3]
   params['a'].dynamics = None

   params['b'].values = [10,20,30]
   params['b'].dynamics = {2:20e3}


To perform the actual simulations, users may utilize the :code:`stochastic\_evol()` routine:

.. code-block:: python

   rho = sos.stochastic_evol(H_fun, params, dt,rho0)

   def H_fun(a,b):
      # H0, H1, H2 are defined elsewhere
      return H0 + H1(a) + H2(b)

The method passes the appropriate combinations of the stochastic parameters 'a' and 'b' to a user-defined function :code:`H_fun` that returns the Hamiltonian / Liouvillian for a given set 
of parameters. The routine then automatically constructs the Fokker--Planck superoperator  :math:`\mathcal{F}` and transforms and extracts the system state to and from the Fokker--Planck space.  Note that the  Quantum Space utilized in the routine can be either the Hilbert space of the quantum system or the Liouville space,
depending on whether a state vector or a density matrix is provided as the input argument :code:`rho0`.  Also, although not shown in the code example above, the method 
also supports inclusion of incoherent dynamics and users may provide a list of collapse operators as a keyword input argument. In this case, the simulation is obviously  always carried out in the Liouville space. 


The :code:`stochastic\_evol()`  routine does not necessarily return the evolved state. If the initial state is omitted, the method returns the propagator in Fokker Planck space. The latter
can be applied to a vectorized Fokker-Planck space density operator. The code below illustrates how users can evolve a state in Fokker Planck space using propagator
caching. Here, users must utilize the :code:`dm2fp()` and :code:`fp2dm()` methods to transform the density matrices of the system to and from the Fokker Planck space. 
Both methods require a series of input argument for the conversion with are described in detail in the syntax reference of the methods. 

.. code-block:: python
   
   # transform dm to fp space
   rhotmp = sos.dm2fp(rho, params.dof, try_ket=False)
   # cache propagator
   U =  sos.stochastic_evol(H_fun, params, dt, method='qutip', space='liouville')
   # apply propagator for N steps 
   store = []
   for i in range(N):
      # evolve
      rhotmp =  U*rhotmp
      # transform dm back from fp space 
      rho =  sos.fp2dm(rhotmp, is_hilbert=False, dims=s.dims, weights=params.tensor_weights())
      store.append(rho)

.. note::
   The :code:`stochastic\_evol()` routine does only work for time-independent Hamiltonians and collapse operators. Time-dependencies of the non-stochastic dynamics are
   currently not supported. 



Syntax Reference
----------------


.. py:currentmodule:: simos.fokker
.. automodule:: simos.fokker
   :members: StochasticLiouvilleParameters, stochastic_evol, StochasticLiouvilleParameters.__init__, StochasticLiouvilleParameters.tensor_values, StochasticLiouvilleParameters.tensor_mixer, StochasticLiouvilleParameters.tensor_mixer, StochasticLiouvilleParameters.dof

.. rubric:: References

.. [IK16] Kuprov, I.  Fokker-Planck formalism in magnetic resonance simulations. J. Magn. Reson., 270, 124-135 (2016)



