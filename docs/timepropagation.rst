.. _TimePropagation:

Time Propagation
================
The time-evolution of quantum systems is a fundamental problem in quantum mechanics. In SimOS, we provide two methods to simulate the time-evolution of quantum systems: time-independent evolution and time-dependent evolution. The time-independent evolution is used to simulate the time-evolution of a quantum system under a time-independent Hamiltonian. 

Time-independent evolution in SimOS
-----------------------------------
The time-independent evolution is implemented in the function :func:`simos.propagation.evol`. The function handles both, the evolution of a quantum system under a time-independent Hamiltonian and under a time-independent Liouvillian. The function is called as follows:

.. py:currentmodule:: simos.propagation
.. automodule:: simos.propagation
   :members: evol

Time-dependent evolution in SimOS
---------------------------------
The time-dependent evolution is implemented in the function :func:`simos.propagation.prop`. The function handles the evolution of a quantum system under a time-dependent Hamiltonian. Furthermore, it can also handle time-dependent Liouvillians. Here, all combinations of time-dependent Hamiltonians and Liouvillians are possible: Either the Hamiltonian can be time-independent with time-dependent collapse operators or the Hamiltonian can be time-dependent with time-independent collapse operators. Handling both time-dependent Hamiltonians and Liouvillians is also possible. 

We assume that the time-depndence takes the form 

.. math::
   H(t) = H_0 + \sum_i^N c_i(t) H_i

where $H_0$ is the time-independent part of the Hamiltonian and $H_i$ are the time-dependent parts of the Hamiltonian. The coefficients $c_i(t)$ are time-dependent modulation functions. 

To specify this, the user needs to specify the static Hamiltonians $H_i$ as a list that is passed to the ``H1`` argument. The time-dependent coefficients $c_i(t)$ are specified as a list of arrays that are passed to the ``carr1`` argument. ``dt`` is the time-step of the simulation. The function :func:`simos.propagation.prop` is called for example as follows:
.. code-block:: python
   rho = simos.prop(H0,dt,rho0,H1=[H1,H2],carr1=[c1,c2])

Analogously, the time-dependent Liouvillians can be specified by passing a list of collapse operators to the c_ops1 argument. 

The :func:`simos.propagation.prop` routine supports multiple so-called engines. The engines are specified by the ``engine`` argument. The engines are used to simulate the time-evolution of the quantum system. For all backends, the engine cpu is supported. When working with the qutip backend, the native qutip propagation can be used by specifying the engine qutip. The input structure of the prop function is very well suited to parallelize the time-propagation by leveraging a GPU. This is achieved via the inclusion of ``parament``. Here, the propagation is handled on the GPU and only a minimal amout of data is transferred between the CPU and the GPU. The GPU engine is specified by the engine parament. The engine parament is only available when the parament package is installed.

.. figure:: /img/parament.jpg
   :width: 500px
   :align: center

.. raw:: html
   <br />&nbsp;<br>

.. py:currentmodule:: simos.propagation
.. automodule:: simos.propagation
   :members: prop


Explicit Simulation of Spatial Dynamics
=======================================
The Stochastic Liouville equation is a powerful tool to simulate the spatial dynamics of quantum systems. In the general form of a Focker-Planck equation a variety of situations such as stochastic rotational diffusion but also magic angle spinning can be simulated under the same umbrella. The key is the realizaion of many trajectories of the quantum system. We use a specialized object, :class:`simos.focker.StochasticLiouvilleParameters`, to simulate the spatial dynamics of quantum systems. The function :func:`simos.focker.stochastic_evol` is then used to carry out the simulation.



.. py:currentmodule:: simos.focker
.. automodule:: simos.focker
   :members: StochasticLiouvilleParameters, stochastic_evol, StochasticLiouvilleParameters.__init__, StochasticLiouvilleParameters.tensor_values, StochasticLiouvilleParameters.tensor_mixer, StochasticLiouvilleParameters.tensor_mixer, StochasticLiouvilleParameters.dof