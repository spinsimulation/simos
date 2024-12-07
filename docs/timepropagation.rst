.. _TimePropagation:

Time Propagation
================

The time-evolution of quantum systems is a fundamental problem in quantum mechanics that can be computationally very demanding,  especially for large Hilbert and Liouville spaces.
SimOS provides two basic routines to simulate the time-evolution of (open) quantum systems:

* the :code:`evol()` routine allows to simulate dynamics under time-independent Hamiltonians or Liouvillians
* the :code:`prop()` routine allows to simulate evolution of a quantum sysmtem under time-dependent Hamiltonian or Liouvillians. Both the coherent as well as the incoherent contributions may be time-dependent.

Time-independent evolution in SimOS
-----------------------------------

In the simplest case, the system Hamiltonian (or Liouvillian, if incoherent dynamics are present) is static, i.e. time-independent. In this case,
the :code:`evol()` routine can be utilized so simulate the time evolution of the quantum system. 
The method calculates the matrix exponential of the Hamiltonian/Liouvillian (i.e. the propagator :math:`U`) and applies it to an initial state as explained in more detail in our :ref:`QM-primer <QMPrimer>`. 
If incoherent dynamics are included, a Lindblad form of the QME is utilized.   

The function is called with the static system-Hamiltonian :code:`H` , the time step (i.e. duration) of the evolution period :code:`t` and, if incoherent dynamics are to be included in the simulation,
a list of collapse operators :code:`c_ops` (keyword-argument).  If an initial state :code:`rho_0` is provided,
the evolved state is returned.  If the initial state is omitted, the function instead returns the propagator of the evolution. The latter is especially useful to calculate the
evolution under the same static Hamiltonian for many time steps, since the costly matrix exponential only has to be evaluated once. 

.. note::
   To perform state rotations on a Bloch sphere, i.e. rotations with well-defined rotation angle and axis in the complex vector space of the quantum system,
   you may use the :code:`simos.propagation.rot` routine instead. It functions identical to :code:`simos.propagation.evol`, however it acccepts a rotation operator 
   and rotation angle as input arguments and is therefore the more convenient interface for simple state rotations. 
   In full analogy to the :code:`simos.propagation.evol` method, state evolution is performed if an initial
   state is being provided while the propagator for the rotation is returned if the initial state is omitted.
   To visualize state rotations, you may also find our `Bloch-sphere visualizer <https://bloch.kherb.io/>`_ helpful. 

Let us consider a simple example, a single electron spin 1/2 that is polarized along an external magnetic field :math:`B_0`.  We want to perform a :math:`\pi/2` pulse (i.e. a 90-degree rotation on the Bloch sphere)
and evolve the state under a static Zeeman interaction with the external magnetic field for a single time step. For the 90-degree rotation, we utilize the :code:`rot()` routine, for the evolution
under the static Zeeman Hamiltonian we use :code:`evol()`.

.. code-block:: python

   # quantum system of a single electron spin 
   S  = {'name': 'S', 'val': 1/2}
   s = sos.System([S]) 
   # the initial state, polarized electron spin 
   rho0  = sos.state(s, 'S[0.5]')
   # rotate state with a pi/2 pulse    
   rho0x = sos.rot(s.Sx, np.pi/2, rho0)  
   # evolve the state
   dt = 1e-6 # time step [s]
   B = 10e-3 # magnetic field [mT]
   HZ = sos.ye*s.Sz # the Zeeman Hamiltonian 
   rhoxevol = sos.evol(HZ, dt, rho0x) 


Time-dependent evolution in SimOS
---------------------------------


Hamiltonian and collapse operators may become time-dependent under the application of time-varying control fields, e.g., magnetic or electric fields or
shaped laser-excitation pulses. In most cases, the dynamic behavior can be efficiently parametrized with a limited set of control fields and collapse operators without any loss of generality.
The time-varying Hamiltonian and collapse operators are formulated as

.. math::
    \hat{H}(t) = \hat{H}_0 + \sum_{i=1}^N c_i(t) \hat{H}_i

and

.. math::
    \hat{L}_k(t) = \hat{L}_{k,0} + \sum_{i=1}^N C_i(t) \hat{L}_{k,i}


with time-independent basis functions  :math:`\hat{H}_i` and :math:`\hat{L}_{k,i}` and control amplitudes :math:`c_i(t)` and :math:`C_i(t)`. 

The :code:`prop()` routine of SimOS accepts the time-independent parts :math:`\hat{H}_0` (:math:`\hat{L}_{k,0}`) as well as lists of time-varying control operators 
:math:`\hat{H}_i` (:math:`\hat{L}_{k,i}`) and their control amplitudes :math:`c_i(t)` (:math:`C_i(t)`) as input arguments.
The time-dependence of the control amplitudes is assumed to be piecewise constant and control amplitudes are specified for series of discretized time-intervals :math:`dt`.
A detailed explanation of all input parameters for :code:`prop()` can be found in the syntax reference of the method. 

.. note::

   The :code:`prop()` routine is only compatible with and thus implemented for our numerical backends and not for our symbolic Sympy backend.
   To propagate quantum systems in a symbolic matter, please utilize the :code:`evol()` routine instead. 



Engines
^^^^^^^

The development of computationally efficient integrator schemes for the numerical propagation of quantum systems is an active field of research.
Besides classical ordinary differential equation (ODE) solvers, Euler-type integrators are common in the field of magnetic resonance simulations.
Here, in full analogy to the static solutions presented in our :ref:`QM-primer <QMPrimer>`, the matrix exponential of the piecewise constant Hamiltonian
or Liouvillian is evaluated separately for each time interval.
The repeated computation of the matrix exponential is computationally expensive.
However, a more efficient implementation can be achieved by using parallel-in-time integrators such as PARAMENT, introduced recently by some of use.
We have further shown  that Euler-type integrators can be easily converted into Magnus-type integrators that benefit from much quicker convergence. 

The :code:`prop()` routine serves as an interface to different so-called **engines** that perform the calculation and are readily selected via the :code:`engine` keyword argument.
You can thus select an integrator scheme that is optimal for your specific problem characteristics such as dimensionality and sparsity of the operators.

.. list-table:: Engines for time propagation.
   :widths: 25 25 50
   :header-rows: 1

   * - Engine
     - Compatible Backends
     - Description
   * - :code:`cpu`
     - All backends
     - CPU-based Euler-type or Magnus-type integrator (selected with keyword argument :code:`Magnus`)
   * - :code:`parament`
     - All backends 
     - GPU-based Magnus integrator. Utilizes the Parament package to perform parallelized time-propagation on a GPU, requires installed Parament and a GPU
   * - :code:`qutip`
     - Only available for qutip backend 
     - Uses native QuTIP :code:`mesolve` method, which is based on a Runge-Kutta solver
   * - :code:`RK45`
     - All backends.
     - Utilizes :code:`scipy.integrate` 4th order Runge-Kutta

.. list-table::


.. note::
   All technical details of our Parament engine are available in the `original publication <https://doi.org/10.3929/ethz-b-000514331>`_. 
   Esentially, Parament leverages the computation power of graphics processing units (GPUs) to perform the integration of all time steps in parallel.



Syntax Reference
----------------

.. py:currentmodule:: simos.propagation
.. automodule:: simos.propagation
   :members: evol, prop
