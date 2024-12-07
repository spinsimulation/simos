.. _NV:

Nitrogen-Vacancy Centers
========================

Color centers in diamond and silicon carbide have been studied for years and demonstrated significant potential in nanoscale sensing of magnetic and electric fields and as long-lived
quantum memories. The nitrogen vacancy center is the most well studied point defect in diamond and provides nanoscale resolution, long coherence
times at room temperature, all-optical spin state initialization and readout as well as coherent manipulation of spin states with resonant microwave fields. 
SimOS offers a submodule :class:`NV` that contains a series of high-level helper functions to simplify simulations of photophysics and nanoscale NMR experiments with NV centers.


System Initialisation
^^^^^^^^^^^^^^^^^^^^^

The :class:`NV` submodule provides an NV-specific class :class:`NVSystem` to initialise the multipartite quantum system of NV centers. A series of keyword arugments 
allow to account for various scenarios.

.. list-table:: Keyword arguments for NV system initialisation. 
   :widths: 25 25 50
   :header-rows: 1

   * - Keyword Argument
     - Default 
     - Description
   * - :code:`optics`
     - :code:`True`
     - If True, electronic structure of the NV is included. 
   * - :code:`orbital`
     - :code:`False`
     - If True, the excited electronic state is initialised as an orbital doublet. 
   * - :code:`nitrogen`
     - :code:`True`
     - If True, coupling to a nuclear nitrogen spin is included.
   * - :code:`synthetic`
     - :code:`True`
     - If True, we assume a synthetic (i.e. non-natural) NV center, formed by ion implanation, and the coupled nitrogen spin has a spin :math:`I= 1/2`. If False, a spin :math:`I = 1` is assumed for the coupled nitrogen.
   * - :code:`further_spins`
     - :code:`[]`
     - Describes further nuclear spin that are coupled to the NV center.


If all keyword-arguments are omitted, the initialised system is representative of an NV center at room temperature and native SimOS code
for its construction would be:

.. code::

    # NV Center Electron Spin 
    S = {'val': 1, 'name':'S', 'type': 'NV-'} 
    # Nuclear Nitrogen Spin
    I = {'val': 1/2, 'name':'I', 'type': 'N'} 
    # Electronic States
    GS = {'val': 0 , 'name':'GS'}
    ES = {'val': 0 , 'name':'ES'}
    SS = {'val': 0 , 'name':'SS'}
    s = sos.System([([(GS, ES), S], SS), N])

The following scheme visualizes the construction of the composite Hilbert space in a stepwise fashion.

.. figure:: /img/NV_hilbert.png
   :width: 650px
   :align: center


|

NV Hamiltonian
^^^^^^^^^^^^^^

The :class:`NVSystem` provides users with a method :code:`field_hamiltonian()` which returns the Hamiltonian of the NV electron spin (:math:`S=1`) in the electronic ground and 
excited state.  
The ground state spin Hamiltonian 

.. math::
  \hat{H}_{\mathrm{GS}} = D_{\mathrm{GS}} \left[ \hat{S}_z^2 - \frac{2}{3} \mathbb{I}_3 \right] +  \gamma_{\mathrm{NV}}\hat{\vec{S}} \cdot \vec{B}  \\
  + d_{\mathrm{GS}}^{\parallel} (E_z + \delta_z) +   d_{\mathrm{GS}}^{\perp}(E_x + \delta_x) (\hat{S}_y^2 - \hat{S}_x^2) + d_{\mathrm{GS}}^{\perp}(E_y + \delta_y)(\hat{S}_x\hat{S}_y + \hat{S}_y\hat{S}_x)

contains the ground state zero-field splitting of magnitude :math:`D_{\mathrm{GS}}` , the Zeeman interaction of the NV electron spin (gyromagnetic ratio :math:`\gamma_{\mathrm{NV}}` ) in an
external magnetic field :math:`\vec{B}` and the stark interaction with the sum of electric  :math:`\vec{E}` and strain :math:`\vec{\delta}` fields. We distinguish parallel :math:`d_{\mathrm{GS}}^{\parallel}` 
and perpendicular :math:`d_{\mathrm{GS}}^{\perp}` electric dipole moments of the NV center ground state.

For the excited state Hamiltonian, we distinguish room- and low-temperature situations. At low temperature the orbital substructure of the excited state is explicitly taken into account and the 
excited state Hamilonian is

.. math::
  \hat{H}_{\mathrm{ES}}^{\mathrm{LT}} = D_{\mathrm{ES}}^{\parallel}\mathbb{I}_2 \otimes \left[ \hat{S}_z^2 - \frac{2}{3} \mathbb{I}_3 \right] + D_{\mathrm{ES}}^{\perp}\left[\hat{\sigma}_z \otimes (\hat{S}_y^2 - \hat{S}_x^2) - \hat{\sigma}_x \otimes (\hat{S}_y\hat{S}_x + \hat{S}_x\hat{S}_y) \right] \\
  - \lambda^{\parallel}\hat{\sigma}_y\otimes\hat{S}_z + \lambda^{\perp} \left[\hat{\sigma}_z\otimes(\hat{S}_x\hat{S}_z + \hat{S}_z\hat{S}_x ) - \hat{\sigma}_x\otimes(\hat{S}_y\hat{S}_z + \hat{S}_z\hat{S}_y ) \right] \\
  +  \gamma_{\mathrm{NV}}\mathbb{I}_2 \otimes (\hat{\vec{S}} \cdot \hat{\vec{B}}) + \gamma_{\lambda} (\hat{\vec{\sigma}}\cdot \vec{B}) \otimes \mathbb{I}_3  \\
  +  d_{\mathrm{ES}}^{\parallel} (E_z + \delta_z) \mathbb{I}_6+ d_{\mathrm{ES}}^{\perp} (E_x + \delta_x)\hat{\sigma}_z\otimes \mathbb{I}_3 +  d_{\mathrm{ES}}^{\perp} (E_y + \delta_y)\hat{\sigma}_x\otimes \mathbb{I}_3


where  :math:`\hat{\sigma}_i` (:math:`i  =  x, y, z`) are standard Paul matrices for the excited state orbital operators, :math:`D_{\mathrm{ES}}^{\parallel}` and :math:`D_{\mathrm{ES}}^{\perp}` are the
parallel and perpendicular of the excited state zero field splitting, :math:`\gamma_{\lambda}` is the  orbital gyromagnetic ratio and :math:`d_{\mathrm{ES}}^{\parallel}` 
and  :math:`d_{\mathrm{ES}}^{\perp}` are the parallel and perpendicular components of the excited state electric dipole moment. At room temperature, the orbital doublet is not explicitly included 
and the excited state Hamiltonian simplifies accordingly. The :code:`field_hamiltonian()` function probes whether the :class:`NVSystem` was initialised with the orbital substructure 
to distinguish low- and room-temperature cases. 

Since all NV-specific constants are automatically sourced from the SimOS library, the user only has to specify magnetic, electric and strain fields to :code:`field_hamiltonian()`. 
The magnetic field is specified with a keyword argument and the method defaults to zero magnetic field if the argument is omitted.  
Since the excited state dipole moment (i.e. :math:`d_{\mathrm{ES}}^{\parallel}` 
and  :math:`d_{\mathrm{ES}}^{\perp}` ) are not known, electric and strain fields are not input as fields but rather as frequencies which represent products
of electric dipole moment and combined electric and strained fields.


Photophysics
^^^^^^^^^^^^

The spin state of the negatively charged NV center can be optically prepared and read out using green light illumination
enabling single-spin detection spanning room temperature to cryogenic temperature. To simulate the spin-dependent photoluminescence dynamics of single NV centers, the electronic
level structure must be taken into account. 

In the simplest case, valid at room temperature, the negatively charged NV center electronic structure may be described by three electronic levels: a triplet ground state, a triplet excited state and a metastable
singlet state. Here, we consider off-resonant excitation of the NV center from the ground to the excited state, which is commonly achieved with green (:math:`\sim` 532 nm) laser light.
The excited state has an optical lifetime of about 16 ns and decays either radiatively by emission of a red photon, or non-radiatively via an intersystem crossing (ISC) to the metastable
singlet state. This ISC is strongly spin dependent, resulting in (i) higher fluorescence intensity for the :math:`m_S = 0` spin state versus the :math:`m_S = 1` spin state and
(ii) polarization of the spin state to :math:`m_S = 0` throughout multiple excitation-decay cycles. 


.. figure:: /img/NV_photoscheme.png
   :width: 650px
   :align: center


The simplest simulation of the NV center photophyiscs is purely incoherent, i.e. it does not take into account any coherent spin dynamics and solely
simuates the effect of incoherent optical excitation and decay events with the simplified room-temperature model. To perform such a simulation, we construct the representative quantum system as:

.. code::

    NV = sos.NV.NVSystem(nitrogen = False)

The :class:`NVSystem` provides a method :code:`transition_operators()` which returns all optical decay rates of the NV center. At room temperature, these can be assumed to be independent of exact temperature and external magnetic
and electric fields. A parameter :code:`beta` allows to indicate the power of the laser illumination relative to the NV center saturation power (saturation for :code:`beta = 1`).

.. code::

    # c_ops_on transition operators if laser is on, c_ops_off if laser is off
    c_ops_on, c_ops_off = NV.transition_operators(beta = 0.2)

We simulate the photoluminescence dynamics during a readout laser pulse for an
NV center in the :math:`m_S = 0` and :math:`m_S = 1` state, respectively, with the code provided below.

.. code::
    
    pl0 = []
    pl1 = []
    U = sos.evol( None, dt , c_ops = c_ops_on )
    rho0 = (NV.Sp[0]).copy()
    rho1 = (NV.GSid*NV.Sp[1]).copy()
    for t in tax :
        pl0.append(sos.expect(NV.ESid, rho_0 ) )
        pl1.append(sos.expect(NV.ESid, rho_1 ) )
    rho0 = sos . applySuperoperator (U , rho0 )
    rho1 = sos . applySuperoperator (U , rho1 )


While this simple and fully incoherent model sufficiently captures NV center photophysics at elevated temperatures, it does not accurately describe effects at cryogenic
and intermediate temperature regimes.  Here, interplay of spin and orbital dynamics in the excited state may result in fast spin relaxation and vanishing ODMR contrast.
Accurate simulations must be performed with a quantum master equation which includes coherent spin dynamics under external magnetic and electric (strain) fields as well 
as incoherent excitation and decay dynamics accounting for the orbital character of the excited state and phonon-induced hopping.
We initialize the representative quantum system and define the system Hamiltonian and collapse operators for a specific temperature :code:`T` and magnetic and electric (strain) fields :code:`Bvec, Evec` using the :code:`transition_operators()` and :code:`field_hamiltonian()` methods
of the NV submodule.

.. code::
    NV= sos.NV.NVSystem(nitrogen = False, orbital = True)
    HGS, HES = NV.field_hamiltonian(Bvec=Bvec,Evec=Evec)
    H = HGS + HES
    c_ops_on, c_ops_off = s.get_transitionoperators(T=300,  beta=0.2, Bvec=Bvec, Evec=Evec)


Again, we simulate the photoluminescence dynamics during a readout laser pulse for an NV center in the :math:`m_S = 0` and :math:`m_S = 1` state, respectively. 

.. code ::

    pl0 = []
    pl1 = []
    U  = sos.evol((HGS+HES), dt, c_ops = c_ops_on)
    rho0 = (NV.GSid*NV.Sp[0]).copy()
    rho1 = (NV.GSid*NV.Sp[1]).copy()
    for t in tax:
        pl0.append(sos.expect(NV.ESid, rho_0))
        pl1.append(sos.expect(NV.ESid, rho_1))  
        rho0 = sos.applySuperoperator(U,rho0)
        rho1 = sos.applySuperoperator(U,rho1)

The figure below illustrates the results for three different temperature, clearly illustrating vanishing ODMR contrast at intermediate temperature regimes where phonon-mediated hopping interferes
with optical spin state preparation and readout. 


.. figure:: /img/NV_photoresults.png
   :width: 650px
   :align: center



.. note::
    For a full code example on photophysics simulations  you can visit our Virtual Lab and follow the `NV.ipynb <https://simos.kherb.io/virtual/lab/index.html?path=examples%2FNV.ipynb>`_ notebook. 



Nanoscale NMR
^^^^^^^^^^^^^


NV centers have been used to perform nanoscale NMR measurements of proximal nuclei both within the diamond lattice and, when within several nm of the surface, 
external nuclei in chemical species.  The NV center is further capable of detecting the free induction decay of nearby nuclear spins using a measurement technique
referred to as weak measurements. Here, the nuclear spin is polarized and rotated with a :math:`\frac{\pi}{2}` pulse to initiate precession.
Repeated weak measurements on the NV center, consisting of optical spin initialization, XY-8 dynamical decoupling frequency encoding,
and optical readout are used to track the nuclear spin precession. 

.. figure:: /img/NV_nmrscheme.png
   :width: 650px
   :align: center


In the simplest case, assuming ideal optical initialization and readout, this measurement protocol can be simulated without considering the photophysics of the NV center. 
The quantum system is then constructed from the NV center's electronic spin :math:`S=1` which is coupled to a single :math:`^{13}` C nuclear spin :math:`S=1/2`. 

.. code::

  C = {"name": "C", "val": 1/2}
  NV = sos.NV.NVSystem(further_spins = [C], optics = False, nitrogen = False)


In a next step, we prepare the initial state that represents the polarized and rotated :code:`^{13}`C spin coupled to a polarized NV center and define the system Hamiltonian
from the native system operators. We work in the rotating frame for the NV center where the hyperfine interaction is truncated to the secular and pseudo-secular contributions 
nd utilize the :code:`^{13}`C gyromagnetic ratio of SimOS. 


.. code::

  # Initial State
  rho0 = sos.state(s, "S[0],C[0.5]")
  rho1 = sos.rot(s.Cx,np.pi/2,rho0)   
  # System Hamiltonian
  H0 = sos.yC13*B0*s.Cz + apara*s.Sz*s.Cz + aperp*s.Sx*s.Cx


Finally, we perform the sensing protocol, consisting of a series of :math:`N` XY-8 blocks that are interspaced by free evolution periods.
Again, we make use of the NV sub-module functionality and utilize the :code:`XY8()` and :code:`meas\_NV` routines to apply the dipolar decoupling sequence and perform NV center readout.
We obtain the free-induction decay of the nuclear spin and, after Fourier transformation, the corresponding spectrum. 

.. code::

  N = 500
  frf = sos.w2f(B0*sos.yC13)
  twait = 0.125/frf
  tau = 1/frf/2
  store = []
  rho = rho1.copy()
  for i in range(N):
      rho = sos.rot(s.Sop_y_red,np.pi/2,rho)
      rho = sos.NV.XY8(H0,tau,s,rho,N=16)
      rho = sos.rot(s.Sop_x_red,np.pi/2,rho)
      meas,rho = sos.NV.meas_NV(rho,s)
      store.append(meas)
      rho = sos.evol(H0,twait,rho)
  Tsample = 16*tau + twait
  foffset = np.round(Tsample*frf)/Tsample

.. figure:: /img/NV_nmrresults.png
   :width: 650px
   :align: center

Syntax Reference
^^^^^^^^^^^^^^^^


.. py:currentmodule:: simos.systems.NV
.. automodule:: simos.systems.NV
   :members: NVSystem, decay_rates, phonon_rates, coupl2geom, auto_pairwise_coupling, gen_rho0, get_NVaxes, lab2NVcoordinates, XY8, meas_NV, exp2cts,normalize_data, Wcp 



