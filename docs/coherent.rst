
.. _Coherent:

Coherent Interactions
=====================

The coherent dynamics of spin ensembles are mostly governed by interactions of spins with external fields
and pairwise spin-spin interactions. Simos provides a series of methods to simplify construction of Hamiltonians
for the most common sources of coherent dynamics of optically adressable spins.


.. figure:: /img/coherent_scheme.png
   :width: 750px
   :align: center

|

Of course you do not have to utilize these functions to introduce coherent interactions - especially for simple,
isotropic interactions (or interactions which are not covered by our helper-functions) it is very easy to just define them "by hand" using the operators of your spin system.
For example, a Zeeman interaction of a spin can simply be introduced as

.. code-block:: python

   S  = {"name": "S", "val": 1/2}
   s = sos.System([S])
   Hz = sos.ye * s.Sz

.. tip::
   SimOS provides you with a set of :ref:`physical constants <Constants>` as well as a 
   complete set of :ref:`gyromagnetic ratios <GyroRatios>`.


General Framework
-----------------

Any coherent spin-spin or spin-spin or spin-field interaction can be described by a Hamiltonian of the form

.. math::
   \hat{H} =  \vec{\hat{S}} \cdot \mathbf{A} \cdot \vec{X}

where :math:`\vec{\hat{S}} = (\hat{S_x}, \hat{S_y}, \hat{S_z})` is a vector of spin operators, :math:`\textbf{A}` is a 3x3 matrix which describes
the spatial part of the interaction and :math:`\vec{X}` is either a spin-operator vector for a second spin, i.e. :math:`\vec{X} \equiv (\hat{I_x}, \hat{I_y}, \hat{I_z})`, or a magnetic
field vector i.e. :math:`\vec{X} \equiv (B_x, B_y, B_z)`.

The interaction matrix :math:`\textbf{A}`  can alternatively be formulated as spherical tensors of rank 0, 1 and 2 whose construction and utilization is 
detailed in our section on :ref:`spherical tensors <Spherical>`. 
Further, we can for each interaction distinguish the laboratory frame of reference (LAB) from a principal axis system (PAS) of reference. In the PAS of the interaction, the rank 0 and rank 2 
contributions, which together constitute the symmetric part of :math:`\mathbf{A}`, are diagonal. The PAS representation of :math:`\mathbf{A}`, in the following referred to as  :math:`\mathbf{a}`,  
is conveniently parametrized via

*  the trace :math:`\overline{a}  = (a_{xx} + a_{yy} + a_{zz}) / 3` for rank 0
*  the elements :math:`a_{xy}, a_{xz}, a_{yz}` for the rank 1 component, and,
*  the anisotropy :math:`\delta = a_{zz} - \overline{a}` and asymmetry :math:`\eta = \frac{a_{yy}-a_{xx}}{\delta}` for the rank 2 component.

Here we assume the ordering

.. math::
   |a_{zz} -\overline{a}| \geq  |a_{xx} -\overline{a}|  \geq |a_{yy} -\overline{a}| 

of the principal components which ensures that :math:`\eta` is always positive and smaller than one. Transformation from  PAS to LAB frame is further parametrized with a set of three euler angles :math:`\alpha, \beta, \gamma` assuming a :math:`zyz` convention which is
illustrated on the right side of the scheme below.


.. figure:: /img/coherent_formalisms.png
   :width: 600px
   :align: center

|


Initialising the Spatial Part  
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To facilitate handling of various formalisms as well as LAB and PAS representations for the spatial component of coherent interactions, we provide the class :class:`AnisotropicCoupling`. 
It can be initialised from LAB or PAS matrix or spherical tensor representations of the spatial interaction component or a subset of parameters :math:`\overline{a}, a_{xy}, a_{xz}, a_{yz}, \delta, \eta`.
If the interaction is initialised from a PAS representation, a set of euler-anngles  :math:`\alpha, \beta, \gamma` can be specified to also generate the LAB representation. When initialising
from a LAB representation the euler angles are automatically determined via matrix diagonalisation of the rank 2 component.  

The code below showcases the  initialisation of a generic anisotropic interaction with rank 0 and 2 components  from trace, anisotropy, asymmetry and euler angles.

.. code-block:: python

   import simos as sos 

   # Interaction parameters
   iso = sos.f2w(1e6) # 6 MHz
   delta = sos.f2w(0.5e6)
   eta = 0.5
   alpha = sos.deg2rad(20)
   beta = sos.deg2rad(40)
   gamma = sos.deg2rad(30)

   # Initialise from parameters
   A = sos.AnisotropicCoupling(iso = iso, delta = delta, eta = eta , euler = [alpha, beta, gamma])

We can now exctract the matrix and spherical tensor representations of PAS and LAB frames as well as the interaction parameters as illustrated below.

.. code-block:: python

   A_matrix_lab = A.mat("lab") # Matrix representation, LAB
   A_matrix_pas = A.mat("pas") # Matrix representation, PAS
   A_spher_lab = A.spher("lab") # Spherical tensor representation, LAB
   A_spher_pas = A.spher("pas") # Spherical tensor representation, PAS
   A_params = A.parameters() #  Interaction parameters can be accessesed as a dictionary 

Alternatively, the same interaction object could be obtained directly from matrix or spherical tensor representation as indicated below. 

.. code-block:: python

   A = sos.AnisotropicCoupling(mat = A_matrix_lab)
   A = sos.AnisotropicCoupling(mat = A_matrix_pas,  euler = [alpha, beta, gamma] )
   A = sos.AnisotropicCoupling(spher = A_spher_lab )
   A = sos.AnisotropicCoupling(spher = A_spher_pas, euler = [alpha, beta, gamma] )



Initialising the Interaction Hamiltonian 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 

The actual interaction Hamiltonian can be initiated with the :code:`interaction_hamiltonian()` method. For example,  a generic spin-spin coupling would be obtained as

.. code-block:: python

   S  = {"name": "S", "val": 1/2}
   I  = {"name": "I", "val": 1/2}
   s = sos.System([S,I])

   A = sos.AnisotropicCoupling(iso = iso, delta = delta, eta = eta , euler = [alpha, beta, gamma])
   Hcoupling = sos.interaction_hamiltonian(s, "S", A, spin2 = "I", frame = "lab", approx = "None")

The keyword argument :code:`approx` supports a series of common approximations  as summarized in the table below.


.. list-table:: Approximations
   :widths: 25 25 50
   :header-rows: 1

   * - Argument
     - :math:`\hat{H}^{\prime}\propto`
     - Assumes ... 
   * - :code:`full` / :code:`none`
     - :math:`\hat{H}`
     - no approximation made
   * - :code:`secular`
     - :math:`\hat{S}_{z}\hat{I}_{z}`
     - high-field approximation (A term of the dipolar alphabet)
   * - :code:`pseudosecular`
     - :math:`\hat{I}_{x}\hat{S}_{x} + \hat{S}_{y}\hat{I}_{y} +  \hat{S}_{z}\hat{I}_{z}`
     - high-field approximation for a pair of spins with the same Larmor frequency (homonuclear case; A+B term of the dipolar alphabet)
   * - :code:`hfi`
     - :math:`\hat{S}_{z}\hat{I}_{x} + \hat{S}_{z}\hat{I}_{y} + \hat{S}_{z}\hat{I}_{z}`
     - selective high-field approximation for the first spin of the pair , commonly applies for hyperfine interactions. 
.. list-table::


For spin-field interactions, the routine is called with a :code:`field` instead of the :code:`spin2` keyword arugment and only no or secular approximations
are applicable. 

.. code-block:: python

   Hfield = sos.interaction_hamiltonian(s, "S", A, spin2 = "I", field = [0,0,1e-3], approx = "None")



Specific Interaction Hamiltonians
---------------------------------

For a series of common interactions, we provide further high-level functions to facilitate generation of the Hamiltonian. All these methods utilize the generic 
:code:`interaction_hamiltonian()` routine in the background.


Zeeman Interaction
^^^^^^^^^^^^^^^^^^

Electron and nuclear spin is associated with a magnetic dipole moment and therefore sensitive to surrounding magnetic fields via the Zeeman interaction.  
The magnetic-field sensitivity of a given spin is defined by the :ref:`gyromagnetic ratio <GyroRatios>` of the nucleus or electron which is proportional
to the magnitude of the spin magnetic dipole moment. The interaction of nuclear or electronic spins with static magnetic fields is described by 
the Zeeman Hamiltonian

.. math::
   \hat{H}_{\mathrm{Z}} =  \vec{B} \cdot \mathbf{\delta} \cdot \vec{\hat{S}}

where :math:`\mathbf{\delta}` is the chemical shift tensor that accounts for the microscopically corrected Zeeman interaction
in an external magnetic field :math:`\vec{B}`. 

In isotropic media, the Hamiltonian simplifies to 

.. math::
   \hat{H}_{\mathrm{Z}}  =  \delta_{iso} \vec{B} \cdot \vec{\hat{S}}

with a scalar, isotropic chemical shift :math:`\delta_{\mathrm{iso}}`. In the simplest case  :math:`\delta_{\mathrm{iso}} = \gamma` , the gyromagnetic
ratio of the bare nucleus or electron.  SimOS provides a method :code:`zeeman_interaction()` to specify generic Zeeman interactions for individual spins.  The interaction strength, 
i.e. :math:`\mathbf{\delta}` can be provided as a scalar (for isotropic interactions), a 3x3 matrix, an instance of the :class:`AnisotropicCoupling` class or a tuple with 3 or 6 entries,
specifying the isotropic chemical shift :math:`\delta_{\mathrm{iso}}`, span :math:`\Omega`, skew :math:`\kappa` and, possibly, three euler angles in :math:`zyz` ordering. 
Span and skew are commonly utilized to parametrize the anisotropic chemical shift and are defined as

.. math::
   
   \Omega = \delta_{xx} - \delta_{zz}

.. math::
   
   \kappa = \frac{3 (\delta_{yy} - \delta_{iso})}{\Omega}

from the principal components of the chemical shift tensor assuming an ordering :math:`\delta_{xx} \geq \delta_{yy} \geq \delta_{zz}` of the principal components following the IUPAC convention for chemical shift.
The code below illustrates initialisation of a zeeman interaction for the isotropic case and for an anisotropic CSA, utilizing initialisation from CSA parameters.

.. code-block:: python


   S  = {"name": "S", "val": 1/2}
   I  = {"name": "I", "val": 1/2}
   s = sos.System([S,I])

   iso = sos.f2w(400e6)
   omega = sos.f2w(2e6)
   kappa = 0.4
   alpha = sos.deg2rad(30)
   beta = sos.deg2rad(20)
   gamma =sos.deg2rad(50)
   Hz_iso = sos.zeeman_interaction(s, "A", iso, [0,0, B0])
   Hz_aniso = sos.zeeman_interaction(s, "A", (iso, omega, kappa, alpha, beta, gamma), [0,0, B0])

Further,  a method :code:`auto_zeeman_interaction()` automatically calculates the isotropic Zeeman interactions for an entire spin system and 
thus requires that the isotopes are specified for all spin members of the quantum system.


RF Pulses and Spin Rotations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The effect of radiofrequency control pulses with well-defined phase and rotation angle is readily incorporated with our :code:`rot()` method which is further
documented in our section on :ref:`time propagation <TimePropagation>`. 


Dipolar Coupling
^^^^^^^^^^^^^^^^

The dipolar coupling between pairs of nuclear or electronic spins is a through-space magnetic interaction of their associated magnetic dipoles. 
The dipolar coupling Hamiltonian

.. math::
   \hat{H}_{\mathrm{dip}} = -\frac{\mu_0 \hbar \gamma_1 \gamma_2}{4 \pi}\frac{1}{r^3} \left(  3(\vec{S_1}\cdot \vec{r_u})(\vec{S_2}\cdot \vec{r_u}) - \vec{S_1}\cdot \vec{S_2}\right)   \\ 
   =  \vec{S_1}^T \cdot \mathbf{D} \cdot \vec{S_2}

can be calculated from the gyromagnetic ratios :math:`\gamma_1, \gamma_2` of the spins and the distance
and orientation of their connecting vector :math:`r` with respect to the external magnetic field. The resulting dipolar coupling tensor
:math:`\mathbf{D}` is apure rank-2 3x3 matrix which may be expressed in terms of the so-called "dipolar alphabet" with components:

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

SimOS provides two methods for purely dipolar couplings. The :code:`dipolar_spatial()` method only returns the spatial part, i.e. the dipolar coupling tensor, of the interaction while 
the :code:`dipolar_coupling()` method returns the complete dipolar coupling Hamiltonian. 


Zero-Field Splitting
^^^^^^^^^^^^^^^^^^^^


To incorporate zero-field splittings (ZFS)

.. math::
   \hat{H}_{\mathrm{ZFS}} =  \vec{\hat{S}} \cdot \mathbf{D} \cdot \vec{\hat{S}}

of electron spins :math:`\geq` 1 SimOS provides a method :code:`zfs_interaction()`. The spatial part :math:`\mathbf{D}` of the of the ZFS interaction can be specified as a 3x3 matrix, an instance of the :class:`AnisotropicCoupling` class
or via the parallel :math:`D` and :math:`E` and anti-parallel components of the ZFS tensor and (optionally) a set of euler angles. In the PAS of the interaction, the ZFS Hamiltonian as a function of :math:`D` and :math:`E` results as

.. math::
   \hat{H}_{\mathrm{ZFS}}=  D \left[ \hat{S}_z ^2 - \frac{1}{3} S(S+1) \right] + E \left[ \hat{S}_x^2 - \hat{S}_y^2 \right]


Quadrupole Interaction
^^^^^^^^^^^^^^^^^^^^^^

Nuclear spins with :math:`I \geq 1`  are quadrupolar and a non-spheric charge distribution of their nucleus translates into a finite electric quadrupole moment :math:`Q`. This quadrupole
moment couples to electric field gradient (EFG) :math:`\mathbf{V}` at the nuclear site, resulting in a nuclear quadrupole interaction

.. math::
   \hat{H}_{\mathrm{Q}} = \frac{e Q}{2I(2I-1)\hbar}\ \vec{I} \cdot \mathbf{V} \cdot \vec{I} .

SimOS provides method :code:`quad_interaction()` to facilitate incorporation of quadrupole interactions. The spatial part of the interaction, i.e.  :math:`\frac{e Q  \mathbf{V}}{2I(2I-1)\hbar}` can either be provided as a 3x3 matrix, as an instance of
the :class:`AnisotropicCoupling` or via the anisotropy :math:`\delta_{\mathrm{Q}}`, the asymmetry :math:`\eta_{\mathrm{Q}}` and, if PAS and LAB frames are not aligned, a set of euler angles. Parameters :math:`\delta_{\mathrm{Q}}` and :math:`\eta_{\mathrm{Q}}` are defined as follows:

.. math::
   \delta_{\mathrm{Q}} =  \frac{eQV_{ZZ}}{2I(2I-1)\hbar} =  \frac{C_{\mathrm{Q}}}{2I(2I-1)}

and 

.. math::
   \eta_{\mathrm{Q}} = \frac{V_{YY}- V_{XX}}{V_{ZZ}}

using an ordering :math:`|V_{ZZ}| \geq |V_{XX}| \geq |V_{YY}|` of the principal components of the EFG tensor and introducing the quadrupolar coupling constant :math:`C_{\mathrm{Q}}`, a common out put parameter of many quantum mechanical modelling packages. 


Syntax Reference
----------------

.. py:currentmodule:: simos.coherent
.. automodule:: simos.coherent
   :members: AnisotropicCoupling, interaction_hamiltonian, zeeman_interaction, auto_zeeman_interaction, dipolar_spatial, dipolar_coupling, zfs_interaction, quad_interaction
