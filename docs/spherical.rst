.. _Spherical:

Spherical Tensor Notation
=========================


.. note::
    Currently, only a minimal documentation on spherical tensor notation and utilization in SimOS is available. 
    We are working on a more extensive version! Thank you for your patience.


Hamiltonian operators can be formulated in spherical tensor notation,  i.e. as products of spatial tensors
:math:`A_{l,k}` and spin tensor operators :math:`\hat{T}_{l,k}` of well-defined rank :math:`l` :

.. math::
    \hat{H} = \sum_l \sum_{k=-1}^l (-1)^k  A_{l,k} \hat{T}_{l,k}.

This notation is especially popular in modern solid-state NMR spectroscopy as it allows to assess 
the effects of rotational dynamics in classical and spin space (i.e. magic-angle spinning or rf irradiation) and resulting
resonance conditions with relative ease. For an extensive discussion of spherical tensor operator formalism 
and applications in ss-NMR, we refer our users to standard NMR literature.

SimOS core functionality does not utilize the spherical tensor formalism, however we provide a series of functions that facilitate usage of this 
notation.


Spatial Tensors 
----------------

The spatial component of common Hamiltonian operators of spin-systems are  3x3 matrices  and thus cartesian second-rank spatial tensors. 
The methods :code:`mat2spher()` and :code:`spher2mat()` allow to transform 3x3 matrices into spherical tensor notation and vice-versa.



Spin-Tensor Operators
---------------------

The method :code:`spherspin()` returns the spin tensor operators for a given spin member of a quantum
system while the :code:`spherbasis()` constructs a complete spin tensor operator basis for a given spin system. 



Syntax Reference
----------------

.. py:currentmodule:: simos.spherical
.. automodule:: simos.spherical
   :members: mat2spher, spher2mat, spherspin, spherbasis
