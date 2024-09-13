#!/usr/bin/env python

from setuptools import setup

longdescription = """
# SimOS
SimOS (**SIM**ulation of **O**ptically-adressable **S**pins) is a library for the simulation of open quantum systems consisting of spins, electronic levels and combinations thereoff. It can simulate the spin dynamics of conventional magnetic or electron paramagnetic resonance, but is further capable to simulate optically adressable spins or purely optical systems. ODMR combines electron spin resonance with optical measurements and is a highly sensitive technique to study systems possessing spin-dependent fluorescence. In our examples section, we feature two prototypical examples of optically adressable spins - the nitrogen vacancy center in diamond and photogenerated spin-correlated radical pairs.

Modelling the dynamics of open quantum system can be a non-trivial task. Their incoherent interaction with the system environment interferes with their coherent time evolution and an accurate simulation of their dynamics requires a quantum master equation. SimOS provides an interface for the construction of operators and super-operators of arbitrarily complex systems of spins and electronic levels and facilitates simulations of their dynamics on various levels of theory. Pythons popular numpy, scipy, qutip and sympy libraries are readily integrated and can be chosen as backends. We try to minimize high-level functions and keep the style of the simulations as close to a pen-and paper style as possible.

Our main focus lies on the QME in Lindblad form, for which we provide various engines for computationally efficient time propagation. In addition spatial dynamics such as rotational diffusion, linear flow or magic angle spinning can be simulated with our Fokker-Planck framework.
"""

setup(name='simos',
      version='0.1.1',
      description='The simos libraray provides tools for spin simulations in Python.',
      long_description=longdescription,
      long_description_content_type='text/markdown',
      author='Konstantin Herb',
      author_email='science@rashbw.de',
      url='https://simos.kherb.io',
      license='all rights reserved',
       packages=['simos','simos.backends','simos.constants','simos.utils','simos.systems'],
      install_requires=['numpy','scipy'],
      extras_require={
            'recommended': ['qutip','sympy','matplotlib','IPython'],
            'full': ['qutip','sympy','matplotlib','numba','IPython','tkinter'],
            'docs': ['sphinx','sphinx_rtd_theme','simos[recommended]'],
            'test': ['pytest','pytest-cov','simos[recommended]']
      },
      classifiers=[
            'Development Status :: 3 - Alpha',
            'Environment :: GPU',
            'Intended Audience :: Science/Research',
            'Intended Audience :: Developers',
            'License :: Other/Proprietary License',
            'Operating System :: Microsoft :: Windows',
            'Operating System :: MacOS',
            'Operating System :: POSIX',
            'Operating System :: Unix',
            'Operating System :: OS Independent',
            'Programming Language :: Python :: 3',
            'Topic :: Education',
            'Topic :: Scientific/Engineering',
            'Topic :: Scientific/Engineering :: Physics',
            'Topic :: Scientific/Engineering :: Chemistry',
      ],
      platforms='any',
      zip_safe=False)
