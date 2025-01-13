#!/usr/bin/env python

from setuptools import setup

longdescription = """
# SimOS
SimOS (**SIM**ulation of **O**ptically-adressable **S**pins) is a versatile library for the simulation of open quantum systems consisting of spins, electronic levels and combinations thereoff. Originally developed for optically detected magnetic resonance (ODMR) experiments, SimOS is also capable of performing general nuclear magnetic resonance (NMR) and electron paramagnetic resonance (EPR/ESR) simulations. Its main focus is to simulate optically adressable spins or purely optical systems. 

SimOS provides an interface for constructing operators and super-operators for complex systems of spins and electronic levels, facilitating simulations of their dynamics across various theoretical levels. The library integrates seamlessly with popular Python libraries such as numpy, scipy, qutip, and sympy, which can be chosen as backends. 

Our primary focus is on the quantum master equation (QME) in Lindblad form, for which we offer various engines for computationally efficient time propagation. Additionally, SimOS supports the simulation of spatial dynamics, including rotational diffusion, linear flow, and magic angle spinning, through our Fokker-Planck framework.

SimOS aims to provide a pen-and-paper style of simulation, minimizing high-level functions to keep the syntax clean and mathematically driven. While SimOS can perform a wide range of simulations, it was designed to be a Python-based alternative to other established simulation tools, offering similar capabilities in a more accessible and flexible environment.
"""

setup(name='simos',
      version='0.2.2',
      description='Spin simulations in Python (NMR, EPR/ESR as well as ODMR).',
      long_description=longdescription,
      long_description_content_type='text/markdown',
      author='Konstantin Herb',
      author_email='science@rashbw.de',
      url='https://simos.kherb.io',
      license='GNU Lesser General Public License v3 or later (LGPLv3+)',
      packages=['simos','simos.backends','simos.constants','simos.utils','simos.systems'],
      install_requires=['numpy>=1.25.0','scipy>=1.12.0'],
      extras_require={
            'recommended': ['qutip>=5.0.0','sympy>=1.11','matplotlib','IPython'],
            'full': ['qutip','sympy','matplotlib','numba>=0.58.0','IPython','tkinter'],
            'docs': ['sphinx','sphinx_rtd_theme','simos[recommended]'],
            'test': ['pytest','pytest-cov','simos[full]']
      },
      classifiers=[
            'Development Status :: 3 - Alpha',
            'Environment :: GPU',
            'Intended Audience :: Science/Research',
            'Intended Audience :: Developers',
            'License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)',
            'Operating System :: Microsoft :: Windows',
            'Operating System :: MacOS',
            'Operating System :: POSIX',
            'Operating System :: Unix',
            'Operating System :: OS Independent',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.9',
            'Programming Language :: Python :: 3.10',
            'Programming Language :: Python :: 3.11',
            'Programming Language :: Python :: 3.12',
            'Topic :: Education',
            'Topic :: Scientific/Engineering',
            'Topic :: Scientific/Engineering :: Physics',
            'Topic :: Scientific/Engineering :: Chemistry',
      ],
      platforms='any',
      zip_safe=False)
