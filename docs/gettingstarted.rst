Getting Started
=======================

What we Do
----------


SimOS is a library for the simulation of open quantum systems consisting of spins, electronic levels and combinations thereoff. 
It can simulate the spin dynamics of conventional magnetic or electron paramagnetic resonance, but is further capable to simulate optically adressable spins or purely optical systems.
ODMR combines electron spin resonance with optical measurements and is a highly sensitive technique to study systems possessing spin-dependent fluorescence. In our examples section, 
we feature two prototypical examples of optically adressable spins - the nitrogen vacancy center in diamond and photogenerated spin-correlated radical pairs.

Modelling the dynamics of open quantum system can be a non-trivial task. Their incoherent interaction with the system environment interferes with
their coherent time evolution and an accurate simulation of their dynamics requires a quantum master equation. SimOS provides an interface for the construction of operators and super-operators of arbitrarily complex systems of spins and electronic levels and facilitates simulations of their dynamics on various levels of theory. Pythons
popular numpy, scipy, qutip and sympy libraries are readily integrated and can be chosen as backends.
We try to minimize high-level functions and keep the style of the simulations as close to a pen-and paper style as possible. 

Our main focus lies on the QME in Lindblad form, for which we provide various engines for computationally efficient time propagation. 
In addition spatial dynamics such as rotational diffusion, linear flow or magic angle spinning can be simulated with our Fokker-Planck framework. 


.. _Installation:

Installation Guide
------------------
The SimOS library is available at PyPI, the Python Package Index. It can be installed via pip, the Python package installer. As the backend structure is modular, packages for the backends have to be installed separately. To install the recommended set of dependencies, run the following command:

.. code-block:: bash

    pip install simos[recommended]

To install the full set of dependencies, run the following command:

.. code-block:: bash

    pip install simos[full]

This will install the recommended set of dependencies, the parament GPU integrator and the numba JIT compiler. Note that the NUMBA support is still experimental and might not work on all systems.

To install the minimal set of dependencies, run the following command:

.. code-block:: bash

    pip install simos

Note, that if you install the minimal set of dependencies, executing the functions with the default backend method "qutip" will raise an error and you have to specify the backend method explicitly.

To install in developer mode (e.g. after cloning the Github repository manually), please refer to the :ref:`developer guide <Developer Guide>`. You will need this type of installation if you want to contribute to the development of SimOS or modify its components.


Python beginners
^^^^^^^^^^^^^^^^	
If you are completely new to Python, we recommend that you install the Anaconda distribution, which includes Python, the Jupyter notebook environment and many scientific libraries. You can download Anaconda from the following link: 

https://www.anaconda.com/products/distribution

We generally recommend that you use a distinct Python environment for SimOS, e.g. a virtual environment or a conda environment. This will prevent conflicts with other Python packages that you might have installed on your system. After installation, open the Anaconda prompt and create a new environment with the following command:

.. code-block:: bash

    conda create --name spin python=3.12

This will create a new environment called "spin" with Python 3.12. To activate the environment, run the following command:

.. code-block:: bash

    conda activate spin

Now you can install SimOS in this environment with the following command:

.. code-block:: bash

    pip install simos[recommended]

This will automatically install the recommended set of dependencies. You can now start the Jupyter notebook environment with the following command:

.. code-block:: bash

    jupyter notebook

If you prefer a non-browser-based environment, you can install e.g. Visual Studio Code, which is a popular code editor with Python support. You can download Visual Studio Code from the following link:

https://code.visualstudio.com/

.. _Virtual:

SimOS Virtual Lab
------------------
Due to the amazing work of the Pyodide team, we are able to provide a virtual lab for SimOS. This virtual lab is a Jupyter notebook environment running in your browser, which allows you to run SimOS without installing it on your local machine. Nonetheless, all code execution is done locally in your browser. This is possible due to the WebAssembly technology, which allows to run Python code in the browser. The virtual lab is available at the following link:

.. raw:: html

    <a href="https://simos.kherb.io/virtual/lab/index.html?path=Welcome.ipynb" target="_blank" style="padding: 0.5em 3em; background-color: grey; color: black; text-decoration:underline; border-radius:0.3em;">Start SimOS Virtual Lab</a><br/>&nbsp;


Please note that the virtual lab is an order of magnitude slower than running SimOS natively on your machine. The initial import can take up to one minute. It is intended for educational purposes and quick testing of the library. For more complex simulations, we recommend to install SimOS on your local machine. 