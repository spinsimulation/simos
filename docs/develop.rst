
Developer Guide
===============
This guide is intended for developers who want to contribute to the development of SimOS 
or modify its components. If you are new to Python, we recommend that you first read the 
:ref:`getting started guide <Getting Started>`. We are happy to welcome new contributors 
to the project at any time and from any field of research where SimOS could contribute. 
If you have any questions, please do not hesitate to contact us. 


Contribution guidelines
-----------------------

Documentation is important! 

Please use expressive and meaningful variable names. 
This makes it easier for all who are using the code and try to understand what's happening. 
In addition, please comment you implementations. When you use constants, approximations 
or a special implementation scheme, please consider including a DOI of a publication from
 which you retrieved your information.

In general, please follow the recommendations given in the Python Enhancement 
Proposals, especially `PEP 8 <https://www.python.org/dev/peps/pep-0008/>`_ and `PEP 257 <https://www.python.org/dev/peps/pep-0257/>`_. 


Install as a Developer
----------------------

clone the git repository to you local machine. Then install the
package with a symlink installation, so that changes to the source files will
be immediately available to other users of the package on our system:

.. code-block:: python

    pip install -e .


Add a Backend
-------------
To add a new backend, you need to create a backend file in the backends folder. Use the existing backends numpy or scipy.sparse as a template. Backend functions that are not preceeded by an underscore are considered public and should be implemented. 
Most crucially, you need to implement a datatype class that inherits from the base class of the new backend datatype. This class should implement the basic operations that are required for the user to interact with the backend (e.g. expm, trace, etc.). 