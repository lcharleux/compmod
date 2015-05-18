.. CompMod documentation master file, created by
   sphinx-quickstart on Thu Nov 27 15:56:27 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

CompMod Documentation
===================================

Installation can be performed in many ways, here a two:
  
* The right way:
  
.. code-block:: bash

   pip install git+https://github.com/lcharleux/compmod.git

Under Windows systems, installing the Python(x,y) bundle allows this procedure to be performed easily.

* If you are contributing to the module, you can just clone the repository:
    
.. code-block:: bash

   git clone https://github.com/lcharleux/compmod.git   

And remember to add the abapy/abapy directory to your ``PYTHONPATH``. For example, the following code can be used under Linux (in ``.bashrc`` or ``.profile``):

.. code-block:: bash

  export PYTHONPATH=$PYTHONPATH:yourpath/compmod


Contents:

.. toctree::
   :maxdepth: 2

   models
   materials
   distributions
   rheology


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

