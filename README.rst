Mathematical Model of Einstein-Szilard Refrigeration Cycle 
==============================

This codebase is the implementation of mathematical modeling of ESRC in Python 3.2. The current version of the model accommodates a range of temperatures, two pressures, and any basis of ammonia composition in the generator unit.

The codebase requires the default Anaconda libraries (numpy, sympy, pandas, matplotlib). The sympy equation solver
is used to solve systems of equations. 

The mixed enthalpy/entropy data is pulled from AspenPlus using Peng-Robinson equation of state at various temperatures and pressures.

3D graphs and contour plots were created both in Python and in MATLAB. We, the hosts, recommend MATLAB for the 3D plots.


Installing
------------------------------

Currently installation involves cloning this repository.

.. code-block:: bash

   [user@local]~% cd ~; git clone https://github.com/spencerhongcornell/varnerthermof18.git

.. code-block:: bash

   [user@local]~% cd ~; git clone git@github.com:spencerhongcornell/varnerthermof18.git

Aftewards, open up *install.py* and adjust settings accordingly.  Then, simply run:

.. code-block:: bash

   [user@local]~% python install.py

and you are good to go with using Squid.

Contributing/Extensions
------------------------------

If you would like to be a collaborator, first contact Spencer Hong (me) either through github or email and request permissions. 

This codebase is governed under the Apache License, so you are free to fork, modify, and distribute this code. 

For future extensions of this project, we recommend the following paths:

1) Investigating non-ideal mixture interactions in the absorber unit by using mixed fugacities model
2) Making the code more efficient by introducing vectorization (currently calculates 1000 parameters in 10 minutes)
3) Rewriting the code as a linear matrix solver problem

Documentation
------------------------------

The main.py file imports the helper functions from helper.py. Main.py and helper.py run calculations at 4 bar whle main1bar.py and helper1bar.py run calculations at 1 bar. The data is stored at /data/txy and /data/mixed folders as CSV files.

