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

   [user@local]~% cd ~; git clone https://github.com/spencerhongcornell/ESRC-Model.git

.. code-block:: bash

   [user@local]~% cd ~; git clone git@github.com:spencerhongcornell/ESRC-Model.git


Contributing/Extensions
------------------------------

If you would like to be a collaborator, first contact Spencer Hong (me) either through github or email and request permissions. 

This codebase is governed under the Apache License, so you are free to fork, modify, and distribute this code. 

For future extensions of this project, we recommend the following paths:

1) Investigating non-ideal mixture interactions in the absorber unit by using mixed fugacities model
2) Making the code more efficient by introducing vectorization (currently calculates 1000 parameters in 10 minutes)
3) Rewriting the code as a linear matrix solver problem
4) Trying different working fluids in the model
5) Building this cycle in ASPEN by modeling the absorber unit as a 2-unit process
6) Incorporating heat exchangers into the model (we assumed no heat exchangers between cold and hot fluids)

Documentation
------------------------------

The main.py file imports the helper functions from helper.py. Main.py and helper.py run calculations at 4 bar whle main1bar.py and helper1bar.py run calculations at 1 bar. The data is stored at /data/txy and /data/mixed folders as CSV files.

main.py will output a CSV file of coefficient of performances from temperature range 370K - 385K and 0 to 1 ammonia composition. The code will check the following things: heat of the evaporator is positive (the fact that this is the refrigeration box), the molar flow rates of stream 3 and 4 are positive, and the heat of the generator is positive (the fact that we are adding heat to provide refrigeration). 

What this means is that all of the COP provided by the model obey these conditions.

To check conservation of energy, return Qevap, Qflas, Qabs, and Qgen and then add them to a sum to see that net change in energy is 0 kJ/s.

NOTE: The "m" in the code refers to molar flow rates, not mass flow rates. The units of these molar flow rates are mol/s. Furthermore, all of the data in the /data/mixed CSV files are in kJ/mol.

.. code-block:: bash

   [user@local]~% cd ~; python main.py

