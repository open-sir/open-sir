Getting Started
===============

Dependencies
************

* Python >= 3.2
* Numpy
* Sklearn
* Scipy

Installation
************
Open-SIR and all its dependencies can be installed from the repository using
``pip``:
::

    git clone https://github.com/open-sir/open-sir.git
    cd open-sir
    pip install .

In order to unintall Open-SIR simply execute:
::

    pip uninstall opensir

Usage example
*************

Command line interface
######################

It's possible to run the model using the CLI:
::

    usage: opensir-cli [-h] [-m {sir,sirx}] [-t TIME] [-f FILE] [-s] [-d DELIMITER]

    Run SIR or SIR-X model given a set of initial conditions and model parameters.

    optional arguments:
      -h, --help            show this help message and exit
      -m {sir,sirx}, --model {sir,sirx}
      -t TIME, --time TIME  Number of days for the simulation
      -f FILE, --file FILE  Name of the input file. If missing, the input file is read from STDIN
      -s, --suppress_header
                            Suppress CSV header
      -d DELIMITER, --delimiter DELIMITER
                            CSV delimiter

The input file is a TOML file with the following format
::

    [initial_conds]
    <cond> = value

    [parameters]
    <param> = value

For example, You can use Open-SIR to create a 6 days prediction of the number
of susceptible (S), infected (I) and removed (R) population.  The initial
conditions represent Ealing data as of 04/04/2020. The parameters provide a
prediction in the hypothetical case that no lockdown would be taking place.
::

    [initial_conds]
    S0 = 341555
    I0 = 445
    R0 = 0

    [parameters]
    alpha = 0.95
    beta = 0.38

Then, it's possible to run the model with ``T=6`` days:
::

    opensir-cli --model sir --time 6 --file input_file.txt

Or reading the input file from STDIN:
::

    cat input_file | opensir-cli --model sir --time 6

The output of opensir-cli is a .csv file with the output of the model.

.. note:: *Note*: On Windows, the CLI must be run from Powershell or any bash 
    shell such as `Git BASH <https://gitforwindows.org/>`_

Python API
##########

You can replicate the predictions of the CLI with the following python script:

.. code-block::

    from models import SIR
    my_sir = SIR() # Initialize an empty SIR model
    params = [0.95, 0.38] # Define model parameters (alpha, beta)
    w0 = [341555, 445, 0] # Define initial conditions (S0, I0, R0)
    my_sir.set_params(p=params, initial_conds=w0) # Set model parameters
    n_days = 6 # Define the amount of days to predict
    my_sir.solve(n_days, n_days+1) # Call model.solve functions
    sol = my_sir.fetch() # Fetch model solution

Try the Jupyter Notebook
########################

Open and run the
`Jupyter Notebook <https://github.com/open-sir/open-sir/blob/master/SIR.ipynb>`_
to:

* Get an overview of the SIR model
* Explore case studies

And learn how the API can be used to:

* Build compartmental models
* Fit parameters to existing data 
* Predict susceptible, infected and removed population
* Calculate confidence intervals of the predictions
