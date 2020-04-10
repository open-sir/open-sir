# open-sir

Open-SIR is an Open Source Python project for modelling pandemics and infectious diseases using Compartmental Models, such as the widely used [Susceptible-Infected-Removed (SIR) model](http://rocs.hu-berlin.de/corona/docs/forecast/model/#classic-sir-dynamics). 
The current stage of the software is *Alpha*.

## Features
- Model the dynamics of infectious diseases
- Parameter fitting
- Calculation of confidence intervals
- CLI for interfacing with non Python environments (Bash, Node.JS, Matlab, etc).

So far, Open-SIR provides an implementation of the SIR model and the novel [SIR-X model, developed by Maier and Dirk](https://science.sciencemag.org/content/early/2020/04/07/science.abb4557.full) from the [Robert Koch Institut](http://rocs.hu-berlin.de/corona/docs/forecast/model/#sir-x-dynamics-outbreaks-with-temporally-increasing-interventions).

### Dependencies

* Python >= 3.2
* Numpy
* Sklearn
* Scipy

## Installation
Open-SIR and all its dependencies can be installed from the repository using
`pip`:

```
git clone https://github.com/open-sir/open-sir.git
cd open-sir
pip install .
```

In order to unintall Open-SIR simply execute:
```
pip uninstall opensir
```

## Usage example

You can use Open-SIR to create a 6 days prediction of the number of susceptible (S), infected (I) and removed (R) population. 
The initial conditions '-i' represent Ealing data as of 04/04/2020. The parameters '-p' provide a prediction in the hypothetical 
case that no lockdown would be taking place.

### Command line interface

It's possible to run the model using the CLI:

```
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
```

The input file is a TOML file with the following format

```
[initial_conds]
<cond> = value

[parameters]
<paran> = value
```

For instance, the following data is valid for running the SIR model with
alpha=0.95, beta=0.38 and initial conditions S0=341555, I0=445 and R0=0

```
[initial_conds]
S0 = 341555
I0 = 445
R0 = 0

[parameters]
alpha = 0.95
beta = 0.38
```

Then, it's possible to run the model with T=6 days:

```
opensir-cli --model sir --time 6 --file input_file.txt
```

Or reading the input file from STDIN:

```
cat input_file | opensir-cli --model sir --time 6
```

The output of opensir-cli is a .csv file with the output of the model.

*Note: On Windows, the CLI must be run from Powershell or any bash shell such as [Git BASH](https://gitforwindows.org/)*

### Python API

You can replicate the predictions of the CLI with the following python script:
```python
from models import SIR
my_sir = SIR() # Initialize an empty SIR model
params = [0.95, 0.38] # Define model parameters
w0 = [341555, 445, 0] # Define initial conditions
my_sir.set_params(p=params, initial_conds=w0) # Set model parameters
n_days = 6 # Define the amount of days to predict
my_sir.solve(n_days, n_days+1) # Call model.solve functions
sol = my_sir.fetch() # Fetch model solution
```

### Try the Jupyter Notebook

Open and run the Jupyter Notebook [SIR.ipynb](SIR.ipynb) to:
* Get an overview of the SIR model
* Explore case studies

And learn how the API can be used to:

* Build compartmental models
* Fit parameters to existing data 
* Predict susceptible, infected and removed population
* Calculate confidence intervals of the predictions

## Development Setup

Open-SIR Development Environment uses
[Pipenv](https://pipenv.pypa.io/en/latest/) to automatically create a virtual
environment and manage python packages. The python packages required by
Open-SIR are listed in the [Pipfile](Pipfile).

After cloning the repository, change the current directory to the repository
via `cd open-sir` and automatically install the environment from the Pipfile
using Pipenv:
```
pipenv install
```
Next, activate the Pipenv shell:
```
pipenv shell
```
You can run the following command to check that the installation succeeded.
```
pipenv run start -i input_file.txt -t 6
```

### Build documentation

Build the Sphinx documentation out of the open-sir Pipenv environment.
```
pipenv run doc
```

### Running the tests

Test the package running the tests out of the open-sir Pipenv environment.
```
pipenv run test
```

The documentation is placed under `_build/html`.

### Coding style tests

This project uses [Pylint](https://www.pylint.org/) and [Black
19.10b0](https://black.readthedocs.io/en/stable/) to ensure consistent coding
practices, and enforced by CircleCI. Non-default parameters are available in
[.pylintrc](.pylintrc).

## Authors

* **[José Álamos](https://github.com/jia200x)** - [RIOT-OS](https://github.com/RIOT-OS)
* **[Felipe Huerta](https://github.com/felipehuerta17)** - [PhD Student](https://www.imperial.ac.uk/people/f.huerta-perez17) at [Imperial College London](https://github.com/ImperialCollegeLondon)
* **[Sebastián Salata](https://github.com/sasalatart)** - Software Engineer - Full Stack

See also the list of [contributors](https://github.com/open-sir/open-sir/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

## Acknowledgements

* [Robert Koch Institut](https://www.rki.de/EN/Home/homepage_node.html) for the clear explanation of SIR and SIR-X models.
