Open-SIR API
============

SIR Model
*********

Most epidemic models share a common approach on modelling the spread of a disease. 
The susceptible-infectious-removed (SIR) model is a simple deterministic compartmental 
model to predict disease spread. An objective population is divided in three groups: 
the susceptible (:math:`S`), the infected (:math:`I`) and the recovered or removed (:math:`R`). 
These quantities enter the model as fractions of the total population :math:`P`.

.. math::
   S = \frac{\text{Number of susceptible individuals}}{\text{Population size}},

.. math::
   I = \frac{\text{Number of infected individuals}}{\text{Population size}},

.. math::
   R = \frac{\text{Number of recovered or removed individuals}}{\text{Population size}},

As a pandemics infects and kills much more quickly than human natural rates of birth and death, 
the population size is assumed constant  except for the individuals that recover or die. 
Hence, :math:`S+I+R=P/P=1`. The pandemics dynamics is modelled  as a system of ordinary differential equations
which governs the rate of change at which the percentage of susceptible, infected 
and recovered/removed individuals in a population evolve.

The number of possible transmissions is proportional to the number of interactions between 
the susceptible and infected compartments, :math:`S \times I`:

.. math::
   \frac{dS}{dt} = -\alpha SI,

Where :math:`\alpha` / :math:`[\text{time}]^{-1}` is the transmission rate of the process which quantifies how 
many of the interactions between susceptible and infected populations yield to new infections per day.

The population of infected individuals will increase with new infections and decrease with recovered or removed people. 

.. math::
   \frac{dI}{dt} = \alpha S I  - \beta I,
   
.. math::
   \frac{dR}{dt} = \beta I,

Where :math:`\beta` is the percentage of the infected population that is removed from the transmission process per day.

The infectious period, :math:`T_I` / :math:`[\text{time}]` , is defined as the reciprocal of the removal rate: 

.. math::
   T_I=\frac{1}{\beta}.

In early stages of the infection, the number of infected people is much lower than the susceptible population. 
Hence, :math:`S \approx 1` making :math:`dI/dt` linear and the system has the analytical solution
:math:`I(t) = I_0 \exp (\alpha - \beta)t`.


.. autoclass:: opensir.models.SIR
   :members:
   :inherited-members:

SIR-X Model
***********

The SIR-X model extends the SIR model adding two parameters:
the quarantine rate :math:`\kappa` and the containement rate 
:math:`\kappa_0`. This extension allows the model to capture the
"decrease" of susceptible population owing containment and quarantine
measures.

.. autoclass:: opensir.models.SIRX
   :members:
   :inherited-members:
