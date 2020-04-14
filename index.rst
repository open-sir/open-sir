.. Open-SIR documentation master file, created by
   sphinx-quickstart on Sun Apr  5 12:09:27 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Open-SIR's documentation!
====================================

Open-SIR is an Open Source Python project for modelling pandemics and
infectious diseases using Compartmental Models, such as the widely used
`Susceptible-Infected-Removed (SIR) model
<http://rocs.hu-berlin.de/corona/docs/forecast/model/#classic-sir-dynamics>`_

The current stage of the software is *Alpha*.

Features
========

- Model the dynamics of infectious diseases
- Parameter fitting
- Calculation of confidence intervals
- CLI for interfacing with non Python environments (Bash, Node.JS, Matlab, etc).

So far, Open-SIR provides an implementation of the SIR model and the novel
`SIR-X model, developed by Maier and Dirk
<https://science.sciencemag.org/content/early/2020/04/07/science.abb4557.full>`_
from the `Robert Koch Institut
<http://rocs.hu-berlin.de/corona/docs/forecast/model/#sir-x-dynamics-outbreaks-with-temporally-increasing-interventions>`_

.. toctree::
   :glob:
   :maxdepth: 2
   :caption: Contents:

   doc/getting-started
   doc/API
   doc/SIR

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
