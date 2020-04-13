""" Python setup configuration file for Open-SIR"""
from setuptools import setup

setup(
    name="opensir",
    version="1.0.0",
    author="Open-SIR community",
    author_email="open.sir.project@gmail.com",
    description=("Open-SIR is an Open Source Python project for modelling"
                 " pandemics and infectious diseases using Compartmental"
                 " Models"),
    license="MIT",
    keywords="SIR SIRX pandemics modelling",
    url="http://github.com/open-sir/open-sir",
    packages=['opensir', 'opensir.models'],
    scripts=["opensir-cli"],
    install_requires=[
        'scipy',
        'sklearn',
        'numpy',
        'toml',
    ],
    long_description=("Open-SIR is an Open Source Python project for"
                      " modelling epidemics like COVID-19"
                      " using Compartmental Models, such as the widely used"
                      " Susceptible-Infected-Removed (SIR model)"
                      " Features:"
                      " Model the dynamics of infectious diseases, parameter fitting,"
                      " calculation of confidence intervals,"
                      " CLI for interfacing with non Python environments such as"
                      " Bash, Node.JS, Matlab, etc."),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Utilities",
        "License :: OSI Approved :: BSD License",
    ],
)
