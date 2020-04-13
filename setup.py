import os
from setuptools import setup

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    file_path = os.path.join(os.path.dirname(__file__), fname)
    return open(file_path, encoding="utf-8").read()


setup(
    name="opensir",
    version="1.0.0",
    author="Open-SIR community",
    author_email="open.sir.project@gmail.com",
    description=(
        "Open-SIR is an Open Source Python project for modelling"
        "pandemics and infectious diseases using Compartmental"
        "Models"
    ),
    license="MIT",
    keywords="SIR SIRX pandemics modelling",
    url="http://github.com/open-sir/open-sir",
    packages=["opensir", "opensir.models"],
    scripts=["opensir-cli"],
    install_requires=["scipy", "sklearn", "numpy", "toml",],
    long_description=read("README.md"),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Utilities",
        "License :: OSI Approved :: BSD License",
    ],
)
