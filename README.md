# Gromacs Job Production/Analysis

This code aims at providing tools to easily generate gromacs simulation protocol and analyse results with parts based on already existing codes.

## Getting Started

The following instructions will get you a copy of the project up and running on your local machine.

### Prerequisites

The whole code is written in Python3 so you will need python3 package on your machine.

In addition to Python3 standard library, the following packages are needed for Production.py:

```
numpy

```
And for Analysis.py:

```
numpy
pandas
matplotlib
scipy
mdanalysis
fatslim
```

which you can install using [pip3](https://pip.pypa.io/en/latest/):

#### On Debian

```
sudo (or su -c) apt-get install python3-pip
```

And then for any system having pip3:

```
pip3 install --upgrade numpy pandas matplotlib scipy mdanalysis fatslim
```

### Installing

To get the code on your machine simply use the following command:

```
git clone https://github.com/fben94/GJB
```

