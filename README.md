<p align="center">
  <img src="https://raw.githubusercontent.com/juanaldayparejo/acs-dist/main/docs/images/acs_logo.png" alt="Atmospheric Chemistry Suite Logo" width="400"/>
</p>

__________

Python package to process and analyse solar occultation measurements made with the Atmospheric Chemistry Suite aboard the ExoMars Trace Gas Orbiter.

This package is currently maintained by [Juan Alday](https://gapt.iaa.es/content/juan-alday).

If interested users are missing key points in the documentation, would appreciate seeing jupyter notebooks for certain purposes, or want to report issues, please do so by contacting us.

## Installation

To install this package, please follow the next steps:

### 1. Installing the latest version of the code from GitHub

The latest version of code has to be downloaded from [Github](https://github.com/juanaldayparejo/acs-dist.git) under a [GNU General Public License v3](LICENSE). To do so, type in the command window:

```bash
git clone https://github.com/juanaldayparejo/acs-dist.git
```

Before installing the library, we recommend users to create and load a new Python [virtual environment](https://docs.python.org/3/library/venv.html) for a clean install:

```bash
python -m venv name_of_virtual_environment/
source name_of_virtual_environment/bin/activate
```

Then move into the package directory:

```bash
cd acs-dist
```

Finally, we need to install the library. Given that this is a highly dynamic package were new additions are frequently introduced, we recommend installing the package but keeping it editable by typing:

```bash
pip install --editable .
```

This will install **acs**, but with the ability to update any changes made to the code. In addition, it will install all the required libraries it depends on.

### 2. Installing the Mars Climate Database (optional but recommended)

Some of the functionality included in this Python package relies on the extraction of vertical profiles of the Martian atmosphere from the [Mars Climate Database](https://www-mars.lmd.jussieu.fr/mars/access.html). 

To enable this functionality, the users must download the database and compile the Fortran functions that enable the extraction of the profiles from the files. The full version of the MCD can be downloaded from this [link](https://www-mars.lmd.jussieu.fr/MCD_pro/mcd_pro.html).

Once downloaded, we need to compile the fortran functions of the MCD using the *acs-dist/acs/mcd/fmcd_gfortran.sh* compilation file stored within the acs-dist distribution. Before compiling it, we must modify with following fields:

- NETCDF: Path to the directory where the NETCDF libraries are compiled.
- wheremcd: Path to the Mars Climate Database.
- version: Version of the MCD that will be used.

Once these fields are modified according to our local settings, we can compile it by writing in the terminal:

```bash
source fmcd_gfortran.sh
```

This will create a file starting with *fmcd* and with a *.so* extension. This will be automatically loaded when importing the **acs** library into our Python programs.


Instead, if it is desirable to disable the MCD functionality, then the user must modify the *acs-dist/acs/__init__.py* file by commenting out the line importing the code from the *mcd* folder.
