<p align="center">
  <img src="https://raw.githubusercontent.com/juanaldayparejo/acs-dist/main/docs/images/acs_logo.png" alt="Atmospheric Chemistry Suite Logo" width="400"/>
</p>

__________

Python package to process and analyse solar occultation measurements made with the Atmospheric Chemistry Suite aboard the ExoMars Trace Gas Orbiter.

This package is currently maintained by [Juan Alday](https://gapt.iaa.es/content/juan-alday).

If interested users are missing key points in the documentation, would appreciate seeing jupyter notebooks for certain purposes, or want to report issues, please do so by contacting us.

## Installation

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



