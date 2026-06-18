from setuptools import setup, find_packages
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text(encoding="utf-8")

setup(
    name="acs",
    version="1.0.0",
    author="Juan Alday",
    description="Python package to process solar occultation data from the Atmospheric Chemistry Suite aboard the ExoMars Trace Gas Orbiter",
    long_description=long_description,
    long_description_content_type="text/markdown",  # important for Markdown rendering
    url="https://github.com/juanaldayparejo/acs-dist.git",
    project_urls={
        "Source": "https://github.com/juanaldayparejo/acs-dist.git",
    },
    #packages=["archnemesis"],
    packages=find_packages(), 
    install_requires=[
      'numpy',
      'matplotlib',
      'numba>=0.57.0',
      'scipy',
      'joblib',
      'h5py',
    ],
    extras_require={
        'docs': ['sphinx', 'sphinx_rtd_theme'],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.10',
)
