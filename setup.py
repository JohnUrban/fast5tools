import sys

from glob import glob
from setuptools import setup


if sys.version_info[0] != 2 or sys.version_info[1] < 7:
    sys.exit("Sorry, only Python 2.7+ supported.")

setup(
    name="fast5tools",
    version="0.4",
    description="Tools for working with fast5 files (nanopore data).",
    author="John Urban",
    author_email="mr.john.urban@gmail.com", 
    url="https://github.com/JohnUrban/fast5tools",
    install_requires=[
        "biopython",
        "h5py",
        "matplotlib",
        "numpy",
    ],
    packages=["fast5tools", "tests"],
    scripts=glob("bin/*"),
    test_suite="tests",
    zip_safe = False,
    include_package_data=True, 
    classifiers=[ 
        'Development Status :: 4 - Beta', 
        'Intended Audience :: Science/Research', 
        'License :: OSI Approved :: GNU General Public License (GPL)', 
        'Topic :: Scientific/Engineering :: Bio-Informatics']
)
