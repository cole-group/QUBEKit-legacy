# QUBEKit: QUantum BEspoke FF toolKit


[QUBEKit](https://blogs.ncl.ac.uk/danielcole/qube-force-field/) is python based force field derivation toolkit that allows users to derive accurate molecular mechanics parameters directly from quantum mechanical calculations. 

Users who have used QUBEKit to derive any new force field parameters should cite the following papers:
* [Biomolecular Force Field Parameterization via Atoms-in-Molecule Electron Density Partitioning](https://pubs.acs.org/doi/abs/10.1021/acs.jctc.6b00027)
* [Harmonic Force Constants for Molecular Mechanics Force Fields via Hessian Matrix Projection](https://pubs.acs.org/doi/10.1021/acs.jctc.7b00785)

To install we recommend cloning the QUBEKit folder into a home directory and running the install script which will append the script to the search path:

    git clone git@github.com:cole-group/QuBeKit.git
    cd QUBEKit
    ./QUBE_install.sh
    
After relaunching the terminal you should now be able to use QUBEKit, try the following command to bring up the help documentation and check the installation (note tab completion should work):

    QUBEKit.py -h

## Requirements:
* [Anaconda3](https://www.anaconda.com/download/)
* [Biochemical and Organic Simulation System (BOSS)](http://zarbi.chem.yale.edu/software.html)
* [OpenMM](http://openmm.org/)
* [Gaussian09](http://gaussian.com/)
* [ONETEP](http://www.onetep.org/)
* [Matlab 2017](https://uk.mathworks.com/products/matlab.html)
### Python modules used:
* numpy
* argparse
* collections
* colorama
* matplotlib

## In Development

QUBEKit should currently be considered a work in progress. While it is stable we are constantly working to improve the code and increase the amount of compatible software. A user tutorial can be found on our Github [wiki](https://github.com/cole-group/QuBeKit/wiki) page. 
