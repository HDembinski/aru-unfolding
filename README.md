What is ARU?
============

ARU is short for Automatic Regularized Unfolding (ARU). It is a non-parametric algorithm to unfold detector effects from one-dimensional data distributions. The unfolded solution minimizes the mean integrated squared error (MISE) of the fit with respect to the data. ARU uses a completely unbinned analysis and a powerful regularization based on the information gain (relative entropy) in the solution. It provides a smooth estimate of the true underlying distribution with an error band that has proper Frequentist coverage.

In contrast other algorihtms, ARU requires no manual tuning. In particular, the regularization parameter is picked automatically such that it minimizes the MISE. This means that ARU can be used as a black box algorithm, if necessary, and in a scenario where data needs to be unfolded automatically.

Installation
============

The installation uses cmake, version 2.8 or higher. Type:

	cmake . # check requirements and generate Makefile
	make    # actual compilation

There is no support for `make install`.

ARU relies on nlopt, SWIG, and ROOT. In addition, the Python
modules numpy and scipy need to be installed.

The required software can be downloaded here:
* cmake         - http://www.cmake.org
* nlopt         - http://ab-initio.mit.edu/wiki/index.php/Main_Page
* ROOT          - http://root.cern.ch/drupal
* SWIG          - http://www.swig.org
* numpy & scipy - http://www.scipy.org
* matplotlib    - http://matplotlib.sourceforge.net (optional)

The installation of matplotlib is optional, but the tests won't
work without it, neither will the automatic plotting of the solution.

Using ARU
=========

In order to use ARU, you need to
- adapt MyKernel.h
- supply a data file

You can then use unfold.py to unfold your data.
Usage: unfold.py <inputFile> <outputFile>

MyKernel.h needs to implement the detector kernel of your experiment,
have a look inside and into src/VKernel.h. At the very least you need to
implement the member function ResolutionPdf(...).

By default, unfold.py expects your raw data in a ROOT file with a single
TNtupleD inside which has the observations stored the first variable. The
output is then also written in a ROOT file. See the documentation inside
unfold.py for the format.

unfold.py can be adapted to your situation. You may want to restrict the range
of the unfolded variable, for example. You can also easily change the way how
unfold.py loads the data and writes the solution.

Tests
=====

You can run some tests to check the installation and ARU's performance.
The tests are located in the tests folder and need to be called with from
the root directory of the installation (where CMakeLists.txt is located).

make_toy_data.py:
Generate a trivial toy data set, ready to be used with unfold.py.

blobel.py:
Unfolds the example found in
V. Blobel, Unfolding methods in high energy physics experiments,
Proceedings of the 1984 CERN School of Computing (1984).
An optional integer argument may be supplied, which is used as the random
seed.

regfit.py:
Unfolds another example, two Gaussian peaks, one very narrow, one wide.
You need to supply two integer arguments. First is the random seed, second
is the number of MC events to unfold.

Licence and citation
====================

The ARU source code is released under GPL v3, see the file COPYRIGHT.
Copyright holder is the Karlsruhe Institute of Technology (KIT).
Author contact is Dr. Hans Dembinski.

If you use ARU in your publication, please cite our paper:

H.P. Dembinski and M. Roth,
An algorithm for automatic unfolding of one-dimensional data distributions,
Nuclear Instruments and Methods in Physics Research A 729 (2013) 410-416.

The latest version of ARU is available from:
https://github.com/HDembinski/aru-unfolding

Older versions of ARU's source code can be found on hepforge:
http://aru.hepforge.org
