#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This file is part of ARU, a software to unfold detector effects
# from data distributions.
#
# Copyright (C) 2011 Karlsruhe Instiute of Technology,
#   contact Hans Dembinski <hans.dembinski@gmail.com>
# 
# ARU is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Foobar is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

import os
import sys
import numpy as np
from math import *
from utils import *
from aru import *
from kernel import MyKernel

def GetData(inputFileName):
  from ROOT import TFile
  ntuple = TFile.Open(inputFileName).GetListOfKeys()[0].ReadObj()
  data = np.empty(ntuple.GetEntries())
  for i in xrange(ntuple.GetEntries()):
    ntuple.GetEntry(i)
    data[i] = ntuple.GetArgs()[0]
  return data

# define input/output, minimal checks
try:
  inputFilename  = sys.argv[1]
  outputFilename = sys.argv[2]
  if not os.path.exists(inputFilename):
    raise SystemExit("Input file %s does not exist."%inputFilename)
  if os.path.exists(outputFilename):
    raise SystemExit("Output file %s does exist, I don't dare to overwrite."%outputFilename)
except IndexError:
  raise SystemExit("Wrong number of arguments. Usage: unfold.py <inputFile> <outputFile>")

# inputFile contains the raw data as a TNtupleD with a single branch
print "reading data..."
data = GetData(inputFilename)

# replace implementation of MyKernel by a proper kernel
kernel = MyKernel()

# x-range to be unfolded, may be set to other values
xmin = np.min(data)
xmax = np.max(data)

# customize number of knots (use a sufficiently large number) and their positions
nKnots = 20
xknots = np.linspace(xmin,xmax,nKnots)

# main program, encapsulated in ARU object
print "unfolding..."
aru = ARU(kernel, ToVector(xknots))
c_coefs, c_cov = aru.Unfold(ToVector(data))
unfolded = aru.GetUnfoldedSpline()
folded = aru.GetFoldedSpline()
c_refCoefs = aru.GetReferenceCoefficients()

coefs = Extract(c_coefs)
cov = Extract(c_cov)

# writing the output file with contents
# - TNtupleD with variables
#   x        : sample value
#   yref     : final reference distribution (for checking)
#   yfitted  : folded distribution (the actual fit to the data)
#   yunfolded: unfolded distribution (the solution)
#   ysigma   : uncertainty of yunfolded
# - TVectorD, coefficients of the solution
# - TMatrixDSym, covariance matrix of the solution
print "Writing result..."

class Sigma(object):
  def __init__(self,unfolded, coefs, cov):
    self.unfolded = unfolded
    self.coefs = coefs
    self.cov = cov
  def __call__(self,x):
    n = self.unfolded.GetSize()
    bs = [ unfolded.Basis(k)(x) for k in xrange(n) ]
    return np.sqrt( np.dot(np.dot(bs,cov),bs) )

xs = np.linspace(xmin,xmax,200)

yref      = np.vectorize(lambda x: unfolded(c_refCoefs,x))(xs)
yfitted   = np.vectorize(lambda x: folded  (c_coefs,x))(xs)
yunfolded = np.vectorize(lambda x: unfolded(c_coefs,x))(xs)
sigma = Sigma(unfolded, coefs, cov)
ysigma    = np.vectorize(sigma)(xs)

from ROOT import TNtupleD, TFile, TMatrixDSym, TVectorD
nt = TNtupleD("nt","","x:yref:yfitted:yunfolded:ysigma")
for i in xrange(len(xs)):
  nt.Fill(xs[i],yref[i],yfitted[i],yunfolded[i],ysigma[i])

rootCoefs = TVectorD(len(coefs))
rootCov = TMatrixDSym(len(coefs))
for i in xrange(len(coefs)):
  rootCoefs[i] = coefs[i]
  for j in xrange(len(coefs)):
    rootCov[i][j] = cov[i][j]
f = TFile.Open(outputFilename,"RECREATE")
nt.Write()
rootCoefs.Write()
rootCov.Write()
f.Close()

try:
  # draw everything, if matplotlib is installed
  from matplotlib import pyplot as plt

  from matplotlib.patches import Polygon
  poly = Polygon(MakeSigmaPatch(unfolded, c_coefs, cov, xknots), closed=True)
  poly.set(ec="b",fc="b",alpha=0.1,fill=True,zorder=0)
  plt.gca().add_patch(poly)

  plt.plot(xs,yref,"y",lw=2,label="$g(x)$",zorder=1)

  plt.plot(xs,yfitted,"r",lw=1,label="$f(y)$",zorder=2)

  plt.plot(xs,yunfolded,"b",lw=2,label="$b(x)$",zorder=2)

  xh, w, werr = GetScaledHistogram(data,bins=20,range=(xmin,xmax))
  plt.errorbar(xh,w,werr,fmt="ko",capsize=0,label="data",zorder=3)

  plt.ylim(ymin=0)
  plt.xlim(xmin,xmax)
  plt.xlabel("x")
  plt.ylabel(r"$N_\mathrm{data} \times p.d.f.$")
  plt.legend(loc="upper left")
  plt.show()

except ImportError:
  print "Install matplotlib to get plots of the solution."
