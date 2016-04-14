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

import numpy as np
from math import *
from scipy.stats import norm

def Extract(rootobj):
  """
  Extracts data from ROOT objects and returns corresponding Numpy arrays.

  Extract can be used with a string argument. In this case, the string is
  interpreted as a filename of a ROOT file. The ROOT file is opened and
  Extract is applied on each object inside the file. The result is returned
  as a dictionary, mapping the key names to the outputs.

  Authors
  -------
  Hans Dembinski <hans.dembinski@gmail.com>
  """

  import ROOT
  import numpy as np
  from aru import vector_double, matrix_double

  if isinstance(rootobj, vector_double):
    a = np.empty(rootobj.size())
    for i in xrange(len(a)):
      a[i] = rootobj[i]
    return a

  if isinstance(rootobj, matrix_double):
    a = np.empty((len(rootobj),len(rootobj[0])))
    for i in xrange(a.shape[0]):
      for j in xrange(a.shape[1]):
        a[i,j] = rootobj[i][j]
    return a

  if isinstance(rootobj, TMatrixDSym):
    a = np.empty((rootobj.GetNcols(),rootobj.GetNrows()))
    for i in xrange(rootobj.GetNcols()):
      for j in xrange(rootobj.GetNrows()):
        a[i,j] = rootobj(i,j)
    return a

  if isinstance(rootobj, ROOT.std.vector("double")):
    a = np.empty(rootobj.size())
    for i in xrange(len(a)):
      a[i] = rootobj[i]
    return a

  if isinstance(rootobj, ROOT.std.vector("unsigned long")):
    a = np.empty(rootobj.size(),dtype=np.uint)
    for i in xrange(len(a)):
      a[i] = rootobj[i]
    return a

  if isinstance(rootobj, ROOT.std.vector("int")):
    a = np.empty(rootobj.size(),dtype=np.int)
    for i in xrange(len(a)):
      a[i] = rootobj[i]
    return a

  if isinstance(rootobj, str):
    f = ROOT.TFile.Open(rootobj)
    d = {}
    for key in f.GetListOfKeys():
      d[key.GetName()] = Extract(key.ReadObj())
    return d

  if isinstance(rootobj, ROOT.TVectorD):
    n = rootobj.GetNoElements()
    a = np.empty(n)
    for i in xrange(n):
      a[i] = rootobj[i]
    return a

  if isinstance(rootobj, ROOT.TMatrixD) or isinstance(rootobj, ROOT.TMatrixDSym):
    a = np.empty((rootobj.GetNcols(),rootobj.GetNrows()))
    for i in xrange(rootobj.GetNcols()):
      for j in xrange(rootobj.GetNrows()):
        a[i,j] = rootobj(i,j)
    return a

  # this has to be come before "if isinstance(rootobj,ROOT.TH1)"
  if isinstance(rootobj,ROOT.TH2):
    nx = rootobj.GetXaxis().GetNbins()
    ny = rootobj.GetYaxis().GetNbins()
    x  = np.empty(nx+1)
    y  = np.empty(ny+1)
    z  = np.empty((nx,ny))
    ez = np.empty((nx,ny))
    for ix in xrange(nx):
      x[ix] = rootobj.GetXaxis().GetBinLowEdge(ix+1)
      for iy in xrange(ny):
        y[iy] = rootobj.GetYaxis().GetBinLowEdge(iy+1)
        z[ix,iy] = rootobj.GetBinContent(ix+1,iy+1)
        ez[ix,iy] = rootobj.GetBinError(ix+1,iy+1)
      x[nx] = rootobj.GetXaxis().GetBinLowEdge(nx+1)
    y[ny] = rootobj.GetYaxis().GetBinLowEdge(ny+1)
    return x,y,z,ez
  
  if isinstance(rootobj,ROOT.TH1):
    n  = rootobj.GetXaxis().GetNbins()
    x  = np.empty(n+1)
    y  = np.empty(n)
    ey = np.empty(n)
    for i in xrange(n):
      x[i]  = rootobj.GetXaxis().GetBinLowEdge(i+1)
      y[i]  = rootobj.GetBinContent(i+1)
      ey[i] = rootobj.GetBinError(i+1)
    x[n] = rootobj.GetXaxis().GetBinLowEdge(n+1)
    return x,y,ey

  raise ValueError("Extract cannot handle object %s yet, please take a look at pyik.rootext and implement it"%(type(rootobj)))

def ToVector(array):
  from aru import vector_double
  n = len(array)
  res = vector_double(n)
  for i in xrange(n):
    res[i] = array[i]
  return res

def centers(x):
  """
  Computes the centers of an array of bin edges.

  Parameters
  ----------
  x: array-like
    A 1-d array containing lower bin edges.

  Returns
  -------
  c: array of dtype float
    Returns the centers of the bins.
  hw: array of dtype float
    Returns the half-width of the bins.

  Examples
  --------
  >>> c,hw = centers([0.0,1.0,2.0])
  >>> print c, hw
  [ 0.5  1.5] [ 0.5  0.5]

  Authors
  -------
  Hans Dembinski <hans.dembinski@gmail.com>
  """

  n = len(x)-1
  c = np.empty(n)
  hw = np.empty(n)
  for i in xrange(len(x)-1):
    hw[i] = 0.5*(x[i+1]-x[i])
    c[i] = x[i] + hw[i]
  return c,hw

class ToyModel:
  def __init__(self,norms,mus,sigmas):
    self.norms  = np.array(norms)
    self.norms /= np.sum(self.norms)
    self.mus    = mus
    self.sigmas = sigmas

  def __call__(self,x):
    result = 0.0
    n = len(self.norms)
    for i in xrange(n):
      result += self.norms[i]*norm(self.mus[i],self.sigmas[i]).pdf(x)
    return result

  def rvs(self, nEvents):
    n = len(self.norms)
    ps = np.zeros(n)
    for i in xrange(n):
      ps[i:] += self.norms[i]
    ps /= np.sum(self.norms)
    counts = np.zeros(n,dtype=int)
    for i in xrange(nEvents):
      for j,p in enumerate(ps):
        uniran = np.random.rand()
        if uniran < p:
          counts[j] += 1
          break
    rvs = np.array([])
    for i in xrange(n):
      rvs = np.append(rvs,norm(self.mus[i],self.sigmas[i]).rvs(counts[i]))
    return rvs

def MinimumMeanSquaredErrorFit(obj):
  from aru import vector_double, matrix_double, MatrixAdd
  nBasis = obj.GetBasisSize()

  foldedCache = np.empty((nBasis, obj.GetUnbinnedSize()))
  foldedMuCache = np.empty((nBasis, obj.GetBinnedSize()))
  for k in xrange(nBasis):
    for i in xrange(obj.GetUnbinnedSize()):
      foldedCache[k,i] = obj.GetFoldedCache(k,i)
    for i in xrange(obj.GetBinnedSize()):
      foldedMuCache[k,i] = obj.GetFoldedMuCache(k,i)

  fmatrix = obj.GetFMatrix()

  # start values
  c_coefs = vector_double(nBasis,1.0)
  c_cov = matrix_double(nBasis,vector_double(nBasis,0.0))
  c_cov2 = matrix_double(nBasis,vector_double(nBasis,0.0))

  def MeanSquaredError(sqrtw):
    w = sqrtw*sqrtw

    obj.SetWeight(w)
    obj.Minimize(c_coefs)
    obj.CovarianceStat(c_cov, c_coefs)
    obj.UpdateReference(c_coefs)

    if w > 0:
      obj.CovarianceReg(c_cov2, c_coefs)
      MatrixAdd(c_cov, c_cov2)

    result = 0.0
    if w != 0.0:
      for i in xrange(nBasis):
        for j in xrange(nBasis):
          result += (c_cov[i][j]+c_coefs[i]*c_coefs[j])*fmatrix[i][j]

      for i in xrange(foldedCache.shape[1]):
        folded = 0.0
        for k in xrange(nBasis):
          folded += c_coefs[k]*foldedCache[k,i]
        result -= 2*folded

      for ib in xrange(foldedMuCache.shape[1]):
        hw = obj.GetHistogramWeight(ib)
        if hw>0:
          foldedMu = 0.0
          for k in xrange(nBasis):
            foldedMu += c_coefs[k]*foldedMuCache[k,ib]
          result -= 2*hw*foldedMu/obj.GetHistogramStep(ib)

    print "iterating (%.5e): current weight = %g"%(result,w)

    return result

  from scipy.optimize import brent
  sqrtw = brent(MeanSquaredError, brack=(0,1), tol=1e-2)
  w = sqrtw*sqrtw
  obj.SetWeight(w)
  obj.Minimize(c_coefs)
  return c_coefs

def GetScaledHistogram(data,bins,range=None):
  ws, xedges = np.histogram(data,bins=bins,range=range)
  if type(bins) is not int:
    bins = len(xedges)-1
  xs = centers(xedges)[0]
  werrs = np.sqrt(ws)
  factor = bins/float(xedges[-1]-xedges[0]) # due to binning
  ws = np.asarray(ws,np.float)
  ws *= factor
  werrs *= factor
  return xs, ws, werrs

def MakeSigmaPatch(basis, coef, cov, xknot):

  nover = 10

  xs = []
  for i in xrange(len(xknot)-1):
    x0 = xknot[i]
    x1 = xknot[i+1]

    dx = (x1-x0)/nover

    for j in xrange(nover-1):
      xs.append(x0+j*dx)
  xs.append(xknot[-1])

  def CalcSigma(basis,cov,x):
    n = basis.GetSize()
    s2 = 0.0
    for l in xrange(n):
      for m in xrange(n):
        s2 += basis.Basis(l)(x)*basis.Basis(m)(x)*cov[l][m]
    if s2 < 0:
      s2 = 0
    return sqrt(s2)

  y1 = [ basis(coef,x) for x in xs ]

  y2 = [ CalcSigma(basis,cov,x) for x in xs ]

  xy = []

  n = len(xs)
  # upper
  for i in xrange(n):
    xy.append([ xs[i], y1[i]+y2[i] ])
  # lower
  for i in xrange(n-1,0,-1):
    xy.append([ xs[i], y1[i]-y2[i] ])

  return np.array(xy)

def JackknifeShifts(obj,c_coefs):
  from aru import TMatrixDSym
  nBasis = obj.GetBasisSize()

  grad0 = Extract(obj.Gradient(c_coefs))
  c_hesse1 = TMatrixDSym(nBasis)
  obj.HesseStat(c_hesse1, c_coefs)
  c_hesse2 = TMatrixDSym(nBasis)
  obj.HesseReg(c_hesse2, c_coefs)
  hesse0 = Extract(c_hesse1) + obj.GetWeight()*Extract(c_hesse2)

  foldedCache = obj.GetFoldedCache
  foldedMuCache = obj.GetFoldedMuCache

  shifts = []
  coefs = Extract(c_coefs)
  for i in xrange(obj.GetUnbinnedSize()):
    folded = 0.0
    for k in xrange(nBasis):
      folded += coefs[k]*foldedCache(k,i)

    grad = np.copy(grad0)
    hesse = np.copy(hesse0)
    for k in xrange(nBasis):
      grad[k] += foldedCache(k,i)/folded
      for l in xrange(k,nBasis):
        delta = foldedCache(k,i)*foldedCache(l,i)/folded**2
        hesse[k,l] -= delta
        if k!=l:
          hesse[l,k] -= delta

    u,s,vt = np.linalg.svd(hesse)
    hesseInv = np.zeros((nBasis,nBasis))
    for k in xrange(nBasis):
      for l in xrange(k,nBasis):
        for m in xrange(nBasis):
          hesseInv[k,l] = vt[k,m]*u[l,m]/s[m]
        if k!=l:
          hesseInv[l,k] = hesseInv[k,l]
    shifts.append(-np.dot(hesseInv,grad))

  for ib in xrange(obj.GetBinnedSize()):
    hw = obj.GetHistogramWeight(ib)

    if hw==0:continue

    foldedMu = 0.0
    for k in xrange(nBasis):
      foldedMu += coefs[k]*foldedMuCache(k,ib)

    grad = np.copy(grad0)
    hesse = np.copy(hesse0)
    for k in xrange(nBasis):
      grad[k] += hw*foldedMuCache(k,i)/foldedMu
      for l in xrange(k,nBasis):
        delta = hw*foldedMuCache(k,i)*foldedMuCache(l,i)/foldedMu**2
        hesse[k,l] -= delta
        if k!=l:
          hesse[l,k] -= delta

    u,s,vt = np.linalg.svd(hesse)
    hesseInv = np.zeros((nBasis,nBasis))
    for k in xrange(nBasis):
      for l in xrange(k,nBasis):
        for m in xrange(nBasis):
          hesseInv[k,l] = vt[k,m]*u[l,m]/s[m]
        if k!=l:
          hesseInv[l,k] = hesseInv[k,l]
    shifts.append(-np.dot(hesseInv,grad))

  return shifts

def GetEmpiricalDensity(kernel,unfolded,data):
  xa = unfolded.GetStart()
  xb = unfolded.GetStop()
  result = 0.0
  for y in data:
    x = kernel.X(y)
    if xa <= x and x < xb:
      result += 1.0/kernel.Efficiency(y)
  return result/(xb-xa)

class ARU(object):
  """
  Class with a simple interface for unfolding data with ARU.

  Parameters
  ----------
  kernel: VKernel
    Instance of the detector kernel, derived from Aru::VKernel.
  xknots: std.vector("double")
    STL vector of doubles containing the knot positions.

  Author
  ------
  Hans Dembinski <hans.dembinski@gmail.com>
  """

  def __init__(self, kernel, xknots):
    from aru import SplineFunction1D,FoldedSplineFunction1D,KnotVector
    self.kernel = kernel
    self.unfolded = SplineFunction1D(KnotVector(ToVector(xknots)))
    self.folded = FoldedSplineFunction1D(kernel, self.unfolded)
    self.referenceCoefs = None
    self.weight = None

  def Unfold(self, data, nEventsForBinning=100):
    """
    Perform the unfolding.

    Parameters
    ----------
    data: std.vector("double")
      STL vector of doubles containing the observations.

    Returns
    -------
    coefs: std.vector("double")
      STL vector of doubles containing the spline coefficients (solution).
    cov: TMatrixDSym
      Covariance matrix of the spline coefficients (solution).
    """

    from aru import ObjectiveFunction, matrix_double, vector_double, MatrixAdd
    nBasis = self.unfolded.GetSize()

    print "creating objective function..."
    objective = ObjectiveFunction(data, self.kernel, self.folded, self.unfolded, 8, nEventsForBinning)

    print "fitting..."
    coefs = MinimumMeanSquaredErrorFit(objective)
    self.refCoefs = vector_double(objective.GetReferenceCoefficients())
    self.weight = objective.GetWeight()
    cov = matrix_double(nBasis,vector_double(nBasis,0.0))
    cov2 = matrix_double(nBasis,vector_double(nBasis,0.0))
    objective.CovarianceStat(cov, coefs)
    objective.CovarianceReg(cov2, coefs)
    MatrixAdd(cov,cov2)
    return coefs, cov

  def GetUnfoldedSpline(self):
    return self.unfolded

  def GetFoldedSpline(self):
    return self.folded

  def GetReferenceCoefficients(self):
    if self.refCoefs is None:
      raise StandardError("no reference computed yet")
    return self.refCoefs

  def GetWeight(self):
    if self.weight is None:
      raise StandardError("no weight computed yet")
    return self.weight
