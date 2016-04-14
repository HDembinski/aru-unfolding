#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import numpy as np
from matplotlib import pyplot as plt
from pyik.mplext import auto_adjust
from utils import *
from aru import *
from kernel import *

print "generate"
np.random.seed(int(sys.argv[1]))

nMC = int(sys.argv[2])

norms = (0.3,0.7)
means = (0.3,0.7)
sigmas = (0.1,0.05)
kernelSigma = 0.1

# ~norms = (1.,20.)
# ~means = (0.4,0.6)
# ~sigmas = (0.05,0.2)
# ~kernelSigma = 0.1

model = ToyModel(norms,means,sigmas)
data = model.rvs(nMC)

kernel = GaussianKernel(kernelSigma)
data += kernelSigma*np.random.randn(nMC)

xknots = np.linspace(0,1,20)

aru = ARU(kernel, ToVector(xknots))
c_coefs, c_cov = aru.Unfold(ToVector(data))
unfolded = aru.GetUnfoldedSpline()
folded = aru.GetFoldedSpline()
c_refCoefs = aru.GetReferenceCoefficients()

modelNorm = np.sum((folded.GetStart()<=data)*(data<=folded.GetStop()))

csum = 0.0
for k in xrange(unfolded.GetSize()):
  csum += c_coefs[k]*unfolded.Basis(k).Integral()

print "integral comparison",modelNorm, csum

# draw everything

xs = np.linspace(0,1,200)

ysmodel    = nMC*model(xs)
ysprior    = np.vectorize(lambda x: unfolded(c_refCoefs,x))(xs)
ysfitted   = np.vectorize(lambda x: folded  (c_coefs,x))(xs)
ysunfolded = np.vectorize(lambda x: unfolded(c_coefs,x))(xs)

from matplotlib.patches import Polygon
poly = Polygon(MakeSigmaPatch(unfolded, c_coefs, Extract(c_cov), xknots), closed=True)
poly.set(ec="b",fc="b",alpha=0.1,fill=True,zorder=0)
plt.gca().add_patch(poly)

plt.plot(xs,ysmodel,"g--",lw=2,label="$t(x)$",zorder=1)

plt.plot(xs,ysprior,"y",lw=3,label="$g(x)$",zorder=1)

plt.plot(xs,ysfitted,"r",lw=1,label="$f(y)$",zorder=3)

plt.plot(xs,ysunfolded,"b",lw=2,label="$b(x)$",zorder=2)

xh, w, werr = GetScaledHistogram(data,bins=20,range=(0,1))
plt.errorbar(xh,w,werr,fmt="ko",capsize=0,label="data",zorder=4)

# ~print "noreg fit"
# ~c_coefs_noreg = std.vector("double")(nBasis,1)
# ~objective.SetWeight(0)
# ~Fit(objective,c_coefs_noreg)
# ~yunreg    = np.vectorize(lambda xx: unfolded(c_coefs_noreg,xx))(x)
# ~plt.plot(x,yunreg,"b:",lw=1,label="$b_\mathrm{w=0}(x)$",zorder=1)

plt.ylim(ymin=0)
plt.xlabel("x")
plt.ylabel(r"$\mathrm{d}N / \mathrm{d} x$")
plt.legend(loc="upper left").get_frame().set_visible(0)

plt.show()
