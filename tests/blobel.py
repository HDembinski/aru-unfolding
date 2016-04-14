#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import numpy as np
from matplotlib import pyplot as plt
from utils import *
from aru import *
from kernel import BlobelKernel

kernel = BlobelKernel()

print "generate"
try:
  np.random.seed(int(sys.argv[1]))
except IndexError:
  np.random.seed(1)

nMC = 5000

def Model(x):
  b = (1.0,10.0,5.0)
  g = (2.0,0.2,0.2)
  c = (0.4,0.8,1.5)
  result = 0.0
  for k in xrange(3):
    result += b[k]*g[k]**2/((x-c[k])**2+g[k]**2)
  return result

rxs = []
while True:
  xy = np.random.rand(2)
  x = xy[0]*2
  y = xy[1]*12
  if y < Model(x):
    rxs.append(x)
    if len(rxs) == nMC:
      break

data = []
for x in rxs:
  y = kernel.Y(x)
  z = np.random.rand()
  if z < kernel.Efficiency(y):
    data.append(y + 0.1*np.random.randn())

xknots = np.linspace(0,2,20)

aru = ARU(kernel, ToVector(xknots))
coefs, cov = aru.Unfold(ToVector(data))
unfolded = aru.GetUnfoldedSpline()
folded = aru.GetFoldedSpline()
refCoefs = aru.GetReferenceCoefficients()

# draw everything

xs = np.linspace(xknots[0],xknots[-1],200)

from scipy.integrate import quad
ymodel    = nMC*Model(xs)/quad(Model,0,2)[0]
yprior    = np.vectorize(lambda x: unfolded(refCoefs,x))(xs)
yfitted   = np.vectorize(lambda x: folded  (coefs,x))(xs[xs<1.8])
yunfolded = np.vectorize(lambda x: unfolded(coefs,x))(xs)

plt.figure(1)

from matplotlib.patches import Polygon
poly = Polygon(MakeSigmaPatch(unfolded, coefs, cov, xknots), closed=True)
poly.set(ec="b",fc="b",alpha=0.1,fill=True,zorder=0)
plt.gca().add_patch(poly)

plt.plot(xs,ymodel,"g--",lw=2,label="$t(x)$",zorder=1)

plt.plot(xs,yprior,"y",lw=3,label="$g(x)$",zorder=1)

plt.plot(xs[xs<1.8],yfitted,"r",lw=1,label="$f(y)$",zorder=3)

plt.plot(xs,yunfolded,"b",lw=2,label="$b(x)$",zorder=2)

xh, w, werr = GetScaledHistogram(data,bins=20,range=(folded.GetStart(),folded.GetStop()))
plt.errorbar(xh,w,werr,fmt="ko",capsize=0,label="data",zorder=4)

plt.ylim(ymin=0)
plt.xlabel("x")
plt.ylabel(r"$\mathrm{d}N / \mathrm{d} x$")
plt.legend(loc="upper left").get_frame().set_visible(0)
plt.show()
