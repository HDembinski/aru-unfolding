#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import numpy as np
from ROOT import *

nt = TNtupleD("nt","","y")

np.random.seed(1)
data = 0.5+0.2*np.random.randn(200)

for y in data:
  nt.Fill(y,0) # second zero works around ROOT bug

f = TFile.Open("data.root","RECREATE")
nt.Write()
f.Close()
