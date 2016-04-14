/*
This file is part of ARU, a software to unfold detector effects
from data distributions.

Copyright (C) 2011 Karlsruhe Instiute of Technology,
  contact Hans Dembinski <hans.dembinski@gmail.com>

ARU is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ARU is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with ARU. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _MyKernel_h_
#define _MyKernel_h_

#include "VKernel.h"
#include <cmath>

class MyKernel : public Aru::VKernel {
private:
  const double fSigma;
public:

  MyKernel() :
    fSigma(0.1)
  {
  }

  double
  ResolutionPdf(const double yObserved, const double yExpected)
    const
  {
    const double z = (yObserved-yExpected)/fSigma;
    const double c = 1.0/std::sqrt(2*M_PI)/fSigma;
    return c*std::exp(-0.5*z*z);
  }

};

#endif

