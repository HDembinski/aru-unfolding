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

#ifndef _BlobelKernel_h_
#define _BlobelKernel_h_

#include "Utils.h"
#include "VKernel.h"
#include <cmath>

class BlobelKernel : public Aru::VKernel {
private:
  const double fSigma;

public:

  BlobelKernel()
    : fSigma(0.1)
  {}

  double
  ResolutionPdf(const double yObserved, const double y)
    const
  {
    const double z = (yObserved-y)/fSigma;
    const double c = 1.0/(std::sqrt(2*M_PI)*fSigma);
    return c*std::exp(-0.5*z*z);
  }

  double
  Efficiency(const double y)
    const
  {
    return 1.0 - 0.5*Aru::Sqr(y-1.0);
  }

  double
  X(const double y)
    const
  {
    return 10.0 - 2.0*std::sqrt(5.0)*std::sqrt(5.0-y);
  }

  double
  Y(const double x)
    const
  {
    return x - 0.05*x*x;
  }

};

#endif

