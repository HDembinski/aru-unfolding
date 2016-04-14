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

#ifndef _Aru_VKernel_h_
#define _Aru_VKernel_h_

namespace Aru {

  class VKernel {
  public:

    double
    operator()(const double yObserved, const double x)
      const
    {
      const double y = Y(x);
      return Efficiency(yObserved)*ResolutionPdf(yObserved, y);
    };

    virtual
    double
    ResolutionPdf(const double /* yObserved */, const double /* yTruth */)
      const = 0;

    // efficiency is assumed to be one unless specified otherwise
    virtual
    double
    Efficiency(const double /* yObserved */)
      const
    {
      return 1.0;
    }

    // x -> y = identity unless specified otherwise
    virtual
    double
    Y(const double x)
      const
    {
      return x;
    }

    // x -> y = identity unless specified otherwise
    virtual
    double
    X(const double y)
      const
    {
      return y;
    }

  };

} // NS Aru

#endif

