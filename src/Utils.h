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

#ifndef _Aru_Utils_h_
#define _Aru_Utils_h_

#include <cmath>
#include <vector>
#include <stdexcept>
#include <string>
#include <sstream>

namespace Aru {

  template<typename ANumber>
  inline ANumber Sqr(ANumber t) { return t*t; }

  template<typename AnArray>
  std::string ToString(unsigned int n, const AnArray& a)
  {
    std::ostringstream msg;
    msg << "[";
    for (unsigned int i=0;i<n;++i)
    {
      msg << " "; msg.precision(3); msg << a[i];
    }
    msg << " ]";
    return msg.str();
  }

  namespace OfflineIntegrator {

    // based on Offline/Utilities/Math/Integrator.h
    // from the Pierre Auger Observatory;
    // LICENSE.Auger applies to everything inside
    // the current namespace

    inline
    double
    PolynomialInterpolation(const unsigned int n,
                            const double px[],
                            const double py[],
                            const double x,
                            double& dy)
    {
      if (!n)
        return 0;

      double minDist = std::fabs(x - px[0]);
      unsigned int k = 0;
      double c[n];
      double d[n];
      c[0] = py[0];
      d[0] = py[0];
      double oldx = px[0];
      for (unsigned int i = 1; i < n; ++i) {
        const double newx = px[i];
        oldx = newx;
        const double dist = std::fabs(x - newx);
        if (dist < minDist) {
          minDist = dist;
          k = i;
        }
        c[i] = d[i] = py[i];
      }

      double y = py[k];

      for (unsigned int m = 1; m < n; ++m) {
        const unsigned int nm = n - m;
        for (unsigned int i = 0; i < nm; ++i) {
          const double xa = px[i];
          const double xb = px[i+m];
          const double dx = xb - xa;
          const double cd = (c[i+1] - d[i]) / dx;
          c[i] = (x-xa)*cd;
          d[i] = (x-xb)*cd;
        }
        dy = (2*k < nm) ? c[k] : d[--k];
        y += dy;
      }

      return y;
    }

    /// average of a function represented with equidistant boxes
    template<class Functor>
    inline
    double
    GetBoxAverage(const Functor& functor, const double start, const double step, const int n)
    {
      double sum = 0;
      for (int i = 0; i < n; ++i) {
        const double x = start + i*step;
        sum += functor(x);
      }
      return sum/n;
    }

    template<class Functor>
    inline
    double
    GetTrapezoidalAverage(const Functor& functor,
                          const double previousApproximation,
                          const double a, const double delta, const int level)
    {
      const int n = 1 << (level - 1);
      const double step = delta/n;
      return 0.5*(previousApproximation +
                  GetBoxAverage(functor, a+0.5*step, step, n));
    }

    template<class Functor>
    double
    Integrate(const Functor& functor, const double a, const double b,
              const double accuracy = 1e-8,
              const size_t order = 5, const size_t maxIterations = 15)
    {
      const double delta = b - a;

      std::vector<double> tInt(maxIterations+1);
      std::vector<double> tStep(maxIterations+1);
      double oldInt = tInt[0] = GetBoxAverage(functor, a, delta, 2);
      tStep[0] = delta;
      size_t i = 1;
      for (; i < order; ++i) {
        oldInt = tInt[i] = GetTrapezoidalAverage(functor, oldInt, a, delta, i);
        tStep[i] = delta/(1 << (2*i));  // h^2
      }

      const double eps = 0.5*accuracy;
      double extrapolation = 0;
      for (; i <= maxIterations; ++i) {
        oldInt = tInt[i] = GetTrapezoidalAverage(functor, oldInt, a, delta, i);
        tStep[i] = delta/(1 << (2*i));  // extrapolate h^2 dependence for h->0
        const int first = i - order + 1;
        double residual;
        extrapolation =
          PolynomialInterpolation(order, &tStep[first], &tInt[first], 0, residual);
        if (std::isnan(delta*extrapolation))
          throw std::runtime_error("got nan during integration");
        if (std::fabs(residual) < eps*std::fabs(extrapolation))
          return delta*extrapolation;
      }

      if (std::fabs(extrapolation) > accuracy)
        throw std::runtime_error("maximal level reached in integration, result still inaccurate");

      return delta*extrapolation;
    }

  } // NS OfflineIntegrator

  using OfflineIntegrator::Integrate;

}

#endif

