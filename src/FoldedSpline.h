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

#ifndef _Aru_FoldedSpline_h_
#define _Aru_FoldedSpline_h_

#include "Utils.h"
#include "Spline.h"
#include "VKernel.h"

#include <limits>
#include <cmath>

namespace Aru {

  namespace FoldedSpline {

    class BasisFunction {

    private:
      class Integrand {
      private:
        const FoldedSpline::BasisFunction& fParent;
        double fY;
      public:
        Integrand(const FoldedSpline::BasisFunction& parent)
          : fParent(parent), fY(0)
        {}

        void SetY(const double y) { fY = y; }

        double
        operator()(const double x)
          const
        {
          return fParent.fKernel(fY,x)*fParent.fBasisSpline(x);
        }
      };

      class Integrand2 {
      private:
        const FoldedSpline::BasisFunction& fParent;
        mutable Integrand fIntegrand;
      public:
        Integrand2(const FoldedSpline::BasisFunction& parent)
          : fParent(parent),
            fIntegrand(parent)
        {}

        double
        operator()(const double y)
          const
        {
          fIntegrand.SetY(y);
          return Integrate(fIntegrand,
            fParent.fBasisSpline.GetStart(),
            fParent.fBasisSpline.GetStop()
          );
        }
      };

      double
      FindRange(const double yref, const double ymax)
        const
      {
        const double eps = 1e-12;

        if (fKernel.ResolutionPdf(ymax,yref) > eps)
          return ymax;

        double yfar = ymax;
        double ynear = yref;
        size_t iter = 0;
        for (;iter<200;++iter)
        {
          const double y = 0.5*(yfar+ynear);
          if (std::fabs(y-ynear) <= eps*std::fabs(ynear))
            break;
          if (fKernel.ResolutionPdf(y,yref) > eps)
            ynear = y;
          else
            yfar = y;
        }
        if (iter == 200)
          throw std::runtime_error("find range did not converge");
        return ynear;
      }

      const VKernel& fKernel;
      const Spline::BasisFunction& fBasisSpline;
      mutable Integrand2 fIntegrand2;
      double fStart;
      double fStop;

    public:
      BasisFunction(const VKernel& kernel,
                    const Spline::Function1D& unfolded,
                    const size_t k)
        : fKernel(kernel),
          fBasisSpline(unfolded.Basis(k)),
          fIntegrand2(*this),
          fStart(0),
          fStop(0)
      {
        fStart = fKernel.Y(fBasisSpline.GetStart());
        fStop  = fKernel.Y(fBasisSpline.GetStop());
        if (fBasisSpline.GetStart() != unfolded.GetStart())
          fStart = FindRange(fStart,fKernel.Y(unfolded.GetStart()));
        if (fBasisSpline.GetStop() != unfolded.GetStop())
          fStop = FindRange(fStop,fKernel.Y(unfolded.GetStop()));
      }

      double
      operator()(const double y)
        const
      {
        if (y < GetStart() || y > GetStop())
          return 0.0;
        return fIntegrand2(y);
      }

      double
      Integral(const double a, const double b)
        const
      {
        double ya = a, yb = b;
        if (ya < GetStart()) ya = GetStart();
        if (yb > GetStop()) yb = GetStop();
        if (ya >= yb)
          return 0.0;
        else
          return Integrate(fIntegrand2,ya,yb);
      }

      double
      Integral()
        const
      {
        return Integrate(fIntegrand2,fStart,fStop);
      }

      double
      GetStart()
        const
      {
        return fStart;
      }

      double
      GetStop()
        const
      {
        return fStop;
      }

    };


    class Function1D {

    private:
      Spline::KnotVector fKnotVector;
      std::vector< FoldedSpline::BasisFunction* > fBasisVector;

    public:
      Function1D(const VKernel& kernel,
                 const Spline::Function1D& unfolded)
        : fKnotVector(unfolded.GetKnotVector().size()),
          fBasisVector(unfolded.GetSize())
      {
        for (size_t i=0,n=fKnotVector.size();i<n;++i)
          fKnotVector[i] = kernel.Y(unfolded.GetKnotVector()[i]);

        for (size_t k=0;k<fBasisVector.size();++k)
          fBasisVector[k] = new FoldedSpline::BasisFunction(kernel,unfolded,k);
      }

      ~Function1D()
      {
        for (size_t k=0;k<fBasisVector.size();++k)
          delete fBasisVector[k];
      }

      inline
      const FoldedSpline::BasisFunction&
      Basis(const size_t k)
        const
      {
        return *fBasisVector[k];
      }

      double
      operator()(const std::vector<double>& coefs, double y)
        const
      {
        double result = 0.0;
        for (size_t k=0,n=fBasisVector.size();k<n;++k)
          result += coefs[k]*(*fBasisVector[k])(y);
        return result;
      }

      size_t
      GetSize()
        const
      {
        return fBasisVector.size();
      }

      double
      GetStart(const size_t k = 0)
        const
      {
        return (*fBasisVector[k]).GetStart();
      }

      double
      GetStop(const size_t k = std::numeric_limits<size_t>::max())
        const
      {
        const size_t n = fBasisVector.size();
        return (*fBasisVector[ k<n ? k : n-1 ]).GetStop();
      }

      const Spline::KnotVector&
      GetKnotVector()
        const
      {
        return fKnotVector;
      }

    };

  } // NS FoldedSpline

} // NS Aru

#endif
