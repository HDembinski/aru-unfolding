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

#ifndef _Aru_Spline_h_
#define _Aru_Spline_h_

#include <cmath>
#include <vector>
#include <string>
#include <stdexcept>
#include <iostream>
#include <cfloat>
#include <algorithm>
#include <limits>

namespace Aru {

  namespace Spline {

    // std::vector<double> of positions with virtual extention at both ends
    class KnotVector : public std::vector<double> {
    public:

      KnotVector() {}

      KnotVector(const size_t n)
        : std::vector<double>(n)
      {
        if (n<2)
          throw std::runtime_error("size of KnotVector < 2");
      }

      KnotVector(const std::vector<double>& data)
        : std::vector<double>(data)
      {
        const size_t n = size();
        if (n<2)
          throw std::runtime_error("size of KnotVector < 2");
        for (size_t i=0;i<n-1;++i)
          if ((*this)[i]>(*this)[i+1])
            throw std::runtime_error("knots are not in ascending order");
      }

      double
      operator()(const int i)
        const
      {
        const size_t n = size();

        // virtual knots
        if (i<0)
        {
          const double dx = (*this)[1]-(*this)[0];
          return (*this)[0]+i*dx;
        }

        if (i>int(n-1))
        {
          const double dx = (*this)[n-1]-(*this)[n-2];
          return (*this)[n-1]+(i-(n-1))*dx;
        }

        // real knots
        return (*this)[i];
      }

      size_t
      Locate(const double x)
        const
      {
        return std::upper_bound(begin(),end(),x)-begin()-1;
      }

    };


    class BasisFunction {
    public:

      BasisFunction()
        : fKnotPtr(NULL)
      {}

      BasisFunction(const KnotVector* xknot, int j)
        : fKnotPtr(xknot),
          fIndex(j-3), // shift to internal index
          fIntegral(Internal(GetStop(),-1))
      {}

      BasisFunction(const BasisFunction& other)
      {
        *this = other;
      }

      BasisFunction&
      operator=(const BasisFunction& other)
      {
        fKnotPtr = other.fKnotPtr;
        fIndex = other.fIndex;
        fIntegral = other.fIntegral;      
        return *this;
      }

      double
      operator()(const double x, const char derivative = 0)
        const
      {
        if (derivative == -1)
          return Integral(x);

        if (x < GetStart())
          return 0.0;

        if (x > GetStop())
          return 0.0;

        return Internal(x, derivative)/fIntegral;
      }

      double
      Integral(const double x)
        const
      {
        if (x <= GetStart())
          return 0.0;

        if (x >= GetStop())
          return 1.0;

        return Internal(x, -1)/fIntegral;
      }

      double
      Integral()
        const
      {
        return 1.0;
      }

      double
      GetStart()
        const
      {
        return (*fKnotPtr)[std::max(fIndex,0)];
      }

      double
      GetStop()
        const
      {
        return (*fKnotPtr)[std::min(fIndex+4,int(fKnotPtr->size())-1)];
      }

    private:

      inline
      void
      ComputePolynomialFactors(double& t0, double& t1, double& t2, double& t3, const size_t i)
        const
      {
        const KnotVector& xknot = *fKnotPtr;

        t0 = 0;
        t1 = 0;
        t2 = 0;
        t3 = 0;

        #define SPLINE_BASIS_POLYNOMAL_FACTORS(ia0,ia1,ia2,ib0,ib1,ib2,ic0,ic1,ic2) \
          t0 += -xknot(fIndex+ia0)*xknot(fIndex+ia1)*xknot(fIndex+ia2)/((xknot(fIndex+ib0)-xknot(fIndex+ic0))*(xknot(fIndex+ib1)-xknot(fIndex+ic1))*(xknot(fIndex+ib2)-xknot(fIndex+ic2))); \
          t1 += (xknot(fIndex+ia0)*xknot(fIndex+ia1)+xknot(fIndex+ia0)*xknot(fIndex+ia2)+xknot(fIndex+ia1)*xknot(fIndex+ia2))/((xknot(fIndex+ib0)-xknot(fIndex+ic0))*(xknot(fIndex+ib1)-xknot(fIndex+ic1))*(xknot(fIndex+ib2)-xknot(fIndex+ic2))); \
          t2 += -(xknot(fIndex+ia0)+xknot(fIndex+ia1)+xknot(fIndex+ia2))/((xknot(fIndex+ib0)-xknot(fIndex+ic0))*(xknot(fIndex+ib1)-xknot(fIndex+ic1))*(xknot(fIndex+ib2)-xknot(fIndex+ic2))); \
          t3 += 1.0/((xknot(fIndex+ib0)-xknot(fIndex+ic0))*(xknot(fIndex+ib1)-xknot(fIndex+ic1))*(xknot(fIndex+ib2)-xknot(fIndex+ic2)));

        switch(i)
        {
          case 0:
            SPLINE_BASIS_POLYNOMAL_FACTORS(0,0,0,1,2,3,0,0,0);
          break;
          case 1:
            SPLINE_BASIS_POLYNOMAL_FACTORS(0,0,2,2,3,1,0,0,2);
            SPLINE_BASIS_POLYNOMAL_FACTORS(0,1,3,3,2,1,0,1,3);
            SPLINE_BASIS_POLYNOMAL_FACTORS(1,1,4,2,3,1,1,1,4);
          break;
          case 2:
            SPLINE_BASIS_POLYNOMAL_FACTORS(0,3,3,3,1,2,0,3,3);
            SPLINE_BASIS_POLYNOMAL_FACTORS(1,4,3,3,1,2,1,4,3);
            SPLINE_BASIS_POLYNOMAL_FACTORS(2,4,4,3,1,2,2,4,4);
          break;
          case 3:
            SPLINE_BASIS_POLYNOMAL_FACTORS(4,4,4,1,2,3,4,4,4);
          break;
        }

        #undef SPLINE_BASIS_POLYNOMAL_FACTORS
      }

      double
      Internal(const double x, const char derivative)
        const
      {
        const KnotVector& xknot = *fKnotPtr;

        const size_t iStart = std::max(fIndex,0);
        const size_t iStop  = std::min(fIndex+4,int(fKnotPtr->size())-1);

        size_t i=iStart;
        for (; i < iStop; ++i)
          if (x <= xknot[i+1])
            break;

        double t0,t1,t2,t3;
        ComputePolynomialFactors(t0,t1,t2,t3,i-fIndex);

        switch(derivative)
        {
          case -1:
          {
            const double a = xknot[i];
            double result =
              x*(t0 + x*(t1/2 + x*(t2/3 + x*t3/4)))-
              a*(t0 + a*(t1/2 + a*(t2/3 + a*t3/4)));
            for (size_t j=iStart;j<i;++j)
            {
              ComputePolynomialFactors(t0,t1,t2,t3,j-fIndex);
              const double a = xknot[j];
              const double b = xknot[j+1];
              result +=
                b*(t0 + b*(t1/2 + b*(t2/3 + b*t3/4)))-
                a*(t0 + a*(t1/2 + a*(t2/3 + a*t3/4)));
            }
            return result;
          }
          case 0:
            return t0 + x*(t1 + x*(t2 + x*t3));
          case 1:
            return t1 + x*(2*t2 + x*3*t3);
          case 2:
            return 2*t2 + x*6*t3;
          case 3:
            return 6*t3;
          default:
            return 0;
        }
      }

      const KnotVector* fKnotPtr;
      int fIndex;
      double fIntegral;
    };


    class Function1D {
    public:

      Function1D() {}

      Function1D(const KnotVector& xknot)
      {
        Configure(xknot);
      }

      void
      Configure(const KnotVector& xknot)
      {
        fKnotVector = xknot;
        const size_t size = xknot.size();

        fBasisVector.resize(size+2);
        for (size_t k=0;k<size+2;++k)
          fBasisVector[k] = BasisFunction(&fKnotVector, k);
      }

      inline
      const BasisFunction&
      Basis(const size_t k)
        const
      {
        return fBasisVector[k];
      }

      inline
      double
      operator()(const std::vector<double>& coefs, const double x, const char derivative=0)
        const
      {
        double result = 0.0;
        for (size_t k=0;k<GetSize();++k)
          result += coefs[k]*fBasisVector[k](x,derivative);
        return result;
      }

      inline
      double
      Integral(const std::vector<double>& coefs, const double x = DBL_MAX)
        const
      {
        if (x <= GetStart())
          return 0.0;

        return (*this)(coefs, x, -1);
      }

      inline
      const KnotVector&
      GetKnotVector()
        const
      {
        return fKnotVector;
      }

      size_t
      GetSize()
        const
      {
        return fBasisVector.size();
      }

      double
      GetStart(const size_t k=0)
        const
      {
        return fBasisVector[k].GetStart();
      }

      double
      GetStop(const size_t k=std::numeric_limits<size_t>::max())
        const
      {
        return fBasisVector[k < fKnotVector.size() ? k : fKnotVector.size()-1].GetStop();
      }

    private:
      KnotVector fKnotVector;
      std::vector<BasisFunction> fBasisVector;
    };

  } // NS Spline

} // NS Aru

#endif
