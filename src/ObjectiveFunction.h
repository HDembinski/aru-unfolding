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

#ifndef _Aru_ObjectiveFunction_h_
#define _Aru_ObjectiveFunction_h_

#include "VKernel.h"
#include "Utils.h"
#include "Spline.h"
#include "FoldedSpline.h"

#include <cmath>
#include <vector>
#include <string>
#include <exception>
#include <iostream>
#include <limits>
#include <algorithm>
#include <stdexcept>

#include <TVectorD.h>
#include <TMatrixDSym.h>
#include <TDecompSVD.h>

#include <nlopt.hpp>

namespace Aru {

  double
  NloptObjectiveFunctionInterface(
    const std::vector<double>& pars,
    std::vector<double>& grad,
    void* meta);

  typedef std::vector< std::vector<double> > Matrix;

  void
  Copy(Matrix& output, const TMatrixDSym& input)
  {
    const size_t n = input.GetNcols();
    output.resize(n,std::vector<double>(n));
    for (size_t i=0;i<n;++i)
      for (size_t j=0;j<n;++j)
        output[i][j] = input[i][j];
  }

  void
  Copy(TMatrixDSym& output, Matrix& input)
  {
    const size_t n = input.size();
    output.ResizeTo(n,n);
    for (size_t i=0;i<n;++i)
      for (size_t j=0;j<n;++j)
        output[i][j] = input[i][j];
  }

  void
  MatrixAdd(Matrix& output, const Matrix& input)
  {
    const size_t n = input.size();
    const size_t m = input[0].size();
    for (size_t i=0;i<n;++i)
      for (size_t j=0;j<m;++j)
        output[i][j] += input[i][j];
  }
}

namespace Aru {

  class ObjectiveFunction {
  private:
    const VKernel& fKernel;
    const FoldedSpline::Function1D& fFolded;
    const Spline::Function1D& fUnfolded;
    std::vector<double> fFoldedIntegral;
    std::vector<double> fYEdges;
    std::vector<size_t> fWs;
    Matrix fFoldedMuCache;
    Matrix fFoldedCache;
    std::vector<double> fReferenceCoefs;
    Matrix fReferenceCovariance;
    size_t fDataSize;
    double fWeight;
    Matrix fFMatrix;

  public:
    ObjectiveFunction(const std::vector<double>& data,
                      const VKernel& kernel,
                      const FoldedSpline::Function1D& folded,
                      const Spline::Function1D& unfolded,
                      const size_t histogramSubSampling = 8,
                      const size_t nEventsForHistogram = 500) :
      fKernel(kernel),
      fFolded(folded),
      fUnfolded(unfolded),
      fFoldedIntegral(folded.GetSize()),
      fYEdges(histogramSubSampling*(folded.GetKnotVector().size()-1)+1),
      fWs(fYEdges.size()-1),
      fFoldedMuCache(folded.GetSize(),std::vector<double>(fYEdges.size()-1)),
      fReferenceCoefs(folded.GetSize(),0),
      fReferenceCovariance(folded.GetSize(),std::vector<double>(folded.GetSize())),
      fDataSize(0),
      fWeight(0),
      fFMatrix(folded.GetSize(),std::vector<double>(folded.GetSize()))
    {
      const size_t nBasis = folded.GetSize();
      const size_t nData = data.size();

      const Spline::KnotVector& yknots = folded.GetKnotVector();
      const size_t nBinned = GetBinnedSize();
      for (size_t ib=0;ib<nBinned+1;++ib)
      {
        const double delta = double(ib % histogramSubSampling)/histogramSubSampling;
        fYEdges[ib] = delta==0.0 ?
          yknots[ib/histogramSubSampling] :
          (1.0-delta)*yknots[ib/histogramSubSampling] + delta*yknots[ib/histogramSubSampling+1];
      }

      std::vector<size_t> index(nData);
      for (size_t i=0;i<nData;++i)
        if (fYEdges.front() <= data[i] && data[i] <= fYEdges.back())
        {
          index[i] = upper_bound(fYEdges.begin(),fYEdges.end(),data[i])-fYEdges.begin()-1;
          ++fWs[index[i]];
          ++fDataSize;
        }

      std::vector<double> unbinned;
      for (size_t i=0;i<nData;++i)
        if (fYEdges.front() <= data[i] && data[i] <= fYEdges.back() && fWs[index[i]] < nEventsForHistogram)
          unbinned.push_back(data[i]);

      for (size_t k=0;k<nBasis;++k)
        fFoldedIntegral[k] = folded.Basis(k).Integral();

      for (size_t ib=0;ib<nBinned;++ib)
        if (fWs[ib] < nEventsForHistogram)
        {
          for (size_t k=0;k<nBasis;++k)
          {
            fWs[ib] = 0;
            fFoldedMuCache[k][ib] = 0.0;
          }
        }
        else
          for (size_t k=0;k<nBasis;++k)
            fFoldedMuCache[k][ib] = folded.Basis(k).Integral(fYEdges[ib],fYEdges[ib+1]);

      const size_t nUnbinned = unbinned.size();
      fFoldedCache.resize(nBasis, std::vector<double>(nUnbinned));
      for (size_t i=0;i<nUnbinned;++i)
        for (size_t k=0;k<nBasis;++k)
          fFoldedCache[k][i] = folded.Basis(k)(unbinned[i]);

      std::cout << "nData " << data.size() << " "
                << "nUnbinned " << nUnbinned << std::endl;

      FMatrixIntegrand t(folded);
      for (size_t i=0;i<nBasis;++i)
        for (size_t j=i;j<nBasis;++j)
      {
        t.Set(i,j);
        const double a = t.GetStart();
        const double b = t.GetStop();
        if (a < b)
          fFMatrix[i][j] = Integrate(t, a, b);
        if (i!=j) fFMatrix[j][i] = fFMatrix[i][j];
      }
    }

    size_t
    GetDataSize()
      const
    {
      return fDataSize;
    }

    size_t
    GetBinnedSize()
      const
    {
      return fYEdges.size()-1;
    }

    size_t
    GetUnbinnedSize()
      const
    {
      return fFoldedCache[0].size();
    }

    size_t
    GetBasisSize()
      const
    {
      return fFoldedCache.size();
    }

    const std::vector<double>&
    GetReferenceCoefficients()
      const
    {
      return fReferenceCoefs;
    }

    void
    SetReferenceCoefficients(const std::vector<double>& coefs)
    {
      if (coefs.size()!=fReferenceCoefs.size())
        throw std::runtime_error("sizes of coefficient vectors don't match");
      fReferenceCoefs = coefs;
    }

    const Matrix&
    GetReferenceCovariance()
      const
    {
      return fReferenceCovariance;
    }

    void
    SetReferenceCovariance(const Matrix& cov)
    {
      fReferenceCovariance = cov;
    }

    void
    SetWeight(const double weight)
    {
      fWeight = weight;
    }

    double
    GetWeight()
      const
    {
      return fWeight;
    }

    double
    GetFoldedMuCache(const size_t k, const size_t i)
      const
    {
      return fFoldedMuCache[k][i];
    }

    double
    GetFoldedCache(const size_t k, const size_t i)
      const
    {
      return fFoldedCache[k][i];
    }

    double
    GetHistogramWeight(const size_t i)
      const
    {
      return fWs[i];
    }

    double
    GetHistogramStep(const size_t i)
      const
    {
      return fYEdges[i+1]-fYEdges[i];
    }

    const Matrix&
    GetFMatrix()
      const
    {
      return fFMatrix;
    }

    void
    ValueStat(double& value, const std::vector<double>& coefs)
      const
    {
      const size_t nBasis = GetBasisSize();

      if (coefs.size() != nBasis)
        throw std::runtime_error("size of coefficient vector and function vector differ");

      value = 0.0;
      for (size_t k=0;k<nBasis;++k)
        value += coefs[k]*fFoldedIntegral[k];

      for (size_t i=0;i<GetUnbinnedSize();++i)
      {
        double folded = 0;
        for (size_t k=0;k<nBasis;++k)
          folded += coefs[k]*fFoldedCache[k][i];
        value -= std::log(folded);
      }

      for (size_t ib=0;ib<GetBinnedSize();++ib)
      {
        if (fWs[ib]==0.0) continue;
        double foldedMu = 0;
        for (size_t k=0;k<nBasis;++k)
          foldedMu += coefs[k]*fFoldedMuCache[k][ib];
        value -= fWs[ib]*std::log(foldedMu);
      }
    }

    void
    GradStat(std::vector<double>& gradient, const std::vector<double>& coefs)
      const
    {
      const size_t nBasis = GetBasisSize();

      if (coefs.size() != nBasis)
        throw std::runtime_error("size of coefficient vector and function vector differ");

      if (gradient.size() != nBasis)
        throw std::runtime_error("size of gradient and function vector differ");

      for (size_t k=0; k<nBasis; ++k)
      {
        gradient[k] = fFoldedIntegral[k];

        for (size_t i=0;i<GetUnbinnedSize();++i)
        {
          double folded = 0.0;
          for (size_t l=0;l<nBasis;++l)
            folded += coefs[l]*fFoldedCache[l][i];
          gradient[k] -= fFoldedCache[k][i]/folded;
        }

        for (size_t ib=0;ib<GetBinnedSize();++ib)
        {
          if (fWs[ib] == 0.0) continue;
          double foldedMu = 0.0;
          for (size_t l=0;l<nBasis;++l)
            foldedMu += coefs[l]*fFoldedMuCache[l][ib];
          gradient[k] -= fWs[ib]*fFoldedMuCache[k][ib]/foldedMu;
        }
      }
    }

    void
    HesseStat(Matrix& hesse, const std::vector<double>& coefs)
      const
    {
      const size_t nBasis = GetBasisSize();

      if (coefs.size() != nBasis)
        throw std::runtime_error("size of coefficient vector and function vector differ");

      if (hesse.size() != nBasis || hesse[0].size() != nBasis)
        throw std::runtime_error("sizes of hesse matrix and function vector differ");

      for (size_t k=0; k<nBasis; ++k)
        for (size_t l=k; l<nBasis; ++l)
      {
        hesse[k][l] = 0;

        for (size_t i=0;i<GetUnbinnedSize();++i)
        {
          double folded = 0.0;
          for (size_t m=0;m<nBasis;++m)
            folded += coefs[m]*fFoldedCache[m][i];
          hesse[k][l] += fFoldedCache[k][i]*fFoldedCache[l][i]/Sqr(folded);
        }

        for (size_t ib=0;ib<GetBinnedSize();++ib)
        {
          if (fWs[ib] == 0.0) continue;
          double foldedMu = 0.0;
          for (size_t m=0;m<nBasis;++m)
            foldedMu += coefs[m]*fFoldedMuCache[m][ib];
          hesse[k][l] += fWs[ib]*fFoldedMuCache[k][ib]*fFoldedMuCache[l][ib]/Sqr(foldedMu);
        }

        if (k!=l) hesse[l][k] = hesse[k][l];
      }
    }

    void
    ValueReg(double& value, const std::vector<double>& coefs)
      const
    {
      const size_t n = GetBasisSize();

      if (coefs.size() != n)
        throw std::runtime_error("size of coefficient vector and function vector differ");

      Integrand t(this, coefs);

      value = Integrate(t, t.GetStart(), t.GetStop());

      for (size_t k=0;k<n;++k)
        value -= coefs[k]*fUnfolded.Basis(k).Integral();
    }

    void
    GradReg(std::vector<double>& gradient, const std::vector<double>& coefs)
      const
    {
      const size_t n = GetBasisSize();

      if (coefs.size() != n)
        throw std::runtime_error("size of coefficient vector and function vector differ");

      if (gradient.size() != n)
        throw std::runtime_error("size of gradient and function vector differ");

      for (size_t l=0;l<n;++l)
      {
        IntegrandGradient t(this, coefs, l);
        gradient[l] = Integrate(t, t.GetStart(), t.GetStop());
      }
    }

    void
    HesseReg(Matrix& hesse, const std::vector<double>& coefs)
      const
    {
      const size_t nBasis = GetBasisSize();

      if (coefs.size() != nBasis)
        throw std::runtime_error("size of coefficient vector and function vector differ");

      if (hesse.size() != nBasis || hesse[0].size() != nBasis)
        throw std::runtime_error("sizes of hesse matrix and function vector differ");

      if (GetWeight()==0.0)
      {
        for (size_t i=0;i<nBasis;++i)
          for (size_t j=0;j<nBasis;++j)
            hesse[i][j] = i==j ? 1.0 : 0.0;
        return;
      }

      for (size_t k=0;k<nBasis;++k)
        for (size_t l=k;l<nBasis;++l)
      {
        const Spline::Function1D& basis = fUnfolded;

        if (basis.GetStop(k) <= basis.GetStart(l) ||
            basis.GetStop(l) <= basis.GetStart(k) )
          hesse[k][l] = 0.0;
        else
        {
          IntegrandHesse t(this, coefs, k, l);
          hesse[k][l] = Integrate(t, t.GetStart(), t.GetStop());

          if (k!=l) hesse[l][k] = hesse[k][l];
        }
      }
    }

    void
    ComputeEigenSystem(TVectorD& s,
                       TMatrixD& m,
                       const std::vector<double>& coefs)
      const
    {
      const size_t n = GetBasisSize();

      TMatrixDSym hesse1(n), hesse2(n);
      Matrix mhesse1(n,std::vector<double>(n)), mhesse2(n,std::vector<double>(n));
      HesseStat(mhesse1, coefs);
      HesseReg(mhesse2, coefs);
      Copy(hesse1,mhesse1);
      Copy(hesse2,mhesse2);

      TVectorD d2;
      const TMatrixD u2 = hesse2.EigenVectors(d2);

      // S-Tmp = D_2^{-1/2} U_2^T H_1 U_2 D_2^{-1/2}
      TMatrixDSym sTmp(n);
      for (size_t i=0;i<n;++i)
        for (size_t o=0;o<n;++o)
      {
        sTmp(i,o) = 0.0;
        for (size_t k=0;k<n;++k)
        {
          double tmpk = 0.0;
          for (size_t m=0;m<n;++m)
            tmpk += hesse1(k,m)*u2(m,o);
          sTmp(i,o) += u2(k,i)*tmpk;
        }
        sTmp(i,o) /= std::sqrt(d2[i]*d2[o]);
      }

      const TMatrixD u1 = sTmp.EigenVectors(s);

      // m = u2 d2^{-1/2} u1
      for (size_t i=0;i<n;++i)
        for (size_t j=0;j<n;++j)
      {
        m(i,j) = 0.0;
        for (size_t k=0;k<n;++k)
          m(i,j) += u2(i,k)*u1(k,j)/std::sqrt(d2[k]);
      }
    }

    void
    CovarianceStat(Matrix& covariance,
                   const std::vector<double>& coefs)
      const
    {
      const size_t n = GetBasisSize();

      //
      // propagation of statistical uncertainty
      //

      TVectorD s(n);
      TMatrixD m(n,n);
      ComputeEigenSystem(s, m, coefs);
      const double weight = GetWeight();

      // var[i,j] = sum_kl d[c]_i/d[s^{1/2} \bar{c1}]_k d[c]_j/d[s^{1/2} \bar{c1}]_l delta(k,l)
      //   because [s^{1/2} \bar{c1}]_i have unit variance
      // var[i,j] = sum_k M_ik M_jk S_kk/(S_kk+w)^2
      for (size_t i=0;i<n;++i)
        for (size_t j=i;j<n;++j)
      {
        covariance[i][j] = 0.0;
        for (size_t k=0;k<n;++k)
          covariance[i][j] += m(i,k)*m(j,k)*s[k]/Sqr(s[k]+weight);
        if (i!=j) covariance[j][i] = covariance[i][j];
      }
    }

    void
    CovarianceReg(Matrix& covariance,
                  const std::vector<double>& coefs)
      const
    {
      const size_t n = GetBasisSize();

      //
      // propagation of statistical uncertainty of reference distribution
      //

      Matrix hesse1(n,std::vector<double>(n));
      Matrix hesse2(n,std::vector<double>(n));
      HesseStat(hesse1, coefs);
      HesseReg(hesse2, coefs);
      const double weight = GetWeight();

      TMatrixDSym hesseInv(n);
      for (size_t i=0;i<n;++i)
        for (size_t j=0;j<n;++j)
          hesseInv(i,j) = hesse1[i][j] + weight*hesse2[i][j];

      TVectorD d;
      const TMatrixD u = hesseInv.EigenVectors(d);
      for (size_t i=0;i<n;++i)
        for (size_t j=0;j<n;++j)
      {
        hesseInv(i,j) = 0.0;
        for (size_t k=0;k<n;++k)
          if (d[k]>1e-10)
            hesseInv(i,j) += u(i,k)*u(j,k)/d[k];
      }

      TMatrixD a(n,n);
      for (size_t i=0;i<n;++i)
        for (size_t j=0;j<n;++j)
      {
        const Spline::Function1D& basis = fUnfolded;

        if (basis.GetStop(i) <= basis.GetStart(j) ||
            basis.GetStop(j) <= basis.GetStart(i) )
          a(i,j) = 0; // no overlap
        else
        {
          IntegrandGradientPropagation t(this,coefs,i,j);
          a(i,j) = Integrate(t, t.GetStart(), t.GetStop());
          a(i,j) *= weight;
        }
      }

      a = hesseInv*a;

      for (size_t i=0;i<n;++i)
        for (size_t j=i;j<n;++j)
      {
        covariance[i][j] = 0.0;
        for (size_t k=0;k<n;++k)
          for (size_t l=0;l<n;++l)
            covariance[i][j] += a(i,k)*a(j,l)*fReferenceCovariance[k][l];
        if (i!=j) covariance[j][i] = covariance[i][j];
      }
    }

    void
    UpdateReference(const std::vector<double>& coefs)
    {
      const Spline::KnotVector& xknots = fUnfolded.GetKnotVector();
      const size_t nBasis = fUnfolded.GetSize();
      const size_t sub = 4;
      const size_t nData = sub*xknots.size()-1;

      Matrix cov(nBasis,std::vector<double>(nBasis));
      CovarianceStat(cov, coefs);

      TMatrixD f(nData,nBasis);
      TMatrixD a(nData,nBasis);
      for (size_t i=0;i<nData;++i)
        for (size_t k=0;k<nBasis;++k)
      {
        const double delta = double(i%sub)/sub;
        const double xi = (i%sub==0 || i==nData-1) ?
          xknots[i/sub] : (1.0-delta)*xknots[i/sub]+delta*xknots[i/sub+1];
        const double yi = fKernel.Y(xi);
        f(i,k) = fFolded.Basis(k)(yi)/fKernel.Efficiency(yi);
        a(i,k) = fUnfolded.Basis(k)(xi);
      }

      TVectorD b(nData);
      TMatrixDSym bvar(nData);

      for (size_t i=0;i<nData;++i)
        for (size_t k=0;k<nBasis;++k)
          b[i] += coefs[k]*f(i,k);

      for (size_t i=0;i<nData;++i)
        for (size_t j=i;j<nData;++j)
      {
        for (size_t k=0;k<nBasis;++k)
          for (size_t l=k;l<nBasis;++l)
            bvar(i,j) += (k==l?1:2)*f(i,k)*f(j,l)*cov[k][l];
        if (i!=j) bvar(j,i) = bvar(i,j);
      }

      TDecompSVD svd(a);
      TMatrixD pinv(nBasis,nData);
      for (size_t i=0;i<nBasis;++i)
        for (size_t j=0;j<nData;++j)
          for (size_t k=0;k<nBasis;++k)
            if (svd.GetSig()[k]>1e-12)
              pinv(i,j) += svd.GetV()(i,k)/svd.GetSig()[k]*svd.GetU()(j,k);

      for (size_t i=0;i<nBasis;++i)
      {
        fReferenceCoefs[i] = 0.0;
        for (size_t j=0;j<nData;++j)
          fReferenceCoefs[i] += pinv(i,j)*b[j];

        if (fReferenceCoefs[i] < 1e-3)
          fReferenceCoefs[i] = 1e-3;
      }

      for (size_t i=0;i<nBasis;++i)
        for (size_t j=i;j<nBasis;++j)
      {
        fReferenceCovariance[i][j] = 0.0;
        for (size_t k=0;k<nData;++k)
          for (size_t l=k;l<nData;++l)
            fReferenceCovariance[i][j] += (k==l?1:2)*pinv(i,k)*pinv(j,l)*bvar(k,l);
        if (i!=j) fReferenceCovariance[j][i] = fReferenceCovariance[i][j];
      }
    }

    void
    Minimize(std::vector<double>& coefs)
      const
    {
      const size_t nBasis = GetBasisSize();
      if (coefs.size() != nBasis)
        throw std::runtime_error("coefs has wrong size");

      const double lower_bound = 1e-2;
      const double accuracy = 1e-5;
      
      std::vector<double> starts(coefs);
      for (size_t k=0;k<nBasis;++k)
        if (starts[k] < lower_bound) starts[k] = 2*lower_bound;
      
      nlopt::opt opt(nlopt::LD_LBFGS,nBasis);
      opt.set_min_objective(NloptObjectiveFunctionInterface,
        const_cast<void*>(static_cast<const void*>(this)));
      opt.set_lower_bounds(lower_bound);
      opt.set_ftol_abs(accuracy);
      opt.set_maxeval(5000);

      try {
        coefs = opt.optimize(starts);
      }
      catch(std::runtime_error& e)
      {
        std::cerr << "warning: LBFGS failed, switching to MMA" << std::endl;

        nlopt::opt opt(nlopt::LD_MMA,nBasis);
        opt.set_min_objective(NloptObjectiveFunctionInterface,
          const_cast<void*>(static_cast<const void*>(this)));
        opt.set_lower_bounds(lower_bound);
        opt.set_ftol_abs(accuracy);
        opt.set_maxeval(5000);

        coefs = opt.optimize(starts);
      }
    }

  private:

    class FMatrixIntegrand {
    private:
      const FoldedSpline::Function1D& fFolded;
      size_t fI;
      size_t fJ;

    public:
      FMatrixIntegrand(const FoldedSpline::Function1D& folded)
        : fFolded(folded), fI(0), fJ(0)
      {
      }

      void
      Set(const size_t i, const size_t j)
      {
        fI = i; fJ = j;
      }

      inline
      double
      operator()(const double y)
        const
      {
        return fFolded.Basis(fI)(y)*fFolded.Basis(fJ)(y);
      }

      double
      GetStart()
        const
      {
        return std::max(fFolded.GetStart(fI),fFolded.GetStart(fJ));
      }

      double
      GetStop()
        const
      {
        return std::min(fFolded.GetStop(fI),fFolded.GetStop(fJ));
      }
    };

    class Integrand {
    private:
      const Spline::Function1D& fBasis;
      const std::vector<double>& fCoef;
      const std::vector<double>& fReferenceCoefs;

    public:
      Integrand(const ObjectiveFunction* parent,
                const std::vector<double>& coefs) :
        fBasis(parent->fUnfolded),
        fCoef(coefs),
        fReferenceCoefs(parent->fReferenceCoefs)
      {
      }

      inline
      double
      operator()(double x)
        const
      {
        const double bx = fBasis(fCoef, x);
        const double hx = fBasis(fReferenceCoefs, x);
        return bx*std::log(bx/hx);
      }

      double
      GetStart() const { return fBasis.GetStart(); }

      double
      GetStop() const { return fBasis.GetStop(); }
    };

    class IntegrandGradient {
    private:
      const Spline::Function1D& fBasis;
      const std::vector<double>& fCoef;
      const std::vector<double>& fReferenceCoefs;
      const Spline::BasisFunction& fBasisSpline;

    public:
      IntegrandGradient(const ObjectiveFunction* parent,
                        const std::vector<double>& coefs,
                        size_t i) :
        fBasis(parent->fUnfolded),
        fCoef(coefs),
        fReferenceCoefs(parent->fReferenceCoefs),
        fBasisSpline(fBasis.Basis(i))
      {
      }

      inline
      double
      operator()(double x) const
      {
        const double bi = fBasisSpline(x);
        const double bx = fBasis(fCoef, x);
        const double hx = fBasis(fReferenceCoefs, x);
        return bi*std::log(bx/hx);
      }

      double
      GetStart() const { return fBasisSpline.GetStart(); }

      double
      GetStop() const { return fBasisSpline.GetStop(); }
    };

    class IntegrandHesse {
    private:
      const Spline::Function1D& fBasis;
      const std::vector<double>& fCoef;
      const Spline::BasisFunction& fBasisSplineI;
      const Spline::BasisFunction& fBasisSplineJ;

    public:
      IntegrandHesse(const ObjectiveFunction* parent,
                     const std::vector<double>& coefs,
                     size_t i, size_t j) :
        fBasis(parent->fUnfolded),
        fCoef(coefs),
        fBasisSplineI(fBasis.Basis(i)),
        fBasisSplineJ(fBasis.Basis(j))
      {
      }

      inline
      double
      operator()(double x) const
      {
        const double bx = fBasis(fCoef, x);
        const double bi = fBasisSplineI(x);
        const double bj = fBasisSplineJ(x);
        return bi*bj/bx;
      }

      inline
      double
      GetStart()
        const
      {
        return std::max(fBasisSplineI.GetStart(),fBasisSplineJ.GetStart());
      }

      inline
      double
      GetStop()
        const
      {
        return std::min(fBasisSplineI.GetStop(),fBasisSplineJ.GetStop());
      }
    };

    class IntegrandGradientPropagation {
    private:
      const Spline::Function1D& fBasis;
      const std::vector<double>& fCoef;
      const std::vector<double>& fReferenceCoefs;
      const Spline::BasisFunction& fBasisSpline;
      const Spline::BasisFunction& fReferenceSpline;

    public:
      IntegrandGradientPropagation(const ObjectiveFunction* parent,
                                   const std::vector<double>& coefs,
                                   size_t i, size_t j) :
        fBasis(parent->fUnfolded),
        fCoef(coefs),
        fReferenceCoefs(parent->fReferenceCoefs),
        fBasisSpline(fBasis.Basis(i)),
        fReferenceSpline(fBasis.Basis(j))
      {
      }

      inline
      double
      operator()(double x) const
      {
        const double b1 = fBasisSpline(x);
        const double p2 = fReferenceSpline(x);
        const double hx = fBasis(fReferenceCoefs, x);
        return b1*p2/hx;
      }

      inline
      double
      GetStart()
        const
      {
        return std::max(fBasisSpline.GetStart(),fReferenceSpline.GetStart());
      }

      inline
      double
      GetStop()
        const
      {
        return std::min(fBasisSpline.GetStop(),fReferenceSpline.GetStop());
      }
    };

  }; // Objective

  double
  NloptObjectiveFunctionInterface(
    const std::vector<double>& pars,
    std::vector<double>& grad,
    void* meta)
  {
    const ObjectiveFunction& obj = *static_cast<const ObjectiveFunction*>(meta);

    double value=0;
    obj.ValueStat(value, pars);
    obj.GradStat(grad, pars);

    const double weight = obj.GetWeight();
    if (weight!=0.0)
    {
      double valueR=0;
      obj.ValueReg(valueR,pars);
      value += weight*valueR;

      std::vector<double> gradR(grad.size());
      obj.GradReg(gradR,pars);
      for (size_t k=0;k<grad.size();++k)
        grad[k] += weight*gradR[k];
    }
    return value;
  }

} // NS Aru

#endif
