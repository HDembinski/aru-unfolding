%module aru
%include "std_vector.i"

namespace std {
  %template(vector_double) vector<double>;
  %template(matrix_double) vector< vector<double> >;
}

%{
#include "Spline.h"
#include "VKernel.h"
#include "FoldedSpline.h"
#include "ObjectiveFunction.h"
%}

%rename(SplineBasisFunction) Aru::Spline::BasisFunction;
%rename(SplineFunction1D) Aru::Spline::Function1D;

%rename(FoldedSplineBasisFunction) Aru::FoldedSpline::BasisFunction;
%rename(FoldedSplineFunction1D) Aru::FoldedSpline::Function1D;

%ignore nlopt_wrapper;
%ignore Aru::Spline::BasisFunction::operator=;

%include "Spline.h"
%include "VKernel.h"
%include "FoldedSpline.h"
%include "ObjectiveFunction.h"
