%module kernel

%{
#include "GaussianKernel.h"
#include "BlobelKernel.h"
#include "MyKernel.h"
%}

%include "VKernel.h"
%include "GaussianKernel.h"
%include "BlobelKernel.h"
%include "MyKernel.h"
