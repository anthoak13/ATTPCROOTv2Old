#ifndef ATSAMPLEMETHODS_H
#define ATSAMPLEMETHODS_H

#include "AtChargeWeighted.h"
#include "AtGaussian.h"
#include "AtSample.h"
#include "AtUniform.h"
#include "AtWeightedGaussian.h"

enum class AtTools::AtSample::SampleMethod { kUniform = 0, kChargeWeighted = 1, kGaussian = 2, kWeightedGaussian = 3 };

#endif // #ifndef ATSAMPLEMETHODS_H
