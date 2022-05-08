#ifndef ATSAMPLEMETHODS_H
#define ATSAMPLEMETHODS_H

#include "AtChargeWeighted.h"
#include "AtGaussian.h"
#include "AtSample.h"
#include "AtUniform.h"
#include "AtWeightedGaussian.h"

/**
 * @brief Collection for sampling AtHits
 *
 * Group of classes for randomly sampling AtHits according to some probability density function (PDF).
 *
 * Provides subclassed interfaces for the case where the hits sampled are independent (AtIndependentSample) and when the
 * the PDF depends on some reference hit (AtSampleFromReference). To add an additional sampling method, at minimum it
 * must inherit AtSample. It should also be added to the SampleMethod enum, and the static factory method
 * AtSample::CreateSampler.
 *
 * @defgroup AtHitSampling Random Sampling
 */

/**
 *
 *
 */
enum class AtTools::AtSample::SampleMethod { kUniform = 0, kChargeWeighted = 1, kGaussian = 2, kWeightedGaussian = 3 };

#endif // #ifndef ATSAMPLEMETHODS_H
