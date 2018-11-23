#ifndef SAMPLER_H
#define SAMPLER_H

#include "spectrum.h"
#include "geometry.h"

const float _pi = 3.141592654f;
const float _2pi = 2.0f*_pi;
const float _over2pi = 1.0f / _2pi;
const float _over4pi = 1.0f / (4.0f * _pi);

namespace Sampler
{
  void sample_hemisphere(Vec3& out, float& pdf);
  void sample_sphere(Vec3& out, float& pdf);
};

#endif
