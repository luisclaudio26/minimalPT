#ifndef SAMPLER_H
#define SAMPLER_H

#include "spectrum.h"
#include "geometry.h"

namespace Sampler
{
  void sample_hemisphere(Vec3& out, float& pdf);
  void sample_sphere(Vec3& out, float& pdf);
};

#endif
