#include "../../include/core/sampler.h"

void Sampler::sample_hemisphere(Vec3& out, float& pdf)
{
  float u1 = (float)rand()/RAND_MAX;
  float u2 = (float)rand()/RAND_MAX;

  const float r = sqrt(1.0f - u1*u1);
  const float phi = _2pi * u2;

  out = Vec3(r*cos(phi), u1, r*sin(phi));
  pdf = _over2pi;
}

void Sampler::sample_sphere(Vec3& out, float& pdf)
{
  // uniformly sample point on a unitary sphere centered at the origin,
  // then scale and translate it so it lies on the surface of this sphere
  float u1 = (float)rand()/RAND_MAX;
  float u2 = (float)rand()/RAND_MAX;

  const float y = 1 - 2*u1;
  const float r = sqrt(std::max(0.0f, 1.0f - y * y));
  const float phi = _2pi * u2;

  out = Vec3(r*cos(phi), y, r*sin(phi));
  pdf = _over4pi;
}
