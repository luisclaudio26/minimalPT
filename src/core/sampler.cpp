#include "../../include/core/sampler.h"

// internal random number engine
static std::mt19937 SAMPLER_MT;

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

void Sampler::uniform_sample_disk(Vec2& out, float& pdf)
{
  std::uniform_real_distribution<float> random_float(0.0f, 1.0f);
  float u1 = random_float(SAMPLER_MT);
  float u2 = random_float(SAMPLER_MT);

  float r = sqrt(u1);
  float theta = _2pi * u2;

  out = Vec2(r*cos(theta), r*sin(theta));
  pdf = _overpi;
}

void Sampler::cosine_weight_sample_hemisphere(Vec3& out, float& pdf)
{
  Vec2 d; float pdf_disk;
  uniform_sample_disk(d, pdf_disk);
  float z = sqrt(std::max(0.0f, 1.0f-d.x*d.x-d.y*d.y));

  out = Vec3(d.x, z, d.y);
  pdf = z * _overpi; // This is the dot product of OUT with the normal [0 1 0]
}
