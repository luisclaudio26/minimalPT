#include "../../include/core/sampler.h"

Vec2 Sampler::get2D()
{
  return Vec2(get1D(), get1D());
}
