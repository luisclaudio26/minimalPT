#ifndef SAMPLER_H
#define SAMPLER_H

#include "geometry.h"

class Sampler
{
public:
  virtual float get1D() = 0;
  Vec2 get2D();
};

#endif
