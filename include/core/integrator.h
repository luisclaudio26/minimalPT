#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "geometry.h"
#include "scene.h"
#include "spectrum.h"

class Integrator
{
public:
  Integrator();

  // film data
  int vRes, hRes;

  //TODO: lock access to color_buffer_cpu, because the
  //GUI thread should not read from color_buffer_cpu
  //while the rendering thread is writing samples to it
  ColorBuffer color_buffer;

  void render(const Scene& scene);
};

#endif
