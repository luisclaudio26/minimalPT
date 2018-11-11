#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "geometry.h"
#include "scene.h"
#include "spectrum.h"

class Integrator
{
private:
  RGBA direct_illumination(const Scene& scene);
  RGBA normal_shading(const Scene& scene,
                      const Ray& primary_ray,
                      const Isect& isect);

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
