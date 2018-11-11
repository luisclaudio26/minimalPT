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

  // the samples buffer which accumulate sample contributions,
  // and its respective per pixel weights, which we use for filtering.
  // frame holds the elementwise division of sample_buffer by weight.
  RadianceBuffer samples;
  std::vector<float> weights;
  ColorBuffer frame;

  void render(const Scene& scene);
};

#endif
