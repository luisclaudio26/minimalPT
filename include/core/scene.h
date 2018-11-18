#ifndef SCENE_H
#define SCENE_H

#include <vector>
#include "intersection.h"
#include "camera.h"
#include "shape.h"

class Scene
{
public:
  // at first will explicitly store light sources, and
  // everything will be as brute-force as it can be.
  // Once we importance sample lightsources to optimize
  // things we we'll perform a preprocess step looking
  // for emissive things
  std::vector<Shape> prims;
  std::vector<int> emissive_prims;
  Camera cam;

  void add_primitive(const Shape& s);

  bool cast_ray(const Ray& r, Isect& isect) const;
  bool cast_shadow_ray(const Ray& r, float t) const;

  Scene() {}
};

#endif
