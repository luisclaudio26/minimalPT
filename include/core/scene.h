#ifndef SCENE_H
#define SCENE_H

#include <vector>
#include "intersection.h"
#include "camera.h"
#include "shape.h"

class Scene
{
public:
  std::vector<Shape> prims;
  std::vector<int> emissive_prims;
  Camera cam;

  //TODO float emissive_area;
  //TODO sample_emissive(); <- importance samples lights according to surface area

  void add_primitive(const Shape& s);

  bool cast_ray(const Ray& r, Isect& isect) const;
  bool cast_shadow_ray(const Ray& r, float t) const;

  Scene() {}
};

#endif
