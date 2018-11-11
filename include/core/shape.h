#ifndef SHAPE_H
#define SHAPE_H

#include "intersection.h"
#include "geometry.h"

class Shape
{
public:
  Shape(const Vec3& o, float r)
    : o(o), r(r), diff_color(0.0f,1.0f,0.0f) { }

  // we'll handle spheres only at first!
  Vec3 o; float r;

  // material will be a simple Lambertian color
  Vec3 diff_color;

  bool intersect(const Ray& r, Isect& tgt) const;
};

#endif
