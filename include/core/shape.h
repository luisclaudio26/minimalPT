#ifndef SHAPE_H
#define SHAPE_H

#include "intersection.h"
#include "geometry.h"

class Shape
{
public:
  Shape(const Vec3& o, float r)
    : o(o), r(r), diff_color(0.0f,1.0f,0.0f), emission(0.0f) { }

  // we'll handle spheres only at first!
  Vec3 o; float r;

  // material will be a simple Lambertian color
  Vec3 diff_color;

  // Emission profile should be a radiance
  // function Le(x,d) that tells us how much power
  // is emitted from this shape surface at point x and
  // direction d. By simply putting a constant, we
  // assume that radiance is uniform and isotropic
  // over the shape's surface (i.e. every point emits
  // the same power at every direction).
  //
  // This means for consistency reasons, when assigning
  // a value we should have in mind the irradiance of
  // each point (i.e. the sum of the radiance in all
  // directions for a given point) and the total power irradiated
  // of the surface (the sum of the irradiances of every point)
  // in order for quantities to make sense. This means that
  // simply assigning "1.0" (W.m⁻².sr⁻¹), for example, will
  // result in a sphere that yields WAY more power to
  // the scene than one might originally have intended to.
  //
  // Finally, keep in mind that integration of radiance for
  // a given point must be cosine weighted to account for
  // Lambert's cosine law.
  Vec3 emission;

  bool intersect(const Ray& r, Isect& tgt) const;
};

#endif
