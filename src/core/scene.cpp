#include "../../include/core/scene.h"
#include <glm/gtx/string_cast.hpp>
#include <cfloat>
#include <cstdio>

bool Scene::cast_ray(const Ray& r, Isect& target) const
{
  float closest = FLT_MAX;
  target.t = -1.0f;

  // loop over primitives, checking for intersections.
  // if we succeeded to intersect the shape with r,
  // update closest intersection
  for(const Shape& s : prims)
  {
    Isect isect;
    if( s.intersect(r,isect) && isect.t < closest )
    {
      target = isect;
      closest = isect.t;
    }
  }

  return target.is_valid();
}

bool Scene::cast_shadow_ray(const Ray& r, float t) const
{
  return false;
}
