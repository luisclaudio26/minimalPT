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

void Scene::add_primitive(const Shape& s)
{
  prims.push_back(s);

  // check whether s is emissive or not and if yes, push its index
  if(s.emission.r > 0.0f || s.emission.g > 0.0f || s.emission.b > 0.0f)
  {
    emissive_prims.push_back( prims.size()-1 );
    em_area += s.area();
  }
}
