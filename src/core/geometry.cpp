#include "../../include/core/geometry.h"
#include <cfloat>
#include <cstdio>

bool AABB::intersect(const Ray& r, float& tmin, float& tmax) const
{
  #define EPS 0.00001f

  //use tavianator's slab method: sucessively clip ray in each axis
  tmin = -FLT_MAX; tmax = FLT_MAX;

  for(int i = 0; i < 3; ++i)
  {
    //if r.d[i] == 0, ray is parallel to the current axis
    if(r.d[i] == 0.0f ) continue;

    float t1, t2, ro = r.o[i];

    //this should avoid the cases where the ray intersects infinitely
    //many points on one of the planes
    if( min[i] == r.o[i] || max[i] == r.o[i]) ro += EPS;

    t1 = (min[i] - ro) / r.d[i];
    t2 = (max[i] - ro) / r.d[i];

    tmin = std::fmax(tmin, std::fmin(t1, t2));
    tmax = std::fmin(tmax, std::fmax(t1, t2));
  }

  //tmax = tmin is a hit right in the corner of the box,
  //which we assume to not to be a hit! TODO: is this a problem?
  return tmax >= tmin && tmax > 0.0f;
}
