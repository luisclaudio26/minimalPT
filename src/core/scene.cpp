#include "../../include/core/scene.h"
#include <glm/gtx/string_cast.hpp>
#include <cfloat>
#include <cstdio>

bool Scene::cast_ray(const Ray& r, Isect& target) const
{
  return rand() % 100 < 50;
}

bool Scene::cast_shadow_ray(const Ray& r, float t) const
{
  return false;
}
