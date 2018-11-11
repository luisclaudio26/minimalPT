#include "../../include/core/shape.h"

bool Shape::intersect(const Ray& ray, Isect& tgt) const
{
  // prefill tgt with invalid intersection data
  tgt.t = -1.0f;

  //a = 1, because r.d is normalized
  Vec3 oc = ray.o - o;
  float b = 2*glm::dot(ray.d, oc);
  float c = glm::dot(oc,oc) - r*r;
  float delta = b*b-4*c;

  // no intersection!
  if(delta < 0.0f) return false;

  // care for negative intersections, which are
  // behind the camera (and thus are invalid)
  float t = -b-sqrt(delta)*0.5f;          // closest intersection
  if( t < 0.0f ) t = -b+sqrt(delta)*0.5f; // if negative, check further intersection

  // if still negative, camera is behind the ball and we should return nothing
  if( t < 0.0f) return false;

  // if we reached this point, t > 0.0f and thus we have a valid intersection.
  // but we are not checking whether we're inside the sphere or not!
  tgt.t = t;
  tgt.normal = glm::normalize(ray(t)-o);
  tgt.shape = this;

  return true;
}
