#include "../../include/core/shape.h"

bool Shape::intersect(const Ray& ray, Isect& tgt) const
{
  // prefill tgt with invalid intersection data
  tgt.t = -1.0f;

  //a = 1, because r.d is normalized
  Vec3 oc = ray.o - o;
  float a = glm::dot(ray.d, ray.d);
  float b = 2*glm::dot(ray.d, oc);
  float c = glm::dot(oc,oc) - r*r;
  float delta = b*b-4*c;

  // no intersection!
  if(delta < 0.0f) return false;

  // care for negative intersections, which are
  // behind the camera (and thus are invalid)
  float t = (-b-sqrt(delta))/(2*a);          // closest intersection
  if( t < 0.0f ) t = (-b+sqrt(delta))/(2*a); // if negative, check further intersection

  // if still negative, camera is behind the ball and we should return nothing
  if( t < 0.0f) return false;

  // if we reached this point, t > 0.0f and thus we have a valid intersection.
  // but we are not checking whether we're inside the sphere or not!
  tgt.t = t;
  tgt.d2 = glm::dot(oc, oc);
  tgt.normal = glm::normalize(ray(t)-this->o);
  tgt.shape = this;

  // compute an arbitrary tangent to the intersection point which will be
  // necessary for a local coordinate system.
  // At least for our spheres, it is proven that it is not possible to compute
  // Tangents without having singularities
  //
  //     https://en.wikipedia.org/wiki/Hairy_ball_theorem
  //
  // TODO: In this first version, our tangent will be simply the projection
  // of the ray onto the normal. thus, this will break if a ray hits the surface
  // perpendicular!!!
  tgt.tangent = glm::normalize( ray.d-glm::dot(ray.d,tgt.normal)*tgt.normal );
  tgt.bitangent = glm::cross(tgt.normal, tgt.tangent);

  return true;
}

RGB Shape::brdf(const Vec3& in, const Vec3& out) const
{
  return RGB(1.0f/(2.0f*3.141592654));
}
