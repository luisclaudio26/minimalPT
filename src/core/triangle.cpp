#include "../../include/core/triangle.h"

#include <iostream>
#include <glm/gtx/string_cast.hpp>

void Triangle::compute_tangents()
{
  Vec2 dUV1, dUV2;
  if(uv[0].x == 0.0f && uv[0].y == 0.0f &&
     uv[1].x == 0.0f && uv[1].y == 0.0f &&
     uv[2].x == 0.0f && uv[2].y == 0.0f)
  {
    dUV1 = Vec2(0.0f, 1.0f);
    dUV2 = Vec2(1.0f, 0.0f);
  }
  else
  {
    dUV1 = uv[1]-uv[0];
    dUV2 = uv[2]-uv[0];
  }

  Vec3 P0P2 = v[2]-v[0], P0P1 = v[1]-v[0];

  float invDet = 1.0f / dUV2.x*dUV1.y - dUV1.x*dUV2.y;
  t[0] = glm::normalize(invDet * (dUV1.y*P0P2 - dUV1.x*P0P1));
  t[1] = glm::normalize(invDet * (dUV2.x*P0P1 - dUV2.y*P0P2));
}

float Triangle::area() const
{
  //TODO: precompute this
  Vec3 v0v1 = this->v[1] - this->v[0];
  Vec3 v0v2 = this->v[2] - this->v[0];

  Vec3 cross = glm::cross(v0v1, v0v2);

  return std::fabs(glm::length(cross)) * 0.5f;
}

void Triangle::intersect(const Ray& ray, Isect& isect, bool bf_cull) const
{
  isect.t = -1.0f;

  // Möller-Trumbore algorithm, as described in
  // https://www.scratchapixel.com/lessons/3d-basic-rendering/
  // ray-tracing-rendering-a-triangle/moller-trumbore-ray-triangle-intersection
  Vec3 v0v1 = this->v[1] - this->v[0];
  Vec3 v0v2 = this->v[2] - this->v[0];
  Vec3 pvec = glm::cross(ray.d, v0v2);
  float det = glm::dot(v0v1, pvec);

  //TODO: backface culling
  if(bf_cull && det < 0.00000001) return;
  else if(!bf_cull && fabs(det) < 0.00000001 ) return;

  float invDet = 1 / det;

  Vec3 tvec = ray.o - this->v[0];
  float u = glm::dot(tvec, pvec) * invDet;
  if (u < 0 || u > 1) return;

  Vec3 qvec = glm::cross(tvec, v0v1);
  float v = glm::dot(ray.d, qvec) * invDet;
  if (v < 0 || u + v > 1) return;

  //compute local frame of reference
  Vec3 normal = glm::normalize(glm::cross(v0v1, v0v2));
  Vec3 bitangent = glm::cross( normal, t[0] );
  Mat3 world2local(t[0], normal, bitangent);

  isect.t = glm::dot(v0v2, qvec) * invDet;
  isect.normal = normal;
  isect.uv = u*uv[1] + v*uv[2] + (1-u-v)*uv[0];
  isect.tri = this;
  isect.local2world = world2local; //glm::inverse(world2local);
}

RGB Triangle::sample_emissive(Vec3& p, Vec3& normal, float& pdf) const
{
  //triangle sampling
  float u1 = (float)rand() / RAND_MAX;
  float u2 = (float)rand() / RAND_MAX;
  float x = 1 - sqrt(u1);
  float y = u2*sqrt(u1);
  float z = 1 - x - y;

  //normal will help us decide whether to discard this sample or not
  //without casting shadow rays
  //TODO: precompute normal!
  normal = glm::normalize(glm::cross(v[1]-v[0], v[2]-v[0]));
  p = x*v[0] + y*v[1] + z*v[2];
  pdf = 1.0f / area();
  return material->emissivity();
}

void Triangle::aabb(AABB& target) const
{
  for(int i = 0; i < 3; ++i)
  {
    target.max[i] = std::fmax(v[0][i], std::fmax(v[1][i], v[2][i]));
    target.min[i] = std::fmin(v[0][i], std::fmin(v[1][i], v[2][i]));
  }
}
