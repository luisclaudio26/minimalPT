#include "../include/integrator/directillumination.h"
#include <glm/gtx/string_cast.hpp>

static RGB sampleOneLight(const Vec3& p, const Vec3& wo,
                            const Scene& scene, const Isect& isect)
{
  //randomly choose light to sample
  int n_lights = scene.emissive.size();
  int l_id = scene.emissive[rand() % n_lights];
  const Triangle *l = &scene.prims[l_id];

  //sample a single point on the surface of l, which is
  //an emissive surface (so we know it won't return zero).
  //TODO: how to actually compute this PDF? 1/area of the triangle
  //of 1/area of the whole scene?
  Vec3 p_light; Vec3 normal_light; float pdf_area;
  RGB emission = l->sample_emissive(p_light, normal_light, pdf_area);

  // discard sample if it comes from a light source which is not
  // facing the point.
  // TODO: I'm not sure this is right, but tests with Cornell box
  // went well.
  //if( glm::dot(normal_light, isect.normal) > 0.0f ) return RGB(0.f);

  //check whether this point is visible. if primitive is occluded,
  //the contribution of this sample is zero and we can skip it
  //TODO: handle the case where no intersection is found? Is it possible?
  Vec3 l2v = p - p_light;
  Vec3 wi = glm::normalize(l2v);
  Ray shadow(p, -wi); Isect shadow_isect;

  scene.intersect( shadow, shadow_isect );
  if( !shadow_isect.is_valid() || shadow_isect.tri != l ) return RGB(0.f);

  //primitive is unoccluded; compute the other terms of the integral
  float cosWiN = glm::dot(-wi, isect.normal);
  RGB brdf = isect.tri->material->sample(-wi, wo, isect.normal, isect.uv);

  //convert probability with respect to area to solid angle
  float r2 = glm::length(l2v); float cosWiNL = glm::dot(shadow_isect.normal, wi);
  float pdf_sa = pdf_area * r2 / cosWiNL;

  //final contribution of this sample
  RGB L = brdf * emission * cosWiN / pdf_sa;

  return L * (float)n_lights;
}

static RGB sampleAllLights(const Vec3& p, const Vec3& wo,
                            const Scene& scene, const Isect& isect)
{
  RGB out(0.0f); int n_samples = 15;

  //loop over all lights
  for(int l_id : scene.emissive)
  {
    const Triangle *l = &scene.prims[l_id];

    //Monte Carlo approximation of the integral over hemisphere
    //of directions
    RGB L(0.0f);
    for(int i = 0; i < n_samples; ++i)
    {
      //sample a single point on the surface of l, which is
      //an emissive surface (so we know it won't return zero).
      //TODO: how to actually compute this PDF? 1/area of the triangle
      //of 1/area of the whole scene?
      Vec3 p_light; Vec3 normal_light; float pdf_area;
      RGB emission = l->sample_emissive(p_light, normal_light, pdf_area);

      // discard sample if it comes from a light source which is not
      // facing the point.
      // TODO: I'm not sure this is right, but tests with Cornell box
      // went well -> turns out we cant, as the bunny scenes shows! why?
      //if( glm::dot(normal_light, isect.normal) > 0.0f ) continue;

      //check whether this point is visible. if primitive is occluded,
      //the contribution of this sample is zero and we can skip it
      //TODO: handle the case where no intersection is found? Is it possible?
      Vec3 l2v = p - p_light;
      Vec3 wi = glm::normalize(l2v);
      Ray shadow(p, -wi); Isect shadow_isect;

      scene.intersect( shadow, shadow_isect );
      if( !shadow_isect.is_valid() || shadow_isect.tri != l ) continue;

      //primitive is unoccluded; compute the other terms of the integral
      float cosWiN = glm::dot(-wi, isect.normal);
      RGB brdf = isect.tri->material->sample(-wi, wo, isect.normal, isect.uv);

      //convert probability with respect to area to solid angle
      float r2 = glm::length(l2v); float cosWiNL = glm::dot(shadow_isect.normal, wi);
      float pdf_sa = pdf_area * r2 / cosWiNL;

      //final contribution of this sample
      L += brdf * emission * cosWiN / pdf_sa;
    }

    out += L  * (1.0f / n_samples);
  }

  return out;
}

RGB DirectIllumination::integrate(const Vec2& uv, const Scene& scene) const
{
  Isect first_isect;
  Ray eye_ray = scene.cam->getRay(uv);

  //if we had no intersection, sample the background
  if( !scene.intersect(eye_ray, first_isect) )
    return scene.bgd->sample(eye_ray);

  //we had an intersection. compute lighting at this point
  RGB L = first_isect.tri->material->emissivity();

  //compute irradiance due to incident
  Vec3 p = eye_ray(first_isect.t) + first_isect.normal*0.000001f;
  L += sampleOneLight(p, -eye_ray.d, scene, first_isect);

  return L;
}
