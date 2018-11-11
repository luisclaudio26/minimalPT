#ifndef MTL_H
#define MTL_H

#include "../core/material.h"
#include "../core/spectrum.h"
#include "../core/texture.h"
#include <glm/gtx/string_cast.hpp>

//Nothing special; this is just a material
//which properly implements the materials
//described in .mtl files using nice BxDFs
class MTL : public Material
{
public:
  int illum_model;

  RGB emission;

  //each component has an associated
  //texture ID, which is used to modulate it
  RGB ambient; Texture::ptr amb_tex;
  RGB diffuse; Texture::ptr diff_tex;
  RGB specular; Texture::ptr spec_tex; float shininess;
  RGB transmittance; float ior; // index of refraction

  Texture::ptr bump_tex, height_tex;

  // 1 == opaque; 0 == fully transparent
  float dissolve; Texture::ptr dissolve_tex;

  //BxDF's
  //Lambertian diff;
  //CookTorrance spec;

  RGB sample(const Vec3& wi, const Vec3& wo,
              const Vec3& normal, const Vec2& uv) const override;

  void sample_BSDF(const Vec2& uv, const Ray& wi, const Isect& isect,
                    Vec3& wo, float& wo_pdf, RGB& brdf) const override;

  bool is_emissive() const override
  {
    return emission.r > 0 || emission.g > 0 || emission.b > 0;
  }

  RGB emissivity() const override { return emission; }

  std::string str() const override
  {
    std::string out("");
    out += std::string("amb: ") + glm::to_string(ambient);
    out += std::string(" diff: ") + glm::to_string(diffuse);
    out += std::string(" spec: ") + glm::to_string(specular);
    return out;
  }
};

#endif
