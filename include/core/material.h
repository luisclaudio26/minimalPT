#ifndef MATERIAL_H
#define MATERIAL_H

#include <memory>
#include <string>
#include "../../3rdparty/json.hpp"
#include "texture.h"
#include "spectrum.h"
#include "geometry.h"
#include "intersection.h"

//TODO: gpurenderer's materials are represented as a union
//here we're gonna use polymorphism, so one of the materials
//will be MTLmaterial, for example. I still don't know the
//virtual methods
class Material
{
public:
  typedef std::shared_ptr<Material> ptr;

  //collection of BRDF/BTDF (BSDF)
  //textures modulating the BSDF

  virtual bool is_emissive() const = 0;
  virtual RGB emissivity() const = 0;

  virtual RGB sample(const Vec3& wi, const Vec3& wo,
                      const Vec3& normal, const Vec2& uv) const = 0;

  //TODO: this should not be virtual. create a vector of vector
  virtual void sample_BSDF(const Vec2& uv, const Ray& wi, const Isect& isect,
                            Vec3& wo, float& wo_pdf, RGB& brdf) const = 0;

  //load material data from JSON file into material/texture vectors.
  //in case we need more than one material to be loaded (in the case of
  //.mtl files, for example, where many different materials are described),
  //we take target vectors
  static void load_from_json(const nlohmann::json& in,
                              std::vector<Material::ptr>& target_mat,
                              std::vector<Texture::ptr>& target_tex);

  virtual std::string str() const = 0;
};

#endif
