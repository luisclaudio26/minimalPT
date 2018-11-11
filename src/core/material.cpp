#include "../../include/core/material.h"
#include "../../include/material/mtl.h"
#include "../../include/texture/imagetexture.h"
#include "../../3rdparty/tiny_obj_loader.h"
#include <glm/gtc/type_ptr.hpp>
#include <fstream>

static std::string fix_path(const std::string& path)
{
  std::string out;
  for(int i = 0; i < path.length(); ++i)
    out += (path[i] == '\\') ? '/' : path[i];
  return out;
}

static std::string basedir_from_path(const std::string& path)
{
  std::size_t found = path.find_last_of('/');
  return path.substr(0, found+1);
}

//-----------------------------------
//-------- FROM MATERIAL.H ----------
//-----------------------------------
void Material::load_from_json(const nlohmann::json& in,
                              std::vector<Material::ptr>& target_mat,
                              std::vector<Texture::ptr>& target_tex)
{
  std::string type = in["type"].get<std::string>();
  nlohmann::json param = in["param"];

  if( type.compare("materialMTL") == 0 )
  {
    std::string path = param["pathMTL"].get<std::string>();
    std::fstream file(path);
    std::map<std::string, int> material_map;
    std::vector<tinyobj::material_t> mats;
    std::string err;

    tinyobj::LoadMtl(&material_map, &mats, &file, &err);

    std::size_t slash = path.find_last_of('/');
    std::string basedir = path.substr(0, slash+1);

    for(auto m = mats.begin(); m != mats.end(); ++m)
    {
      MTL *new_m = new MTL;
      new_m->illum_model = m->illum;
      new_m->emission = glm::make_vec3(m->emission);
      new_m->ambient = glm::make_vec3(m->ambient);
      new_m->diffuse = glm::make_vec3(m->diffuse);
      new_m->specular = glm::make_vec3(m->specular);
      new_m->transmittance = glm::make_vec3(m->transmittance);
      new_m->dissolve = m->dissolve;
      new_m->ior = m->ior;
      new_m->shininess = m->shininess;

      //load texture (if they exist)
      //TODO: this assumes texture is an image, but .obj
      //file supports procedural textures also
      if(!m->ambient_texname.empty())
      {
        Texture::ptr im_amb(new ImageTexture(basedir+fix_path(m->ambient_texname)));
        new_m->amb_tex = im_amb;
        target_tex.push_back(im_amb);
      }

      if(!m->diffuse_texname.empty())
      {
        Texture::ptr im_diff(new ImageTexture(basedir+fix_path(m->diffuse_texname)));
        new_m->diff_tex = im_diff;
        target_tex.push_back(im_diff);
      }

      if(!m->specular_texname.empty())
      {
        Texture::ptr im_spec(new ImageTexture(basedir+fix_path(m->specular_texname)));
        new_m->spec_tex = im_spec;
        target_tex.push_back(im_spec);
      }

      if(!m->bump_texname.empty())
      {
        Texture::ptr im_bump(new ImageTexture(basedir+fix_path(m->bump_texname)));
        new_m->bump_tex = im_bump;
        target_tex.push_back(im_bump);
      }

      if(!m->displacement_texname.empty())
      {
        Texture::ptr im_disp(new ImageTexture(basedir+fix_path(m->displacement_texname)));
        new_m->height_tex = im_disp;
        target_tex.push_back(im_disp);
      }

      if(!m->alpha_texname.empty())
      {
        Texture::ptr im_alpha(new ImageTexture(basedir+fix_path(m->alpha_texname)));
        new_m->dissolve_tex = im_alpha;
        target_tex.push_back(im_alpha);
      }

      target_mat.push_back( Material::ptr(new_m) );
    }
  }

}
