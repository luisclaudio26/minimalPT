#include "../../include/core/sceneloader.h"
#include "../../include/core/triangle.h"
#include "../../include/core/material.h"
#include "../../include/core/texture.h"


bool SceneLoader::load_scene_from_json(const nlohmann::json& in)
{
  this->camera = Camera::load_from_json( in["camera"] );
  this->bgd = Background::load_from_json( in["background"] );

  nlohmann::json shapes = in["geometry"];
  for( auto shape : shapes )
    this->shapes.push_back( Shape::load_from_json(shape) );

  //log stuff
  printf("-----------------------\n");
  printf("Camera -> %s\n", this->camera->str().c_str());
  printf("Background -> %s\n", this->bgd->str().c_str());
  for( auto s : this->shapes ) printf("Shape -> %s\n", s->str().c_str());
  printf("-----------------------\n");

  return true;
}

void SceneLoader::generate_scene(Scene& target)
{
  //TODO: after fixing everything converning material/texture
  //loading, this method will become way more simple; we'll
  //simply loop over all the shapes and copy the pointers to
  //a single vector. In the end, materials/textures will be
  //pointer by pointers inside the primitives AND inside the
  //vectors. We need this to correctly deallocate the resources
  //(or maybe not, as shared_ptr's will do this automatically)
  //and maybe apply some operation to all materials/textures
  //at once

  //generate primitives and store all materials and textures
  //into target scene
  for( auto s : shapes )
  {
    s->generate_primitives( target.prims );
    for( auto m : s->materials ) target.materials.push_back(m);
    for( auto t : s->textures ) target.textures.push_back(t);
  }

  target.cam = this->camera;
  target.bgd = this->bgd;
}
