#ifndef MESHOBJ_H
#define MESHOBJ_H

#include "../../3rdparty/tiny_obj_loader.h"
#include "../core/shape.h"
#include "../core/texture.h"
#include "../core/triangle.h"


class MeshOBJ : public Shape
{
public:
  std::vector<Triangle> tris;

  void load_geometry_data(const std::vector<tinyobj::shape_t>& shapes,
                          const tinyobj::attrib_t& attrib);

  void generate_primitives(std::vector<Triangle>& target) const override;

  std::string str() const override;
};

#endif
