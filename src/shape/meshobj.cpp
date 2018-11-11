#include "../../include/shape/meshobj.h"
#include "../../include/material/mtl.h"
#include "../../include/texture/imagetexture.h"
#include "../../include/util/log.h"
#include <cstdio>
#include <glm/gtx/string_cast.hpp>

//----------------------------------
//--------- FROM MESHOBJ.H ---------
//----------------------------------
void MeshOBJ::generate_primitives(std::vector<Triangle>& target) const
{
  //store triangles in world space
  for(auto t : tris)
  {
    target.push_back( t );
    target.back().transform(this->model2world);
  }
}

void MeshOBJ::load_geometry_data(const std::vector<tinyobj::shape_t>& shapes,
                                  const tinyobj::attrib_t& attrib)
{
  this->tris.clear();

  //load triangles for each shape
  for(auto s = shapes.begin(); s != shapes.end(); ++s)
  {
    unsigned int attrib_offset = 0;

    for(int f_id = 0; f_id < s->mesh.num_face_vertices.size(); ++f_id)
    {
      int f = s->mesh.num_face_vertices[f_id];
      Triangle face;

      //load vertices for this triangle
      //TODO: WILL FAIL IF MESH IS NOT MADE OF TRIANGLES!
      for(int v_id = 0; v_id < f; v_id++)
      {
        tinyobj::index_t v = s->mesh.indices[attrib_offset + v_id];

        float vx = attrib.vertices[3*v.vertex_index + 0];
        float vy = attrib.vertices[3*v.vertex_index + 1];
        float vz = attrib.vertices[3*v.vertex_index + 2];
        face.v[v_id] = Vec3( vx, vy, vz );

        if( !attrib.texcoords.empty() )
        {
          float tx = attrib.texcoords[2*v.texcoord_index + 0];
          float ty = attrib.texcoords[2*v.texcoord_index + 1];
          face.uv[v_id] = Vec2( tx, ty );
        }
      }

      //compute tangent and bitangent to this triangle
      face.compute_tangents();

      //load pointer to material. all shapes share the same materials!
      //WARNING: this won't work unless .obj files with no associated .mtl
      //file set all material IDs to zero!!! also, this is why we need
      //to load materials before
      int m_id = s->mesh.material_ids[f_id];
      if(m_id < 0 || m_id > this->materials.size()) printf("m_id = %d\n", m_id);

      face.material = this->materials[s->mesh.material_ids[f_id]];

      this->tris.push_back( face );
      attrib_offset += f;
    }
  }
}

std::string MeshOBJ::str() const
{
  std::string out("");
  out += std::to_string(tris.size()) + std::string(" triangles\n");

  for( auto m : materials )
    out += std::string("\t") + m->str() + std::string("\n");

  for( auto t : textures )
    out += std::string("\t") + t->str() + std::string("\n");

  return out;
}
