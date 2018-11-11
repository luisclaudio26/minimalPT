#ifndef JSONHELPER_H
#define JSONHELPER_H

#include "../core/geometry.h"
#include "../../3rdparty/json.hpp"

namespace JSONHelper
{
  Vec4 load_vec4(const nlohmann::json& in);
  Vec3 load_vec3(const nlohmann::json& in);
  Vec2 load_vec2(const nlohmann::json& in);
  Mat4 load_transformation(const nlohmann::json& in);
}

#endif
