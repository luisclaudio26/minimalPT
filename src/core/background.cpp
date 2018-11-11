#include "../../include/core/background.h"
#include "../../include/background/solid.h"
#include "../../include/util/jsonhelper.h"


Background::ptr Background::load_from_json(const nlohmann::json& in)
{
  std::string type = in["type"].get<std::string>();
  nlohmann::json param = in["param"];

  if( type.compare("solid") == 0 )
  {
    RGB color = JSONHelper::load_vec3( param["color"] );
    return Background::ptr( new Solid(color) );
  }
  else {} //TODO: log error

  return NULL;
}
