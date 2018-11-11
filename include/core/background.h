#ifndef BACKGROUND_H
#define BACKGROUND_H

#include "spectrum.h"
#include "geometry.h"
#include "../../3rdparty/json.hpp"
#include <memory>
#include <string>

class Background
{
public:
  typedef std::shared_ptr<Background> ptr;

  virtual RGB sample(const Ray& r) const = 0;

  static Background::ptr load_from_json(const nlohmann::json& in);
  virtual std::string str() const = 0;
};

#endif
