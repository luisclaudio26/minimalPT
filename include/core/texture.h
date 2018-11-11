#ifndef TEXTURE_H
#define TEXTURE_H

#include <memory>
#include <string>
#include "spectrum.h"
#include "geometry.h"

class Texture
{
private:
public:
  typedef std::shared_ptr<Texture> ptr;

  virtual RGB sample(const Vec2& uv) const = 0;
  virtual std::string str() const = 0;
  virtual ~Texture() {}
};

#endif
