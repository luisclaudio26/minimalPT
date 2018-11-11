#ifndef SOLID_H
#define SOLID_H

#include "../core/background.h"

class Solid : public Background
{
private:
  RGB color;

public:
  Solid(const RGB& color) : color(color) {}

  RGB sample(const Ray& r) const override { return this->color; }
  std::string str() const override
  {
      std::string out("solid ");
      out += std::to_string(color.x) + std::string(", ");
      out += std::to_string(color.y) + std::string(", ");
      out += std::to_string(color.z);
      return out;
  }
};

#endif
