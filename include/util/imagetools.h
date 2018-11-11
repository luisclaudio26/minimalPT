#ifndef IMAGETOOLS_H
#define IMAGETOOLS_H

#include <glm/glm.hpp>
#include "../core/spectrum.h"

namespace ImageTools
{
  typedef struct {
    unsigned char r, g, b;
  } RGBuchar;

  RGBuchar rgb_float_to_uchar(RGB p);
  void rgb_write_to_file(const char* filename, int w, int h, RGBuchar* data);  
}

#endif
