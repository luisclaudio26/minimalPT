#include "../../include/util/imagetools.h"
#include <cstdio>
#include <cmath>

namespace ImageTools
{
  RGBuchar rgb_float_to_uchar(RGB p)
  {
  	RGBuchar out;
  	out.r = (unsigned char)std::fmin(255.0f, 255.0f * p.r);
  	out.g = (unsigned char)std::fmin(255.0f, 255.0f * p.g);
  	out.b = (unsigned char)std::fmin(255.0f, 255.0f * p.b);
  	return out;
  }

  void rgb_write_to_file(const char* filename, int w, int h, RGBuchar* data)
  {
  	FILE* out = fopen(filename, "wb");
  	fprintf(out, "P6\n%d %d\n255\n", w, h);
  	fwrite((const void*)data, sizeof(unsigned char), 3 * h*w, out);
  	fclose(out);
  }
}
