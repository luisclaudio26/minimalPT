#include "../../include/texture/imagetexture.h"
#include "../../include/util/log.h"
#include <glm/gtx/wrap.hpp>
#include <glm/gtx/string_cast.hpp>
#include <cstdlib>
#include <cstdio>

static FREE_IMAGE_FORMAT fif_from_str(const std::string& ext)
{
  if(ext.compare("png") == 0 || ext.compare("PNG") == 0) return FIF_PNG;
  else if(ext.compare("jpg") == 0 || ext.compare("JPG") == 0) return FIF_JPEG;
  else return FIF_UNKNOWN;
}

//--------------------------------------
//-------- FROM IMAGETEXTURE.H ---------
//--------------------------------------
RGB ImageTexture::sample(const Vec2& uv) const
{
  //wrap texture coordinates.
  //Vec2 uv_ = glm::repeat(uv); -> BROKEN, for some reason
  Vec2 uv_ = uv - glm::floor(uv); //this should take the fractional part of uv

  int i = (uv_.y == 1.0f) ? h-1 : uv_.y*h;
  int j = (uv_.x == 1.0f) ? w-1 : uv_.x*w;

  return RGB(data[i*w+j]);
}

ImageTexture::ImageTexture(const std::string& path)
{
  this->w = this->h = -1;
  this->data = NULL; this->img = NULL;

  size_t n = path.find_last_of('.');
  std::string ext = path.substr(n+1, path.npos);
  printf("Loading tex from %s\n", path.c_str());

  FREE_IMAGE_FORMAT fmt = fif_from_str(ext);
  if(fmt == FIF_UNKNOWN) ERROR("Unrecognized image format while loading texture")

  FIBITMAP* img_raw = FreeImage_Load(fmt, path.c_str(), 0);
  this->img = FreeImage_ConvertToRGBAF(img_raw);
  FreeImage_Unload(img_raw);

  this->w = FreeImage_GetWidth(this->img);
  this->h = FreeImage_GetHeight(this->img);

  //16 bytes (4 floats) per pixel, in RGBA order,
  //so we can cast it directly to RGBA
  this->data = (RGBA*)FreeImage_GetBits(this->img);
}

ImageTexture::~ImageTexture()
{
  FreeImage_Unload(this->img);
}
