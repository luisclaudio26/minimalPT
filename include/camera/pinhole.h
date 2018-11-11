#ifndef PINHOLE_H
#define PINHOLE_H

#include "../core/camera.h"
#include "../core/geometry.h"

class Pinhole : public Camera
{
private:
  Vec3 origin;
  Vec3 x, y, z; //camera coordinate system

  //TODO: fix units! better if in mm
  float ar, fw, fh, fd; //aspect ratio, film width, height and distance
  Vec3 film_bl, film_ur; //film bottom-left and upper-right corners

  Mat3 cam2world, world2cam;

public:
  Pinhole(const Vec3& origin, const Vec3 up, const Vec3& look_at,
          float film_width, float aspect_ratio, float focal_dist,
          int hRes, int vRes);

  Ray getRay(const Vec2& uv) const override;

  std::string str() const override
  {
    std::string out("");
    out += "Aspect ratio " + std::to_string(ar);
    out += ", film width " + std::to_string(fw);
    return out;
  }
};

#endif
