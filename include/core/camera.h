#ifndef CAMERA_H
#define CAMERA_H

#include "geometry.h"

class Camera
{
public:
  Camera();
  Camera(const Vec3& origin, const Vec3& up, const Vec3& look_at,
          float film_width, float aspect_ratio, float focal_dist);

  void compute_parameters(const Vec3& origin,
                          const Vec3& up,
                          const Vec3& look_at,
                          float film_width,
                          float aspect_ratio,
                          float focal_dist);

  //camera coordinate system
  Vec3 origin;
  Vec3 x, y, z;
  Mat4 cam2world, world2cam;

  //things measured in mm
  float ar, fw, fh, fd; //aspect ratio, film width, height and focal distance

  // film bottom-left and upper-right corners, computed from the above
  // parameters so we can easily compute rays
  Vec3 film_bl, film_ur;

  // get primary ray for sample (u,v) in [0,1]². for the simple case of a pinhole
  // camera, the primary can only have one direction; once we have a lens, the
  // same (u,v) will receive contribution from many rays coming from different
  // points of the lens
	Ray get_primary_ray(const Vec2& uv) const;

  void sample_lens(Vec3& pos_world, Vec2& pos_lens, float& pdf) const;
  Vec2 deposit_sample(const Vec2& pos_lens, const Vec3& dir_world) const;
};

#endif
