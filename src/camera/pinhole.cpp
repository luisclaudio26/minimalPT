#include "../../include/camera/pinhole.h"

Pinhole::Pinhole(const Vec3& origin, const Vec3 up, const Vec3& look_at,
                  float film_width, float aspect_ratio, float focal_dist,
                  int hRes, int vRes)
{
  this->hRes = hRes; this->vRes = vRes;
  this->ar = aspect_ratio;
  this->fw = film_width;
  this->fd = focal_dist;
  this->origin = origin;

  z = glm::normalize( origin - look_at );
  y = glm::normalize( up - glm::dot(up, z)*z );
  x = glm::normalize( glm::cross(y, z) );

  cam2world = Mat3(x, y, z);
  world2cam = glm::inverse(cam2world);

  //compute film corners in the camera's coordinate system
  fh = fw/ar;
  film_bl = Vec3( -fw/2, -fh/2, -fd );
  film_ur = Vec3( fw/2, fh/2, -fd );
}

Ray Pinhole::getRay(const Vec2& uv) const
{
  //u and v are assumed to be in range [-1,1] so we map it to [0,1]
  Vec2 uv_ = (uv + Vec2(1.0f)) * 0.5f;

  //this is not exactly screen space, as we have a 3D vector; it
  //is just to make calculations easier in the next step
  Vec3 p_screenspace = Vec3(uv_.x, uv_.y, -fd);

  //this product should be element-wise!
  Vec3 p_camspace = film_bl + p_screenspace * (film_ur - film_bl);

  //return ray in world space
  return Ray( origin, cam2world*glm::normalize(p_camspace) );
}
