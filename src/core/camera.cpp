#include "../../include/core/camera.h"

void Camera::compute_parameters(const Vec3& origin,
                                const Vec3& up,
                                const Vec3& look_at,
                                float film_width,
                                float aspect_ratio,
                                float focal_dist)
{
  this->ar = aspect_ratio;
  this->fw = film_width;
  this->fd = focal_dist;
  this->origin = origin;

  z = glm::normalize( origin - look_at );
  y = glm::normalize( up - glm::dot(up, z)*z );
  x = glm::normalize( glm::cross(y, z) );

  cam2world = Mat3(x, y, z);
  world2cam = glm::inverse(cam2world);

  // compute film corners in the camera's coordinate system.
  // even though things are in millimeters, the ratios cancel
  // out the units, so we should have no problems.
  fh = fw/ar;
  film_bl = glm::normalize( Vec3(-fw/2, -fh/2, -fd) );
  film_ur = glm::normalize( Vec3( fw/2,  fh/2, -fd) );
}

Camera::Camera(const Vec3& origin, const Vec3& up, const Vec3& look_at,
                  float film_width, float aspect_ratio, float focal_dist)
{
  compute_parameters(origin, up, look_at, film_width, aspect_ratio, focal_dist);
}

Camera::Camera()
{
  compute_parameters( Vec3(0.0f, 0.0f, 0.0f),
                      Vec3(0.0f, 1.0f, 0.0f),
                      Vec3(0.0f, 0.0f, -1.0f),
                      35.0f,
                      4.0/3.0f,
                      35.0f);
}

Ray Camera::get_primary_ray(const Vec2& uv) const
{
  //this is not exactly screen space, as we have a 3D vector; it
  //is just to make calculations easier in the next step
  Vec3 p_screenspace = Vec3(uv.x, uv.y, -fd);

  //this product should be element-wise!
  Vec3 p_camspace = film_bl + p_screenspace * (film_ur - film_bl);

  //return ray in world space
  return Ray( origin, cam2world*glm::normalize(p_camspace) );
}
