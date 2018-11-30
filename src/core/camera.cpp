#include "../../include/core/camera.h"
#include "../../include/core/sampler.h"

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
  film_bl = Vec3(-fw/2, -fh/2, 0.0f);
  film_ur = Vec3( fw/2,  fh/2, 0.0f);
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
  // PINHOLE MODEL
  /*
  //this is not exactly screen space, as we have a 3D vector; it
  //is just to make calculations easier in the next step
  Vec3 p_screenspace = Vec3(uv.x, uv.y, 0.0f);

  //this product should be element-wise!
  Vec3 p_camspace = film_bl + p_screenspace * (film_ur - film_bl);

  //return ray in world space
  return Ray( origin, cam2world*glm::normalize(p_camspace) );
  */

  // THIS LENS MODEL
  // film sample in camera space
  Vec3 p_camspace( film_bl+Vec3(uv,0.0f)*(film_ur-film_bl) );

  // sample point on the lens, assumed to be at a distance fd from the
  // the film plane. Given the f-stop we can compute the lens radius (in mm)
  const float f_stop = 128.0f;
  const float lens_diameter = fd / f_stop;

  // uniform_sample_disk() is on [-0.5,0.5]². We thus multiply it by the
  // diameter to map it from [-r,r]² (where r is the lens radius)
  Vec2 q_lensspace; float pdf;
  Sampler::uniform_sample_disk(q_lensspace, pdf);
  Vec3 q_camspace(q_lensspace * lens_diameter, -fd);

  // TODO: compute ray direction using thin lens equation.
  // we flip x and y just to deinvert the inverted image
  Vec3 d = glm::normalize(q_camspace-p_camspace);
  d.x *= -1.0f; d.y *= -1.0f;

  //return ray in world space
  // TODO: THIS IS NOT STRICTLY CORRECT. Rays should originate not on the origin
  // but on their actual position on the film. This works because the film is
  // too small compared to the dimension of the objects on the scene.
  return Ray(origin, cam2world * d );
}
