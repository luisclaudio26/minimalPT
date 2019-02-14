#include "../../include/core/camera.h"
#include "../../include/core/sampler.h"

#include <cstdio>

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

  // matrix construction is column-wise!
  cam2world = Mat4( Vec4(x, 0.0f),
                    Vec4(y, 0.0f),
                    Vec4(z, 0.0f),
                    Vec4(origin, 1.0f));
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

Ray Camera::get_primary_ray(const Vec2& uv_) const
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
  // flip axes in order to deinvert image
  Vec2 uv(1.0f-uv_.x, 1.0f-uv_.y);

  // sample point on the lens, assumed to be at a distance fd from the
  // the film plane. Given the f-stop we can compute the lens radius (in mm)
  const float f_number = 64.0f;
  const float lens_diameter = fd / f_number;
  const float focus_point = 500.0f; //in mm

  // uniform_sample_disk() is on [-0.5,0.5]². We thus multiply it by the
  // diameter to map it from [-r,r]² (where r is the lens radius)
  Vec2 q_lensspace; float pdf;
  Sampler::uniform_sample_disk(q_lensspace, pdf);
  Vec3 q_camspace(q_lensspace * lens_diameter, 0.0f);

  // film sample in camera space
  float img_dist = 1.0f / (1.0f/fd - 1.0f/focus_point);
  Vec3 p_camspace( film_bl+Vec3(uv,0.0f)*(film_ur-film_bl) );
  p_camspace.z = img_dist;

  // ray direction is simply [0 0 0] - p_camspace (lens center is at origin)
  Vec3 d = glm::normalize(-p_camspace);
  Ray main_ray(p_camspace, d);

  float t = -(img_dist+focus_point)/d.z;
  Vec3 cp = main_ray( t );

  Vec3 d_out = glm::normalize(cp - q_camspace);

  //return ray in world space
  // TODO: THIS IS NOT STRICTLY CORRECT. Rays should originate not on the origin
  // but on their actual position on the film. This works because the film is
  // too small compared to the dimension of the objects on the scene.
  #define mm_to_m(x) (x*0.001f)

  Vec4 o_worldspace( cam2world * Vec4(mm_to_m(q_camspace),1.0f) );
  Vec4 d_worldspace( cam2world * Vec4(d_out,0.0f) );

  return Ray( Vec3(o_worldspace), Vec3(d_worldspace) );
}

Ray Camera::get_primary_ray(const Vec2& uv_, Vec3& normal) const
{
  // This is a thin lens model, so the lens is infinitely thin and thus
  // the normal is the same everywhere. realistic lens would not be like this
  normal = -z;
  return get_primary_ray(uv_);
}

void Camera::sample_lens(Vec3& pos_world, Vec2& pos_lens, float& pdf) const
{
  // TODO: precompute and store this
  // sample point on the lens, assumed to be at a distance fd from the
  // the film plane. Given the f-stop we can compute the lens radius (in mm)
  const float f_number = 64.0f;
  const float lens_diameter = fd / f_number;
  const float focus_point = 500.0f; //in mm
  #define mm_to_m(x) (x*0.001f)

  // camera lens space -> camera space -> world space
  Sampler::uniform_sample_disk(pos_lens, pdf);
  Vec3 pos_camera(pos_lens * lens_diameter, 0.0f);
  pos_world = Vec3( cam2world * Vec4(mm_to_m(pos_camera), 1.0f) );
}

Vec2 Camera::deposit_sample(const Vec2& pos_lens, const Vec3& dir_world) const
{
  const float f_number = 64.0f;
  const float lens_diameter = fd / f_number;
  const float focus_point = 500.0f; //in mm
  const float img_dist = 1.0f / (1.0f/fd - 1.0f/focus_point);

  // get dir_world and pos_lens both in camera space
  Vec3 dir_cam = Vec3(world2cam * Vec4(dir_world, 0.0f));
  Vec3 pos_cam = Vec3(pos_lens * lens_diameter, 0.0f);

  // get intersection of the ray (pos_cam, dir_cam) with the focal plane
  float t = -(fd + pos_cam.z)/dir_cam.z;
  Vec3 p = pos_cam + t*dir_cam;

  // get position of the projected sample in camera space
  // (those are elementwise product/division operations) by
  // pinhole-projecting the point p onto the image plane
  Vec3 sample_cam = (img_dist * p) / -focus_point;

  // the point p must be a linear interpolation between
  // film_bl and film_ur (independent in each axis)
  Vec3 a( (sample_cam-film_bl)/(film_ur-film_bl) );

  return Vec2(a.x, a.y);
}
