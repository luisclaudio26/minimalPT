#ifndef SHAPE_H
#define SHAPE_H

#include "intersection.h"
#include "geometry.h"
#include "spectrum.h"

typedef enum {
  LAMBERTIAN, DELTA, GLASS
} MaterialType;

class Shape
{
private:
    float a;

public:
  Shape(const Vec3& o, float r)
    : o(o), r(r), diff_color(0.0f), emission(0.0f),
      type(LAMBERTIAN), eta(1.0f), a(4.0f*PI*r*r) { }

  // we'll handle spheres only at first!
  Vec3 o; float r;

  // material will be a simple Lambertian color
  RGB diff_color; MaterialType type; float eta;
  RGB brdf(const Vec3& in, const Vec3& out, const Vec3& p) const;
  Vec3 sample_brdf(const Vec3& p, const Vec3& in, float& pdf_solidangle) const;
  float pdf_brdf(const Vec3& in, const Vec3& out, const Vec3& p) const;

  // geometry related stuff
  inline float area() const { return a; }
  bool intersect(const Ray& r, Isect& tgt) const;
  void sample_surface(Vec3& point, Vec3& normal, float& pdf_area) const;

  // compute coordinate system around an specific normal.
  // we do this by simply taking a perpendicular vector on the plane
  // defined by the normal and the xaxis. If they're parallel, we use
  // the z axis. It returns a matrix that maps the local coordinate
  // system to the world coordinate system (because the glm::mat3
  // constructor acts column-wise, thus each column of the matrix is
  // one of the basis vectors).
  glm::mat3 get_local_coordinate_system(const Vec3& normal) const;

  // Emission profile should be a radiance
  // function Le(x,d) that tells us how much power
  // is emitted from this shape surface at point x and
  // direction d. By simply putting a constant, we
  // assume that radiance is uniform and isotropic
  // over the shape's surface (i.e. every point emits
  // the same power at every direction).
  //
  // This means for consistency reasons, when assigning
  // a value we should have in mind the irradiance of
  // each point (i.e. the sum of the radiance in all
  // directions for a given point) and the total power irradiated
  // of the surface (the sum of the irradiances of every point)
  // in order for quantities to make sense. This means that
  // simply assigning "1.0" (W.m⁻².sr⁻¹), for example, will
  // result in a sphere that yields WAY more power to
  // the scene than one might originally have intended to.
  //
  // Finally, keep in mind that integration of radiance for
  // a given point must be cosine weighted to account for
  // Lambert's cosine law.
  Vec3 emission; //in W.m⁻².sr⁻¹

  // sample this shape as an area light source.
  // TODO: include flags for diffuse/directional and more things related to
  // the lightfield emitted by the shape.
  void sample_as_light(Vec3& p, float& pdf_p,
                        Vec3& dir, float &pdf_dir, Vec3& n) const;
};

#endif
