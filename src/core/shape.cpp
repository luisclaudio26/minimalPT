#include "../../include/core/shape.h"
#include "../../include/core/sampler.h"
#include <cstdio>

bool Shape::intersect(const Ray& ray, Isect& tgt) const
{
  // prefill tgt with invalid intersection data
  tgt.t = -1.0f;

  //a = 1, because r.d is normalized
  Vec3 oc = ray.o - o;
  float a = glm::dot(ray.d, ray.d);
  float b = 2*glm::dot(ray.d, oc);
  float c = glm::dot(oc,oc) - r*r;
  float delta = b*b-4*c;

  // no intersection!
  if(delta < 0.0f) return false;

  // care for negative intersections, which are
  // behind the camera (and thus are invalid)
  float t = (-b-sqrt(delta))/(2*a);          // closest intersection
  if( t < 0.0f ) t = (-b+sqrt(delta))/(2*a); // if negative, check further intersection

  // if still negative, the ball is behind the camera and we should return nothing.
  if(t < 0.0f) return false;

  // if we reached this point, t > 0.0f and thus we have a valid intersection.
  // but we are not checking whether we're inside the sphere or not!
  tgt.t = t;
  tgt.shape = this;
  tgt.normal = glm::normalize(ray(t) - o);

  Vec3 ot = ray.o - ray(t);
  tgt.d2 = glm::dot(ot, ot);

  // flip normals if the origin of the ray is inside the sphere:
  // ||oc|| < r  =>  <oc,oc> < r²  =>  <oc,oc> - r² < 0   =>   c < 0
  bool inside_sphere = c < 0.0f;
  if( inside_sphere ) tgt.normal = -tgt.normal;

  // compute an arbitrary tangent to the intersection point which will be
  // necessary for a local coordinate system.
  // At least for our spheres, it is proven that it is not possible to compute
  // Tangents without having singularities
  //
  //     https://en.wikipedia.org/wiki/Hairy_ball_theorem
  //
  // TODO: In this first version, our tangent will be simply the projection
  // of the ray onto the normal. thus, this will break if a ray hits the surface
  // perpendicular!!!
  tgt.tangent = glm::normalize( ray.d-glm::dot(ray.d,tgt.normal)*tgt.normal );
  tgt.bitangent = glm::cross(tgt.normal, tgt.tangent);

  return true;
}

RGB Shape::brdf(const Vec3& in, const Vec3& out, const Vec3& p) const
{
  Vec3 normal = glm::normalize(p-o);

  switch(type)
  {
    case GLASS:
    {
      bool outer_surface = glm::dot(in, normal) > 0.0f;
      Vec3 refracted = outer_surface ? glm::refract(-in,normal,eta) : glm::refract(-in,-normal,1.0f/eta);
      return glm::length(out-refracted) < 0.0001f ? RGB(1.0f) : RGB(0.0f);
    };

    case DELTA:
    {
      Vec3 reflected = glm::reflect(-in, normal);
      return glm::length(out-reflected) < 0.0001f ? RGB(1.0f) : RGB(0.0f);
    };

    default:
    case LAMBERTIAN:
    {
      // 1/pi because the reflectance integral is cosine-weighted, and the integral
      // of a cosine-weighted hemisphere is pi
      const float k = 1.0f/(3.141592654);
      return k * diff_color;
    }
  };
}

Vec3 Shape::sample_brdf(const Vec3& p, const Vec3& in, float& pdf_solidangle) const
{
  // normal is needed for all materials
  // QUESTION: how to flip normals in the case of sampling a BTDF
  // in the inner surface of the sphere and we need to send rays to the
  // outter one? Also we need this to invert the eta value
  Vec3 normal = glm::normalize(p-o);

  // select BRDF behavior
  switch(type)
  {
    case GLASS:
    {
      pdf_solidangle = 1.0f;

      // if IN is pointing in the same direction of the normal, this
      // means that intersection is occuring in the outer surface (remember
      // that our BRDF setting flips the incoming direction). Thus, if IN
      // points to the opposite side, this means that intersetion is ocurring
      // in the inner surface
      bool outer_surface = glm::dot(in, normal) > 0.0f;
      return outer_surface ? glm::refract(-in,normal,eta) : glm::refract(-in,-normal,1.0f/eta);
    }

    case DELTA:
    {
      pdf_solidangle = 1.0f;
      return glm::reflect(-in, normal);
    }

    default:
    case LAMBERTIAN:
    {
      Vec3 w_;
      Sampler::cosine_weight_sample_hemisphere(w_, pdf_solidangle);
      //Sampler::sample_hemisphere(w_, pdf_solidangle);
      glm::mat3 local2world = get_local_coordinate_system(normal);
      return local2world * w_;
    };
  }
}

float Shape::pdf_brdf(const Vec3& in, const Vec3& out, const Vec3& p) const
{
  // needed for all materials
  Vec3 normal = glm::normalize(p-o);

  switch(type)
  {
    case GLASS:
    {
      bool outer_surface = glm::dot(in, normal) > 0.0f;
      Vec3 refracted = outer_surface ? glm::refract(-in,normal, eta) : glm::refract(-in,-normal,1.0f/eta);
      return glm::length(out-refracted) < 0.0001f ? 1.0f : 0.0f;
    }

    case DELTA:
    {
      Vec3 reflected = glm::reflect(-in, normal);
      return glm::length(out-reflected) < 0.0001f ? 1.0f : 0.0f;
    }
    default:
    case LAMBERTIAN:
    {
      // TODO: check this. if we uniformly sample the hemisphere I think
      // this is correct, but if we cosine-weight it this might be wrong
      //return _overpi;
      return glm::dot(out, normal) * _overpi;
    }
  }
}

glm::mat3 Shape::get_local_coordinate_system(const Vec3& normal) const
{
  // get vector perpendicular to the normal and contained on the plane
  // defined by the normal and the x axis.
  // TODO: Problems occur if normal = axis!
  Vec3 axis(1.0f, 0.0f, 0.0f);
  float NdotA = glm::dot(axis, normal);

  Vec3 tangent = glm::normalize(axis - NdotA * normal);
  Vec3 bitangent = glm::cross(normal, tangent);

  // TODO: check this
  return glm::mat3(bitangent, normal, tangent);
}

void Shape::sample_surface(Vec3& point, Vec3& normal, float& pdf_area) const
{
  // uniformly sample point on a unitary sphere centered at the origin,
  // then scale and translate it so it lies on the surface of this sphere
  float u1 = (float)rand()/RAND_MAX;
  float u2 = (float)rand()/RAND_MAX;
  float y = 1 - 2*u1;
  float r = sqrt(std::max(0.0f, 1.0f - y * y));
  float phi = 2.0f*3.141592654f * u2;
  Vec3 q(r*cos(phi), y, r*sin(phi));

  normal = q;
  point = q*this->r + this->o;
  pdf_area = 1.0f / (4.0f * 3.141592654f * this->r * this->r);
}

void Shape::sample_as_light(Vec3& p, float& pdf_p,
                            Vec3& dir, float &pdf_dir,
                            Vec3& n) const
{
  // pick uniformly distributed point on sphere
  sample_surface(p, n, pdf_p);

  // pick uniformly distributed point on hemisphere of directions around p.
  // everything's uniform because we assume it to be a DIFFUSE light source!
  Vec3 dir_ts;
  Sampler::sample_hemisphere(dir_ts, pdf_dir);

  // get direction in world coordinates
  glm::mat3 local2world = get_local_coordinate_system(n);
  dir = local2world * dir_ts;
}
