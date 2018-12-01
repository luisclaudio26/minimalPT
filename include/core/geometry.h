#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <glm/glm.hpp>

typedef glm::vec3 Vec3;
typedef glm::vec4 Vec4;
typedef glm::vec2 Vec2;
typedef glm::mat4 Mat4;
typedef glm::mat3 Mat3;

const float PI = 3.141592654;
const float OVER_PI = 0.318309886;

//--------------------------------------------------------
//---------------- "High level" structures ---------------
//--------------------------------------------------------
class Ray
{
public:
  Vec3 o; Vec3 d;

  Ray(): o(Vec3(0,0,0)), d(Vec3(0,0,0)) {  }
  Ray(const Vec3& o, const Vec3& d) : o(o), d(d) { }

  Vec3 operator()(float t) const { return o + d*t; }
  Ray operator-() const { return Ray(o, -d); }
  Ray operator=(const Ray& r) { o = r.o; d = r.d; return *this; }
  Ray operator*(const Mat4& T) const
  {
    Ray out;
    out.o = Vec3( T * Vec4(o, 1.0f) );

    //TODO: allow pre-computed/uniform matrices!
    //QUESTION: WTF is that?
    Mat3 Tit = glm::transpose( glm::inverse( Mat3(T) ) );
    out.d = Vec3( Tit * d );

    return out;
  }
};

class AABB
{
public:
  Vec3 min, max;

  bool intersect(const Ray& r, float& tmin, float& tmax) const;
};

#endif
