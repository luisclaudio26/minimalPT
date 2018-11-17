#ifndef INTERSECTION_H
#define INTERSECTION_H

#include "geometry.h"

//forward declare Shape to avoid circularity, as it
//needs to know Isect for its intersect() method
class Shape;

class Isect
{
public:
  const Shape* shape; //DO NOT DELETE THIS POINTER!!!
  float t; float d2;  //t of the intersection, distance squared
  Vec3 normal, tangent, bitangent;


  //TODO: are t = INF intersections valid?
  bool is_valid() const { return t > 0 && t < FLT_MAX; }
};

#endif
