#ifndef KDTREE_H
#define KDTREE_H

#include "../core/geometry.h"
#include "../core/triangle.h"

class KdNode
{
public:
  AABB aabb;
  int axis; float split;
  KdNode *left, *right;

  //TODO: we're mixing Leaf and internal node
  //data. split this (union? inheritance?)
  std::vector<int> prims_ids;

  KdNode() : left(NULL), right(NULL) {}
  KdNode(const std::vector<AABB>& aabbs,
          const std::vector<int>& prims_ids,
          const AABB& aabb, int recursion_depth = 0);

  void make_leaf(const std::vector<int>& prims_ids, const AABB& aabb);

  bool should_split(const std::vector<AABB>& aabbs,
                    const std::vector<int>& prims_ids,
                    const AABB& aabb, int& axis, float& t_split);

  bool is_leaf() const;
  bool intersect(const Ray& r, const std::vector<Triangle>& prims,
                  float tmin, float tmax, Isect& target) const;
};

class KdTree
{
private:
  KdNode root;

public:
  void build(const std::vector<Triangle>& prims, AABB& bounds);
  bool intersect(const Ray& r, const std::vector<Triangle>& prims,
                  Isect& isect) const;

  bool intersect(const Ray& r, const std::vector<Triangle>& prims,
                  float t) const;
};

#endif
