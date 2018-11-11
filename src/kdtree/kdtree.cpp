#include "../../include/kdtree/kdtree.h"
#include <glm/gtx/string_cast.hpp>
#include <algorithm>
#include <stack>

static const int MAX_PRIM = 15;
static const int MAX_DEPTH = 10;

static void compute_aabb(const std::vector<AABB>& aabbs,
                          const std::vector<int>& prims_ids, AABB& target)
{
  target.min = Vec3(FLT_MAX);
  target.max = Vec3(-FLT_MAX);

  for(auto id : prims_ids)
  {
    AABB p_aabb = aabbs[id]; //prims[id]->aabb(p_aabb);
    for(int i = 0; i < 3; ++i)
    {
      //TODO: is it necessary to take the minimum of
      //p_aabb's min and max corners? I don't think so
      target.min[i] = std::fmin(p_aabb.min[i], target.min[i]);
      target.max[i] = std::fmax(p_aabb.max[i], target.max[i]);
    }
  }

}

struct Edge
{
  float t; //position of this edge
  int prim; //to which primitive it belongs
  enum {START,END} type; //this edge is the right (start) or
                          //left (end) edge of the primitive

  Edge(float t, int prim, bool start) : t(t), prim(prim)
  {
    type = start ? START : END;
  }

  bool operator<(const Edge& rhs)
  {
    return t < rhs.t || (t == rhs.t && (int)type < (int)rhs.type);
  }
};

//-----------------------------------
//---------- FROM KDTREE.H ----------
//-----------------------------------
//DEPRECATED!
bool KdNode::intersect(const Ray& r, const std::vector<Triangle>& prims,
                        float tmin, float tmax, Isect& target) const
{
  if( is_leaf() )
  {
    bool hit = false;
    for(auto id : prims_ids)
    {
      Triangle t = prims[id];
      Isect p_isect; t.intersect(r, p_isect);

      if( p_isect.is_valid() && p_isect.t < target.t)
      {
        target = p_isect;
        hit = true;
      }
    }
    return hit;
  }

  //not a leaf
  //if ray is parallel to the splitting plane, there's no problem:
  //tsplit will be +INF, so it will be >= tmax (and won't be <= tmin),
  //thus we'll explore the NEAR box only, which is correct.
  //Problems may occur if the ray's origin is exactly on the splitting
  //plane, in which case we simply shift it a bit to the right.
  float split = this->split; if(split == r.o[this->axis]) split += 0.000001f;
  float tsplit = (split - r.o[this->axis])/r.d[this->axis];

  //define NEAR and FAR boxes.
  //this is defined by checking at which side of the
  //splitting plane the origin is. If the ray is exactly
  //on the origin, the direction defines which plane is near
  KdNode *near, *far;
  if( r.o[this->axis] < split ||
      (r.o[this->axis] == split && r.d[this->axis] < 0) )
  {
    near = this->left;
    far = this->right;
  }
  else
  {
    near = this->right;
    far = this->left;
  }

  //decide whether we must traverse NEAR only, FAR only or both
  //The case where both boxes should be visited (tmin < tsplit
  //&& tsplit < tmax) is equivalent to !near_only and !far_only.
  //Technically we should require that tsplit > 0 for far_only
  //to be true, but adding this clause won't change anything.
  //
  //Border cases, where tmin = tmax, are correctly also correctly
  //covered by this logic
  //
  //There's a subtlety here: everytime tsplit < 0, both near_only
  //and far_only are TRUE. However, as near_only is the first condition
  //tested, it returns true and we test only the nearest box, which is
  //CORRECT
  bool near_only = tsplit >= tmax || tsplit <= 0;
  bool far_only = tsplit <= tmin;

  if( near_only ) return near->intersect(r, prims, tmin, tmax, target);
  else if( far_only ) return far->intersect(r, prims, tmin, tmax, target);
  else
  {
    //if we found an intersection in the near
    //box, stop looking
    if( near->intersect(r, prims,tmin, tsplit,  target) ) return true;
    else return far->intersect(r, prims, tsplit, tmax, target);
  }
}

struct TraversalNode
{
  const KdNode* node;
  float tmin, tmax;

  TraversalNode(const KdNode* node, float tmin, float tmax)
    : node(node), tmin(tmin), tmax(tmax) {}
};

bool KdTree::intersect(const Ray& r, const std::vector<Triangle>& prims,
                        Isect& isect) const
{
  /*
  isect.t = FLT_MAX;
  float tmin, tmax;
  if( root.aabb.intersect(r, tmin, tmax) )
    return root.intersect(r, prims, tmin, tmax, isect);
  else return false; //we missed the biggest box; stop!
  */

  isect.t = FLT_MAX;
  float full_tmin, full_tmax;
  if( !root.aabb.intersect(r, full_tmin, full_tmax) ) return false;

  bool hit = false;
  std::stack<TraversalNode> stack;
  stack.push( TraversalNode(&root, full_tmin, full_tmax) );

  while( !stack.empty() )
  {
    //get node on top of stack
    TraversalNode tn = stack.top(); stack.pop();
    const KdNode* cur = tn.node;
    float tmin = tn.tmin, tmax = tn.tmax;

    //we can exit here because in this case we can guarantee that no
    //intersetion will happen closer than isect (because they will forcely
    //happen beyond tmin) 
    if( isect.t < tmin ) break;

    if( cur->is_leaf() )
    {
      for(auto id : cur->prims_ids)
      {
        const Triangle& t = prims[id];
        Isect p_isect; t.intersect(r, p_isect);

        if( p_isect.is_valid() && p_isect.t < isect.t)
        {
          isect = p_isect;
          hit = true;
        }
      }

      //we can't exit now because this intersection may be outside the
      //bounding box and thus another primitive inside one of the boxes
      //behind might contain a closer intersection
      //-- if(hit) break;
    }
    else
    {
      //if ray is parallel to the splitting plane, there's no problem:
      //tsplit will be +INF, so it will be >= tmax (and won't be <= tmin),
      //thus we'll explore the NEAR box only, which is correct.
      //Problems may occur if the ray's origin is exactly on the splitting
      //plane, in which case we simply shift it a bit to the right.
      float split = cur->split; if(split == r.o[cur->axis]) split += 0.0001f;
      float tsplit = (split - r.o[cur->axis])/r.d[cur->axis];

      //define NEAR and FAR boxes.
      //this is defined by checking at which side of the
      //splitting plane the origin is. If the ray is exactly
      //on the origin, the direction defines which plane is near
      KdNode *near, *far;
      if( r.o[cur->axis] < split ||
          (r.o[cur->axis] == split && r.d[cur->axis] < 0) )
      {
        near = cur->left;
        far = cur->right;
      }
      else
      {
        near = cur->right;
        far = cur->left;
      }

      //decide whether we must traverse NEAR only, FAR only or both
      //The case where both boxes should be visited (tmin < tsplit
      //&& tsplit < tmax) is equivalent to !near_only and !far_only.
      //Technically we should require that tsplit > 0 for far_only
      //to be true, but adding this clause won't change anything.
      //
      //Border cases, where tmin = tmax, are correctly also correctly
      //covered by this logic
      //
      //There's a subtlety here: everytime tsplit < 0, both near_only
      //and far_only are TRUE. However, as near_only is the first condition
      //tested, it returns true and we test only the nearest box, which is
      //CORRECT
      bool near_only = tsplit >= tmax || tsplit <= 0;
      bool far_only = tsplit <= tmin;

      if( near_only ) stack.push( TraversalNode(near, tmin, tmax) );
      else if( far_only ) stack.push( TraversalNode(far, tmin, tmax) );
      else
      {
        stack.push( TraversalNode(far, tsplit, tmax) );
        stack.push( TraversalNode(near, tmin, tsplit) );
      }
    }

  }

  return hit;
}

bool KdTree::intersect(const Ray& r, const std::vector<Triangle>& prims,
                        float t) const
{
    //identical to the regular traversal, except
    //it returns immediately if any primitive closer
    //than t was found and we do no backface culling
    //This could be merged with the first Intersect
    //method, but it would introduce needless, extra branches
    float full_tmin, full_tmax;
    if( !root.aabb.intersect(r, full_tmin, full_tmax) ) return false;

    std::stack<TraversalNode> stack;
    stack.push( TraversalNode(&root, full_tmin, full_tmax) );

    while( !stack.empty() )
    {
      //get node on top of stack
      TraversalNode tn = stack.top(); stack.pop();
      const KdNode* cur = tn.node;
      float tmin = tn.tmin, tmax = tn.tmax;

      if( cur->is_leaf() )
      {
        for(auto id : cur->prims_ids)
        {
          const Triangle& tri = prims[id];
          Isect p_isect; tri.intersect(r, p_isect, false);

          //intersection in point closer then r(t); this means
          //that the point r(t) is occluded!
          if( p_isect.is_valid() && p_isect.t < t ) return true;
        }
      }
      else
      {
        float split = cur->split; if(split == r.o[cur->axis]) split += 0.0001f;
        float tsplit = (split - r.o[cur->axis])/r.d[cur->axis];

        KdNode *near, *far;
        if( r.o[cur->axis] < split ||
            (r.o[cur->axis] == split && r.d[cur->axis] < 0) )
        {
          near = cur->left;
          far = cur->right;
        }
        else
        {
          near = cur->right;
          far = cur->left;
        }

        bool near_only = tsplit >= tmax || tsplit <= 0;
        bool far_only = tsplit <= tmin;

        if( near_only ) stack.push( TraversalNode(near, tmin, tmax) );
        else if( far_only ) stack.push( TraversalNode(far, tmin, tmax) );
        else
        {
          stack.push( TraversalNode(far, tsplit, tmax) );
          stack.push( TraversalNode(near, tmin, tsplit) );
        }
      }
    }

    //if we reached this point, no intersection closer than
    //r(t) was found
    return false;
}

bool KdNode::is_leaf() const { return left == NULL && right == NULL; }

bool KdNode::should_split(const std::vector<AABB>& aabbs,
                          const std::vector<int>& prims_ids,
                          const AABB& aabb, int& axis, float& t_split)
{
  //--------- Surface Area Heuristic ----------
  const float ISECT_COST = 70.0f;
  const float TRAV_COST = 1.0f;
  const float EMPTY_BONUS = 0.8f; //when the split leaves one of the children
                                  //completely empty, reduce the cost of this
                                  //split in EMPTY_BONUS percent

  int n_prims = prims_ids.size();
  float leaf_cost = n_prims * ISECT_COST;

  //1. compute area of the full bounding box
  Vec3 box_extent = aabb.max - aabb.min;
  float SA = 2.0f*(box_extent[0]*box_extent[1] +
                    box_extent[0]*box_extent[2] +
                    box_extent[1]*box_extent[2]);
  float invSA = 1.0f / SA;

  //2. for each AXIS
  float best_split_cost = FLT_MAX;
  float best_split_point;
  int best_split_axis = -1;

  for(int axis = 0; axis < 3; ++axis)
  {
    //2.1. put bounding boxes edges inside vector and sort it
    std::vector<Edge> edges;
    for(auto id : prims_ids)
    {
      AABB p_aabb = aabbs[id]; //prims[id]->aabb(p_aabb);
      edges.push_back( Edge(p_aabb.min[axis], id, true) );
      edges.push_back( Edge(p_aabb.max[axis], id, false) );
    }
    std::sort(edges.begin(), edges.end());

    //2.2. loop over the edges:
    int n_right = n_prims, n_left = 0;
    float split_cost = FLT_MAX, split_point;
    for(auto e : edges)
    {
      if(e.type == Edge::END) n_left++;

      //it will happen that split points will be outside
      //the bounding boxes boundaries, when primitive is
      //shared between two different boxes. we can't split
      //in a point outside the box!
      if( e.t > aabb.min[axis] && e.t < aabb.max[axis] )
      {
        //2.2.1 compute area of the children nodes
        int oa0 = (axis+1)%3, oa1 = (axis+2)%3; //the other axes besides
                                                //the current one

        float SL = 2.0f*(box_extent[oa0]*box_extent[oa1] +
                        (e.t - aabb.min[axis]) *
                        (box_extent[oa0]+box_extent[oa1]));

        float SR = 2.0f*(box_extent[oa0]*box_extent[oa1] +
                        (aabb.max[axis] - e.t) *
                        (box_extent[oa0]+box_extent[oa1]));

        float pL = SL*invSA, pR = SR*invSA;

        //2.2.2 compute cost of split
        float eb = (n_left == 0 || n_right == 0) ? EMPTY_BONUS : 1.0f;
        float cur_split_cost = TRAV_COST + (pL*n_left+pR*n_right)*ISECT_COST*eb;

        //2.2.3 check whether it is less than the best cost
        if( cur_split_cost < split_cost )
        {
          split_cost = cur_split_cost;
          split_point = e.t;
        }
      }

      if(e.type == Edge::START) n_right--;
    }

    //check whether splitting in this axis is better than the previous
    if( best_split_axis == -1 || split_cost < best_split_cost )
    {
      best_split_cost = split_cost;
      best_split_axis = axis;
      best_split_point = split_point;
    }
  }

  //3. check whether we should split or make a leafs
  t_split = best_split_point;
  axis = best_split_axis;
  return best_split_cost < leaf_cost;
}

void KdNode::make_leaf(const std::vector<int>& prims_ids, const AABB& aabb)
{
  this->left = this->right = NULL;
  this->aabb = aabb;
  this->axis = -1; this->split = 0.0f;
  this->prims_ids = std::move( prims_ids );
}

KdNode::KdNode(const std::vector<AABB>& aabbs,
                const std::vector<int>& prims_ids,
                const AABB& aabb, int recursion_depth)
{
  //stop criterion
  if(prims_ids.size() < MAX_PRIM || recursion_depth >= MAX_DEPTH)
    make_leaf( prims_ids, aabb );
  else
  {
    //choose splitting point
    float split; int axis;
    if( !should_split(aabbs, prims_ids, aabb, axis, split) )
      make_leaf(prims_ids, aabb);
    else
    {
      //classify primitives into left/right
      std::vector<int> left, right;
      for(auto id : prims_ids)
      {
        AABB p_aabb = aabbs[id];

        if(p_aabb.min[axis] <= split) left.push_back( id );
        if(p_aabb.max[axis] >= split) right.push_back( id );
      }

      //set internal node data
      this->split = split;
      this->axis = axis;
      this->aabb = aabb;

      //compute new bounding boxes
      AABB left_aabb = aabb;
      left_aabb.max[axis] = split;

      AABB right_aabb = aabb;
      right_aabb.min[axis] = split;

      //recursively build children
      this->left = new KdNode(aabbs, left, left_aabb, recursion_depth+1);
      this->right = new KdNode(aabbs, right, right_aabb, recursion_depth+1);
    }
  }

}

void KdTree::build(const std::vector<Triangle>& prims, AABB& bounds)
{
  //this vector will contain indices to all primitives
  //we want to store inside the tree
  std::vector<int> prims_ids;
  for(int i = 0; i < prims.size(); ++i) prims_ids.push_back(i);

  //precompute bounding boxes
  std::vector<AABB> aabbs; aabbs.reserve(prims.size());
  for(int i = 0; i < prims.size(); ++i)
    prims[i].aabb(aabbs[i]);

  //bounding box for the whole scene
  AABB aabb; compute_aabb(aabbs, prims_ids, aabb);
  bounds = aabb;

  root = KdNode(aabbs, prims_ids, aabb);
}
