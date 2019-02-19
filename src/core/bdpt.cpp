#include "../../include/core/integrator.h"
#include "../../include/core/sampler.h"
#include <cstdio>
#include <chrono>
#include <thread>
#include <vector>

typedef struct
{
  Isect isect;     // intersection info
  Vec3 pos;        // position of this vertex
  Ray last;        // the ray used to find this vertex
  float pdf_fwd;   // PDF (surface area) of being selected by the last vertex
  float pdf_bwd;   // PDF (surface area) of being selected by the next vertex
  bool valid;      // is this vertex valid or not?
} Vertex;

static inline float geometric_coupling(const Vertex& v1, const Vertex& v2)
{
  Vec3 r_ = v2.pos - v1.pos;
  Vec3 r = glm::normalize(r_);
  float d2 = glm::dot(r_, r_);

  float N1_r = glm::dot(v1.isect.normal, r);
  float N2_r = glm::dot(v2.isect.normal, -r);

  float out = (N1_r * N2_r) / d2;
  return out;
}

static inline RGB three_point_brdf(const Vertex& p_last,
                                    const Vertex& p_current,
                                    const Vertex& p_next)
{
  // TODO: what errors does this introduce? why are they so "correct"
  // and of random nature, as one would expect from numerical errors
  // more or less?
  //Vec3 wi = glm::normalize(p_last.pos - p_current.pos);
  //Vec3 wo = glm::normalize(p_next.pos - p_current.pos);

  Vec3 wi = -p_current.last.d;
  Vec3 wo = p_next.last.d;

  return p_current.isect.shape->brdf(wi, wo, p_current.pos);
}

static RGB connect_paths(int n_cam_vertices, int n_light_vertices,
                          const std::vector<Vertex>& vertices,
                          const Scene& scene, float& path_pdf)
{
  // UPDATE: the semantics of cam_vertex/light_vertex are now:
  // n_cam_vertices and n_light_vertices, meaning how many vertices
  // from each path we pretend to use (starting from the vertex on the lens
  // or on the light source). Thus, vertex_s stores the last cam vertex and
  // vertex_t stores the first light vertex.
  // TODO: throw error if vertex_s = -1, as we do not allow full light paths
  const int vertex_t = vertices.size() - n_light_vertices;
  const int vertex_s = n_cam_vertices - 1;

  // preset path_pdf to 1.0 to avoid NaN samples when returning 0
  path_pdf = 1.0f;

  // discard path if it is not complete!
  // TODO: we could still use it to compute paths of smaller lenghts!
  const Vertex &v_l = vertices[vertex_t];
  if( n_light_vertices > 0 && !v_l.valid ) return RGB(0.0f);

  const Vertex &v_c = vertices[vertex_s];
  if( !v_c.valid ) return RGB(0.0f);

  // -----------------------------------------
  // ----------- Visibility check ------------
  // -----------------------------------------
  // <editor-fold> Visibility check at connection
  // cast shadow ray to check whether vertices see each other
  // TODO: CRASHES WHEN CONNECTION HAPPENS BETWEEN LENS SAMPLE AND A LIGHT PATH VERTEX
  // TODO: o problema das quinas que não pegam o color bleed direito é exatamente
  // aqui na conexão, na estratégia s2t2!
  if( n_light_vertices > 0 )
  {
    // IDEA: This vertex offset position could be moved inside the shadow raycasting routing
    Vec3 ray_o = v_c.pos + v_c.isect.normal*(v_c.isect.shape->type == GLASS ? -0.0001f : 0.0001f);
    Vec3 ray_d = glm::normalize(v_l.pos - ray_o);
    Ray shadow_ray(ray_o, ray_d); Isect shadow_isect;
    scene.cast_ray(shadow_ray, shadow_isect);

    // two points do not see each other if they reach a surface from opposite sides.
    // this will happen everytime a light source is behind a thin wall: the shadow
    // ray will successfully connect the two paths, as intersection will happen at
    // the very same position, but the do not actually see each other!
    bool same_side = glm::dot(ray_d, v_l.isect.normal) < 0.0f;

    // TODO: the numerical problem near the corners is still happening, with rays
    // cast from the near the corner being missed
    // TODO: also due to numerical problems, setting distance >= 0.00001f makes
    // the "disks" appears just like the pathtracer implementation. this means that
    // the direct lighting is correct!
    // -> I have no idea what is going on here!
    Vec3 shadow_p = shadow_ray(shadow_isect.t);
    float dist = glm::distance(shadow_p, v_l.pos);

    if( !same_side || dist >= 0.001f ) return RGB(0.0f);
  }
  // </editor-fold>

  // ----------------------------------------------
  // ----------- Throughput computation -----------
  // ----------------------------------------------
  // <editor-fold> New throughput for MIS
  // 0. Create three pointers p_last, p_current and p_next which will
  //    store the working vertices. This is useful to make the loops
  //    easier to write, without hacking around to account for the light
  //    and camera vertices. p_last = NULL, p_current = lens, p_next depends
  //    on the number of camera vertices used (if we have only one camera vertex,
  //    the next one should be the first vertex of the light path)
  int id_next = (n_cam_vertices == 1) ? vertex_t : 1;

  const Vertex* p_last = NULL;
  const Vertex* p_current = &vertices[0];
  const Vertex* p_next = &vertices[id_next];

  RGB tp(1.0f);
  path_pdf *= vertices[0].pdf_fwd;

  // the actual halting criteria is: p_next = last vertex OR n_light_vertices = 0
  // and p_next = last camera vertex. Negate everything and you arrive in the
  // following expression.
  while( p_next != &vertices.back()
        && (n_light_vertices != 0 || p_next != &v_c) )
  {
    // advance pointers. next_id will take care of "jumping" between camera and
    // light vertices
    id_next = (id_next == vertex_s) ? vertex_t : id_next+1;
    p_last = p_current;
    p_current = p_next;
    p_next = &vertices[id_next];

    // Compute BRDF(p_last -> p_current -> p_next) and GCT(p_current -> p_next).
    // The use of BRDF and GCT caches will make this more efficient while keeping
    // code readable. This is certainly more expensive then storing things and
    // designing a more complex algorithm to handle things, but the goal is to
    // keep code easy to understand.
    RGB brdf = three_point_brdf(*p_last, *p_current, *p_next);
    float G = geometric_coupling(*p_current, *p_next);

    tp *= G * brdf;
    path_pdf *= p_current->pdf_fwd;
  }

  // PDF of the last vertex
  path_pdf *= p_next->pdf_fwd;

  // </editor-fold>

  // ----------------------------------------------
  // ----------- MIS weight computation -----------
  // ----------------------------------------------
  // <editor-fold> MIS weight


  // </editor-fold>

  RGB emission = n_light_vertices == 0 ?
                  vertices[vertex_s].isect.shape->emission
                    : vertices.back().isect.shape->emission;

  return tp * emission;
}

RGB Integrator::bd_path(const Scene& scene,
                        const Ray& primary_ray,
                        const Vec3& lens_normal,
                        const Isect& isect,
                        int path_length)
{
  // direct path from light source
  if(path_length == 2) return isect.shape->emission;

  // vertex storage
  std::vector<Vertex> vertices;
  vertices.resize(2*path_length);
  for( auto& v : vertices )
  {
    v.valid = false;
    v.pdf_fwd = 1.0f;
    v.pdf_bwd = -1.0f;
  }

  // ---------------------------------------------------
  // ------------------ Camera path --------------------
  // ---------------------------------------------------
  // <editor-fold> Camera path construction
  Vertex &lens_v = vertices[0];
  lens_v.valid = true;
  lens_v.pos = primary_ray.o;
  lens_v.isect.normal = lens_normal;
  lens_v.pdf_fwd = 1.0f; // TODO: place actual lens sample PDF!
  lens_v.pdf_bwd = -1.0f; // TODO: BRDF sample from first_v

  // we start by filling the second vertex cell (as the first EDGE is fixed due
  // to lens delta specular behavior).
  //
  // Its PDF is the only one which doesn't cut
  // out terms when we divide compute L(X)/p(x), thus we need the explicit
  // conversion from solid angle to surface area.
  // UPDATE: actually it does cancel out because there's also a geometric coupling
  // term that comes from the measurement equation. the jacobian of the PDF cuts
  // two cosine terms and leave one behind (cosine between normal at sensor and
  // the sampled direction).
  Vertex &first_v = vertices[1];
  first_v.valid = true;
  first_v.pos = primary_ray(isect.t);
  first_v.last = primary_ray;
  first_v.isect = isect;
  first_v.pdf_fwd = 1.0f; // the lens sample completely defines where the second
                          // vertex we fall i.e. p(x2) = p(x2 | x1) = 1.0 (??)
  first_v.pdf_bwd = -1.0f; //TODO: BRDF sample from vertices[2]

  // keep sampling BRDF and building camera subpath
  int cp_idx = 1;
  for(int i = 2; i < path_length; ++i)
  {
    Vertex &cur_v = vertices[cp_idx];

    // sample BRDF for outgoing direction
    float dir_pdf;
    Vec3 out_dir = cur_v.isect.shape->sample_brdf(cur_v.pos, -cur_v.last.d, dir_pdf);

    // cast ray in the sampled direction
    Vec3 ray_o = cur_v.pos + cur_v.isect.normal*(cur_v.isect.shape->type == GLASS ? -0.0001f : 0.0001f);
    Ray cur_to_next(ray_o, out_dir); Isect next_isect;
    if( !scene.cast_ray(cur_to_next, next_isect) ) break; //TODO: Ray escaped. Do something?

    // store new vertex
    cp_idx += 1;

    Vertex &next_v = vertices[cp_idx];
    next_v.valid = true;
    next_v.last = cur_to_next;
    next_v.isect = next_isect;
    next_v.pos = cur_to_next(next_isect.t);
    next_v.pdf_fwd = dir_pdf * glm::dot(next_isect.normal, -out_dir) / next_isect.d2;

    // backward probability: this is the PDF of the CURRENT vertex being chosen
    // by NEXT, while the forward probability is the PDF of the NEXT vertex being
    // chosen by the CURRENT one - all in a "BRDF sampling" sense.
    Vertex &last_v = vertices[cp_idx-2];
    float dir_pdf_bwd = cur_v.isect.shape->pdf_brdf(cur_to_next.d, -cur_v.last.d, cur_v.pos);

    // cur_v.isect.d2 => distance between last_v and cur_v
    // cur_v.last.d => direction of the ray linking last_v to cur_v (in this order)
    last_v.pdf_bwd = dir_pdf_bwd * glm::dot(last_v.isect.normal, cur_v.last.d) / cur_v.isect.d2;
  }

  // The work done on the loop does not computed backward probabilities for the
  // last two vertices on the camera path:
  // P( last vertex ) = 1.0f / TotalEmissiveArea, if it is a light source; 0 otherwise
  Vertex &last_cam = vertices[path_length-1];
  Vertex &last_but_one = vertices[path_length-2];

  if( !last_cam.isect.shape || last_cam.isect.shape->emission == RGB(0.0f) )
    last_cam.pdf_bwd = last_but_one.pdf_bwd = 0.0f;
  else
  {
    last_cam.pdf_bwd = 1.0f / scene.emissive_area();

    // TODO this assumes all lights are diffuse and sampled using uniform
    // hemisphere sampling! that's an awful lot of hard coded stuff
    const float pdf_uniform_hemisphere = _over2pi;
    last_but_one.pdf_bwd = pdf_uniform_hemisphere * glm::dot(last_but_one.isect.normal, last_but_one.last.d) / last_but_one.isect.d2;
  }

  // </editor-fold>

  // -------------------------------------------------------------
  // ----------------------- Light path --------------------------
  // -------------------------------------------------------------
  //<editor-fold> Light path construction
  // randomly pick a lightsource and sample a point on its surface
  // TODO: this will be biased if we have more than one light source!
  // must consider the case where there's more than one light source,
  // where the probability of picking a point is 1.0f / A, where A
  // is the total surface area of the light sources (and we importance
  // sample the area, also).
  const int l_idx = rand() % scene.emissive_prims.size();
  const Shape& l = scene.prims[ scene.emissive_prims[l_idx] ];
  Vec3 l_p; float l_p_pdf; // we're sampling radiance, thus a point
  Vec3 l_d; float l_d_pdf; // on the surface and a direction
  Vec3 l_n;
  l.sample_as_light(l_p, l_p_pdf, l_d, l_d_pdf, l_n);

  // unlike the lens sample, we'll need to store at least partial info on the
  // light sample because we need it explicitly for the radiance computation step.
  // Check the calculations on notebook to understand why this is also the only
  // PDF in light subpath construction that is expressed in surface area and not
  // solid angle (spoiler alert: all the remaining vertices will have their jacobians
  // cancelled out when we compute L(X)/p(x)).
  Vertex &last = vertices.back();
  last.valid = true;
  last.pdf_fwd = l_p_pdf;
  last.pdf_bwd = -1.0f;
  last.pos = l_p;
  last.isect.normal = l_n;
  last.isect.shape = &l;

  // cast ray from this sample and find the first intersection
  Ray l_ray(l_p + 0.0001f*l_n, l_d); Isect l_isect;

  //TODO: DO SOMETHING IF WE MISSED THIS RAY. Camera subpath has the same issue.
  // for the specific case of the light path having size 1, we could simply skip
  // and work with the light sample only
  if( !scene.cast_ray(l_ray, l_isect) ) return RGB(0.0f);

  // fill the last but one vertex
  Vertex &last_v = vertices[vertices.size()-2];
  last_v.valid = true;
  last_v.last = l_ray;
  last_v.isect = l_isect;
  last_v.pos = l_ray(l_isect.t);
  last_v.pdf_fwd = l_d_pdf * glm::dot(l_isect.normal, -l_d) / l_isect.d2;
  last_v.pdf_bwd = -1.0f;

  // build light subpath by sampling BRDFs.
  // TODO: maybe we won't need to go up to path_length but rather path_length-1,
  // as we are desconsidering "pure" light paths (this would imply having to
  // compute raster coordinates and stuff, which would need lots of modifications
  // on the code structure).
  int lp_idx = vertices.size()-2;
  for(int i = 2; i < path_length; ++i)
  {
    // get current vertex, its position and sample BRDF
    Vertex &cur_v = vertices[lp_idx];
    Vertex &last_v = vertices[lp_idx+1];

    float pdf_dir;
    Vec3 out_dir = cur_v.isect.shape->sample_brdf(cur_v.pos, -cur_v.last.d, pdf_dir);

    // cast ray in this direction to get next vertex
    Vec3 ray_o = cur_v.pos + cur_v.isect.normal * (cur_v.isect.shape->type == GLASS ? -0.0001f : 0.0001f);
    Ray next_ray(ray_o, out_dir); Isect next_isect;
    if( !scene.cast_ray(next_ray, next_isect) ) break; //TODO: ray escaped. do something?

    // store vertex info
    lp_idx -= 1;

    Vertex &next_v = vertices[lp_idx];
    next_v.valid = true;
    next_v.last = next_ray;
    next_v.isect = next_isect;
    next_v.pdf_fwd = pdf_dir * glm::dot(next_isect.normal, -out_dir) / next_isect.d2;
    next_v.pos = next_ray(next_isect.t);

    // from the second vertex on we can start computing the backward PDFs
    float dir_pdf_bwd = cur_v.isect.shape->pdf_brdf(next_v.last.d, -cur_v.last.d, cur_v.pos);
    last_v.pdf_bwd = dir_pdf_bwd * glm::dot(last_v.isect.normal, cur_v.last.d) / cur_v.isect.d2;
  }
  // </editor-fold>

  // The way we coded, the last two vertices will have their backward PDFs set to
  // -1. This is not a problem, as we never use the last vertices of a light path
  // because we need camera subpaths to be of length at least 2.
  // TODO: This is something we could optimize!

  // -------------------------------------------
  // ------------ Path connections -------------
  // -------------------------------------------
  // <editor-fold> Comments on MIS
  // TODO: the optimal way of doing this is to combine all paths using the
  // balance/power heuristics instead of simply averaging everything (which might
  // be causing the darkened look in refractive materials).
  // In order to do this, we can:
  // 1) compute all paths of all lenghts, just like the code already does
  // 2) if p_i(X) of a certain path X is zero, skip it (this will cause 0/0 problems)
  // 3) if not, accumulate the running sum of the p_i(X)'s.
  //
  //
  // I'M NOT SURE THAT THE POWER HEURISTICS USED IN PATHTRACING IS 100% CORRECT!
  // Specifically, the way probabilities are computed seems odd; we should be able
  // to compute the probability of choosing a given path by any of the sampling
  // strategies: thus, in the case of BDPT, if we have 5 path sampling strategies,
  // for example, when sampling one path we would need 5 PDFs? Check this!
  //
  // UPDATE: indeed, one must be able to pick a strategy (light sampling, let's
  // say), and not only sample a path with it, but also evaluate the probability
  // of picking a given path using this strategy.
  //
  // In order to correctly compute the weights for MIS, we need FORWARD and
  // REVERSE PDFs after incrementally building a path: we incrementally build a
  // path, accumulating the PDFs (in surface area) of each vertex; then we start
  // from the last vertex and compute the PDFs (in surface area) of picking each
  // vertex in reverse direction. Once we do this, we can compute the PDF os
  // sampling this very same path using any strategy. For example:
  //
  // Camera path:
  //      x0 -> x1 -> x2 -> x3 -> x4 -> x5
  // Light path:
  //      y5 <- y4 <- y3 <- y2 <- y1 <- y0
  //
  // Path X built using <3,3> connection strategy:
  //     x0 -> x1 -> x2 --- y2 <- y1 <- y0
  //
  // In a general sense, the PDF os picking this particular path is always:
  //     p(x0).p(x1).p(x2).p(y2).p(y1).p(y0)
  //
  // What changes between strategies, though, is the PDF os picking a given vertex
  // based on the previous choice. For example, the <3,3> strategy samples the
  // path in a way that the PDF of X is:
  //
  //     p(X) = p(x0).p(x1|x0).p(x2|x1) . p(y0).p(y1|y0).p(y2|y1)
  //
  // The <2,4> connection, however, would sample it in a way that p(X) would be:
  //
  //    p(X) = p(x0).p(x1|x0) . p(y0).p(y1|y0).p(y2|y1).p(x2|y2)
  //
  // Notice that we need to compute the PDF of picking x2 when in the vertex y2.
  // We need to check for visibility only once (the shadow ray cast in the first
  // part of the connection asserts that), then we need to lookup the BRDF at y2
  // to know what is the probability of picking x2 from then.
  //
  // Also, the <4,2> connection would be:
  //
  //    p(X) = p(x0).p(x1|x0).p(x2|x1).p(y2|x2) . p(y0).p(y1|y0)
  //
  // Which means we need to store p(y2|x2) and p(x2|y2). Indeed, the <5,1> strategy
  // would also require p(y1|y0) (and conversely, the <1,5> strategy would require
  // p(x0|x1)). Thus, once we build the camera and the light paths, we should
  // store all the probabilities in forward and reverse directions. The final
  // weight for balance heuristic needs all the probabilities after building the
  // path using a given connection strategy.
  // </editor-fold>

  RGB acc(0.0f); int n_strategies = 0;
  for(int c = 2; c <= path_length; ++c)
  {
    float path_pdf;
    RGB path_rad = connect_paths(c, path_length - c, vertices, scene, path_pdf);
    acc += path_rad * (1.0f / path_pdf);

    n_strategies += 1;
  }
  return acc * (1.0f / n_strategies);

  /*
  float path_pdf;
  RGB acc = connect_paths(2, 3, vertices, scene, path_pdf);
  return acc * (1.0f / path_pdf);
  */
}

RGB Integrator::bdpt(const Scene& scene,
                      const Ray& primary_ray,
                      const Vec3& lens_normal,
                      const Isect& isect)
{
  /*
  RGB out(0.0f);
  for(int i = 2; i <= 7; ++i)
    out += bd_path(scene, primary_ray, lens_normal, isect, i);
  return out;
  */

  return bd_path(scene, primary_ray, lens_normal, isect, 3);
}

// -----------------------------------------------------------------------------
// -------------------------- COMMENTS AND LEGACY ------------------------------
// -----------------------------------------------------------------------------
// QUESTION: how does one correctly weights those paths? remember that the
// importance sampling Monte Carlo integral is given by:
//
//  1   N   f(X_i)
// --- SUM --------
//  N   i   p(X_i)
//
// Which means we need to divide by the correct number of samples. This worked
// in the pathtracer above because each sample took only one path of each length;
// thus, f(X_i) was something of the form:
//
// f(X_i) = P_1(A_i) + P_2(B_i) + ... + P_n(C_i)
//
// where each P_i is an integral over paths of length i.
// When accumulated into the sample buffer and divided by N, we had exactly N
// samples for each integral, thus the estimative was right. In this case, our
// bd_path() function will collect many more paths as a single "sample".
//
// ANSWER: the naïve approach is to simply average the path contributions; the
// strong of BDPT, however, is to multiple importance sample them.
//
// Overall, bidirectional pathtracing is the "ultimate" multiple importance
// sampling technique; while in regular pathtracing we combine two paths only
// (one where the last vertex comes from BRDF sampling and another from light
// sampling), here we combine up to N different paths, which cover a variety
// of hard sampling cases (such as specular materials, delta light sources,
// caustic paths, hard indirect lighting paths, etc.)

// throughput at path junction
/*
RGB brdf_v_c = v_c.isect.shape->brdf(-v_c.last.d, v2l, v_c.pos);
float cosVC = glm::dot(v_c.isect.normal, v2l);
if( v_c.isect.shape->type == GLASS ) cosVC *= -1.0f;

RGB brdf_v_l = v_l.isect.shape->brdf(-v2l, -v_l.last.d, v_l.pos);
float cosVL = glm::dot(v_l.isect.normal, -v_l.last.d);
*/

/*
float cosVL = glm::dot(v_l.isect.normal, -v_l.last.d);
if( v_l.isect.shape->type == GLASS ) cosVL *= -1.0f;
*/

// throughput for the rest of the light path.
//
// we compute this throughput in the OPPOSITE sense then what is computed
// in the light subpath construction!
//
// also, there's no need for shadow ray the last vertex as it is guaranteed
// to be visible by construction.
/*
RGB tp_lp(1.0f); Ray from_cam = -shadow_ray;
for(int i = vertices.size()-2; i > vertex_t; --i)
{
  Vertex &cur = vertices[i];

  // update throughput
  RGB brdf = cur.isect.shape->brdf(from_cam.d, -cur.last.d, cur.pos);
  float cosTerm = glm::dot(cur.isect.normal, -cur.last.d);
  if( cur.isect.shape->type == GLASS ) cosTerm *= -1.0f;

  tp_lp *= cosTerm * brdf;

  // update variables for next loop
  from_cam = cur.last;
}
*/

// TODO: IT LOOKS LIKE MANY V_L'S ARE FALLING CLOSE TO THE SOURCE LIGHT, GIVING
// THE HARD SHADOW LOOK! Although it does not look right, conceptually it appears
// to be ok: half of the incoming light on all points come from the indirect light
// reflected on the ceiling.
// Does this mean, then, that using only one strategy for BDPT (for example,
// s = 2 and t = 2), will give us wrong results, and we need all strategies
// in order to have the whole thing correct? Pathtracing uses only BRDF sampling
// and light sampling, meaning that some paths will be way more rare than others
// (which will be compensated by the low PDF), so some strategies always leads
// to correct results whereas other do not??
//
// A PDF é que vai garantir que as diversas estratégias convergem pro mesmo canto! alguns
// caminhos para uma estratégia são muito mais prováveis que para outra; a PDF deve ponderar
// inversamente essas coisas de forma que elas convirjam para o mesmo resultado.
//
// a potência da fonte de luz pequena é a chave! ela emite pouca energia, por isso que a cena com as probabilidades "certas" é tão escura.
// a distância é a chave na hora de computar a PDF. os caminhos de luz que batem na
// parede verde/vermelha são mais distantes que os que batem na parede próxima: eles
// deviam então ser mais "improváveis" -> PDF menor, daí peso maior pra eles.
//
// no fim, o throughput dos caminhos não deve estar tão errado assim (se é que está,
// em primeiro lugar), a PDF é que deveria garantir que todo mundo é igual num
// sentido assintótico

// BIZU MÁXIMO!!!
/*
static int n_ceil = 0;
static int n = 0;
n++; float w = 1.0f;
if( l_isect.shape == &scene.prims[1] )
{
  n_ceil++;
  w = 0.1f;
}
printf("\rCeil hit rate %f", (float)n_ceil/n*100);
*/

// RAYS THAT DO NOT REACH THIS PART HAVE THE WRONG(?) PDFS! WHY?!
// This also points to the fact that most rays in the ceiling are failing to
// connect after! maybe because of the light source itself occluding it?
// !!!!!!! THE PROBLEM IS NOT THE PDF ITSELF
// Well, the problem seems not to be the PDF itself, but all points that can
// directly see the light source end up with a high contribution.
// This means that we can just terminate the ray computation here, as we were
// doing before.

/*
static float acc_pdf = 0.0;
static int n = 0;
acc_pdf += vertices[5].pdf * vertices[4].pdf; n++;
printf("\ravg pdf %f", acc_pdf/n);
*/

// TODO: numerical stability/precision on the method used to find the roots
// on the sphere intersection method might cause some flaws in corner regions (?),
// as well as "curved walls" (?).
// At corners, the lower the threshold in (dist < ...) section, the bigger the
// problem near the corners!
// -> it seems that numerical stability isn't really a concern here. altering
// thresholds doesn't seem to chance things. I honestly don't know what's going on

// ----------------- TEST LIGHT PATH ------------------
// 1. Link the first two vertices of the camera subpath with the last one. This
// should be equivalent to light sampling.
/*
RGB tp(1.0f); float path_pdf = 1.0f;
Vertex &v_last = vertices[vertices.size()-1];
Vertex &v = vertices[1];

// try to connect paths
Vec3 v2l_ = v_last.pos - v.pos;
float d2 = glm::dot(v2l_, v2l_);
Vec3 v2l = glm::normalize(v2l_);

// this seems to work, but may cause problems with the light sources themselves
Ray shadow_ray(v.pos, v2l); Isect shadow_isect;
scene.cast_ray(shadow_ray, shadow_isect);
float vis = glm::length(shadow_ray(shadow_isect.t)-v_last.pos) < 0.00001f ? 1.0f : 0.0f;

// vertices actually see each other. compute contribution.
RGB brdf = v.isect.shape->brdf(-v.last.d, v2l, v.pos);

float cosTerm = glm::dot(v.isect.normal, v2l);
if( v.isect.shape->type == GLASS ) cosTerm *= -1.0f;

      tp *= brdf * cosTerm;
path_pdf *= v_last.pdf * d2 / glm::dot(v_last.isect.normal, -v2l);

return vis*(v_last.isect.shape->emission * tp) * (1.0f/path_pdf);
*/

// !!!!! THE UNDERLYING ASSUMPTION IN THESE PATH CONSTRUCTIONS IS THAT PDF's
// ARE EXPRESSED IN SOLID ANGLE, NOT SURFACE AREA! THUS WE CANT CONVERT THESE
//last_v.pdf = l_d_pdf; // SHOULD WE EXPECT THE SAME RESULTS USING PDF IN
                        // SURFACE AREA AND SOLID ANGLE?! YES, BUT WE NEED THE
                        // PROPER CHANGE OF VARIABLES. the way we build the paths
                        // is expecting pdfs in solid angle, not surface area.
                        // if we wanna do everything in surface area pdfs, we'll
                        // need some extra cosine terms!!!
// once we pick a lens sample for a given pixel, there's probability 1 that we
// will choose a given ray (because of geometric optics). all the subsequente
// probabilities need to be expressed in surface area, thus we onvert pdfs in
// solid angle to it (this conversion is implicit in our path construction).
// the pdf in solid angle is all we have

// THERE ARE LOTS OF POINTS FALLING TOO CLOSE! THUS L_ISECT.D2 IS ALMOST ZERO
// AND, AS OF THE TERM 1.0/PDF, IT EXPLODES

// ----- TEST: build full path from the camera vertices -----
// compute path throughput and pdf
// Seems correct overall, but must take care of the escaping rays, which are
// doing some random thing!
/*
RGB tp(1.0f); float path_pdf = 1.0f;

for(int i = 1; i < vertices.size()-1; ++i)
{
  Vertex &cur = vertices[i];
  Vertex &next = vertices[i+1];

  // do not compute any contribution if this path has not the proper length
  if( !cur.valid || !next.valid ) tp = RGB(0.0f);

  // BRDF and cosine terms
  RGB brdf = cur.isect.shape->brdf(-cur.last.d, next.last.d, cur.pos);

  float cosTerm = glm::dot(cur.isect.normal, next.last.d);
  if( cur.isect.shape->type == GLASS ) cosTerm *= -1.0f;
  tp *= brdf * cosTerm;

  // path pdf is the product of the pdf of sampling each vertex.
  // TODO: include probability of the first vertex! (lens sample)
  path_pdf *= cur.pdf;
}

Vertex &v = vertices[path_length-1];
if( v.valid )
{
  RGB emission = v.isect.shape->emission;
  path_pdf *= v.pdf;

  return (tp * emission) * (1.0f / path_pdf);
}
else return RGB(0.0f);
*/
