#include "../../include/core/integrator.h"
#include "../../include/core/sampler.h"
#include <cstdio>
#include <chrono>
#include <thread>
#include <vector>

// THE FIRST TESTS SEEM TO BE OK OVERALL AND CONVERGENCE SPEED FOR THE HARD TEST
// SCENES IS REASONABLE. HOWEVER, IMAGES ARE A BIT DARKER THEN THEY SHOULD, MAINLY
// ON THE SPECULAR MATERIALS. MAYBE SOME MISSING COSINES OR DIVISIONS BY PI?

typedef struct
{
  Isect isect; // intersection info
  Vec3 pos;    // position of this vertex
  Ray last;    // the ray used to find this vertex
  float pdf;   // PDF (solid angle) of being selected by the last vertex
  bool valid;  // is this vertex valid or not?
} Vertex;

static RGB connect_paths(int cam_vertex, int light_vertex, const std::vector<Vertex>& vertices,
                          const Scene& scene, float& path_pdf)
{
  const int vertex_t = vertices.size() - 1 - light_vertex;
  const int vertex_s = cam_vertex;

  const Vertex &v_l = vertices[vertex_t];
  const Vertex &v_c = vertices[vertex_s];

  // preset path_pdf to 1.0 to avoid NaN samples when returning 0
  path_pdf = 1.0f;

  // discard path if it is not complete!
  // TODO: we could still use it to compute paths of smaller lenghts!
  if( !v_l.valid ) return RGB(0.0f); //return RGB(1.0f, 1.0f, 0.0f);
  if( !v_c.valid ) return RGB(0.0f); //return RGB(0.0f, 1.0f, 1.0f);

  // -----------------------------------------
  // ----------- Visibility check ------------
  // -----------------------------------------
  // cast shadow ray to check whether vertices see each other
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
  Vec3 shadow_p = shadow_ray(shadow_isect.t);
  float dist = glm::distance(shadow_p, v_l.pos);
  if( !same_side || dist >= 0.0001f ) return RGB(0.0f);

  // -----------------------------------------------
  // ---------- throughput - LIGHT PATH ------------
  // -----------------------------------------------
  // check the notebook to understand why this starts with a cosine term, but it
  // has something to do with regrouping the geometric terms. this is only the
  // case if vertex_t > 0 (pure camera paths do not need this)
  RGB tp_lp(1.0f);

  // this first term is only present if we have an actual light path (not only
  // using the light sample alone)
  if( light_vertex > 0 )
    tp_lp *= glm::dot(vertices.back().isect.normal, vertices[vertices.size()-2].last.d);

  for(int i = vertices.size()-2; i > vertex_t; --i)
  {
    const Vertex &v = vertices[i];
    Vec3 out_dir = vertices[i-1].last.d;

    RGB brdf = v.isect.shape->brdf(-v.last.d, out_dir, v.pos);

    float cosV = glm::dot(v.isect.normal, out_dir);
    if( v.isect.shape->type == GLASS ) cosV *= -1.0f;

    tp_lp *= brdf * cosV;
  }

  // ------------------------------------------------
  // ---------- throughput - CAMERA PATH ------------
  // ------------------------------------------------
  RGB tp_cp(1.0f);
  for(int i = 1; i < vertex_s; ++i)
  {
    const Vertex &v = vertices[i];
    Vec3 out_dir = vertices[i+1].last.d;

    RGB brdf = v.isect.shape->brdf(-v.last.d, out_dir, v.pos);

    float cosV = glm::dot(v.isect.normal, out_dir);
    if( v.isect.shape->type == GLASS ) cosV *= -1.0f;

    tp_cp *= brdf * cosV;
  }

  // ---------------------------------------
  // ------------ Full path PDF ------------
  // ---------------------------------------
  // run through the whole path accumulating its pdf
  float pdf_lp = 1.0f, pdf_cp = 1.0f;

  for(int i = vertices.size()-1; i >= vertex_t; --i)
    pdf_lp *= vertices[i].pdf;
  for(int i = 0; i <= vertex_s; ++i)
    pdf_cp *= vertices[i].pdf;

  path_pdf = pdf_lp * pdf_cp;

  //-----------------------------------------------
  // ---------- throughput - CONNECTION -----------
  //-----------------------------------------------
  // There's no PDF for the junction! the probability of picking this specific
  // path depends only on the probability of selecting the vertices. there's no
  // meaning in computing PDF for the junction, as once we build the subpaths,
  // we have chance 1 of connecting them, let's say (?)
  Vec3 v2l_ = v_l.pos - v_c.pos;
  Vec3 v2l = glm::normalize(v2l_);

  // geometric coupling term and BRDFs for the connection point
  // TODO: care for the GLASS BRDFs here!
  float G_cosVC = glm::max(0.0f, glm::dot(v_c.isect.normal, +v2l));
  if( v_c.isect.shape->type == GLASS ) G_cosVC *= -1.0f;

  float G_cosVL = glm::max(0.0f, glm::dot(v_l.isect.normal, -v2l));
  if( v_l.isect.shape->type == GLASS ) G_cosVL *= -1.0f;

  float G_r2 = glm::dot(v2l_, v2l_);
  float G = (G_cosVC * G_cosVL) / G_r2;

  RGB brdfVL = v_l.isect.shape->brdf(-v2l, -v_l.last.d, v_l.pos);
  RGB brdfVC = v_c.isect.shape->brdf(+v2l, -v_c.last.d, v_c.pos);

  // full path throughput
  RGB tp(1.0f);
  if( light_vertex == 0 )
    tp = tp_cp * (brdfVC * G);
  else
    tp = tp_cp * (brdfVC * G * brdfVL) * tp_lp;

  return tp * vertices.back().isect.shape->emission;
}

RGB Integrator::bd_path(const Scene& scene,
                        const Ray& primary_ray,
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
    v.pdf = 1.0f;
  }

  // ---------------------------------------------------
  // ------------------ Camera path --------------------
  // ---------------------------------------------------
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
  first_v.pdf = 1.0f;

  // keep sampling BRDF and building camera subpath
  int cp_idx = 1;
  for(int i = 2; i < path_length; ++i)
  {
    Vertex &v = vertices[cp_idx];

    // sample BRDF for outgoing direction
    float dir_pdf;
    Vec3 out_dir = v.isect.shape->sample_brdf(v.pos, -v.last.d, dir_pdf);

    // cast ray in the sampled direction
    Vec3 ray_o = v.pos + v.isect.normal*(v.isect.shape->type == GLASS ? -0.0001f : 0.0001f);
    Ray to_next_point(ray_o, out_dir); Isect next_isect;
    if( !scene.cast_ray(to_next_point, next_isect) ) break; //TODO: Ray escaped. Do something?

    // store vertex. PDFs are stored in solid angle, as the
    // jacobian is cancelled when computing L(x)/pdf(x).
    cp_idx += 1;

    Vertex &next_v = vertices[cp_idx];
    next_v.valid = true;
    next_v.pdf = dir_pdf;
    next_v.last = to_next_point;
    next_v.isect = next_isect;
    next_v.pos = to_next_point(next_isect.t);
  }

  // -------------------------------------------------------------
  // ----------------------- Light path --------------------------
  // -------------------------------------------------------------
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
  last.pdf = l_p_pdf;
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
  last_v.pdf = l_d_pdf;

  // build light subpath by sampling BRDFs.
  // TODO: maybe we won't need to go up to path_length but rather path_length-1,
  // as we are desconsidering "pure" light paths (this would imply having to
  // compute raster coordinates and stuff, which would need lots of modifications
  // on the code structure).
  int lp_idx = vertices.size()-2;
  for(int i = 2; i < path_length; ++i)
  {
    // get current vertex, its position and sample BRDF
    Vertex &v = vertices[lp_idx];

    float pdf_dir;
    Vec3 out_dir = v.isect.shape->sample_brdf(v.pos, -v.last.d, pdf_dir);

    // cast ray in this direction to get next vertex
    Vec3 ray_o = v.pos + v.isect.normal * (v.isect.shape->type == GLASS ? -0.0001f : 0.0001f);
    Ray next_ray(ray_o, out_dir); Isect next_isect;
    if( !scene.cast_ray(next_ray, next_isect) ) break; //TODO: ray escaped. do something?

    // store vertex info
    lp_idx -= 1;

    Vertex &next_v = vertices[lp_idx];
    next_v.valid = true;
    next_v.last = next_ray;
    next_v.isect = next_isect;
    next_v.pdf = pdf_dir;
    next_v.pos = next_ray(next_isect.t);
  }

  // -------------------------------------------
  // ------------ Path connections -------------
  // -------------------------------------------
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

  RGB acc(0.0f); int n_strategies = 0;
  for(int c = 1; c < path_length-1; ++c)
  {
    int l = path_length - 2 - c;

    float path_pdf;
    RGB path_rad = connect_paths(c, l, vertices, scene, path_pdf);

    acc += path_rad * (1.0f / path_pdf);
    n_strategies += 1;
  }

  return acc * (1.0f / n_strategies);
}

RGB Integrator::bdpt(const Scene& scene,
                      const Ray& primary_ray,
                      const Isect& isect)
{
  RGB out(0.0f);
  for(int i = 3; i < 7; ++i) out += bd_path(scene, primary_ray, isect, i);
  return out;
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
