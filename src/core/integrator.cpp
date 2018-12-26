#include "../../include/core/integrator.h"
#include "../../include/core/sampler.h"
#include <cstdio>
#include <chrono>
#include <thread>
#include <vector>

using namespace std::chrono;

// --------------------------------------------------------------------
// ------------------------- Testing stuff ----------------------------
// --------------------------------------------------------------------

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



// ----------------------------------------------------------

// --------------------------------------------------------------------

// --------------------------------------------------------------------
// ------------------ Pathtracing helper functions --------------------
// --------------------------------------------------------------------
static RGB sample_light(const Ray& primary_ray,
                          const Isect& isect,
                          const Scene& scene,
                          float& pdf_solidangle)
{
  Vec3 p(primary_ray(isect.t) + isect.normal*(isect.shape->type == GLASS ? -0.0001f : 0.0001f));
  Vec3 d = -primary_ray.d;

  // randomly pick a lightsource...
  const int idx = rand() % scene.emissive_prims.size();
  const Shape& l = scene.prims[ scene.emissive_prims[idx] ];

  // ... compute a Monte Carlo estimative for the exitant radiance on the point
  // p in the direction of the input ray
  RGB exitant_radiance(0.0f);

  // sample a point on the surface
  Vec3 q, n;       // sampled point q and normal at q
  float pdf_area;  // probability of choosing q on the surface of l
  l.sample_surface(q, n, pdf_area);

  // compute radiance coming from q in direction w = q - p
  // TODO: in this first version, emission is always isotropic and uniform,
  // so every direction at every point on the surface will have the same radiance
  Vec3 w = q-p;
  Vec3 w_n = glm::normalize(w);
  float r2 = glm::dot(w,w);
  RGB L_qw = l.emission;

  // cast shadow ray to compute visibility; we trace a ray r that simply
  // returns the t of the closest intersection; if r(t) != q, then q is
  // occluded.
  // TODO: an specific visibility method may be less complex than the usual
  // (and expensive) cast_ray() method.
  //Ray shadow_ray(p + 0.0001f*isect.normal, w_n);
  Ray shadow_ray(p, w_n);
  Isect isect_shadow;
  scene.cast_ray(shadow_ray, isect_shadow);
  float v = glm::length(shadow_ray(isect_shadow.t)-q) < 0.00001f ? 1.0f : 0.0f;

  // return a single the contribution f(x) and the pdf p(x) of choosing this point.
  // We clamp the PDF to zero in case glm::dot(n, -w_n) is negative because even
  // the radiance of this path is zero (because of v), it will break the MIS computation
  // which depends on this value to be zero.
  //
  // TODO: It seems that some radiance paths have negative PDF (which is expected
  // when the sampled point is in the opposite direction of the sphere) but with
  // visibility != 0.0. This is directly related to the way we compute v, as the
  // error occurs with -w_n.n always being less then 0.0001.
  RGB rad = v * glm::dot(w_n,isect.normal) * isect.shape->brdf(w_n, -d, p) * L_qw;
  pdf_solidangle = pdf_area * r2 / glm::dot(n, -w_n);

  return rad;
}

static RGB sample_brdf(const Ray& primary_ray,
                        const Isect& isect,
                        const Scene& scene,
                        float& pdf_solidangle)
{
  RGB rad_out(0.0f);

  // intersection point
  Vec3 p(primary_ray(isect.t) + isect.normal*(isect.shape->type == GLASS ? -0.0001f : 0.0001f));

  // sample BRDF direction
  Vec3 brdf_dir = isect.shape->sample_brdf(p, -primary_ray.d, pdf_solidangle);

  // cast ray, check whether it hits a light source
  Ray brdf_ray(p, brdf_dir); Isect isect_brdf;
  if( scene.cast_ray(brdf_ray, isect_brdf) )
  {
    float cosNW = glm::dot(isect.normal, brdf_dir);
    if(isect.shape->type == GLASS) cosNW *= -1.0f;

    //QUESTION: wouldn't it be necessary to consider the cosine term of the radiance going out?
    rad_out = isect_brdf.shape->emission
                * isect.shape->brdf(brdf_dir, -primary_ray.d, p)
                * cosNW;
  }
  else
    rad_out = RGB(0.0f); // TODO: environment contribution

  return rad_out;
}

// --------------------------------------------------------------------
// ------------------ From integrator.h -------------------------------
// --------------------------------------------------------------------
Integrator::Integrator()
  : vRes(600), hRes(800)
{
  // initialize buffers
  int n = vRes*hRes;
  samples.insert(samples.begin(), n, RGB(0.0f,0.0f,0.0f));
  weights.insert(weights.begin(), n, 0.0f);
    frame.insert(  frame.begin(), n, RGBA(0.1f,0.1f,0.1f,1.0f));

  // compute pixel area from film info. TODO: how to have access to this??
  /*
  float film_width = 35.0f; //in mm
  float aspect_ratio = 4.0f/3.0f;
  float film_height = film_width / aspect_ratio;
  pixel_area = (film_width/hRes) * (film_height/vRes);

  printf("Pixel area: %f mm²\n", pixel_area);
  */
}

RGB Integrator::bd_path(const Scene& scene,
                        const Ray& primary_ray,
                        const Isect& isect,
                        int path_length)
{
  // ------------------------------
  // ---- Auxiliary structures ----
  // ------------------------------
  typedef struct
  {
    Isect isect; // intersection info
    Vec3 pos;    // position of this vertex
    Ray last;    // the ray that upon intersection returned this vertex
    float pdf;   // PDF (solid angle) of being selected by the last vertex
    bool valid;  // is this vertex valid or not?
  } Vertex;
  // ------------------------------

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

  // ------------------ Camera path --------------------
  // we start by filling the second vertex cell (as the first EDGE is fixed due
  // to lens delta specular behavior).
  int cp_idx = 1;

  Vertex &first_v = vertices[cp_idx];
  first_v.valid = true;
  first_v.pos = primary_ray(isect.t)
                  + isect.normal*(isect.shape->type == GLASS ? -0.0001f : 0.0001f);
  first_v.last = primary_ray;
  first_v.isect = isect;
  first_v.pdf = 1.0f; // TODO: review this. should be pdf of the lens sample (?)

  // keep sampling BRDF and building camera subpath
  Isect cp_isect = isect;
  Ray last_ray = primary_ray;

  for(int i = 2; i < path_length; ++i)
  {
    // last intersection point
    Vec3 p = vertices[cp_idx].pos;

    // sample BRDF for outgoing direction
    float dir_pdf;
    Vec3 out_dir = cp_isect.shape->sample_brdf(p, -last_ray.d, dir_pdf);

    // cast ray in the sampled direction
    Ray to_next_point(p, out_dir); Isect next_isect;
    if( !scene.cast_ray(to_next_point, next_isect) ) break; //TODO: Ray escaped. Do something?

    // store vertex
    cp_idx += 1;

    Vertex &v = vertices[cp_idx];
    v.valid = true;
    v.pdf = dir_pdf;
    v.last = to_next_point;
    v.isect = next_isect;
    v.pos = to_next_point(next_isect.t)
              + next_isect.normal*(next_isect.shape->type == GLASS ? -0.0001f : 0.0001f);
  }
  // ----------------------- Light path --------------------------
  // randomly pick a lightsource and sample a point on its surface
  const int l_idx = rand() % scene.emissive_prims.size();
  const Shape& l = scene.prims[ scene.emissive_prims[l_idx] ];
  Vec3 l_p, l_n;   // sampled point q and normal at q
  float l_pdf;     // probability of choosing q on the surface of l
                   // TODO: this will be biased if we have more than one light source!
                   // must consider the case where there's more than one light source,
                   // where the probability of picking a point is 1.0f / A, where A
                   // is the total surface area of the light sources (and we importance
                   // sample the area, also).
  l.sample_surface(l_p, l_n, l_pdf);

  // unlike the lens sample, we'll need to store at least partial info on the
  // light sample because we need it explicitly for the radiance computation step.
  Vertex &last = vertices[vertices.size()-1];
  last.valid = true;
  last.pdf = l_pdf; // PDF in respect to surface area. must convert it after!
  last.pos = l_p;
  last.isect.normal = l_n;
  last.isect.shape = &l;

  // cast ray from this sample and find the first intersection
  Vec3 l_d = l_n; // TODO: actually sample hemisphere around l_p!!!
  Ray l_ray( l_p + 0.0001f*l_n, l_d ); Isect l_isect;

  //TODO: DO SOMETHING IF WE MISSED THIS RAY. Camera subpath has the same issue.
  if( !scene.cast_ray(l_ray, l_isect) ) return RGB(1.0f, 0.0f, 0.0f);

  // fill the last but one vertex
  int lp_idx = vertices.size()-2;

  Vertex &last_v = vertices[lp_idx];
  last_v.valid = true;
  last_v.pdf = 1.0f; //TODO: because we used the normal direction
  last_v.last = l_ray;
  last_v.isect = l_isect;
  last_v.pos = l_ray(l_isect.t)
                + l_isect.normal*(l_isect.shape->type == GLASS ? -0.0001f : 0.0001f);

  // build light subpath by sampling BRDFs
  for(int i = 2; i < path_length; ++i)
  {
    // get current vertex, its position and sample BRDF
    Vertex &cur = vertices[lp_idx];
    Vec3 p = cur.pos;

    float pdf_dir;
    Vec3 out_dir = cur.isect.shape->sample_brdf(p, -cur.last.d, pdf_dir);

    // cast ray in this direction to get next vertex
    Ray next_ray(p, out_dir); Isect next_isect;
    if( !scene.cast_ray(next_ray, next_isect) ) break; //TODO: ray escaped. do something?

    // convert pdf from solid angle to surface area units
    float pdf_surface = pdf_dir * glm::dot(next_isect.normal, -out_dir) / next_isect.d2;

    // store vertex info
    lp_idx -= 1;

    Vertex &next = vertices[lp_idx];
    next.valid = true;
    next.last = next_ray;
    next.isect = next_isect;
    next.pos = next_ray(next_isect.t)
                + next_isect.normal*(next_isect.shape->type == GLASS ? -0.0001f : 0.0001f);
    next.pdf = pdf_surface;
  }

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

  // 2. Link the first two vertices of the camera subpath with the first vertex
  // of the light subpath. This should give a 5 vertex path (when path_length = 3).

  // discard path if it is not complete!
  const int vertex_l = vertices.size()/2+1;
  Vertex &v_l = vertices[vertex_l];
  Vertex &v_c = vertices[1];

  if( !v_l.valid || !v_c.valid ) return RGB(1.0f, 0.0f, 0.0f);

  // try to connect paths
  Vec3 v2l_ = v_l.pos - v_c.pos;
  float d2 = glm::dot(v2l_, v2l_);
  Vec3 v2l = glm::normalize(v2l_);

  // this seems to work, but may cause problems with the light sources themselves
  Ray shadow_ray(v_c.pos, v2l); Isect shadow_isect;
  scene.cast_ray(shadow_ray, shadow_isect);
  float vis = glm::length(shadow_ray(shadow_isect.t)-v_l.pos) < 0.001f ? 1.0f : 0.0f;

  // vertices actually see each other. compute contribution.
  // throughput at path junction
  RGB brdf_v_c = v_c.isect.shape->brdf(-v_c.last.d, v2l, v_c.pos);
  float cosVC = glm::dot(v_c.isect.normal, v2l);
  if( v_c.isect.shape->type == GLASS ) cosVC *= -1.0f;

  RGB brdf_v_l = v_l.isect.shape->brdf(-v2l, -v_l.last.d, v_l.pos);
  float cosVL = glm::dot(v_l.isect.normal, -v_l.last.d);
  if( v_l.isect.shape->type == GLASS ) cosVL *= -1.0f;

  // pdf of path junction
  // TODO: this seems to be the problem!!
  float pdf_junction_angle = v_c.isect.shape->pdf_brdf(-v_c.last.d, v2l, v_c.pos);
  float pdf_junction = pdf_junction_angle * glm::dot(v_l.isect.normal, -v2l) / shadow_isect.d2;

  // throughput for the rest of the light path.
  //
  // we compute this throughput in the OPPOSITE sense then what is computed
  // in the light subpath construction!
  //
  // also, there's no need for shadow ray the last vertex as it is guaranteed
  // to be visible by construction.
  RGB tp_lp(1.0f); Ray from_cam = -shadow_ray;
  for(int i = vertex_l+1; i < vertices.size()-1; ++i)
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

  // run through the lightpath accumulating its pdf
  float pdf_lp = 1.0f;
  for(int i = vertices.size()-1; i >= vertex_l; --i)
    pdf_lp *= vertices[i].pdf;

  // throughput and pdf of the whole path. the final version of
  // this must also compute the throughput and pdf of the camera subpath!
  float pdf = pdf_lp;

  RGB tp = (brdf_v_c*cosVC) * (brdf_v_l*cosVL) * tp_lp;
  RGB e = vertices[vertices.size()-1].isect.shape->emission;

  return vis * e * tp * (1.0f / pdf);

  // ----------------------------------------------------

  // --------------------------------------------------------------
  // TODO: build full path by linking a prefix of the camera path with a suffix
  //       of the light path. cast shadow ray to check for occlusions.

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


  //return RGB(0.1f);
}

RGB Integrator::bdpt(const Scene& scene,
                      const Ray& primary_ray,
                      const Isect& isect)
{
  return bd_path(scene, primary_ray, isect, 3);
}

RGB Integrator::camera_path(const Scene& scene,
                            const Ray& primary_ray,
                            const Isect& isect, int path_length)
{
  // we want to know the radiance L exiting from p = primary_ray(isect.t) and
  // going to the origin of primary_ray, along the direction primary_ray.d.
  RGB path_rad(0.0f, 0.0f, 0.0f);

  // the first "component" is the possible emission from p
  // TODO: this is included in the Direct Illumination computation!
  //path_rad += isect.shape->emission;
  if(path_length == 1) return isect.shape->emission;

  // now we need to compute all the light arriving at p, from all directions; or,
  // speaking in terms of surface, we want to know the contribution of each point
  // q in the scene. For this, we need to know the radiance L(p <- q), i.e. the
  // radiance that goes from q to p, but also we need to account for the cosine
  // weighting of radiance, the BRDF attenuation and also express the probability
  // p(w) of picking a direction in terms of the probability p(q) of picking that
  // specific point.

  // TODO: pdf of primary_ray. For a pinhole camera, this is always 1 (because
  // for a single sample on the film, there's only one ray passing through the
  // aperture).
  float path_pdf = 1.0f;
  RGB throughput(1.0f);
  Vec3 p = primary_ray(isect.t) + isect.normal * (isect.shape->type == GLASS ? -0.0001f : 0.0001f);
  Isect isect_p = isect;
  Ray o_to_p = primary_ray;

  // this loop importance samples the BRDF
  // TODO: next event estimation, which will make this loop compute all the subpaths
  // up to path_length. only doing so we can russian roulette the path termination
  for(int i = 0; i < path_length-2; ++i)
  {
    // pick point by uniformly sampling hemisphere around p and tracing a ray
    // from p. Remember that we are computing a integral over the whole surface
    // of the scene, thus we could in principle simply pick random points on the
    // scene to compute the Monte Carlo estimative, but this would result in lots
    // of points which would be occluded and thus have zero contribution. By
    // tracing a ray and getting the closest point, we are "guaranteed" to get a
    // point that will contribute, and the probability in terms of solid angle
    // of picking it will be simply the probability of picking that direction. The
    // same probability in terms of surface area would be something related to
    // the total area visible from p with the usual cosine weighting.
    //
    // The last bounce will simply add the direct illumination from the last
    // surface to the last but one and skip.
    float pdf_angle;
    Vec3 w = isect_p.shape->sample_brdf(p, -o_to_p.d, pdf_angle);

    // If our ray doesn't hit any surface, we simply return a sample with
    // with contribution zero.
    // TODO: escaping rays must receive the contribution of environment lighting
    Ray p_to_q(p, w); Isect isect_q;
    if( !scene.cast_ray(p_to_q, isect_q) ) return RGB(0.0f);

    // update throughput with BRDF * cos and path PDF
    // OLD: throughput *= isect_p.shape->brdf(w, -o_to_p.d, p) * glm::dot(w, isect_p.normal);
    //
    // QUESTION:
    // This hotfix kinda solves the problem of refraction, although (1) I'm not sure
    // of the physical meaning of this cosine term for refraction - although PBR includes
    // it, (2) still there's a problem with the LIGHT SOURCE image seen on the glass
    // sphere (it is getting the diffuse color of the emissive primitive),
    // (3) the mirror ball is not correctly reflected, we do not see the reflex on
    // the refracted image and (4) there's no sign of caustic paths.
    //
    // Probably there's some issue with direct illumination for paths inside the
    // sphere.
    //
    // Observations:
    //  - Problem (2) happens only for paths of length 4. Groundtruth reveals that
    //    with 1 bounce only it should be possible to create the light source highlight,
    //    as with one bounce only we can reach the other side of the sphere and compute
    //    the direct illumination contribution.
    //float cosNW = isect.shape->type == GLASS ? glm::dot(-isect_p.normal, w) : glm::dot(isect_p.normal, w);
    float cosNW = glm::dot(isect_p.normal, w);
    if( isect_p.shape->type == GLASS ) cosNW *= -1.0f;

    throughput *= isect_p.shape->brdf(w, -o_to_p.d, p) * cosNW;
    path_pdf *= pdf_angle;

    // update variables to recursively compute next light bounce
    //p = p_to_q(isect_q.t) + 0.0001f*isect_q.normal;
    p = p_to_q(isect_q.t) + isect_q.normal * (isect_q.shape->type == GLASS ? -0.0001f : 0.0001f);
    isect_p = isect_q;
    o_to_p = p_to_q;
  }

  // the last iteration importance samples direct lighting by combining
  // BRDF sampling and light sampling using multiple importance sampling.
  // sample BRDF
  float brdf_pdf;
  RGB di_brdf = sample_brdf(o_to_p, isect_p, scene, brdf_pdf);
  RGB rad_brdf = di_brdf * throughput;
  float pdf_brdf = path_pdf * brdf_pdf;

  // if material has specular properties, it is in general useless to sample
  // the light sources, as this will return paths with pdf zero which will NaN
  // the output. if this is the case, simply return the BRDF sampling path.

  /* ---- TEST ----
  if( isect_p.shape->type == GLASS || isect_p.shape->type == DELTA )
    return rad_brdf * (1.0f / pdf_brdf);
     -------------- */

  // sample light sources. reaching this point means that the material is
  // glossy/diffuse (non-delta)
  float light_pdf;
  RGB di_ls = sample_light(o_to_p, isect_p, scene, light_pdf);
  RGB rad_ls = di_ls * throughput;
  float pdf_ls = path_pdf * light_pdf;

  // --- TEST ---
  return rad_ls * (1.0f / pdf_ls);
  // ------------

  // Power heuristic for multiple importance sampling
  float over_sum_pdfs = 1.0f / (pdf_ls*pdf_ls + pdf_brdf*pdf_brdf);
  float w_light = (pdf_ls*pdf_ls) * over_sum_pdfs;
  float w_brdf = (pdf_brdf*pdf_brdf) * over_sum_pdfs;

  // final contribution of this radiance path
  return rad_ls*w_light/pdf_ls + rad_brdf*w_brdf/pdf_brdf;
}

RGB Integrator::pathtracer(const Scene& scene,
                            const Ray& primary_ray,
                            const Isect& isect)
{
  // TODO: russian rouletting. the it is done below is
  // underestimating the total radiance, thus the image
  // is always darker than it should
  const int max_length = 5;

  /*
  RGB rad(0.0f);
  for(int i = 0; i <= max_length; ++i)
    rad += camera_path(scene, primary_ray, isect, i);
  return rad;
  */

  return camera_path(scene, primary_ray, isect, 3);
}

RGB Integrator::normal_shading(const Scene& scene,
                                const Ray& primary_ray,
                                const Isect& isect)
{
  return (isect.normal+1.0f)*0.5f;
}

RGB Integrator::radiance_measurement(const Scene& scene,
                                      const Ray& primary_ray,
                                      const Isect& isect)
{
  // this is a measurement of the radiance arriving at the origin
  // of the ray, thus we do not consider the cosine-weighting which
  // would be associated with the normal of sensor surface. This will
  // be accounted for in the render() function. Here we only account
  // for the cosine-weighting factor from the emiting surface.
  //
  // Despite irradiance changing with distance, radiance does not! (if we do not
  // have a participating medium). This may be explained by the fact that a single
  // beam of photons with the same orientation and origin (which is what radiance
  // describes) has no reason to "disappear" while travelling through vacuum.
  // However, if we imagine a non-colimated beam of photons irradiated from the
  // same point, but with slightly different orientations, measuring irradiance
  // at different distances will result in different values as the same bundle
  // rays arrives in a greater area and thus power is more spread.
  float w = glm::dot(isect.normal, primary_ray.d);
  return w * isect.shape->emission;
}

RGB Integrator::direct_illumination_surfacearea(const Scene& scene,
                                                const Ray& primary_ray,
                                                const Isect& isect)
{
  // Recall that estimating direct illumination amounts to approximate the
  // light transport equation by computing irradiance coming only from emissive
  // surfaces; this means that we can approximate it as a sum of integrals for
  // each light surface:
  //
  //        L(p,d) ~ emission(p,d) +  sum INT_J² L(p,-w).brdf(p,d,w).cos(w,N) dw
  //                                   j
  // Where each term in the summation is an integral over the directions on the
  // solid angle projected by the j-th light source on the hemisphere H² over p,
  // estimated as:
  //
  //            1       L(p,-W).brdf(p,d,W).cos(W,N)
  //           --- SUM ------------------------------
  //            N   W              p(W)
  //
  // As it is hard to "guess" the direction which will hit a given surface on the
  // space, it is easier to simply sample a point q on the lightsource and then
  // compute the probability p(w) given the probability p(q) of choosing q, where
  // w is the normalized vector q-p:
  //
  //            1       L(p,-W).brdf(p,d,W).cos(W,N).V(p,q)
  //           --- SUM -------------------------------------
  //            N   q              p(q).   r²
  //                                    ---------
  //                                     cos(q,W)
  //
  // TODO: why r²/cos instead of cos/r²??
  //
  // Notice that now we have to take care of occlusion, as the sampled point may
  // not be visible from p. In this case, we simply set its contribution to ZERO
  // carry on with calculations, which is expressed in the visibility term V(p,q)
  // (in practice, this is calculated by casting a shadow ray). Although this is
  // not the most optimized way of doing things, as we could (somehow) importance
  // sample the visibility component, it seems to be the simplest way that works.
  // PBRT sets its pdf to zero, which effectively discards the sample, but I'm not
  // sure right now on how to treat the pdf of the non-zero samples: if we simply
  // consider it to be 1/(area of the WHOLE light source), we're not giving the
  // sample the right weight, as in practice we are sampling only the unoccluded
  // region and thus we should weight by the reciproval of the unoccluded, projected
  // area. Doing so would OVERESTIMATE the contribution of partially occluded
  // lightsources, thus biasing the image to be brighter than it should.
  RGB irradiance(0.0f);
  Vec3 p(primary_ray(isect.t));
  Vec3 d = -primary_ray.d;

  // for each light source on the scene...
  for(int idx : scene.emissive_prims)
  {
    const Shape& l = scene.prims[idx];

    // ... compute a Monte Carlo estimative for its irradiance on the point p
    RGB irr_light(0.0f);
    const int n_samples = 1;
    for(int i = 0; i < n_samples; ++i)
    {
      // sample a point on the surface
      Vec3 q, n;       // sampled point q and normal at q
      float pdf_area;  // probability of choosing q on the surface of l
      l.sample_surface(q, n, pdf_area);

      // compute radiance coming from q in direction w = q - p
      // TODO: in this first version, emission is always isotropic and uniform,
      // so every direction at every point on the surface will have the same radiance
      Vec3 w = q-p;
      Vec3 w_n = glm::normalize(w);
      float r2 = glm::dot(w,w);
      RGB L_qw = l.emission;

      // cast shadow ray to compute visibility; we trace a ray r that simply
      // returns the t of the closest intersection; if r(t) != q, then q is
      // occluded.
      // TODO: an specific visibility method may be less complex than the usual
      // (and expensive) cast_ray() method.
      Ray shadow_ray(p + 0.0001f*isect.normal, w_n); Isect isect_shadow;
      scene.cast_ray(shadow_ray, isect_shadow);
      float v = glm::length(shadow_ray(isect_shadow.t)-q) < 0.0001f ? 1.0f : 0.0f;

      // accumulate f(x)/p(x)
      RGB rad = v * glm::dot(w_n,isect.normal) * isect.shape->brdf(w_n, -d, p) * L_qw;
      float pdf_angle = pdf_area * r2 / glm::dot(n, -w_n);
      irr_light += rad * (1.0f/pdf_angle);
    }

    // compute final estimative for direct light contribution of l
    irradiance += irr_light * (1.0f/n_samples);
  }

  // TODO: I'm taking the emission just to make sure I'm not double counting
  // emissions in pathtracer.
  return irradiance + isect.shape->emission;
}

RGB Integrator::direct_illumination_solidangle(const Scene& scene,
                                                const Ray& primary_ray,
                                                const Isect& isect)
{
  // direct illumination computes the radiance at the intersection point p and
  // direction d = -primary_ray.d. we consider now not only the emission term (as in
  // the radiance_measurement() shader), but also the contribution of the light
  // arriving at the whole hemisphere and exiting in direction d (Veach's thesis,
  // p. 85). This can be expressed as integral of the hemisphere of directions
  // over p:
  //
  //        L(p,d) = emission(p,d) +  INT_H² L(p,-w).brdf(p,d,w).cos(w,N) dw
  //
  // Which can be numerically estimated as:
  //
  //                                1         De(p,-W).brdf(p,d,W) cos(W,N)
  //     L(p,d) ~ emission(p,d) +  --- . SUM -----------------------------
  //                                N     W              p(W)
  //
  // Where De(p,-W) is the DIRECT EMISSION of a putative emissive surface that
  // emits light to P in direction -W.
  //
  // This is a Monte Carlo estimator with IMPORTANCE SAMPLING, thus it is more
  // general then the usual estimative for uniformly distributed samples, where
  // the 1/p(W) can be taken out of the summation as the measure of the
  // integration domain H².
  //
  // Evaluating this integral as is would gives us lots of Zero contributions,
  // as it is hard to "guess" the direction W that hits a light source. Thus,
  // we'll choose POINTS l in the light sources and then evaluate the probability
  // p(l) of choosing this point.
  //
  // Obviously incoming light has many different sources, but direct illumination
  // approximates this by considering ONLY light that comes from light sources.
  // Thus.
  const int n_samples = 1;
  const float over_n_samples = 1.0f / n_samples;

  //glm::mat3 local2world(isect.bitangent, isect.normal, isect.tangent);
  glm::mat3 local2world = isect.shape->get_local_coordinate_system(isect.normal);
  RGB out(0.0f, 0.0f, 0.0f);

  for(int i = 0; i < n_samples; ++i)
  {
    RGB Lsample(0.0f, 0.0f, 0.0f);

    // uniform sample hemisphere around Normal
    float u1 = (float)rand()/RAND_MAX;
    float u2 = (float)rand()/RAND_MAX;
    const float r = sqrt(1.0f - u1*u1);
    const float phi = _2pi * u2;
    Vec3 w = local2world*Vec3( r*cos(phi), u1, r*sin(phi) );

    // send ray in direction w and check if there's emission
    Vec3 o = primary_ray(isect.t)+0.001f*isect.normal;
    Ray rW(o, w);
    Isect isectW;

    if( scene.cast_ray(rW, isectW) )
    {
        // outgoing ray hit some surface; capture its emission
        // (remember that a surface's emission is expressed in
        // in W.m⁻².sr⁻¹!!!)
        Lsample = isectW.shape->emission
                  * glm::dot(isect.normal,w)
                  * isect.shape->brdf(isect.normal, w, o);
    }

    // accumulate, weighting it by 1/p(w). as we're uniformly sampling
    // the hemisphere, p(w) = 1/2pi, thus we multiply it by 2pi.
    out += Lsample * _2pi;
  }

  return out * over_n_samples + isect.shape->emission;
}

RGBA Integrator::camera_response_curve(const RGB& irradiance) const
{
  const float k = 1.0f;
  const Vec3 gamma(1.0f/2.2f);
  return RGBA( glm::min(RGB(1.0f), glm::pow(k*irradiance,gamma) ), 1.0f);
}

void Integrator::render(const Scene& scene)
{
  // loop over pixels in film, sampling them
  // on the center. trace primary ray and return
  // black/white if we had an intersection or not.
  //
  // Here we need to keep in mind that we are sampling
  // IRRADIANCE at film positions, and our goal is
  // to integrate it all to know the total POWER arriving
  // at a given pixel. The Monte Carlo strategy then is
  // to sample irradiance at different points in the
  // film pixel, compute its mean and multiply by the
  // pixel size, which is the Monte Carlo estimative for
  // the integral of irradiance over area, which is POWER.
  //
  // Notice that there's a subtlety here: to compute irradiance
  // at a given point, we need to integrate radiance over
  // the hemisphere of directions, weighting   //return (isect.normal + 1.0f)*0.5f;
  //return isect.normal;by the cosine (Lambert's
  // cosine law). A realistic camera model would do this by sampling
  // many different directions for the same point, and each one of these
  // directions would generate a ray that passes through the lens and goes
  // to the scene (or, speaking in Area terms - not Solid Angle, we would sample
  // the lens). Given that our camera is a pinhole, the radiance
  // function is a DELTA, being zero everywhere but at the direction
  // linking the sample position and the camera pinhole. Thus,
  // its integral (i.e. the irradiance at the point) is simply
  // the radiance arriving from this ray, weighted by the cosine.
  //
  // It is very important to correctly compute these quantities,
  // and one should never forget the DOMAIN MEASURE term of the
  // Monte Carlo integral, otherwise the quantities will be
  // meaningless. For a film where all pixels have the same size
  // this would not be a problem and the final image would "make sense",
  // but the quantities itself would be underestimated (our final image is
  // biased). This would pose a problem if we wanted to apply a
  // camera response curve to simulate a given camera behavior.
  //
  // Moreover, realistic films do not measure power, but actual
  // energy (which is just power integrated over time). Even more realistic
  // films would consider dependence on wavelength (?).
  //
  // Plenoptic camera simulation would require a film that records estimatives
  // for RADIANCE, not irradiance, by discretizing the hemisphere of directions.
  for(int j = 0; j < vRes; ++j)
    for(int i = 0; i < hRes; ++i)
    {
      RGB irradiance_sample(0.0f, 0.0f, 0.0f);

      // uniformly sample pixel (i,j) and map sample coordinates
      // to [0,1]² with j-axis flipping (film coordinate system)
      float e1 = (float)rand()/RAND_MAX, e2 = (float)rand()/RAND_MAX;
      Vec2 uv( (i+e1)/hRes, ((vRes-1)-j+e2)/vRes );

      // get primary ray and first intersection,
      // then invoke shader to compute the sample value
      Ray primary_ray = scene.cam.get_primary_ray(uv);

      Isect isect; RGB rad(0.0f, 0.0f, 0.0f);
      if( scene.cast_ray(primary_ray, isect) )
        //rad = pathtracer(scene, primary_ray, isect);
        rad = bdpt(scene, primary_ray, isect);

      // cosine-weight radiance measure coming from emission_measure().
      // as explained above, for a pinhole camera, this is our irradiance sample.
      // Normal of the film is the look_at direction (-z)
      irradiance_sample = rad * glm::dot(-scene.cam.z, primary_ray.d);

      // sample splatting
      int sample_add = j*hRes+i;
      samples[sample_add] += irradiance_sample;
      weights[sample_add] += 1.0f;
    }

  // reconstruct image (power measurements).
  // samples/weights is an estimative of the mean value of the
  // irradiance at the film. We multiply this by its size to
  // obtain the estimative for the irradiance integral over the
  // surface (power).
  //
  // UPDATE: working with power measurements is HARD, as the
  // pixel area is too small and thus things round off to zero.
  // We them employ a rough estimate (a simple average) of irradiance
  // over the pixel area to convert to pixel values.
  for(int i = 0; i < vRes*hRes; ++i)
  {
    RGB avg_irradiance = samples[i]/weights[i];
    frame[i] = camera_response_curve(avg_irradiance);
  }
}

// ------------------------------------
// ------------ Work later ------------
// ------------------------------------
void Integrator::light_path(const Scene& scene, int path_length)
{
  // pick light source
  // TODO: do not use rand()!
  // TODO: importance sample total emitted power
  const int idx = rand() % scene.emissive_prims.size();
  const Shape& l = scene.prims[ scene.emissive_prims[idx] ];

  // sample light source surface
  Vec3 light_sample, light_normal; float light_pdf;
  l.sample_surface(light_sample, light_normal, light_pdf);

  // TODO: ...bounce light path throughout the scene...

  // if path_length is 1, trace ray directly to the camera lens
  Vec3 lens_sample; Vec2 lens_sample_uv; float lens_pdf;
  scene.cam.sample_lens(lens_sample, lens_sample_uv, lens_pdf);

  // trace ray from the last bounce to the lens sample position.
  // remember that this implies (1) checking visibility, (2) computing
  // PDF of tracing this ray based on the BRDF of the last intersection
  // point and (3) an update of the path throughput. this is the final
  // radiance we want to return
  Vec3 dir = glm::normalize(lens_sample-light_sample);
  Ray shadow_ray(light_sample, dir);
  Isect vis_check;
  RGB rad = scene.cast_ray(shadow_ray, vis_check) ? RGB(0.0f) : RGB(1.0f);

  // deposit sample on camera film
  Vec2 sample_film = scene.cam.deposit_sample(lens_sample_uv, -dir);
  int i = (int)(hRes*sample_film.x), j = (int)(vRes*sample_film.y);

  int add = j*hRes + i;
  samples[add] += RGB(1.0f);
  weights[add] += 1.0f;
}

void Integrator::light_tracer(const Scene& scene)
{
  // once we have the radiance arriving at the given lens sample, at
  // the given direction, the camera should be able to "deposit" this
  // sample on the film somehow
  light_path(scene, 1);

  // TODO: do not do this at every iteration!!!!
  for(int i = 0; i < vRes*hRes; ++i)
  {
    RGB avg_irradiance = samples[i]/weights[i];
    frame[i] = camera_response_curve(avg_irradiance);
  }
}
