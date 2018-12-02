#include "../../include/core/integrator.h"
#include "../../include/core/sampler.h"
#include <cstdio>
#include <chrono>
#include <thread>

using namespace std::chrono;

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
  //Vec3 p = primary_ray(isect.t) + 0.0001f*isect.normal;
  Vec3 p = primary_ray(isect.t) + isect.normal * (isect.shape->type == GLASS ? -0.0001f : 0.0001f);
  Isect isect_p = isect;
  Ray o_to_p = primary_ray;

  // this loop importance samples the BRDF
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
  if( isect_p.shape->type == GLASS || isect_p.shape->type == DELTA )
    return rad_brdf * (1.0f / pdf_brdf);

  // sample light sources. reaching this point means that the material is
  // glossy/diffuse (non-delta)
  float light_pdf;
  RGB di_ls = sample_light(o_to_p, isect_p, scene, light_pdf);
  RGB rad_ls = di_ls * throughput;
  float pdf_ls = path_pdf * light_pdf;

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

  RGB rad(0.0f);
  for(int i = 0; i <= max_length; ++i) rad += camera_path(scene, primary_ray, isect, i);
  //rad = camera_path(scene, primary_ray, isect, 3);

  return rad;
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
  const float _2pi = 2.0f * 3.141592654f;
  const float over2pi = 1.0f / _2pi;
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
    const float phi = 2.0f * 3.1416f * u2;
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
  /*
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
        rad = pathtracer(scene, primary_ray, isect);
        //rad = (1.0f+primary_ray.d)*0.5f;

      // cosine-weight radiance measure coming from emission_measure().
      // as explained above, for a pinhole camera, this is our irradiance sample.
      // Normal of the film is the look_at direction (-z)
      irradiance_sample = rad * glm::dot(-scene.cam.z, primary_ray.d);

      // sample splatting
      int sample_add = j*hRes+i;
      samples[sample_add] += irradiance_sample;
      weights[sample_add] += 1.0f;
    }
    */

  // ---- Multithreaded version ----
  const int blocks_x = 4, blocks_y = 2;
  const int patch_x = hRes/blocks_x, patch_y = vRes/blocks_y;
  std::vector<std::thread> render_jobs;

  for(int y = 0; y < blocks_y; ++y)
    for(int x = 0; x < blocks_x; ++x)
    {
      // select patch for thread (i,j) and spawn render job
      int job_x = x*patch_x, job_y = y*patch_y;

      render_jobs.push_back(std::thread([&](int start_x, int start_y) {
        for(int j = start_y; j < start_y+patch_y; ++j)
          for(int i = start_x; i < start_x+patch_x; ++i)
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
              rad = pathtracer(scene, primary_ray, isect);

            // cosine-weight radiance measure coming from emission_measure().
            // as explained above, for a pinhole camera, this is our irradiance sample.
            // Normal of the film is the look_at direction (-z)
            irradiance_sample = rad * glm::dot(-scene.cam.z, primary_ray.d);

            // sample splatting
            int sample_add = j*hRes+i;
            samples[sample_add] += irradiance_sample;
            weights[sample_add] += 1.0f;
          }
      }, job_x, job_y));
    }

  for(auto& rj : render_jobs) rj.join();

  // -------------------------------

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
