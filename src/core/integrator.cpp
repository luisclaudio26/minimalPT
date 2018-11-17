#include "../../include/core/integrator.h"
#include <cstdio>
#include <chrono>
using namespace std::chrono;

/*
  high_resolution_clock::time_point tS = high_resolution_clock::now();
  high_resolution_clock::time_point tE = high_resolution_clock::now();
  duration<double> time_span = duration_cast<duration<double>>(tE - tS);
*/

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

RGB Integrator::direct_illumination(const Scene& scene,
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

  glm::mat3 local2world = glm::mat3(isect.bitangent, isect.normal, isect.tangent);
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
    Ray rW( primary_ray(isect.t)+0.001f*isect.normal, w );
    Isect isectW;

    if( scene.cast_ray(rW, isectW) )
    {
        // outgoing ray hit some surface; capture its emission
        // (remember that a surface's emission is expressed in
        // in W.m⁻².sr⁻¹!!!)
        Lsample = isectW.shape->emission
                  * glm::dot(isect.normal,w) * isect.shape->brdf(isect.normal,w);
    }


    // accumulate, weighting it by 1/p(w). as we're uniformly sampling
    // the hemisphere, p(w) = 1/2pi, thus we multiply it by 2pi.
    out += Lsample * _2pi;
  }

  return out * over_n_samples + isect.shape->emission;
}

RGBA Integrator::camera_response_curve(const RGB& irradiance) const
{
  float k = 1.0f;
  return RGBA( glm::min(RGB(1.0f),k*irradiance), 1.0f);
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
  // the hemisphere of directions, weighting by the cosine (Lambert's
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
        rad = direct_illumination(scene, primary_ray, isect);

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
