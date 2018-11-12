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
  // for the cosine-weighting from the emiting surface.
  const float PI4 = 12.566370614f;
  float dist = isect.d2; // TODO: CHECK THIS
  float w = glm::dot(isect.normal, primary_ray.d);

  return RGB(dist);
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
  // films would consider dependence on wavelength (?)
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
        //rad = radiance_measurement(scene, primary_ray, isect);
        rad = RGB(sqrt(isect.d2));

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
    //frame[i] = camera_response_curve(avg_irradiance);
    frame[i] = RGBA(avg_irradiance, 1.0f);
  }
}
