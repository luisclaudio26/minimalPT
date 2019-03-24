#include "../../include/core/threadpool.h"

void Threadpool::hold()
{
  // sinalize threads to stop working
  halt_request = true;

  // wait until all jobs have halted.
  // is this correct?
  while( halted_jobs < workers.size() )
  {
    // do nothing
  }
}

void Threadpool::resume()
{
  std::unique_lock<std::mutex> lock(mtx);

  halt_request = false;
  halted_jobs = 0;
  cv.notify_all();
}

Threadpool::Threadpool(int n, Integrator& integrator, const Scene& scene)
  : halt_request(true), halted_jobs(0), integrator(integrator), scene(scene)
{
  for(int i = 0; i < n; ++i)
    workers.push_back( std::thread(&Threadpool::thread_callback, std::ref(*this)) );
}

void Threadpool::thread_callback()
{
  std::unique_lock<std::mutex> lock(Threadpool::mtx);

  // get first job
  jobs_mtx.lock();
  Job job = jobs.front();
  jobs.pop();
  jobs_mtx.unlock();

  for(int i = 0; i < 300; ++i)
  {
    // render this patch
    // TODO: first version will be slower because of rand()'s thread contention
    integrator.render_patch(scene, job.i, job.j, job.w, job.h);

    // check if halting was requested before getting next job
    if( halt_request )
    {
      // TODO: what if this thread pre-empted between
      // halted_jobs incrementation and cv() wait?
      halted_jobs++;
      cv.wait(lock);
    }

    // reinsert updated job on the queue and loop
    jobs_mtx.lock();
      job.spp++;
      jobs.push(job);
      job = jobs.front();
      jobs.pop();
    jobs_mtx.unlock();
  }

  return;
}
