#include "../../include/core/threadpool.h"

Threadpool::Threadpool(int n, Integrator& integrator, const Scene& scene)
  : halt_request(true),
    halted_jobs(0),
    terminated(0),
    n_threads(n),
    integrator(integrator),
    scene(scene)
{ }

void Threadpool::spawn()
{
  halt_request = false;
  halted_jobs = 0;

  for(int i = 0; i < n_threads; ++i)
    workers.push_back( std::thread(&Threadpool::thread_callback, std::ref(*this)) );
}

void Threadpool::hold()
{
  // sinalize threads to stop working
  halt_request = true;

  // wait until all jobs have halted (or terminated)
  while( halted_jobs < workers.size() && terminated < workers.size())
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

void Threadpool::thread_callback()
{
  const int SPP = 50;

  // get first job
  Job job;
  jobs_mtx.lock();
    job = jobs.front();
    jobs.pop();
  jobs_mtx.unlock();

  for(;;)
  {
    // render this patch
    // TODO: first version will be slower because of rand()'s thread contention
    integrator.render_patch(scene, job.i, job.j, job.w, job.h);

    // check if halting was requested before getting next job
    if( halt_request )
    {
      // TODO: what if this thread pre-empted between
      // halted_jobs incrementation and cv() wait?
      std::unique_lock<std::mutex> lock(Threadpool::mtx);
      halted_jobs++;
      cv.wait(lock);
    }

    // reinsert updated job if patch has not enough samples.
    // if job queue is empty, release mutex and bail out to terminate thread;
    // otherwise, pick a new job.
    job.spp++;
    jobs_mtx.lock();
      if(job.spp < SPP) jobs.push(job);

      if( jobs.empty() )
      {
        jobs_mtx.unlock();
        break;
      }
      else
      {
        job = jobs.front();
        jobs.pop();
      }
    jobs_mtx.unlock();
  }

  std::cout<<"Thread "<<std::this_thread::get_id()<<" done"<<std::endl;
  terminated++;

  return;
}
