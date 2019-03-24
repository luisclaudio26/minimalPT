#ifndef THREADPOOL_H
#define THREADPOOL_H

#include <thread>
#include <mutex>
#include <vector>
#include <queue>
#include <chrono>
#include <atomic>
#include <condition_variable>
#include "./integrator.h"
#include <iostream>

class Threadpool
{
public:
  typedef struct {
    int i, j;
    int w, h;
    int spp;
  } Job;

  Threadpool(int n);

  std::queue<Job> jobs;
  std::vector<std::thread> workers;
  void thread_callback();

  std::mutex mtx;
  std::condition_variable cv;
  std::atomic<int> halted_jobs;
  bool halt_request;

  void hold();
  void resume();
};

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

  //std::cout<<std::endl;
}

void Threadpool::resume()
{
  std::unique_lock<std::mutex> lock(mtx);

  halt_request = false;
  halted_jobs = 0;
  cv.notify_all();
}

Threadpool::Threadpool(int n)
  : halt_request(true), halted_jobs(0)
{
  for(int i = 0; i < n; ++i)
    workers.push_back( std::thread(&Threadpool::thread_callback, std::ref(*this)) );
}

void Threadpool::thread_callback()
{
  std::unique_lock<std::mutex> lock(Threadpool::mtx);

  for(int i = 0; i < 400; ++i)
  {
    if( halt_request )
    {
      // TODO: what if this thread pre-empted between
      // halted_jobs incrementation and cv() wait?
      halted_jobs++;
      cv.wait(lock);
    }

    std::cout<<std::this_thread::get_id()<<" says hi"<<std::endl;
    std::this_thread::sleep_for(std::chrono::seconds(1));
  }

  return;
}

#endif
