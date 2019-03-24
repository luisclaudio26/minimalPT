#ifndef THREADPOOL_H
#define THREADPOOL_H

#include <thread>
#include <mutex>
#include <vector>
#include <queue>
#include <chrono>
#include <atomic>
#include <condition_variable>
#include <iostream>

#include "./integrator.h"

class Threadpool
{
public:
  typedef struct {
    int i, j;
    int w, h;
    int spp;
  } Job;

  Threadpool(int n, Integrator& integrator, const Scene& scene);

  // resources
  int n_threads;
  Integrator& integrator; const Scene& scene;
  std::mutex jobs_mtx; std::queue<Job> jobs;
  std::vector<std::thread> workers;
  void thread_callback();

  // thread management
  std::mutex mtx;
  std::condition_variable cv;
  std::atomic<int> halted_jobs;
  std::atomic<int> terminated;
  bool halt_request;

  // external interface. hold() is expected to lock
  // until all threads have halted, so it is safe
  // to do anything on the color buffer after calling it
  void spawn();
  void hold();
  void resume();
};

#endif
