#include "parallelism.hpp"

#include "eta.hpp"

WorkerPool::WorkerPool(size_t num_workers) {
  workers.resize(num_workers);
  abort_flag = false;
  count = 0;
  for ( auto & worker : workers ) {
    worker = std::thread([this](){worker_routine();});
  }
}

WorkerPool::~WorkerPool() {
  std::unique_lock<std::mutex> task_lock(task_mut);
  abort_flag = true;
  task_cv.notify_all();
  // release the mutex
  task_lock.unlock();
  // wait for all workers to exit
  for ( auto & worker : workers ) {
    worker.join();
  }
}

void WorkerPool::push_back(Task task) {
  // if there are ZERO workers, just execute the task and return.
  if ( workers.empty() ) {
    task(); return;
  } // else...
  // first increase the task count
  std::unique_lock<std::mutex> count_lock(count_mut);
  ++count;
  count_lock.unlock();
  // then add task to list
  std::unique_lock<std::mutex> task_lock(task_mut);
  tasks.push_back(task);
  // and wake up a potential sleeping worker
  task_cv.notify_one();
}

void WorkerPool::sync() {
  std::unique_lock<std::mutex> count_lock(count_mut);
  count_cv.wait(count_lock, [this](){return count == 0;});
}

void WorkerPool::sync_timed() {
  std::unique_lock<std::mutex> count_lock(count_mut);
  EtaEstimator eta(count);
  while ( count > 0 ) {
    count_cv.wait(count_lock);
    // keep track of the duration of the next count change
    eta.update();
		std::cout << "\033[2K\rETA: " << eta << std::flush;
  }
  std::cout << std::endl;
}

void WorkerPool::worker_routine() {
  while ( true ) {
    // try to lock the task mutex
    std::unique_lock<std::mutex> task_lock(task_mut);
    // task_mut is locked: try to get a new task
    std::list<Task> my_tasks;
    if ( !tasks.empty() ) {
      // move the first element of tasks to my_tasks
      my_tasks.splice(my_tasks.begin(), tasks, tasks.begin());
    } else if ( abort_flag ) {
      break; // breaks the while ( true ) loop
    } else { // tasks is emty, but abort_flag is false: wait for more tasks
      task_cv.wait(task_lock);
    }
    // first release the mutex
    task_lock.unlock();
    // then execute acquired tasks in my_tasks
    for ( auto & task : my_tasks ) {
      task();
      // signal that a task has been completed
      std::unique_lock<std::mutex> count_lock(count_mut);
      --count;
      count_cv.notify_all(); // sync() might be waiting
    }
  } // while ( true ) loop
}
