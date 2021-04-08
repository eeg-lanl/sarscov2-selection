#ifndef PARALLELISM_HPP_
#define PARALLELISM_HPP_

#include <thread>
#include <mutex>
#include <condition_variable>
#include <list>
#include <vector>
#include <functional>


class WorkerPool {
public:
  WorkerPool(size_t num_workers);
  ~WorkerPool();
  /// disable copy contructor
  WorkerPool(const WorkerPool & wp) = delete;
  WorkerPool& operator=(const WorkerPool & wp) = delete;
  typedef std::function<void(void)> Task;
  void push_back(Task task);
  void sync(); // silent sync
  void sync_timed(); // print estimated remaining time
  // TODO: timed sync
  template<class InputIterator>
  void execute_range(InputIterator first, InputIterator last);
private:
  std::vector<std::thread> workers;
  std::list<Task> tasks;
  bool abort_flag;
  std::mutex task_mut;
  std::condition_variable task_cv;
  int count;
  std::mutex count_mut;
  std::condition_variable count_cv;
  void worker_routine();
};

template<class InputIterator>
void WorkerPool::execute_range(InputIterator first, InputIterator last) {
  for ( auto it = first; it != last; ++it ) {
    push_back(*it);
  }
  sync();
}


#endif
