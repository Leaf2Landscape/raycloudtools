// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Glen Eaton
//
// This file provides a generic, thread-safe, bounded queue implementation
// using the C++ standard library synchronization primitives. It is a crucial
// component for the producer-consumer pattern, enabling backpressure to prevent
// excessive memory consumption when the producer is faster than the consumers.

#ifndef RAYLIB_RAYVOXEL_RAYLASTHREADSAFEQUEUE_H
#define RAYLIB_RAYVOXEL_RAYLASTHREADSAFEQUEUE_H

#include <condition_variable>
#include <mutex>
#include <queue>

namespace ray
{
  template <typename T>
  class ThreadSafeQueue
  {
  public:
    /// @brief Constructs a bounded queue with a fixed maximum size.
    /// @param max_size The maximum number of items the queue can hold.
    ThreadSafeQueue(size_t max_size) : max_size_(max_size), done_(false) {}

    /// @brief Pushes an item onto the queue. If the queue is full, this call
    ///        will block until a consumer makes space.
    /// @param item The item to push.
    void push(T&& item)
    {
      std::unique_lock<std::mutex> lock(mutex_);
      // Wait until the queue is not full.
      while (queue_.size() >= max_size_ && !done_) {
        cond_not_full_.wait(lock);
      }

      if (done_) {
        return;
      }

      queue_.push(std::move(item));
      cond_not_empty_.notify_one();
    }

    /// @brief Pops an item from the queue. If the queue is empty, this call
    ///        will block until a producer adds an item or signals completion.
    /// @param item The variable to populate with the popped item.
    /// @return True if an item was successfully popped, false if the queue is
    ///         empty and the producer is finished.
    bool pop(T& item)
    {
      std::unique_lock<std::mutex> lock(mutex_);
      // Wait until the queue is not empty or the producer is done.
      while (queue_.empty() && !done_) {
        cond_not_empty_.wait(lock);
      }

      if (queue_.empty() && done_) {
        return false;
      }

      item = std::move(queue_.front());
      queue_.pop();
      cond_not_full_.notify_one();
      return true;
    }

    /// @brief Signals to all waiting threads that production is complete.
    ///        This will unblock any waiting producers or consumers.
    void notify_done()
    {
      std::unique_lock<std::mutex> lock(mutex_);
      done_ = true;
      // Wake up all threads that might be waiting.
      cond_not_empty_.notify_all();
      cond_not_full_.notify_all();
    }

  private:
    std::queue<T> queue_;
    std::mutex mutex_;
    std::condition_variable cond_not_full_;
    std::condition_variable cond_not_empty_;
    size_t max_size_;
    bool done_;
  };

} // namespace ray

#endif // RAYLIB_RAYVOXEL_RAYLASTHREADSAFEQUEUE_H
