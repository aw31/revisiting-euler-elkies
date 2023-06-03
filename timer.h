#pragma once

#include <chrono>
#include <cstdio>

class Timer {
 private:
  std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
  std::chrono::time_point<std::chrono::high_resolution_clock> last_time;

 public:
  Timer() {
    start_time = std::chrono::high_resolution_clock::now();
    last_time = start_time;
  }

  void log_task(const std::string& task = "") {
    auto current_time = std::chrono::high_resolution_clock::now();

    auto total_time =
        std::chrono::duration<double>(current_time - start_time).count();
    auto task_time =
        std::chrono::duration<double>(current_time - last_time).count();

    printf("\n");
    printf("=== %s ===\n", task.c_str());
    printf("Time: %.3fs\n", task_time);
    printf("Total: %.3fs\n", total_time);

    last_time = current_time;
  }
};