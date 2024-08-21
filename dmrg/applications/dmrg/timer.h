#pragma once

#include <chrono>

inline std::chrono::time_point<std::chrono::high_resolution_clock> start_time;

void startTimer() {
    start_time = std::chrono::high_resolution_clock::now();
}

double stopTimer() {
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    return elapsed.count();
}
