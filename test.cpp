#include <iostream>
#include <vector>

#include <cmath>
#include <omp.h>
#include <chrono>
using namespace std;
int main() {
    
    int num_particles = 10000000;
    int a =12;
    int b = 10;
    vector<int> v(num_particles);
    chrono::time_point<std::chrono::high_resolution_clock> start_time = std::chrono::high_resolution_clock::now();
    #pragma omp parallel
    {
      int thread_id = omp_get_thread_num();
      // Number of threads in the current team
      int nthreads = omp_get_num_threads();

    #pragma omp critical
      {
         std::cout << "Hello world, I'm thread " << thread_id << " out of " << nthreads << " total threads. " << std::endl; 
      }
    #pragma omp for
        for (int i = 0; i < num_particles; i++) {
            v[i] = a + b;
        }    
    }
    chrono::time_point<std::chrono::high_resolution_clock> end_time = std::chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end_time - start_time;
    cout << " Elapsed time is : " << elapsed.count() << " " << std::endl;
        
    return 0;
}