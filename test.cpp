#include <iostream>
#include <vector>
#include <random>
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
    #pragma omp parallel for
        for (int i = 0; i < num_particles; i++) {
            v[i] = a + b;
        }    
    
    chrono::time_point<std::chrono::high_resolution_clock> end_time = std::chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end_time - start_time;
    cout << " Elapsed time is : " << elapsed.count() << " " << std::endl;
        
    return 0;
}