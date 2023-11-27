#include <iostream>
#include <algorithm>
#include <omp.h>
#define ARRAY_SIZE 100000000
#define ARRAY_VALUE 1231

using namespace std;
int main()
{
    
    int *arr = new int[ARRAY_SIZE];
    std::fill_n(arr, ARRAY_SIZE, ARRAY_VALUE);
    chrono::time_point<std::chrono::high_resolution_clock> start_time = std::chrono::high_resolution_clock::now();
    #pragma omp parallel for
    for(int i = 0; i < ARRAY_SIZE; i++)
    {
        arr[i] = arr[i] / arr[i] + arr[i] / 5 - 14;
    }
    chrono::time_point<std::chrono::high_resolution_clock> end_time = std::chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end_time - start_time;
    cout << " Elapsed time is : " << elapsed.count() << " " << std::endl;
    return 0;
}