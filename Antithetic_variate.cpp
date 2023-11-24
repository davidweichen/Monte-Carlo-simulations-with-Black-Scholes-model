#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <chrono>

double mc_options_AV(double K, double S0, double r, double sigma, double T, 
int Nsteps, int Nrep) {
    std::vector<std::vector<double>> spot_prices1(Nrep / 2, std::vector<double>(1 + Nsteps));
    std::vector<std::vector<double>> spot_prices2(Nrep / 2, std::vector<double>(1 + Nsteps));

    double dt = T / Nsteps;
    double nudt = (r - 0.5 * sigma * sigma) * dt;
    double sidt = sigma * sqrt(dt);
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> dist(0, 1);
    
    for (int i = 0; i < Nrep / 2; i++) {
        spot_prices1[i][0] = S0;
        spot_prices2[i][0] = S0;
        for (int j = 0; j < Nsteps; j++) {
            double epsilon = dist(gen);
            spot_prices1[i][j + 1] = spot_prices1[i][j] * exp(nudt + sidt * epsilon);
            spot_prices2[i][j + 1] = spot_prices2[i][j] * exp(nudt - sidt * epsilon);
        }
    }
    
        double C1 = 0.0, C2 = 0.0;
        for (int i = 0; i < Nrep / 2; i++) {
            C1 += std::max(spot_prices1[i][Nsteps] - K, 0.0);
            C2 += std::max(spot_prices2[i][Nsteps] - K, 0.0);
        }
        double C = (C1 + C2) / (2.0 * (Nrep / 2));
        double C0 = exp(-r * T) * C;
        return C0;
 
}

int main() {
    double K = 110.0;
    double S0 = 100.0;
    double r = 0.05;
    double sigma = 0.2;
    double T = 1.0;
    int Nsteps = 52;
    std::vector<int> num_particles{10000, 100000, 1000000, 10000000};
    for(int n: num_particles){
    
    std::chrono::time_point<std::chrono::high_resolution_clock> start_time = std::chrono::high_resolution_clock::now();
    double optionPrice = mc_options_AV(K, S0, r, sigma, T, Nsteps, n);
    std::chrono::time_point<std::chrono::high_resolution_clock> end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    std::cout << " Elapsed time is : " << elapsed.count() << " " << std::endl;
    std::cout << "Option Price: " << optionPrice << std::endl;
    }
    return 0;
}
