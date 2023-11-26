#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <omp.h>
#include <chrono>
using namespace std;


double generate_spot_prices(int num_particles, int num_weeks, double strike_price ,double spot_price, double risk_free_rate, double volatility) {
    // Create a vector to store the spot prices at each time step.
    #pragma omp parallel
    {
    vector<vector<double>> spot_prices(num_particles, vector<double>(num_weeks + 1));

    // Generate a random normal variable.
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<double> dist(0.0, 1.0);
    double dt = 1.0 / num_weeks;
    double C = 0.0;
    // Simulate the spot price at each time step in parallel.
#pragma omp for schedule(dynamic)
    for (int t = 0; t < num_particles; t++) {
        
        spot_prices[t][0] = spot_price;

        // Calculate the spot price at the current time step.
#pragma omp critical
        for (int i = 0; i < num_weeks; i++) {
            spot_prices[t][i + 1] = spot_prices[t][i] * exp((risk_free_rate - 0.5 * volatility * volatility) * dt + volatility * sqrt(dt) * dist(gen));
        }
    }

  #pragma omp for schedule(dynamic)
        for (int i = 0; i < num_particles; i++) {
            C += max(spot_prices[i][num_weeks] - strike_price, 0.0);
        }
        // Average the discounted payoffs to get the call option price.
        C /= num_particles * exp(-risk_free_rate * num_weeks * dt);
        
        return C;
    }
}


int main() {
    // Generate a sample of spot prices at each time step.
    vector<int> num_particles{10000, 100000, 1000000, 10000000};
    int num_weeks = 52;
    double spot_price = 100.0;
    double strike_price = 110.0;
    double risk_free_rate = 0.05;
    double volatility = 0.2;
    for(int n: num_particles){
        chrono::time_point<std::chrono::high_resolution_clock> start_time = std::chrono::high_resolution_clock::now();
        double call_option_price = generate_spot_prices(n, num_weeks, strike_price ,spot_price, risk_free_rate, volatility);
        chrono::time_point<std::chrono::high_resolution_clock> end_time = std::chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed = end_time - start_time;
        cout << " Elapsed time is : " << elapsed.count() << " " << std::endl;
        // Print the call option price.
        cout << "Call option price: " << call_option_price << endl;
    }

    return 0;
}