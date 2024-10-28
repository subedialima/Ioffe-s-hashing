#include <iostream>       
#include <vector>         
#include <random>         
#include <stdexcept>      
#include <cmath>          
#include <limits>         
#include <utility>        

using namespace std;

vector<pair<int, int>> IoffeHash(vector<double> *sketch, int hashLength, vector<unsigned int> *seeds) {
    // Check for invalid input
    if (!sketch || sketch->empty()) {
        throw std::invalid_argument("IoffeHash received an empty sketch vector.");
    }
    if (!seeds || seeds->empty()) {
        throw std::invalid_argument("IoffeHash received an empty seeds vector.");
    }

    vector<pair<int, int>> hash(hashLength);  // To store pairs of (k', t_k) as integers
    
    // Initialize random number generators using seeds
    mt19937 generator;
    gamma_distribution<double> gamma_dist(2.0, 1.0); // Gamma(2, 1) distribution
    uniform_real_distribution<double> uniform_dist(0.0, 1.0); // Uniform(0, 1) distribution

    // Loop to generate hash of the specified length
    for (int i = 0; i < hashLength; i++) {
        // Set seed for reproducibility
        generator.seed(seeds->at(i % seeds->size()));

        // Variables to track minimum a_k
        double min_a = std::numeric_limits<double>::max();
        int min_index = -1;
        double min_tk = 0.0;

        for (size_t k = 0; k < sketch->size(); ++k) {
            if (sketch->at(k) > 0) {
                // Sample random variables
                double r_k = gamma_dist(generator);
                double c_k = gamma_dist(generator);
                double beta_k = uniform_dist(generator);

                // Calculate t_k and apply the floor operation
                double t_k = floor((log(sketch->at(k)) / r_k) + beta_k);
                
                // Calculate y_k
                double y_k = exp(r_k * (t_k - beta_k));

                // Calculate a_k
                double a_k = c_k / (y_k * exp(r_k));

                // Update if a_k is the smallest found so far
                if (a_k < min_a) {
                    min_a = a_k;
                    min_index = k;
                    min_tk = t_k;  // Store t_k
                }
            }
        }

        // Store the index (k') and the floored integer value of t_k
        hash[i] = make_pair(min_index, static_cast<int>(min_tk));
    }

    return hash;
}

int main() {
    // Define multiple sketch vectors
    vector<vector<double>> sketches = {
        {1, 3, 4, 5, 6, 7, 8, 9, 10, 4, 2, 4, 6, 8, 10, 12, 14, 16},
        {1, 3, 4, 5, 6, 7, 8, 9, 10, 4},
        {2, 4, 6, 8, 10, 12, 14, 16}
    };

    // Define the hash length
    int hashLength = 10;

    // Initialize seeds for random generator
    vector<unsigned int> seeds = {1234, 5678, 91011, 1213, 1415, 100, 200, 300, 700, 600};

    try {
        // Loop through each sketch vector
        for (size_t i = 0; i < sketches.size(); ++i) {
            cout << "Hash Results for Sketch " << i + 1 << ":\n";

            // Call IoffeHash function for the current sketch vector
            vector<pair<int, int>> hashResult = IoffeHash(&sketches[i], hashLength, &seeds);

            // Output the result for the current sketch
            for (const auto &pair : hashResult) {
                cout << "(Index: " << pair.first << ", t_k: " << pair.second << ")\n";
            }
            cout << endl;
        }
    } catch (const std::exception &e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }

    return 0;
}
