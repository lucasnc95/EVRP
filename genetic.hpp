#ifndef GENETIC_HPP
#define GENETIC_HPP

#include <vector>
#include <cstdlib>
#include <algorithm>
#include "EVRP.hpp"

// GA parameters
static const int POP_SIZE = 50;
static const int NUM_GENERATIONS = 100;
static const double CROSSOVER_RATE = 0.8;
static const double MUTATION_RATE = 0.02;

extern int BITS_PER_GENE;
extern std::vector<std::vector<bool>> population;
extern std::vector<double> fitness_vals;

void initialize_genetic(int num_customers);
void run_genetic();
int binary_to_decimal(const std::vector<bool>& chrom, int start);
double decode_and_evaluate(const std::vector<bool>& chrom);
void print_best_solution();

#endif // GENETIC_HPP