// genetic.hpp
#ifndef GENETIC_HPP
#define GENETIC_HPP

#include <vector>
#include "EVRP.hpp"

// Declarar variáveis globais de EVRP
extern double **distances;
extern int MIN_VEHICLES;
extern int NUM_OF_CUSTOMERS;


// GA parameters
static const int POP_SIZE = 500;
static const int NUM_GENERATIONS = 800;
static const double CROSSOVER_RATE = 0.8;
static const double MUTATION_RATE = 0.01;

extern int BITS_PER_GENE;
extern std::vector<std::vector<bool>> population;
extern std::vector<double> fitness_vals;

void initialize_genetic(int num_customers);
void run_genetic();
double decode_and_evaluate(const std::vector<bool>& chrom);
int binary_to_decimal(const std::vector<bool>& chrom, int start);
void print_best_solution();

#endif // GENETIC_HPP