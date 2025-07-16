// genetic.hpp
#ifndef GENETIC_HPP
#define GENETIC_HPP

#include <vector>
#include "EVRP.hpp"

// Variáveis do EVRP definidas em EVRP.cpp

extern double **distances;
extern int MIN_VEHICLES;
extern int NUM_OF_CUSTOMERS;

// Parâmetros do GA
static const int POP_SIZE = 200;
static const int NUM_GENERATIONS = 500;
static const double CROSSOVER_RATE = 0.8;
static const double MUTATION_RATE = 0.005;

typedef std::vector<bool> Chromosome;

void initialize_genetic(int num_customers);
void run_genetic();
double decode_and_evaluate(const Chromosome &chrom);
int binary_to_decimal(const Chromosome &chrom, int start, int len);
void print_best_solution();

#endif // GENETIC_HPP