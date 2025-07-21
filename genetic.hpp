// genetic.hpp
#ifndef GENETIC_HPP
#define GENETIC_HPP

#include <vector>
#include <utility>
#include "EVRP.hpp"

// GA parameters
static const int POP_SIZE = 100;
static const int NUM_GENERATIONS = 300;
static const double CROSSOVER_RATE = 0.7;
static const double MUTATION_RATE = 0.01;

using Chromosome = std::vector<int>;

extern std::vector<Chromosome> population;
extern std::vector<double> fitness_vals;
extern std::vector<std::vector<int>> best_routes;
extern double **distances;

void initialize_genetic(int num_customers);
double decode_and_evaluate(const Chromosome &chrom);
Chromosome tournament_selection();
std::pair<Chromosome,Chromosome> crossover(const Chromosome &a, const Chromosome &b);
void mutate(Chromosome &c);
void run_genetic();
void print_best_solution();
void print_solution_for_k(int rawK);

#endif // GENETIC_HPP