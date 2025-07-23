#ifndef GENETIC_HPP
#define GENETIC_HPP

#include <vector>
#include <utility>

// --- Parâmetros do GA ---
static const int POP_SIZE         = 100;
static const int NUM_GENERATIONS  = 500;
static const double CROSSOVER_RATE = 0.8;
static const double MUTATION_RATE  = 0.02;
static const int LOCAL_EVERY       = 20;

// --- Controle de veículos alvo ---
extern int TARGET_VEHICLES;      // número desejado de veículos
extern double VEHICLE_PENALTY;   // penalidade por ter menos veículos

// --- Integração EVRP ---
extern double **distances;       // matriz de distâncias
extern int    NUM_OF_CUSTOMERS;  // número de clientes
extern int    DEPOT;             // índice do depósito

// --- Genoma & Estado ---
using Gene       = std::pair<int,bool>;    // <prioridade, retorna_ao_depot?>
using Chromosome = std::vector<Gene>;
extern std::vector<Chromosome>       population;
extern std::vector<double>           fitness_vals;
extern std::vector<std::vector<int>> best_routes;

// --- Interface ---
void set_target_vehicles(int x, double penalty);
void initialize_genetic(int num_customers);
void run_genetic();
void print_best_solution();

#endif // GENETIC_HPP
