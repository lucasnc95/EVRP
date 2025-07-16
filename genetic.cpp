// genetic.cpp
#include "genetic.hpp"
#include "EVRP.hpp"
#include <iostream>
#include <random>
#include <cmath>
#include <limits>
#include <algorithm>

int BITS_PER_GENE;
std::vector<std::vector<bool>> population;
std::vector<double> fitness_vals;
static std::vector<std::vector<int>> best_routes;
static std::mt19937 rng(std::random_device{}());

// Calcula a distância de uma rota completa (0->...->0) sem efeitos colaterais
static double evaluate_distance(const std::vector<int>& route) {
    double total = 0.0;
    for (size_t i = 0; i + 1 < route.size(); ++i) {
        total += distances[route[i]][route[i+1]];
    }
    return total;
}

void initialize_genetic(int num_customers) {
    BITS_PER_GENE = std::ceil(std::log2(num_customers + 1));
    population.resize(POP_SIZE);
    fitness_vals.resize(POP_SIZE);

    for (int i = 0; i < POP_SIZE; ++i) {
        std::vector<bool> chrom(BITS_PER_GENE * num_customers);
        for (size_t j = 0; j < chrom.size(); ++j) {
            chrom[j] = std::uniform_int_distribution<>(0, 1)(rng);
        }
        population[i] = chrom;
        fitness_vals[i] = decode_and_evaluate(chrom);
    }
}

int binary_to_decimal(const std::vector<bool>& chrom, int start) {
    int val = 0;
    for (int i = 0; i < BITS_PER_GENE; ++i) {
        val = (val << 1) | chrom[start + i];
    }
    return val;
}

double decode_and_evaluate(const std::vector<bool>& chrom) {
    int num_genes = chrom.size() / BITS_PER_GENE;
    std::vector<std::pair<int, int>> order;
    for (int g = 0; g < num_genes; ++g) {
        int val = binary_to_decimal(chrom, g * BITS_PER_GENE);
        order.emplace_back(val, g + 1); // cliente IDs de 1 a N
    }
    std::sort(order.begin(), order.end());

    double best_total = std::numeric_limits<double>::infinity();
    std::vector<std::vector<int>> best_split;

    for (int k = MIN_VEHICLES; k <= NUM_OF_CUSTOMERS; ++k) {
        std::vector<std::vector<int>> routes(k);
        for (size_t i = 0; i < order.size(); ++i) {
            routes[i % k].push_back(order[i].second);
        }
        double total = 0.0;
        for (auto& r : routes) {
            std::vector<int> full_route;
            full_route.reserve(r.size() + 2);
            full_route.push_back(0);
            full_route.insert(full_route.end(), r.begin(), r.end());
            full_route.push_back(0);
            total += evaluate_distance(full_route);
        }
        if (total < best_total) {
            best_total = total;
            best_split = routes;
        }
    }

    best_routes = best_split;
    return best_total;
}

std::vector<bool> tournament_selection() {
    int k = 3;
    double best_fit = std::numeric_limits<double>::infinity();
    int best_index = 0;
    for (int i = 0; i < k; ++i) {
        int idx = std::uniform_int_distribution<>(0, POP_SIZE - 1)(rng);
        if (fitness_vals[idx] < best_fit) {
            best_fit = fitness_vals[idx];
            best_index = idx;
        }
    }
    return population[best_index];
}

std::pair<std::vector<bool>, std::vector<bool>> crossover(const std::vector<bool>& a, const std::vector<bool>& b) {
    std::uniform_real_distribution<> dist(0.0, 1.0);
    if (dist(rng) > CROSSOVER_RATE) return {a, b};

    int n = a.size();
    int pt = std::uniform_int_distribution<>(1, n - 1)(rng);
    std::vector<bool> c1 = a, c2 = b;
    for (int i = pt; i < n; ++i) {
        bool tmp = c1[i];
        c1[i] = c2[i];
        c2[i] = tmp;
    }
    return {c1, c2};
}

void mutate(std::vector<bool>& c) {
    std::uniform_real_distribution<> dist(0.0, 1.0);
    for (size_t i = 0; i < c.size(); ++i) {
        if (dist(rng) < MUTATION_RATE) c[i] = !c[i];
    }
}

void run_genetic() {
    for (int gen = 0; gen < NUM_GENERATIONS; ++gen) {
        std::vector<std::vector<bool>> next_pop;
        std::vector<double> next_fit;
        int best_idx = std::min_element(fitness_vals.begin(), fitness_vals.end()) - fitness_vals.begin();
        next_pop.push_back(population[best_idx]);
        next_fit.push_back(fitness_vals[best_idx]);

        while ((int)next_pop.size() < POP_SIZE) {
            auto p1 = tournament_selection();
            auto p2 = tournament_selection();
            auto [c1, c2] = crossover(p1, p2);
            mutate(c1); mutate(c2);
            next_pop.push_back(c1);
            next_fit.push_back(decode_and_evaluate(c1));
            if ((int)next_pop.size() < POP_SIZE) {
                next_pop.push_back(c2);
                next_fit.push_back(decode_and_evaluate(c2));
            }
        }

        population = std::move(next_pop);
        fitness_vals = std::move(next_fit);
    }

    int best_idx = std::min_element(fitness_vals.begin(), fitness_vals.end()) - fitness_vals.begin();
    extern double current_best;
    current_best = fitness_vals[best_idx];
    decode_and_evaluate(population[best_idx]); // atualiza best_routes
}

void print_best_solution() {
    std::cout << "Best solution:" << std::endl;
    for (size_t i = 0; i < best_routes.size(); ++i) {
        std::cout << "Vehicle " << i + 1 << ": 0";
        for (int c : best_routes[i]) std::cout << " -> " << c;
        std::cout << " -> 0" << std::endl;
    }
}
