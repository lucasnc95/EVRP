// genetic.cpp
#include "genetic.hpp"
#include "EVRP.hpp"
#include <iostream>
#include <random>
#include <cmath>
#include <limits>
#include <algorithm>

// Bits para codificar K e prioridades
static int BITS_VEH;
static int BITS_PER_GENE;

// População e fitness
auto population = std::vector<Chromosome>();
auto fitness_vals = std::vector<double>();
// Melhor divisão
static std::vector<std::vector<int>> best_routes;
// RNG
static std::mt19937 rng(std::random_device{}());

// Avalia distância de uma rota
static double evaluate_distance(const std::vector<int> &route) {
    double total = 0.0;
    for (size_t i = 0; i + 1 < route.size(); ++i)
        total += distances[route[i]][route[i+1]];
    return total;
}

// DP para particionar em mv rotas
static std::pair<double,std::vector<std::vector<int>>> splitRoutes(
    const std::vector<int> &o, int mv) {
    int n = o.size();
    std::vector<double> best(n+1, std::numeric_limits<double>::infinity());
    std::vector<int> prev(n+1, -1);
    best[0] = 0.0;
    for (int j = 1; j <= n; ++j) {
        for (int i = 0; i < j; ++i) {
            double cost = distances[DEPOT][o[i]];
            for (int k = i; k < j-1; ++k)
                cost += distances[o[k]][o[k+1]];
            cost += distances[o[j-1]][DEPOT];
            double tot = best[i] + cost;
            if (tot < best[j]) {
                best[j] = tot;
                prev[j] = i;
            }
        }
    }
    std::vector<int> cuts;
    for (int idx = n; idx > 0; idx = prev[idx])
        cuts.push_back(idx);
    std::reverse(cuts.begin(), cuts.end());
    std::vector<std::vector<int>> R;
    int start = 0;
    for (int cut : cuts) {
        R.emplace_back(o.begin() + start, o.begin() + cut);
        start = cut;
    }
    while ((int)R.size() < mv) {
        int wi = 0;
        double wcost = -1;
        for (int i = 0; i < (int)R.size(); ++i) {
            std::vector<int> full = {DEPOT};
            full.insert(full.end(), R[i].begin(), R[i].end());
            full.push_back(DEPOT);
            double c = evaluate_distance(full);
            if (c > wcost) { wcost = c; wi = i; }
        }
        auto r = R[wi];
        R.erase(R.begin() + wi);
        int mid = r.size() / 2;
        R.insert(R.begin() + wi, std::vector<int>(r.begin() + mid, r.end()));
        R.insert(R.begin() + wi, std::vector<int>(r.begin(), r.begin() + mid));
    }
    double total = 0.0;
    for (auto &r : R) {
        std::vector<int> full = {DEPOT};
        full.insert(full.end(), r.begin(), r.end());
        full.push_back(DEPOT);
        total += evaluate_distance(full);
    }
    return {total, R};
}

void initialize_genetic(int num_customers) {
    int range = NUM_OF_CUSTOMERS - MIN_VEHICLES + 1;
    BITS_VEH = std::ceil(std::log2(range));
    BITS_PER_GENE = std::ceil(std::log2(num_customers + 1));
    int chrom_len = BITS_VEH + BITS_PER_GENE * num_customers;
    population.assign(POP_SIZE, Chromosome(chrom_len));
    fitness_vals.assign(POP_SIZE, 0.0);
    std::uniform_int_distribution<> bit01(0, 1);
    for (int i = 0; i < POP_SIZE; ++i) {
        Chromosome &c = population[i];
        for (int j = 0; j < chrom_len; ++j)
            c[j] = bit01(rng);
        fitness_vals[i] = decode_and_evaluate(c);
    }
}

int binary_to_decimal(const Chromosome &chrom, int start, int len) {
    int v = 0;
    for (int i = 0; i < len; ++i)
        v = (v << 1) | chrom[start + i];
    return v;
}

double decode_and_evaluate(const Chromosome &chrom) {
    int nc = (chrom.size() - BITS_VEH) / BITS_PER_GENE;
    int range = NUM_OF_CUSTOMERS - MIN_VEHICLES + 1;
    int rawK = binary_to_decimal(chrom, 0, BITS_VEH);
    int K = (rawK % range) + MIN_VEHICLES;
    std::vector<std::pair<int, int>> order;
    for (int g = 0; g < nc; ++g) {
        int val = binary_to_decimal(chrom, BITS_VEH + g * BITS_PER_GENE, BITS_PER_GENE);
        order.emplace_back(val, g + 1);
    }
    std::sort(order.begin(), order.end());
    std::vector<int> seq;
    for (auto &p : order)
        seq.push_back(p.second);
    auto [total, routes] = splitRoutes(seq, K);
    best_routes = routes;
    return total;
}

Chromosome tournament_selection() {
    int tsize = 3;
    double best = std::numeric_limits<double>::infinity();
    Chromosome win;
    for (int i = 0; i < tsize; ++i) {
        int idx = std::uniform_int_distribution<>(0, POP_SIZE - 1)(rng);
        if (fitness_vals[idx] < best) {
            best = fitness_vals[idx];
            win = population[idx];
        }
    }
    return win;
}

std::pair<Chromosome, Chromosome> crossover(const Chromosome &a, const Chromosome &b) {
    if (std::uniform_real_distribution<>(0,1)(rng) > CROSSOVER_RATE) return {a, b};
    int n = a.size();
    int pt = std::uniform_int_distribution<>(1, n - 1)(rng);
    Chromosome c1 = a, c2 = b;
    for (int i = pt; i < n; ++i) {
        bool tmp = c1[i]; c1[i] = c2[i]; c2[i] = tmp;
    }
    return {c1, c2};
}

void mutate(Chromosome &c) {
    for (size_t i = 0; i < c.size(); ++i) {
        if (std::uniform_real_distribution<>(0,1)(rng) < MUTATION_RATE) c[i] = !c[i];
    }
}

void run_genetic() {
    for (int gen = 0; gen < NUM_GENERATIONS; ++gen) {
        std::vector<Chromosome> np; std::vector<double> nf;
        int bi = std::min_element(fitness_vals.begin(), fitness_vals.end()) - fitness_vals.begin();
        np.push_back(population[bi]); nf.push_back(fitness_vals[bi]);
        while ((int)np.size() < POP_SIZE) {
            auto p1 = tournament_selection(), p2 = tournament_selection();
            auto [c1, c2] = crossover(p1, p2);
            mutate(c1); mutate(c2);
            nf.push_back(decode_and_evaluate(c1)); np.push_back(c1);
            if ((int)np.size() < POP_SIZE) { nf.push_back(decode_and_evaluate(c2)); np.push_back(c2); }
        }
        population.swap(np); fitness_vals.swap(nf);
    }
    int bi = std::min_element(fitness_vals.begin(), fitness_vals.end()) - fitness_vals.begin();
    extern double current_best; current_best = fitness_vals[bi];
    decode_and_evaluate(population[bi]);
}

void print_best_solution() {
    std::cout << "Best solution:" << std::endl;
    for (size_t i = 0; i < best_routes.size(); ++i) {
        std::cout << "Vehicle " << i+1 << ": 0";
        for (int c : best_routes[i]) std::cout << " -> " << c;
        std::cout << " -> 0" << std::endl;
    }
}
