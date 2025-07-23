#include "genetic.hpp"
#include "EVRP.hpp"
#include <random>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <cmath>

// ------------------------
// Variáveis globais
// ------------------------
std::vector<Chromosome> population;
std::vector<double> fitness_vals;
int TARGET_VEHICLES = 5;
double VEHICLE_PENALTY = 10000.0;

// Melhor solução global encontrada
double global_best_fitness = std::numeric_limits<double>::max();
double global_best_real_distance = 0.0;
std::vector<std::vector<int>> global_best_routes;

static std::mt19937 rng(std::random_device{}());

// ------------------------
// Funções auxiliares
// ------------------------
static std::vector<int> decode_full_route(const Chromosome &C) {
    int n = NUM_OF_CUSTOMERS;
    std::vector<std::pair<int, int>> ord;
    ord.reserve(n);
    for (int i = 0; i < n; ++i) ord.emplace_back(C[i].first, i + 1);
    std::sort(ord.begin(), ord.end());

    std::vector<int> route;
    route.reserve(n + n + 2);
    route.push_back(DEPOT);
    for (int i = 0; i < n; ++i) {
        route.push_back(ord[i].second);
        if (C[i].second && i + 1 < n) {
            route.push_back(DEPOT);
        }
    }
    if (route.back() != DEPOT) route.push_back(DEPOT);
    return route;
}

static int count_vehicles(const std::vector<int> &full_route) {
    int zeros = 0;
    for (int x : full_route) if (x == DEPOT) zeros++;
    return std::max(1, zeros - 1);
}

// ------------------------
// API de configuração
// ------------------------
void set_target_vehicles(int x, double penalty) {
    TARGET_VEHICLES = x;
    VEHICLE_PENALTY = penalty;
}

// ------------------------
// Avaliação de rota
// ------------------------
static double evaluate_distance(const std::vector<int>& route) {
    double tot = 0;
    for (size_t i = 0; i + 1 < route.size(); ++i)
        tot += distances[route[i]][route[i + 1]];
    return tot;
}

static void two_opt(std::vector<int>& r) {
    bool improved;
    int n = r.size();
    do {
        improved = false;
        for (int i = 1; i < n - 2; ++i) {
            for (int k = i + 1; k < n - 1; ++k) {
                double delta = 
                    distances[r[i - 1]][r[k]] +
                    distances[r[i]][r[k + 1]] -
                    (distances[r[i - 1]][r[i]] + distances[r[k]][r[k + 1]]);
                if (delta < -1e-8) {
                    std::reverse(r.begin() + i, r.begin() + k + 1);
                    improved = true;
                }
            }
        }
    } while (improved);
}

static bool swap_improve(std::vector<int>& r) {
    int n = r.size();
    for (int i = 1; i < n - 2; ++i) {
        for (int j = i + 1; j < n - 1; ++j) {
            std::swap(r[i], r[j]);
            double newc = evaluate_distance(r);
            std::swap(r[i], r[j]);
            if (newc + 1e-8 < evaluate_distance(r)) {
                std::swap(r[i], r[j]);
                return true;
            }
        }
    }
    return false;
}

// Retorna: (fitness, distância_real, rotas)
static std::tuple<double, double, std::vector<std::vector<int>>>
decode_and_localsearch(const Chromosome &C) {
    auto full = decode_full_route(C);
    std::vector<std::vector<int>> routes;
    routes.emplace_back();
    for (size_t i = 1; i + 1 < full.size(); ++i) {
        if (full[i] == DEPOT) {
            routes.emplace_back();
        } else {
            routes.back().push_back(full[i]);
        }
    }

    double total_dist = 0.0;
    for (auto &r : routes) {
        std::vector<int> seg = {DEPOT};
        seg.insert(seg.end(), r.begin(), r.end());
        seg.push_back(DEPOT);

        two_opt(seg);
        while (swap_improve(seg));

        for (size_t i = 0; i + 1 < seg.size(); ++i)
            total_dist += distances[seg[i]][seg[i + 1]];

        r.assign(seg.begin() + 1, seg.end() - 1);
    }

    int num_vehicles = count_vehicles(full);
    double fitness = total_dist;
    
    // CORREÇÃO APLICADA AQUI: Removido o termo adicional '1.0 +'
    if (TARGET_VEHICLES > 0 && num_vehicles < TARGET_VEHICLES) {
        int deficit = TARGET_VEHICLES - num_vehicles;
        fitness += deficit * VEHICLE_PENALTY;
    }

    return {fitness, total_dist, routes};
}

// ------------------------
// Inicialização
// ------------------------
void initialize_genetic(int num_customers) {
    population.resize(POP_SIZE);
    fitness_vals.resize(POP_SIZE);
    std::uniform_int_distribution<> pd(0, num_customers * 10);
    std::bernoulli_distribution bd(0.3);

    for (int i = 0; i < POP_SIZE; ++i) {
        Chromosome C(num_customers);
        for (int j = 0; j < num_customers; ++j)
            C[j] = {pd(rng), bd(rng)};
        
        auto [f, real_dist, routes] = decode_and_localsearch(C);
        population[i] = std::move(C);
        fitness_vals[i] = f;
        
        // Atualiza melhor solução global
        if (f < global_best_fitness) {
            global_best_fitness = f;
            global_best_real_distance = real_dist;
            global_best_routes = routes;
        }
    }
}

// ------------------------
// Operadores genéticos
// ------------------------
Chromosome tournament_selection() {
    std::uniform_int_distribution<> id(0, POP_SIZE - 1);
    int best = id(rng);
    double bf = fitness_vals[best];
    for (int k = 1; k < 3; ++k) {
        int r = id(rng);
        if (fitness_vals[r] < bf) {
            bf = fitness_vals[r];
            best = r;
        }
    }
    return population[best];
}

std::pair<Chromosome, Chromosome>
crossover(const Chromosome &A, const Chromosome &B) {
    std::uniform_real_distribution<> rd(0, 1);
    if (rd(rng) > CROSSOVER_RATE) return {A, B};
    
    int n = A.size();
    int pt = std::uniform_int_distribution<>(1, n - 1)(rng);
    Chromosome c1 = A, c2 = B;
    for (int i = pt; i < n; ++i)
        std::swap(c1[i], c2[i]);
    return {c1, c2};
}

void mutate(Chromosome &C) {
    std::uniform_real_distribution<> rd(0, 1);
    std::uniform_int_distribution<> pd(0, C.size() * 10);
    for (auto &g : C) {
        if (rd(rng) < MUTATION_RATE) g.first = pd(rng);
        if (rd(rng) < MUTATION_RATE) g.second = !g.second;
    }
}

// ------------------------
// Loop principal
// ------------------------
// void run_genetic() {
//     for (int gen = 1; gen <= NUM_GENERATIONS; ++gen) {
//         std::vector<Chromosome> NP;
//         std::vector<double> NF;
//         NP.reserve(POP_SIZE);
//         NF.reserve(POP_SIZE);

//         // Preserva o melhor indivíduo
//         int bi = std::min_element(fitness_vals.begin(), fitness_vals.end()) - fitness_vals.begin();
//         NP.push_back(population[bi]);
//         NF.push_back(fitness_vals[bi]);

//         while ((int)NP.size() < POP_SIZE) {
//             auto p1 = tournament_selection();
//             auto p2 = tournament_selection();
//             auto [c1, c2] = crossover(p1, p2);
//             mutate(c1);
//             mutate(c2);

//             auto eval = [](Chromosome &C) -> std::tuple<double, double, std::vector<std::vector<int>>> {
//                 return decode_and_localsearch(C);
//             };

//             auto [f1, rd1, rt1] = eval(c1);
//             NP.push_back(c1);
//             NF.push_back(f1);

//             // Atualiza melhor global
//             if (f1 < global_best_fitness) {
//                 global_best_fitness = f1;
//                 global_best_real_distance = rd1;
//                 global_best_routes = rt1;
//             }

//             if ((int)NP.size() < POP_SIZE) {
//                 auto [f2, rd2, rt2] = eval(c2);
//                 NP.push_back(c2);
//                 NF.push_back(f2);
                
//                 if (f2 < global_best_fitness) {
//                     global_best_fitness = f2;
//                     global_best_real_distance = rd2;
//                     global_best_routes = rt2;
//                 }
//             }
//         }

//         population.swap(NP);
//         fitness_vals.swap(NF);

//         // Debug periódico
//         if (gen % LOCAL_EVERY == 0) {
//             std::cout << "Gen " << gen
//                       << " best fitness=" << global_best_fitness
//                       << "  (#rotas=" << global_best_routes.size() << ")"
//                       << "  dist=" << global_best_real_distance << "\n";

//             for (size_t r = 0; r < global_best_routes.size(); ++r) {
//                 std::cout << " R" << (r + 1) << ": 0";
//                 for (int cust : global_best_routes[r]) {
//                     std::cout << " -> " << cust;
//                 }
//                 std::cout << " -> 0\n";
//             }
//             std::cout << "--------------------------------\n";
//        }
//     }
// }



void run_genetic() {
    // Inicializa o best global
    global_best_fitness       = std::numeric_limits<double>::infinity();
    global_best_real_distance = std::numeric_limits<double>::infinity();
    global_best_routes.clear();

    for (int gen = 1; gen <= NUM_GENERATIONS; ++gen) {
        std::vector<Chromosome> NP;
        std::vector<double>     NF;
        NP.reserve(POP_SIZE);
        NF.reserve(POP_SIZE);

        // Elitismo: copia o melhor atual para a nova população
        int bi = std::min_element(fitness_vals.begin(), fitness_vals.end())
               - fitness_vals.begin();
        NP.push_back(population[bi]);
        NF.push_back(fitness_vals[bi]);

        // Gera o restante da nova população
        while ((int)NP.size() < POP_SIZE) {
            auto p1 = tournament_selection();
            auto p2 = tournament_selection();
            auto [c1, c2] = crossover(p1, p2);
            mutate(c1);
            mutate(c2);

            // Função auxiliar: avalia e atualiza melhor global
            auto eval_and_update = [&](Chromosome &C) {
                double fit, reald;
                std::vector<std::vector<int>> routes;
                std::tie(fit, reald, routes) = decode_and_localsearch(C);
                // atualiza melhor global conforme fitness (com penalidade)
                if (fit < global_best_fitness) {
                    global_best_fitness       = fit;
                    global_best_real_distance = reald;
                    global_best_routes        = routes;
                }
                return fit;
            };

            double f1 = eval_and_update(c1);
            NP.push_back(c1);
            NF.push_back(f1);

            if ((int)NP.size() < POP_SIZE) {
                double f2 = eval_and_update(c2);
                NP.push_back(c2);
                NF.push_back(f2);
            }
        }

        // Substitui população
        population.swap(NP);
        fitness_vals.swap(NF);
    }

    // Ao fim de todas as gerações, imprime o resultado real
    std::cout << "Best real distance = "<< global_best_real_distance<< "   (#vehicles = " << global_best_routes.size() << ")" << std::endl;
}



// ------------------------
// Impressão da solução
// ------------------------
void print_best_solution() {
    std::ofstream out("solution.txt");
    if (!out) {
        std::cerr << "Erro ao abrir solution.txt\n";
        return;
    }

    for (size_t i = 0; i < global_best_routes.size(); ++i) {
        out << "Vehicle " << i + 1 << ": 0";
        for (int c : global_best_routes[i]) {
            out << " -> " << c;
        }
        out << " -> 0\n";
    }
    out << "Total distance: " << global_best_real_distance << "\n";
    std::cout << "Total distance: " << global_best_real_distance << std::endl;
}