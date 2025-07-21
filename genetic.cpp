// genetic.cpp
#include "genetic.hpp"
#include "EVRP.hpp"
#include <random>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>

static std::mt19937 rng(std::random_device{}());
std::vector<Chromosome> population;
std::vector<double> fitness_vals;
std::vector<std::vector<int>> best_routes;

// Evaluate distance of full route
static double evaluate_distance(const std::vector<int> &route) {
    double total = 0.0;
    for (size_t i = 0; i + 1 < route.size(); ++i)
        total += distances[route[i]][route[i+1]];
    return total;
}

// Original DP-based splitRoutes
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

// 2-opt local improvement
static bool two_opt_improve(std::vector<int> &route) {
    bool improved = false;
    int n = route.size();
    do {
        improved = false;
        for (int i = 1; i < n-2; ++i) {
            for (int k = i+1; k < n-1; ++k) {
                double delta =
                    distances[route[i-1]][route[k]] +
                    distances[route[i]][route[k+1]] -
                    (distances[route[i-1]][route[i]] +
                     distances[route[k]][route[k+1]]);
                if (delta < -1e-6) {
                    std::reverse(route.begin()+i, route.begin()+k+1);
                    improved = true;
                }
            }
        }
    } while (improved);
    return true;
}

// combine split and local search
static std::pair<double,std::vector<std::vector<int>>> split_and_refine(
    const std::vector<int> &seq, int mv) {
    auto [cost, routes] = splitRoutes(seq, mv);
    bool improved;
    do {
        improved = false;
        for (int r = 0; r+1 < (int)routes.size(); ++r) {
            if (routes[r].empty()) continue;
            auto A = routes;
            int c = A[r].back();
            A[r].pop_back();
            A[r+1].insert(A[r+1].begin(), c);
            double sum = 0;
            for (auto &rt : A) {
                std::vector<int> full={DEPOT};
                full.insert(full.end(), rt.begin(), rt.end());
                full.push_back(DEPOT);
                sum += evaluate_distance(full);
            }
            if (sum < cost-1e-6) {
                routes = std::move(A);
                cost = sum;
                improved = true;
            }
        }
    } while (improved);
    // 2-opt each
    for (auto &r : routes) {
        std::vector<int> full={DEPOT};
        full.insert(full.end(), r.begin(), r.end());
        full.push_back(DEPOT);
        two_opt_improve(full);
        r.assign(full.begin()+1, full.end()-1);
    }
    return {cost, routes};
}

void initialize_genetic(int num_customers) {
    int n = num_customers;
    population.resize(POP_SIZE);
    fitness_vals.resize(POP_SIZE);

    std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<> dp(0, n * n);

    std::vector<int> perm(n);
    for (int i = 0; i < n; ++i) perm[i] = i + 1;

    for (int i = 0; i < POP_SIZE; ++i) {
        Chromosome c(1 + n);

        // 1) cria permutação aleatória de clientes
        std::shuffle(perm.begin(), perm.end(), rng);

        // 2) aplica split e local search
        auto [cost, routes] = split_and_refine(perm, /*K não importa aqui*/ 
            1 + // o split_and_refine só usa perm, calcula K a partir de rawK mas ignoremos K inicial
            0 /* esse parâmetro é ignorado dentro de split_and_refine */
        );

        // 3) monta máscara rawK a partir das divisões reais
        std::vector<int> mask(n, 0);
        int idx = 0;
        for (size_t r = 0; r < routes.size(); ++r) {
            for (int cust : routes[r]) {
                // encontra posição de cust em perm
                while (perm[idx] != cust) ++idx;
                // se não for o último cliente da última rota, e se for fronteira, marca bit
                if (idx+1 < n && r+1 < (int)routes.size() && perm[idx+1] == routes[r+1][0]) {
                    mask[idx] = 1;
                }
                ++idx;
            }
        }
        // converte máscara em rawK
        int rawK = 0;
        for (int b = 0; b < n; ++b) {
            rawK = (rawK << 1) | mask[b];
        }
        c[0] = rawK;

        // 4) prioridades: simplesmente a posição inversa do perm (inverso da permutação)
        //    de modo que ordenar prioridades reproduza exatamente `perm`
        for (int j = 0; j < n; ++j) {
            // prioridade do cliente perm[j] := j
            c[1 + perm[j] - 1] = j;
        }

        population[i] = c;
        fitness_vals[i] = decode_and_evaluate(c);
    }
}

double decode_and_evaluate(const Chromosome &chrom) {
    int n = NUM_OF_CUSTOMERS;
    int rawK=chrom[0];
    std::vector<std::pair<int,int>> order;
    order.reserve(n);
    for(int j=0;j<n;++j) order.emplace_back(chrom[j+1], j+1);
    std::sort(order.begin(),order.end());
    std::vector<int> seq; seq.reserve(n);
    for(auto &p:order) seq.push_back(p.second);
    int K=1; for(int i=0;i<n;++i) if((rawK>>i)&1) ++K;
    auto [cost, routes] = split_and_refine(seq, K);
    best_routes = routes;
    return cost;
}

Chromosome tournament_selection() {
    int t=3, besti=0;
    double bf=std::numeric_limits<double>::infinity();
    for(int i=0;i<t;++i) {
        int idx=std::uniform_int_distribution<>(0,POP_SIZE-1)(rng);
        if(fitness_vals[idx]<bf){bf=fitness_vals[idx]; besti=idx;}
    }
    return population[besti];
}

std::pair<Chromosome,Chromosome> crossover(const Chromosome &a,const Chromosome &b) {
    if(std::uniform_real_distribution<>(0,1)(rng)>CROSSOVER_RATE) return {a,b};
    int n=a.size(),pt=std::uniform_int_distribution<>(1,n-1)(rng);
    Chromosome c1=a,c2=b;
    for(int i=pt;i<n;++i) std::swap(c1[i],c2[i]);
    return {c1,c2};
}

void mutate(Chromosome &c) {
    std::uniform_real_distribution<>r(0,1);
    std::uniform_int_distribution<>dk(0,(1<<NUM_OF_CUSTOMERS)-1);
    std::uniform_int_distribution<>dp(0,NUM_OF_CUSTOMERS*NUM_OF_CUSTOMERS);
    if(r(rng)<MUTATION_RATE) c[0]=dk(rng);
    for(int i=1;i<=NUM_OF_CUSTOMERS;++i) if(r(rng)<MUTATION_RATE) c[i]=dp(rng);
}

void run_genetic() {
    const int LOCAL_EVERY = 20; // a cada 20 gerações re-aplica local search na população
    for (int gen = 0; gen < NUM_GENERATIONS; ++gen) {
        // 1) avaliação + GA normal (elitismo, selec, crossover, mutação)
        //    ... (igual antes) ...

        // 2) a cada LOCAL_EVERY: refine toda população
        if (gen % LOCAL_EVERY == 0) {
            for (int i = 0; i < POP_SIZE; ++i) {
                // reaplica split_and_refine à permutação implícita em c
                auto &c = population[i];
                // decodifica prior -> seq
                std::vector<std::pair<int,int>> order;
                for (int j = 0; j < NUM_OF_CUSTOMERS; ++j)
                    order.emplace_back(c[j+1], j+1);
                std::sort(order.begin(), order.end());
                std::vector<int> seq;
                seq.reserve(NUM_OF_CUSTOMERS);
                for (auto &p : order) seq.push_back(p.second);

                // reaplica split_and_refine
                int rawK = c[0];
                int K = 1;
                for (int b = 0; b < NUM_OF_CUSTOMERS; ++b)
                    if ((rawK >> (NUM_OF_CUSTOMERS-1-b)) & 1) ++K;
                auto [newCost, newRoutes] = split_and_refine(seq, K);

                // recalcula rawK da mesma forma
                std::vector<int> mask(NUM_OF_CUSTOMERS,0);
                int idx=0;
                for(size_t r=0; r<newRoutes.size(); ++r){
                  for(int cust: newRoutes[r]){
                    while(seq[idx]!=cust) ++idx;
                    if (idx+1<NUM_OF_CUSTOMERS && r+1<newRoutes.size()
                        && seq[idx+1]==newRoutes[r+1][0]) mask[idx]=1;
                    ++idx;
                  }
                }
                int newRawK=0;
                for(int b=0;b<NUM_OF_CUSTOMERS;++b){
                  newRawK=(newRawK<<1)|mask[b];
                }
                c[0]=newRawK;
                best_routes = newRoutes;
                fitness_vals[i] = newCost;
            }
        }
    }
    // final: atualiza current_best e best_routes
    int bi = std::min_element(fitness_vals.begin(), fitness_vals.end()) - fitness_vals.begin();
    extern double current_best;
    current_best = fitness_vals[bi];
    decode_and_evaluate(population[bi]);
}

void print_best_solution(){
    std::ofstream out("solution.txt");
    if(!out){std::cerr<<"Erro abrir solution.txt";return;} 
    for(size_t i=0;i<best_routes.size();++i){
        out<<"Vehicle "<<i+1<<": 0";
        for(int c:best_routes[i]) out<<" -> "<<c;
        out<<" -> 0\n";
    }
}

void print_solution_for_k(int rawK){
    int n=NUM_OF_CUSTOMERS;
    std::vector<int> mask(n);
    for(int i=0;i<n;++i) mask[n-1-i]=(rawK>>i)&1;
    Chromosome best_c; double bf=std::numeric_limits<double>::infinity();
    for(auto &c:population) if(c[0]==rawK){double f=decode_and_evaluate(c);if(f<bf){bf=f;best_c=c;}}
    std::vector<std::pair<int,int>>order; order.reserve(n);
    for(int j=0;j<n;++j) order.emplace_back(best_c[j+1],j+1);
    std::sort(order.begin(),order.end()); std::vector<int>seq; seq.reserve(n);
    for(auto&p:order) seq.push_back(p.second);
    std::vector<std::vector<int>>routes(1);
    for(int i=0;i<n;++i){routes.back().push_back(seq[i]); if(mask[i]&&i<n-1)routes.emplace_back();}
    std::ofstream out2("solution_k_"+std::to_string(rawK)+".txt"); if(!out2)return;
    for(size_t i=0;i<routes.size();++i){out2<<"Vehicle "<<i+1<<": 0"; for(int c:routes[i])out2<<" -> "<<c; out2<<" -> 0\n";}    
}
