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

double evaluate_distance(const std::vector<int> &route) {
    double total = 0.0;
    for (size_t i = 0; i + 1 < route.size(); ++i)
        total += distances[route[i]][route[i+1]];
    return total;
}

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
    // reconstruct cuts
    std::vector<int> cuts;
    for (int idx = n; idx > 0; idx = prev[idx])
        cuts.push_back(idx);
    std::reverse(cuts.begin(), cuts.end());
    // build routes
    std::vector<std::vector<int>> R;
    int start = 0;
    for (int cut : cuts) {
        R.emplace_back(o.begin() + start, o.begin() + cut);
        start = cut;
    }
    // ensure mv routes by splitting largest
    while ((int)R.size() < mv) {
        int wi = 0, maxsz = 0;
        for (int i = 0; i < (int)R.size(); ++i) {
            if ((int)R[i].size() > maxsz) {
                maxsz = R[i].size(); wi = i;
            }
        }
        auto r = R[wi];
        R.erase(R.begin() + wi);
        int mid = r.size() / 2;
        R.insert(R.begin() + wi, std::vector<int>(r.begin(), r.begin() + mid));
        R.insert(R.begin() + wi + 1, std::vector<int>(r.begin() + mid, r.end()));
    }
    // compute total cost
    double total = 0.0;
    for (auto &r : R) {
        std::vector<int> full = {DEPOT};
        full.insert(full.end(), r.begin(), r.end());
        full.push_back(DEPOT);
        total += evaluate_distance(full);
    }
    return {total, R};
}

auto two_opt = [](std::vector<int> &route) {
    bool improved = false;
    int n = route.size();
    do {
        improved = false;
        for (int i = 1; i < n - 2; ++i) {
            for (int k = i + 1; k < n - 1; ++k) {
                double delta =
                    distances[route[i-1]][route[k]] +
                    distances[route[i]][route[k+1]] -
                    (distances[route[i-1]][route[i]] +
                     distances[route[k]][route[k+1]]);
                if (delta < -1e-6) {
                    std::reverse(route.begin() + i, route.begin() + k + 1);
                    improved = true;
                }
            }
        }
    } while (improved);
};




static std::pair<double, std::vector<std::vector<int>>> split_and_refine(
    const std::vector<int> &seq) {
    int n = seq.size();
    // custo original sem split
    std::vector<int> full = {DEPOT}; full.insert(full.end(), seq.begin(), seq.end()); full.push_back(DEPOT);
    double best_cost = evaluate_distance(full);
    std::vector<std::vector<int>> best_routes_loc(1, seq);
    // testa mv=2..n
    for (int mv = 2; mv <= n; ++mv) {
        auto [cst, R] = splitRoutes(seq, mv);
        if (cst < best_cost - 1e-6) {
            best_cost = cst;
            best_routes_loc = std::move(R);
        }
    }
    // aplica 2-opt
    for (auto &r: best_routes_loc) {
        std::vector<int> tmp = {DEPOT}; tmp.insert(tmp.end(), r.begin(), r.end()); tmp.push_back(DEPOT);
        two_opt(tmp);
        r.assign(tmp.begin()+1, tmp.end()-1);
    }
    return {best_cost, best_routes_loc};
}


double decode_and_evaluate(const Chromosome &chrom) {
    int n = NUM_OF_CUSTOMERS;
    // 1) Reconstroi a sequência de clientes a partir dos genes de prioridade
    std::vector<std::pair<int,int>> order;
    order.reserve(n);
    for (int j = 0; j < n; ++j)
        order.emplace_back(chrom[j+1], j+1);
    std::sort(order.begin(), order.end());
    std::vector<int> seq;
    seq.reserve(n);
    for (auto &p : order)
        seq.push_back(p.second);

    // 2) Aplica split_and_refine sem passar K — ele testará todas as partições
    auto [cost, routes] = split_and_refine(seq);

    // 3) Atualiza best_routes para impressão
    best_routes = routes;

    // 4) Retorna o custo refinado
    return cost;
}


void initialize_genetic(int num_customers) {
    int n = num_customers;
    population.resize(POP_SIZE);
    fitness_vals.resize(POP_SIZE);
    std::uniform_int_distribution<> dp(0, n*n);
    for (int i = 0; i < POP_SIZE; ++i) {
        Chromosome c(1+n);
        // prioridades aleatórias
        for (int j = 1; j <= n; ++j)
            c[j] = dp(rng);
        // reconstrói sequência
        std::vector<std::pair<int,int>> ord;
        for (int j = 0; j < n; ++j) ord.emplace_back(c[j+1], j+1);
        std::sort(ord.begin(), ord.end());
        std::vector<int> seq;
        for (auto &p: ord) seq.push_back(p.second);
        // refinamento
        auto [cost, R] = split_and_refine(seq);
        best_routes = R;
        fitness_vals[i] = cost;
        // constrói mask
        std::vector<int> mask(n, 0);
        int idx = 0;
        for (size_t r = 0; r < R.size(); ++r) {
            for (int cust : R[r]) {
                while (seq[idx] != cust) ++idx;
                if (idx+1 < n && r+1 < R.size() && seq[idx+1] == R[r+1][0])
                    mask[idx] = 1;
                ++idx;
            }
        }
        int rawK = 0;
        for (int b = 0; b < n; ++b) rawK = (rawK<<1) | mask[b];
        // garante ao menos um corte
        if (rawK == 0) rawK = 1 << (n/2);
        c[0] = rawK;
        population[i] = c;
        if (i < 3)
            std::cout << "Init#"<<i<<" cost="<<cost<<" routes="<<R.size() <<" rawK="<<rawK<<"\n";
    }
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
    const int LOCAL_EVERY = 20;
    for (int gen = 1; gen <= NUM_GENERATIONS; ++gen) {
        // Elitismo e geração de nova pop
        std::vector<Chromosome> next_pop;
        std::vector<double> next_fit;
        int bi = std::min_element(fitness_vals.begin(), fitness_vals.end()) - fitness_vals.begin();
        next_pop.push_back(population[bi]);
        next_fit.push_back(fitness_vals[bi]);
        while ((int)next_pop.size() < POP_SIZE) {
            Chromosome p1 = tournament_selection();
            Chromosome p2 = tournament_selection();
            auto [c1, c2] = crossover(p1,p2);
            mutate(c1); mutate(c2);
            next_fit.push_back(decode_and_evaluate(c1)); next_pop.push_back(c1);
            if ((int)next_pop.size() < POP_SIZE) {
                next_fit.push_back(decode_and_evaluate(c2)); next_pop.push_back(c2);
            }
        }
        population.swap(next_pop);
        fitness_vals.swap(next_fit);
        // refinamento periódico
        if (gen % LOCAL_EVERY == 0) {
            for (int i = 0; i < POP_SIZE; ++i) {
                Chromosome &c = population[i];
                std::vector<std::pair<int,int>> ord;
                for (int j=0;j<NUM_OF_CUSTOMERS;++j) ord.emplace_back(c[j+1], j+1);
                std::sort(ord.begin(), ord.end());
                std::vector<int> seq;
                for (auto &p:ord) seq.push_back(p.second);
                auto [cost, R] = split_and_refine(seq);
                best_routes = R;
                fitness_vals[i] = cost;
                // rebuild rawK
                std::vector<int> mask(NUM_OF_CUSTOMERS,0);
                int idx=0;
                for (size_t r=0;r<R.size();++r) for (int cust:R[r]){
                    while(seq[idx]!=cust) ++idx;
                    if(idx+1<NUM_OF_CUSTOMERS && r+1<R.size() && seq[idx+1]==R[r+1][0]) mask[idx]=1;
                    ++idx;
                }
                int rawK=0;
                for(int b=0;b<NUM_OF_CUSTOMERS;++b) rawK=(rawK<<1)|mask[b];
                if(rawK==0) rawK=1<<(NUM_OF_CUSTOMERS/2);
                c[0]=rawK;
            }
        }
    }
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
