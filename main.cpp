// main.cpp
#include <iostream>
#include <stdlib.h>
#include <limits.h>
#include <algorithm>

#include "EVRP.hpp"
#include "stats.hpp"
#include "genetic.hpp"

using namespace std;

/* inicializa uma execução */
void start_run(int r) {
    srand(r);
    init_evals();
    init_current_best();
    cout << "Run: " << r << " with random seed " << r << endl;
}

/* registra os resultados da execução */
void end_run(int r) {
    get_mean(r-1, get_current_best());
    cout << "End of run " << r
         << " with best solution quality " << get_current_best()
         << " total evaluations: " << get_evals() << endl << endl;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        cout << "Usage: ./main <instance_file.evrp>" << endl;
        return 1;
    }

    // Lê e monta o problema
    problem_instance = argv[1];
    read_problem(problem_instance);

    // Abre estatísticas
    open_stats();

    for (int run = 1; run <= MAX_TRIALS; ++run) {
        start_run(run);

        // Inicializa e roda o GA
        initialize_genetic(NUM_OF_CUSTOMERS);
        set_target_vehicles( 3, 10000 );
        run_genetic();

        // Imprime a melhor rota encontrada nesta execução
        cout << "\nBest route for run " << run << ":\n";
        print_best_solution();

        // Captura o rawK do melhor indivíduo para reproduzir aquela divisão
      int bestIdx = min_element(fitness_vals.begin(), fitness_vals.end())
                      - fitness_vals.begin();
        int bestK = population[bestIdx][0].first;  // <-- extrai .first do pair

        cout << "\nRoute for rawK = " << bestK << ":\n";
        print_best_solution();
        end_run(run);
    }

    // Finaliza estatísticas e libera tudo
    close_stats();
    free_stats();
    free_EVRP();
    return 0;
}
