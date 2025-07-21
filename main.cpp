#include <iostream>
#include <stdlib.h>
#include <limits.h>
#include <algorithm>

#include "EVRP.hpp"
#include "stats.hpp"
#include "genetic.hpp"
#include <vector>

using namespace std;

/*initializes a run for your heuristic*/
void start_run(int r) {
    srand(r);
    init_evals();
    init_current_best();
    cout << "Run: " << r << " with random seed " << r << endl;
}

/*gets an observation of the run for your heuristic*/
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

    problem_instance = argv[1];
    read_problem(problem_instance);
    open_stats();

    for (int run = 1; run <= MAX_TRIALS; run++) {
        start_run(run);

        initialize_genetic(NUM_OF_CUSTOMERS);
        run_genetic();

        // Imprime rota ótima encontrada
        cout << "Best route for run " << run << ":" << endl;
        print_best_solution();

        // Imprime rota para o mesmo K do melhor indivíduo
        int bestIdx = min_element(fitness_vals.begin(), fitness_vals.end()) - fitness_vals.begin();
        int bestK = population[bestIdx][0];
        cout << "Route for rawK = " << bestK << ":" << endl;
        print_solution_for_k(bestK);
        cout << "Route for rawK = " << 3 << ":" << endl;
      //  print_solution_for_k(3);
        end_run(run);
    }

    close_stats();
    free_stats();
    free_EVRP();
    return 0;
}
