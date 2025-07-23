#include <iostream>
#include <stdlib.h>
#include <limits.h>
#include <algorithm>
#include <chrono>

#include "EVRP.hpp"
#include "stats.hpp"
#include "genetic.hpp"
#include <vector>

using namespace std;
using Clock = std::chrono::high_resolution_clock;

/* inicializa uma execução */
void start_run(int r) {
    srand(r);
    init_evals();
    init_current_best();
    cout << "Run: " << r << " with random seed " << r << endl;
}

/* finaliza uma execução */
void end_run(int r) {
    // armazena o best daquela run
    get_mean(r-1, get_current_best());
    cout << "End of run " << r
         << " with best solution quality " << get_current_best()
         << " total evaluations: " << get_evals() << endl;
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

        // começa a contagem de tempo
        auto t0 = Clock::now();

        initialize_genetic(NUM_OF_CUSTOMERS);
        run_genetic();

        // fim da contagem
        auto t1 = Clock::now();
        double elapsed_s = chrono::duration<double>(t1 - t0).count();

        // imprime apenas o melhor fitness e o tempo
        double best_fitness = get_current_best();
        cout << "Best fitness for run " << run << ": "
             << best_fitness << endl;
        cout << "Time for run " << run << ": "
             << elapsed_s << " s" << endl << endl;

        end_run(run);
    }

    close_stats();
    free_stats();
    free_EVRP();
    return 0;
}
